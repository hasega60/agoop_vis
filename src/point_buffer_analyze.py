import pandas as pd
import geopandas as gpd
import psycopg2
import jismesh.utils as ju
from shapely.geometry import LineString as LineString_shapely
from geojson import Point, LineString, Feature, FeatureCollection
from scipy.spatial import ConvexHull
from tqdm import tqdm
from datetime import datetime
import math
import h3
import os


# 入力csvの配置
#input_path = "../data/city_from_merge_8220_week3"
input_path = "../data/chiba/202004"
buffer_path = "../data/shape/station_2019_buffer_tokyo_JREast.gpkg"
# 出力先
output_dir = "../result/"


# オプションデータがあるか
has_option = False


semiMajorAxis = 6378137.0  # 赤道半径
flattening = 1 / 298.257223563  # 扁平率
e_2 = flattening * (2 - flattening)
trip_interval_second = 10  # tripに記録する秒数の間隔
degree = math.pi / 180
line_color_rgb = (255, 120, 120) #geojsonのlineRGB


def distance_latlon(x1, y1, x2, y2):
    coslat = math.cos((y1 + y2) / 2 * degree)
    w2 = 1 / (1 - e_2 * (1 - coslat * coslat))
    dx = (x1 - x2) * coslat
    dy = (y1 - y2) * w2 * (1 - e_2)
    return math.sqrt((dx * dx + dy * dy) * w2) * semiMajorAxis * degree


if __name__ == '__main__':

    # 読み込むファイル一覧
    csvfiles = []
    for f in os.listdir(input_path):
        if not os.path.isdir(input_path + "/" + f):
            if os.path.splitext(f)[1] == ".csv":
                csvfiles.append(input_path + "/" + f)


    file = csvfiles[0]
    gdf_buffer = gpd.read_file(buffer_path)
    gdf_buffer["station_id"] = gdf_buffer.index
    buffer_cols = list(gdf_buffer.columns)
    pbar = tqdm(csvfiles)
    df_city_stay_h3_group_out, df_mesh_group_out,df_did_convex_out_all = None, None, None
    for file in pbar:
        pbar.set_description(f"Processing..{os.path.basename(file)}")
        file_day = os.path.splitext(os.path.basename(file))[0].split("_")[2]
        df_all = pd.read_csv(file)
        if has_option:
            df_all = df_all[["dailyid", "year", "month", "day", "dayofweek", "hour", "minute", "latitude", "longitude",
                             "os","logtype_subcategory","accuracy","speed","estimated_speed_flag","course",
                             "prefcode", "citycode", "home_prefcode","home_citycode",
                             "workplace_prefcode","workplace_citycode","transportation_type", "gender"]]
        else:
            df_all = df_all[
                ["dailyid", "year", "month", "day", "dayofweek", "hour", "minute", "latitude", "longitude",
                 "os", "logtype_subcategory", "accuracy", "speed", "estimated_speed_flag", "course",
                 "prefcode", "citycode"]]
        #なぜか日付がすべて1になるので，dayofweekで代用
        df_time = df_all[["year", "month", "day", "hour", "minute"]]
        #df_time = df_all[["year", "month", "dayofweek", "hour", "minute"]]
        df_time=df_time.rename(columns={"dayofweek":"day"})
        df_all["time"] = pd.to_datetime(df_time)
        # daylyid別に，時間でソート
        df_all = df_all.sort_values(['dailyid', 'time'])
        points_out, lines_out = [], []

        for buffer in gdf_buffer.values:
            buffer_geom = buffer[buffer_cols.index("geometry")]
            station_id = buffer[buffer_cols.index("station_id")]
            station_name = buffer[buffer_cols.index("N05_011")]

            min_lon, min_lat, max_lon, max_lat = buffer_geom.envelope.bounds
            # 先にEnvelopeで抽出した後に，バッファで空間検索
            df_point_in_buffer = df_all[(df_all["latitude"].between(min_lat,max_lat))&(df_all["longitude"].between(min_lon,max_lon))]
            gdf_point_in_buffer = gpd.GeoDataFrame(df_point_in_buffer, geometry=gpd.points_from_xy(df_point_in_buffer["longitude"], df_point_in_buffer["latitude"]))
            gdf_point_in_buffer.set_crs(epsg=4612)

            gdf_point_in_buffer = gdf_point_in_buffer[gdf_point_in_buffer.within(buffer_geom)]

            # dailyid別にラインを作成
            did_in_buffer_list = gdf_point_in_buffer["dailyid"].unique()
            gdf_point_in_buffer["station_id"]=int(station_id)
            gdf_point_in_buffer["station_name"]=station_name
            df_point_in_buffer_m = gdf_point_in_buffer[["station_id","station_name"]]
            df_did_in_buffer = df_all[df_all["dailyid"].isin(did_in_buffer_list)]
            df_did_in_buffer = pd.merge(df_did_in_buffer, df_point_in_buffer_m, how="left", left_index=True, right_index=True)
            df_did_cols = list(df_did_in_buffer.columns)
            # 二点以上，1分以上滞在しているtripを抽出
            for did in did_in_buffer_list:
                time_before, lat_before, lon_before = -1, -1, -1
                df_did = df_did_in_buffer[df_did_in_buffer["dailyid"]==did]

                trip_count = 0  # daylyidのトリップ番号 一旦外に出た時対応
                dist = 0        # 移動距離
                stay_time = 0   # 滞在時間
                first_time = -1  # 進入時間
                add_point = False
                path = []
                for v in df_did.values:
                    buffer_station_id = v[df_did_cols.index("station_id")]
                    if pd.isna(buffer_station_id) and time_before==-1:
                        # buffer侵入前
                        continue

                    time_value= v[df_did_cols.index("time")]
                    time = int(time_value.timestamp())
                    lat = v[df_did_cols.index("latitude")]
                    lon = v[df_did_cols.index("longitude")]

                    if time_before == -1:
                        time_before=time
                        lat_before=lat
                        lon_before=lon
                        first_time=time
                        continue

                    delta_t = time - time_before
                    if delta_t == 0:
                        # 時間が経過していない
                        continue

                    if pd.isna(buffer_station_id):
                        # バッファの外に出た
                        # buffer内に二点以上記録されていたら,点を追加
                        if add_point:
                            points_out.append(
                                [station_id, station_name, did, trip_count, lat_before, lon_before, time_before, dist,
                                 stay_time])
                            path.append([lon_before, lat_before])
                            # ラインに追加
                            linestring = LineString_shapely(path)
                            lines_out.append(
                                [station_id, station_name, did, trip_count, linestring, first_time, dist, stay_time])

                            path = []
                            trip_count += 1
                            add_point = False

                        # 初期に戻す
                        first_time,time_before, lat_before, lon_before = -1, -1, -1, -1
                        dist,stay_time = 0,0
                        continue

                    else:
                        # 時間経過があるなら追加
                        add_point=True
                        points_out.append(
                            [station_id, station_name, did, trip_count, lat_before, lon_before, time_before, dist,
                             stay_time])
                        path.append([lon_before, lat_before])



                    time_before = time
                    dist += distance_latlon(lon_before, lat_before, lon, lat)
                    lat_before = lat
                    lon_before = lon
                    stay_time = time - first_time

                if len(path)>1:
                    linestring = LineString_shapely(path)
                    lines_out.append([station_id, station_name, did, trip_count,linestring, first_time, dist, stay_time])


        df_points_out = pd.DataFrame(points_out, columns=["station_id", "station_name", "dailyid", "trip_count", "lat",
                                                          "lon", "timestamp", "dist","stay_second"])

        path = f"{output_dir}/point_buffer_analyze/{file_day}_point_in_buffer.csv"
        df_points_out.to_csv(path, index=False)

        df_lines_out = pd.DataFrame(lines_out, columns=["station_id", "station_name", "dailyid", "trip_count", "geometry",
                                                          "first_timestamp", "dist", "stay_second"])
        gdf_lines_out = gpd.GeoDataFrame(df_lines_out, geometry=df_lines_out["geometry"])
        gdf_lines_out = gdf_lines_out.set_crs(epsg=4612)

        path = f"{output_dir}/point_buffer_analyze/{file_day}_line_in_buffer.gpkg"
        gdf_lines_out.to_file(path, layer="line_in_buffer", driver="GPKG")

