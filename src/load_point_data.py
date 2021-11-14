import pandas as pd
import geopandas as gpd
import psycopg2
import jismesh.utils as ju
from shapely.geometry import Polygon
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
# 出力先
output_dir = "../result//service_area"

# アニメーション設定項目
# 集計期間
#start_datetime='2019-05-01 6:00:00'
#end_datetime='2019-05-06 0:00:00'


start_datetime=None
end_datetime=None

#city = "tsukuba"
city = "chiba_all"

# 距離分布集計用
hist_step = 500
max_hist_step = hist_step * 100

max_speed = 200

# 自宅判定時間
home_hour = 4

city_dict = {
    "tsukuba": 8220,
    "tsuchiura": 8203,
    "kashiwa": 12217
}
# オプションデータがあるか
has_option = False

#集計するODを集計するメッシュレベル
meshlevel = 4
h3_index_level = 9

#元データポイントの出力
export_org_point = False
# 出力プリント文の表示
print_export_string= False

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



def dataframe_to_trip(df:pd.DataFrame, trip_id_col, lat_col, lon_col, time_col,altitude_col=None, max_speed=150, add_col=None):
    def _add_output(path, trip, trip_count, trip_origin_destination, trip_out, id_out, id_out_org, add_value_out=[], add_value=None):
        od = path[0]
        d = path[-1]
        od.extend(d)
        trip_origin_destination.append(od)
        trip_out.append(path_to_trip_geojson(path))
        id_out.append(f"{trip}_{trip_count}")
        id_out_org.append(trip)
        if add_value is not None:
            add_value_out.append(add_value)

        return trip_origin_destination, trip_out, id_out, id_out_org, add_value_out

    trip_out, id_out, id_out_org, add_value_out= [],[],[],[]
    trip_origin_destination=[]
    cols = list(df.columns)
    df_trip_id = df[trip_id_col]
    trip_id = df_trip_id.unique()

    for trip in trip_id:
        points = df[df[trip_id_col]==trip].values
        path = []
        trip_count=0
        time_before, lat_before, lon_before = -1, -1, -1

        # 任意追加属性
        add_value = None
        if add_col is not None:
            add_value = points[0][cols.index(add_col)]

        for point in points:
            time_jst= int(point[cols.index(time_col)].timestamp())
            """
            time_unix = int(point[cols.index(time_col)].timestamp())
            # unix_timeだと標準時でtrip表示されてしまうので，強引にjst変換
            utc_offset = datetime.fromtimestamp(time_unix) - datetime.utcfromtimestamp(time_unix)
            time_jst = time_unix + utc_offset.seconds
            """
            lat = point[cols.index(lat_col)]
            lon = point[cols.index(lon_col)]
            elevation = 0
            if altitude_col is not None:
                # 高さを考慮するとき
                elevation = point[cols.index(altitude_col)]

            t = time_jst

            if time_before != -1:
                delta_t = t - time_before
                distance = distance_latlon(lon_before, lat_before, lon, lat)
                speed = 0
                if delta_t > 0:
                    speed = (distance/1000) / (delta_t/3600)
                lat_before = lat
                lon_before = lon
                time_before = t
                if delta_t < 600 and speed < max_speed and speed > 0:
                    # 10分以内かつ二点間の速度が max_speed km/h以下，移動しているならpathとして追加
                    path.append([lon, lat, elevation, t])
                else:
                    # それ以外なら別のpathにする
                    if len(path) > 1:
                        trip_origin_destination, trip_out, id_out, id_out_org, add_value_out=_add_output(path, trip, trip_count,
                                                                              trip_origin_destination,
                                                                              trip_out,
                                                                              id_out, id_out_org,add_value_out,add_value)
                        trip_count += 1
                        path=[]
            else:
                path.append([lon, lat, elevation, t])
                time_before = t
                lat_before = lat
                lon_before = lon

        if len(path) > 1:
            trip_origin_destination, trip_out, id_out, id_out_org, add_value_out = _add_output(path, trip, trip_count,
                                                                                trip_origin_destination,
                                                                                trip_out,
                                                                                id_out, id_out_org,add_value_out,add_value)

    return trip_out, id_out, id_out_org, trip_origin_destination, add_value_out

def rgb2html(rgb):
    R = rgb[0]
    G = rgb[1]
    B = rgb[2]
    color_code = '#{}{}{}'.format(hex(R), hex(G), hex(B))
    return color_code.replace('0x', '')


def path_to_trip_geojson(path):
    feas = []
    color = rgb2html(line_color_rgb)
    feas.append((LineString(path), {"color": color, }))
    return [Feature(geometry=g, properties=p) for g, p in feas]

def output_trip_csv(trip_feature_list:[], user_cols_dict={}, path="trips.csv"):
    """
    :param trip_feature_list: [[{"geometry": {"coordinates": [[139.723445, 35.747456, 71.705, 1597278916], ...], "type": "LineString"}, "properties": {"color": "#ff7878"}, "type": "Feature"}]]
    :param user_cols_dict: out_listに任意の属性がある場合に追加
    :param path: 出力パス
    :return:
    """

    df_out = pd.DataFrame(trip_feature_list, columns=["_geojson"])
    if len(user_cols_dict) > 0:
        for k, v in user_cols_dict.items():
            df_out[k] = v

    df_out.to_csv(path, index=False)
    if print_export_string:
        print(f"export trip animation data:{path}")

if __name__ == '__main__':
    try:
        # 読み込むファイル一覧
        csvfiles = []
        for f in os.listdir(input_path):
            if not os.path.isdir(input_path + "/" + f):
                if os.path.splitext(f)[1] == ".csv":
                    csvfiles.append(input_path + "/" + f)

        file = csvfiles[0]
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
            df_all = df_all.sort_values(['dailyid', 'time'])


            #Androidだけ対象にする
            #df_all = df_all[df_all["os"] == "Android"]
            # 都市でフィルタ (dictになければすべて)
            if city in (city_dict.keys()):
                if has_option:
                    df_city = df_all[df_all["home_citycode"]==city_dict[city]]
                else:
                    df_city = df_all[df_all["citycode"]==city_dict[city]]
            else:
                df_city = df_all

            # ひっかからなかったらskip
            if len(df_city) == 0:
                continue

            if start_datetime is not None:
                df_city = df_city[df_city['time'] >= datetime.strptime(start_datetime, '%Y-%m-%d %H:%M:%S')]

            if end_datetime is not None:
                df_city = df_city[df_city['time'] <= datetime.strptime(end_datetime, '%Y-%m-%d %H:%M:%S')]
            # 都市別ポイント抽出
            if export_org_point:
                path = f"{output_dir}/{city}/point.csv"
                df_city.to_csv(path)
                if print_export_string:
                   print(f"export city point data:{path}")


            # h3_indexを付与
            def get_h3_index(row):
               return h3.geo_to_h3(row.latitude, row.longitude, h3_index_level)

            def get_meshcode(row):
               return ju.to_meshcode(row.latitude, row.longitude, meshlevel)

            # TODO
            # 自宅場所の抽出
            df_city_home = df_city.copy()

            # 規定時間以降に，最初に記録された場所のメッシュにする
            df_city_home["mesh_code"] = df_city_home.apply(get_meshcode, axis=1)
            df_city_home_f = df_city_home[df_city_home["hour"]>=home_hour]
            df_city_home_f = df_city_home_f.groupby(["dailyid"], as_index=False).first()

            df_city_home_f["home_flg"] = True

            df_city_home_f = df_city_home_f[["dailyid", "mesh_code", "home_flg"]]
            df_city_home_m = df_city_home.merge(df_city_home_f, how="left", on=["dailyid", "mesh_code"])

            df_city_home_f = df_city_home_f[["dailyid", "mesh_code"]]

            from_mesh_list = df_city_home_f["mesh_code"].unique()
            df_city_home_f=df_city_home_f.rename(columns={"mesh_code":"home_mesh_code"})
            mesh_str_list=[]
            for v in df_city_home_f["home_mesh_code"].values:
                mesh_str_list.append(f"m{v}")
            df_city_home_f["home_mesh_code"] = mesh_str_list
            df_city_home_m= df_city_home_m.merge(df_city_home_f, how="left", on=["dailyid"])
            df_city = df_city_home_m

            #dailyidごとのconvex hull, 最遠点を作成
            out_did_convex, convex_geos=[],[]
            did_list = df_city_home_f["dailyid"].unique()
            for did in did_list:
                df_did = df_city[df_city["dailyid"]==did].reset_index()
                df_did_home = df_did[df_did["home_flg"] == True]
                from_mesh_code = df_did_home["mesh_code"].values[0]
                lat_home, lon_home = ju.to_meshpoint(from_mesh_code, 0.5, 0.5)

                # convexhull
                points = df_did[["longitude", "latitude"]].values
                # 最遠点
                farthest, dist_list = [None, None], []
                dist = 0
                polygon = None
                try:
                    hull = ConvexHull(points)
                    df_did_hull = df_did[df_did.index.isin(hull.vertices)]
                    if len(df_did_hull)==0:
                        continue

                    for lon, lat in df_did_hull[["longitude", "latitude"]].values:
                        dist_temp = distance_latlon(lon_home, lat_home, lon, lat)
                        dist_list.append(dist_temp)
                        if dist < dist_temp:
                            farthest = [lon, lat]
                            dist = dist_temp

                    polygon = Polygon(hull.points).convex_hull

                except:
                    continue


                out_did_convex.append([
                    did,
                    from_mesh_code,
                    lon_home,
                    lat_home,
                    farthest[0],
                    farthest[1],
                    dist,
                    dist_list])
                convex_geos.append(polygon)


            df_did_convex = pd.DataFrame(out_did_convex,columns=["dailyid","from_mesh_code", "from_lon", "from_lat", "to_lon","to_lat",
                                                                 "distance_far","convex_points_dist"])
            g_convex_geos = gpd.GeoSeries(convex_geos)
            gdf_pol = gpd.GeoDataFrame(g_convex_geos, columns=["geometry"])

            gdf_did_convex = gpd.GeoDataFrame(df_did_convex)
            gdf_did_convex["geometry"] = convex_geos
            gdf_did_convex.set_crs(epsg=4612)
            gdf_did_convex_out = gdf_did_convex[["dailyid","from_mesh_code", "from_lon", "from_lat", "to_lon","to_lat","distance_far","geometry"]]
            gdf_did_convex_out.to_file(f"{output_dir}//{file_day}_dailyid_convex.gpkg",driver="GPKG")
            df_did_convex_out = gdf_did_convex_out[["dailyid","from_mesh_code", "from_lon", "from_lat", "to_lon","to_lat","distance_far"]]

            mesh_str_list=[]
            for v in df_did_convex_out["from_mesh_code"].values:
                mesh_str_list.append(f"m{v}")

            df_did_convex_out["from_mesh_str"] = mesh_str_list
            df_did_convex_out["day"]=file_day

            if df_did_convex_out_all is None:
                df_did_convex_out_all = df_did_convex_out
            else:
                df_did_convex_out_all = pd.concat([df_did_convex_out_all, df_did_convex_out])
            #df_did_convex_out.to_csv(f"{output_dir}//{file_day}_mesh_distance_far.csv", index=False)

            # メッシュごとのconvex_hull polygon, 距離分布用distance_listを作成する
            mesh_convex_list, dist_list_from_mesh = [],[]
            user_count=0
            for mesh in from_mesh_list:
                gdf_mesh_convex = gdf_did_convex[gdf_did_convex["from_mesh_code"]==mesh]
                user_count = len(gdf_mesh_convex)
                if user_count>=2:
                    p_list = gdf_mesh_convex["geometry"].values
                    points=[]
                    for p_m in p_list:
                        coords = p_m.boundary.coords
                        for lon, lat in coords:
                            points.append([lon, lat])

                    hull = ConvexHull(points)
                    polygon = Polygon(hull.points).convex_hull
                elif user_count==1:
                    polygon = gdf_mesh_convex["geometry"].values[0]
                else:
                    continue

                mesh_convex_list.append([mesh,user_count, polygon])

                did_convex_points_dist = gdf_mesh_convex[["dailyid","convex_points_dist"]].values
                for did,dist_l in did_convex_points_dist:
                    for dist in dist_l:
                        dist_list_from_mesh.append([mesh, did, dist])

            df_mesh_convex_dist = pd.DataFrame(dist_list_from_mesh, columns=["mesh_code","dailyID","convex_points_dist"])
            df_mesh_convex_dist.to_csv(f"{output_dir}//{file_day}_mesh_distance.csv", index=False)


            df_mesh_convex_hull = pd.DataFrame(mesh_convex_list, columns=["mesh_code","user_count", "geometry"])
            gdf_mesh_convex_hull = gpd.GeoDataFrame(df_mesh_convex_hull, geometry=df_mesh_convex_hull["geometry"])
            gdf_mesh_convex_hull.set_crs(epsg=4612)


            gdf_mesh_convex_hull.to_file(f"{output_dir}//{file_day}_mesh_convex.gpkg", driver="GPKG")

            # 滞在地点の抽出
            """
            df_city_stay = df_city[(df_city["logtype_subcategory"] == "timer")|
                                   (df_city["logtype_subcategory"] == "􏸛􏷅􏷅􏸁􏸧􏸛􏸘􏸛􏷅􏷅􏸁􏸧􏸛􏸘arrival")|
                                   (df_city["speed"] <= 2)]
            """
            df_city_stay = df_city[(df_city["logtype_subcategory"] == "timer")|
                                   (df_city["logtype_subcategory"] == "􏸛􏷅􏷅􏸁􏸧􏸛􏸘􏸛􏷅􏷅􏸁􏸧􏸛􏸘arrival")]

            #30分単位のカテゴリ付与
            def get_hhb(row):
                return int((row["hour"] * 60 + row["minute"])/30)

            df_city_stay["half_hour_bin"] = df_city_stay.apply(get_hhb, axis=1)

            # 30分単位でのtrip重複除去
            df_city_stay_group = df_city_stay.groupby(["dailyid","dayofweek", "half_hour_bin"], as_index=False).first()
            path = f"{output_dir}/{city}/stay_hhb.csv"
            df_city_stay_group.to_csv(path)
            if print_export_string:
                print(f"export stay point data:{path}")

            pbar.set_description(f"h3_index...{os.path.basename(file)}")


            df_city_stay_group["h3_index"] = df_city_stay_group.apply(get_h3_index, axis=1)
            df_city_stay_h3_group = df_city_stay_group.groupby(["h3_index", "dayofweek", "half_hour_bin"], as_index=False).count()
            df_city_stay_h3_group = df_city_stay_h3_group[["h3_index", "dayofweek", "half_hour_bin", "dailyid"]]
            df_city_stay_h3_group = df_city_stay_h3_group.rename(columns={"dailyid": "count"})

            if df_city_stay_h3_group_out is None:
                df_city_stay_h3_group_out = df_city_stay_h3_group
            else:
                df_city_stay_h3_group_out = df_city_stay_h3_group_out.append(df_city_stay_h3_group)

            if print_export_string:
                print(f"export stay h3 data:{path}")

            # trip animationデータの抽出
            pbar.set_description(f"creating trip...{os.path.basename(file)}")
            # 停止時に累積されるポイント（timer）は除去
            df_city_trip = df_city[df_city["logtype_subcategory"]!="timer"]
            trip_list, id_list, id_out_org, trip_origin_destination, add_value_out = dataframe_to_trip(df_city_trip, trip_id_col="dailyid",
                                          lat_col="latitude", lon_col="longitude",time_col="time", max_speed=max_speed, add_col="home_mesh_code")

            id_col_dict={"dailyid_c": id_list, "dailyid": id_out_org, "home_mesh_code":add_value_out}
            output_trip_csv(trip_feature_list=trip_list, user_cols_dict=id_col_dict, path=f"{output_dir}/{city}/trips_{file_day}.csv")
            df_od = pd.DataFrame(trip_origin_destination, columns=["lon_from", "lat_from", "elevation_from", "time_from",
                                                                   "lon_to", "lat_to", "elevation_to", "time_to"])

            for k, v in id_col_dict.items():
                df_od[k] = v

            # 時間属性を付与
            df_od["hour"] = df_od['time_from'].apply(lambda x: datetime.fromtimestamp(x).hour)
            df_od_grouping = df_od[["dailyid","lon_from", "lat_from","lon_to", "lat_to","hour"]]


            def dist(row):
                return distance_latlon(row["lon_from"], row["lat_from"],row["lon_to"], row["lat_to"])

            def hist_rank(row):
                if row["distance"]>max_hist_step:
                    return max_hist_step+hist_step

                return int(row["distance"]/hist_step + 1)*hist_step

            # mesh集計
            def meshcode_from(row):
                return ju.to_meshcode(row["lat_from"], row["lon_from"], meshlevel)

            def meshcode_to(row):
                return ju.to_meshcode(row["lat_to"], row["lon_to"], meshlevel)


            # 距離分布用にODの距離を付与して出力
            df_od_grouping['mesh_from'] = df_od_grouping.apply(meshcode_from, axis=1)

            df_od_grouping["distance"] = df_od_grouping.apply(dist, axis=1)
            df_od_grouping["hist_rank"] = df_od_grouping.apply(hist_rank, axis=1)

            df_od_grouping_hist = pd.pivot_table(df_od_grouping, index='mesh_from', columns='hist_rank',  aggfunc=len)



            df_od_grouping['mesh_to'] = df_od_grouping.apply(meshcode_to, axis=1)

            # 起終点メッシュで集計
            pbar.set_description(f"creating mesh...{os.path.basename(file)}")
            df_mesh_group = df_od_grouping.groupby(["mesh_from", "mesh_to"], as_index=False).count()
            df_mesh_group = df_mesh_group[["mesh_from", "mesh_to", "dailyid"]]
            if df_mesh_group_out is None:
                df_mesh_group_out = df_mesh_group
            else:
                df_mesh_group_out = df_mesh_group_out.append(df_mesh_group)

            # hour, 起終点メッシュで集計
            """
            df_group = df_od_grouping.groupby(["mesh_from", "mesh_to", "hour"], as_index=False).count()
            df_group = df_group[["mesh_from", "mesh_to", "hour", "dailyid"]]
            """

            if print_export_string:
                print(f"export mesh trip OD data:{path}")



        df_did_convex_out_all.to_csv(f"{output_dir}//mesh_distance_far.csv", index=False)

        df_mesh_group_out = df_mesh_group_out.groupby(["mesh_from", "mesh_to"], as_index=False).mean()

        df_mesh_group_out = df_mesh_group_out.rename(columns={'dailyid': 'count'})
        coords_from = ju.to_meshpoint(df_mesh_group_out['mesh_from'], 0.5, 0.5)
        coords_to = ju.to_meshpoint(df_mesh_group_out['mesh_to'], 0.5, 0.5)

        df_mesh_group_out['lat_mesh_from'] = df_mesh_group_out['mesh_from'].apply(lambda x: ju.to_meshpoint(x, 0.5, 0.5)[0])
        df_mesh_group_out['lon_mesh_from'] = df_mesh_group_out['mesh_from'].apply(lambda x: ju.to_meshpoint(x, 0.5, 0.5)[1])
        df_mesh_group_out['lat_mesh_to'] = df_mesh_group_out['mesh_to'].apply(lambda x: ju.to_meshpoint(x, 0.5, 0.5)[0])
        df_mesh_group_out['lon_mesh_to'] = df_mesh_group_out['mesh_to'].apply(lambda x: ju.to_meshpoint(x, 0.5, 0.5)[1])

        # export
        path = f"{output_dir}/{city}/trip_end_mesh{meshlevel}.csv"
        df_mesh_group_out.to_csv(path, index=False)

    except psycopg2.Error as e:
        print("NG Copy error! ")
        print(e.pgerror)
