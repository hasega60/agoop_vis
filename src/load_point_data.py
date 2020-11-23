import pandas as pd
import psycopg2
import jismesh.utils as ju
import numpy as np
from geojson import Point, LineString, Feature, FeatureCollection
from tqdm import tqdm
from datetime import datetime
import math
import h3

# 入力csvの配置
data_dir = "../data"
# 出力先
output_dir = "../result"

# アニメーション設定項目
# 集計期間
start_datetime='2019-10-03 6:00:00'
end_datetime='2019-10-04 2:00:00'

# 対象外とする
max_speed = 150

# 対象都市と市町村コード
city = "tsukuba"

city_dict = {
    "tsukuba": 8220,
    "tsuchiura": 8203,
}

#集計するODを集計するメッシュレベル
meshlevel = 3

#元データポイントの出力
export_org_point = True

semiMajorAxis = 6378137.0  # 赤道半径
flattening = 1 / 298.257223563  # 扁平率
e_2 = flattening * (2 - flattening)
degree = math.pi / 180
line_color_rgb = (255, 120, 120) #geojsonのlineRGB


def distance_latlon(x1, y1, x2, y2):
    coslat = math.cos((y1 + y2) / 2 * degree)
    w2 = 1 / (1 - e_2 * (1 - coslat * coslat))
    dx = (x1 - x2) * coslat
    dy = (y1 - y2) * w2 * (1 - e_2)
    return math.sqrt((dx * dx + dy * dy) * w2) * semiMajorAxis * degree


def dataframe_to_trip(df:pd.DataFrame, trip_id_col, lat_col, lon_col, time_col, max_speed=150,  out_altitude=False):

    def _add_output(path, trip, trip_count, trip_origin_destination, trip_out, id_out, id_out_org):
        od = path[0]
        d = path[-1]
        od.extend(d)
        trip_origin_destination.append(od)
        trip_out.append(path_to_trip_geojson(path))
        id_out.append(f"{trip}_{trip_count}")
        id_out_org.append(trip)
        return trip_origin_destination, trip_out, id_out, id_out_org

    trip_out, id_out, id_out_org= [],[],[]
    trip_origin_destination=[]
    cols = list(df.columns)
    df_trip_id = df[trip_id_col]
    trip_id = df_trip_id.unique()
    for trip in trip_id:
        df_trip = df.copy()
        df_trip = df_trip[df_trip[trip_id_col]==trip]

        values = df_trip.values
        path = []
        trip_count=0
        time_before, lat_before, lon_before = -1, -1, -1
        for point in values:
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
            if out_altitude:
                # 高さを考慮するとき
                elevation = point.elevation
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
                        trip_origin_destination, trip_out, id_out, id_out_org=_add_output(path, trip, trip_count,
                                                                              trip_origin_destination,
                                                                              trip_out,
                                                                              id_out, id_out_org)
                        trip_count += 1
                        path=[]
            else:
                path.append([lon, lat, elevation, t])
                time_before = t
                lat_before = lat
                lon_before = lon

        if len(path) > 1:
            trip_origin_destination, trip_out, id_out, id_out_org = _add_output(path, trip, trip_count,
                                                                                trip_origin_destination,
                                                                                trip_out,
                                                                                id_out, id_out_org)

    return trip_out, id_out, id_out_org, trip_origin_destination

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
    :param user_cols_dict: 任意の属性を追加　keyが属性名，valueが属性値
    :param path: 出力パス
    :return:
    """

    df_out = pd.DataFrame(trip_feature_list, columns=["_geojson"])
    if len(user_cols_dict) > 0:
        for k, v in user_cols_dict.items():
            df_out[k] = v

    df_out.to_csv(path, index=False)
    print(f"export trip animation data:{path}")

if __name__ == '__main__':
    try:
        df_all = pd.read_csv(f"{data_dir}/test_08.csv")
        df_all = df_all[["dailyid", "year", "month", "day", "dayofweek", "hour", "minute", "latitude", "longitude",
                         "os","logtype_subcategory","accuracy","speed","estimated_speed_flag","course",
                         "prefcode", "citycode", "home_prefcode","home_citycode",
                         "workplace_prefcode","workplace_citycode","transportation_type", "gender"]]


        df_time = df_all[["year", "month", "day", "hour", "minute"]]
        df_all["time"] = pd.to_datetime(df_time)
        df_all = df_all.sort_values(['dailyid', 'time'])

        # データ抽出部分
        df_city = df_all[df_all["home_citycode"]==city_dict[city]]
        df_city = df_city[df_city["os"]=="Android"]
        df_city = df_city[(df_city['time'] > datetime.strptime(start_datetime, '%Y-%m-%d %H:%M:%S'))
                          & (df_city['time'] < datetime.strptime(end_datetime, '%Y-%m-%d %H:%M:%S'))]
        #df_city = df_city[(df_city["hour"] >= start_hour)&(df_city["hour"] <= end_hour)]

        # 都市別ポイント抽出
        if export_org_point:
            path = f"{output_dir}/point_{city}.csv"
            df_city.to_csv(path)
            print(f"export city point data:{path}")

        # trip animationデータの抽出
        # 停止時に累積されるポイント（timer）は除去
        df_city_trip = df_city[df_city["logtype_subcategory"]!="timer"]
        trip_list, id_list, id_out_org, trip_origin_destination = dataframe_to_trip(df_city_trip, trip_id_col="dailyid",
                                      lat_col="latitude", lon_col="longitude",time_col="time", max_speed=max_speed)

        id_col_dict={"dailyid_c": id_list, "dailyid": id_out_org}
        output_trip_csv(trip_feature_list=trip_list, user_cols_dict=id_col_dict, path=f"{output_dir}/trips_{city}.csv")
        df_od = pd.DataFrame(trip_origin_destination, columns=["lon_from", "lat_from", "elevation_from", "time_from",
                                                               "lon_to", "lat_to", "elevation_to", "time_to"])

        for k, v in id_col_dict.items():
            df_od[k] = v

        # 時間属性を付与
        df_od["hour"] = df_od['time_from'].apply(lambda x: datetime.fromtimestamp(x).hour)

        df_od_grouping = df_od[["dailyid","lon_from", "lat_from","lon_to", "lat_to","hour"]]

        # mesh集計
        def meshcode_from(row):
            return ju.to_meshcode(row["lat_from"], row["lon_from"], meshlevel)

        def meshcode_to(row):
            return ju.to_meshcode(row["lat_to"], row["lon_to"], meshlevel)

        df_od_grouping['mesh_from'] = df_od_grouping.apply(meshcode_from, axis=1)
        df_od_grouping['mesh_to'] = df_od_grouping.apply(meshcode_to, axis=1)

        # 起終点メッシュで集計
        df_group = df_od_grouping.groupby(["mesh_from", "mesh_to"], as_index=False).count()
        df_group = df_group[["mesh_from", "mesh_to", "dailyid"]]

        # hour, 起終点メッシュで集計
        """
        df_group = df_od_grouping.groupby(["mesh_from", "mesh_to", "hour"], as_index=False).count()
        df_group = df_group[["mesh_from", "mesh_to", "hour", "dailyid"]]
        """

        df_group = df_group.rename(columns={'dailyid': 'count'})

        coords_from = ju.to_meshpoint(df_group['mesh_from'], 0.5, 0.5)
        coords_to = ju.to_meshpoint(df_group['mesh_to'], 0.5, 0.5)

        df_group['lat_mesh_from'] = df_group['mesh_from'].apply(lambda x :ju.to_meshpoint(x, 0.5, 0.5)[0])
        df_group['lon_mesh_from'] = df_group['mesh_from'].apply(lambda x: ju.to_meshpoint(x, 0.5, 0.5)[1])
        df_group['lat_mesh_to'] = df_group['mesh_to'].apply(lambda x :ju.to_meshpoint(x, 0.5, 0.5)[0])
        df_group['lat_mesh_to'] = df_group['mesh_to'].apply(lambda x: ju.to_meshpoint(x, 0.5, 0.5)[1])

        # export
        path=f"{output_dir}/trip_end_mesh{meshlevel}_{city}.csv"
        df_group.to_csv(path, index=False)
        print(f"export mesh trip OD data:{path}")

        """
        # h3 indexでポイントを集計
        h3_index_level=9
        df_hour = df_all[["year", "month", "day", "hour"]]

        text = lambda x: f't{x}'
        df_h3 = df_all[["dailyid", "latitude", "longitude","transportation_type"]]
        df_h3["time"] = pd.to_datetime(df_hour)
        df_h3["transportation_type"] = df_h3["transportation_type"].map(text)

        h3_index = lambda x: h3.geo_to_h3(x.latitude, x.longitude, h3_index_level)
        df_h3["h3_index"] = df_h3[["latitude", "longitude"]].apply(h3_index, axis=1)
        df_group_h3 = df_h3.groupby(["h3_index", "time"], as_index=False).count()
        df_group_h3 = df_group_h3.rename(columns={'dailyid': 'count'})

        # export
        df_group_h3.to_csv(f"test_h3_{h3_index_level}.csv", index=False)
        user_order = 0
        dailyid = ""
        orders ,trip_ends = [], []
        for r in tqdm(df_all.values):
            # 順番つけて
            if dailyid != r[0]:
                user_order = 0
                dailyid = r[0]
                trip_ends.append(r)
            else:
                if user_order == 0:
                    trip_ends.append(r)

                user_order += 1

            orders.append(user_order)

        df_tripend = pd.DataFrame(trip_ends, columns=df_all.columns)
        df_tripend = df_tripend[["dailyid", "time", "logtype_subcategory", "transportation_type", "latitude", "longitude"]]
        trips = []
        r_before = None
        for r in tqdm(df_tripend.values):
            if r_before is None:
                r_before = r
                continue
            if r_before[0] == r[0] and r_before[1] != r[1]:
                time_to, lat_to, lon_to = r[1], r[-2], r[-1]
                r_before = np.append(r_before,[time_to, lat_to, lon_to])
                trips.append(r_before)

            r_before = r


        df_all["user_order"] = orders

        # 属性絞って全て出力
        #df_export = df_all[["dailyid", "time","user_order", "logtype_subcategory","transportation_type", "latitude", "longitude"]]
        #df_export.to_csv("test_export.csv", index=False)

        # logtype_subcategory属性でトリップ抽出
        df_trip_logtype = df_all[(df_all["logtype_subcategory"] == "move")|(df_all["logtype_subcategory"] == "arrival")|
                         (df_all["logtype_subcategory"] == "departure")|(df_all["logtype_subcategory"] == "timer")]
        df_trip_logtype = df_trip_logtype[["dailyid", "time", "user_order", "logtype_subcategory", "transportation_type", "latitude", "longitude"]]

        # export
        df_trip_logtype.to_csv("test_trip_logtype.csv", index=False)

        df_trip = pd.DataFrame(trips, columns=["dailyid", "time_from", "logtype_subcategory", "transportation_type",
                                               "lat_from", "lon_from", "time_to", "lat_to", "lon_to"])
        text = lambda x: f't{x}'
        df_trip["transportation_type"] = df_trip["transportation_type"].map(text)

        # export
        df_trip.to_csv("test_trip_end.csv", index=False)

        # mesh集計
        meshlevel = 4
        df_trip['mesh_from'] = ju.to_meshcode(df_trip.lat_from, df_trip.lon_from, meshlevel)
        df_trip['mesh_to'] = ju.to_meshcode(df_trip.lat_to, df_trip.lon_to, meshlevel)

        df_group = df_trip.groupby(["mesh_from", "mesh_to", "transportation_type"], as_index=False).count()
        df_group = df_group[["mesh_from", "mesh_to", "transportation_type", "time_from"]]
        df_group = df_group.rename(columns={'time_from': 'count'})

        coords_from = ju.to_meshpoint(df_group['mesh_from'], 0.5, 0.5)
        coords_to = ju.to_meshpoint(df_group['mesh_to'], 0.5, 0.5)

        df_group['lat_mesh_from'] = ju.to_meshpoint(df_group['mesh_from'], 0.5, 0.5)[0]
        df_group['lon_mesh_from'] = ju.to_meshpoint(df_group['mesh_from'], 0.5, 0.5)[1]
        df_group['lat_mesh_to'] = ju.to_meshpoint(df_group['mesh_to'], 0.5, 0.5)[0]
        df_group['lon_mesh_to'] = ju.to_meshpoint(df_group['mesh_to'], 0.5, 0.5)[1]

        # export
        df_group.to_csv(f"test_trip_end_mesh{meshlevel}.csv", index=False)
        """

    except psycopg2.Error as e:
        print("NG Copy error! ")
        print(e.pgerror)
