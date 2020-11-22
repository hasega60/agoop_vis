import pandas as pd
import psycopg2
import jismesh.utils as ju
import numpy as np
from geojson import Point, LineString, Feature, FeatureCollection
from tqdm import tqdm
from datetime import datetime
import math
import h3

data_dir = "../data/mesh"
output_dir = "../result/mesh"
mesh_level = 4
mesh_level_digit = [4,6,8,9,11,12]

if __name__ == '__main__':

    try:
        df_all = pd.read_csv(f"{data_dir}/MDP_20191003.csv")
        df_all["mesh_code"] = df_all["meshid"].apply(lambda x: int(str(x)[:mesh_level_digit[mesh_level-1]]))

        # mesh, 時間帯単位で集計
        df_hour = df_all[df_all["hour"] != "day"]
        # 時間でソート
        df_time = df_hour[["year", "month", "day", "hour"]]
        df_time["minute"]=0
        df_hour["time"] = pd.to_datetime(df_time)
        df_hour = df_hour.sort_values(['mesh_code', 'time'])

        df_hour_group=df_hour[["mesh_code","hour","population"]]
        df_group = df_hour_group.groupby(["mesh_code", "hour"], as_index=False).sum()
        df_group['lat'] = df_group['mesh_code'].apply(lambda x :ju.to_meshpoint(x, 0.5, 0.5)[0])
        df_group['lon'] = df_group['mesh_code'].apply(lambda x: ju.to_meshpoint(x, 0.5, 0.5)[1])
        df_group.to_csv(f"{output_dir}/trip_end_mesh{mesh_level}_hour.csv", index=False)

        # mesh, dayflag(平日休日)で集計
        df_day = df_all[(df_all["hour"] == "day") & ((df_all["dayflag"] == 0)|(df_all["dayflag"] == 1))]
        df_day_group=df_day[["mesh_code", "dayflag", "population"]]
        df_group = df_day_group.groupby(["mesh_code", "dayflag"], as_index=False).sum()


        df_weekend = df_group[df_group["dayflag"]==0]
        df_weekday = df_group[df_group["dayflag"]==1]
        df_mesh_daytype=pd.merge(df_weekday, df_weekend, on="mesh_code", how="outer").fillna(0)


        df_group['lat'] = df_group['mesh_code'].apply(lambda x :ju.to_meshpoint(x, 0.5, 0.5)[0])
        df_group['lon'] = df_group['mesh_code'].apply(lambda x: ju.to_meshpoint(x, 0.5, 0.5)[1])
        df_group.to_csv(f"{output_dir}/trip_end_mesh{mesh_level}_hour.csv", index=False)

    except psycopg2.Error as e:
        print("NG Copy error! ")
        print(e.pgerror)
    データ読み込み
    指定メッシュコードを取得
    時間帯で滞在人数を合計
    一時間単位の平均を取得
    代表点座標取得
    csv出力


