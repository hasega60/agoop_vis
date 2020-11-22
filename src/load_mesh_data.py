import pandas as pd
import psycopg2
import jismesh.utils as ju
import numpy as np
from geojson import Point, LineString, Feature, FeatureCollection
from tqdm import tqdm
from datetime import datetime
import math
import h3


if __name__ == '__main__':

    data_dir = "../data/mesh"
    output_dir= "../result/mesh"
    mesh_level = 5

    try:
        df_all = pd.read_csv(f"{data_dir}/MDP_20191003.csv")
        df_all = df_all[df_all["hour"] != "day"]




    except psycopg2.Error as e:
        print("NG Copy error! ")
        print(e.pgerror)

    データ読み込み
    指定メッシュコードを取得
    時間帯で滞在人数を合計
    一時間単位の平均を取得
    代表点座標取得
    csv出力


