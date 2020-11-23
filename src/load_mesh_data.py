import pandas as pd
import jismesh.utils as ju
from tqdm import tqdm

data_dir = "../data/mesh"
output_dir = "../result/mesh"
mesh_level = 3
mesh_level_digit = [4,6,8,9,11,12]
#フィルタ対象メッシュ（一次メッシュ）
filter_mesh=[5339,5340,5439,5440]

output_hour_agg = False # 時間帯別集計結果を出力するか


def add_mesh_latlon_col(df:pd.DataFrame, mesh_col_index=0):
    #TODO 11桁のメッシュ番号の場合は独自実装が必要 4次メッシュを5等分したもの 10桁目:緯度3/3600，11桁：軽度4.5/3600
    out_list, lat_list, lon_list=[],[],[]
    for v in tqdm(df.values):
        try:
            # meshcodeに規格外のものが入っているケースがあるので，それは取り除く
            coords=ju.to_meshpoint(v[mesh_col_index], 0.5, 0.5)
            lat_list.append(coords[0])
            lon_list.append(coords[1])
            out_list.append(v)
        except:
            continue

    print(f"true mesh code:{len(out_list)} / {len(df)}")

    df_out = pd.DataFrame(out_list, columns=df.columns)
    df_out["lat"] = lat_list
    df_out["lon"] = lon_list
    return df_out

if __name__ == '__main__':
    df_all = pd.read_csv(f"{data_dir}/MDP_20191003.csv")
    df_all["mesh_code"] = df_all["meshid"].apply(lambda x: int(str(x)[:mesh_level_digit[mesh_level-1]]))
    df_all["mesh1"] = df_all["meshid"].apply(lambda x: int(str(x)[:mesh_level_digit[0]]))
    #一次メッシュでフィルタ
    df_all = df_all[df_all["mesh1"].isin(filter_mesh)]

    # mesh, 時間帯単位で集計----------------------------------
    if output_hour_agg:
        df_hour = df_all[df_all["hour"] != "day"]
        # 時間でソート
        df_time = df_hour[["year", "month", "day", "hour"]]
        df_time["minute"]=0
        df_hour["time"] = pd.to_datetime(df_time)
        df_hour = df_hour.sort_values(['mesh_code', 'time'])

        df_hour_group=df_hour[["mesh_code","hour","population"]]
        df_group = df_hour_group.groupby(["mesh_code", "hour"], as_index=False).sum()

        # meshcodeに規格外のものが入っているケースがあるので，それは取り除く
        df_out=add_mesh_latlon_col(df_group, mesh_col_index=0)

        path = f"{output_dir}/trip_end_mesh{mesh_level}_hour.csv"
        df_out.to_csv(path, index=False)
        print(f"export poplation per mesh-hour:{path}")

    # mesh, 日付で集計（）----------------------------------
    df_day = df_all[(df_all["hour"] == "day")]
    df_day_group=df_day[["mesh_code","population"]]
    df_group = df_day_group.groupby(["mesh_code"], as_index=False).sum()
    #df_weekend = df_group[df_group["dayflag"]==0]
    #df_weekday = df_group[df_group["dayflag"]==1]
    #df_mesh_daytype=pd.merge(df_weekday, df_weekend, on="mesh_code", how="outer").fillna(0)

    # meshcodeに規格外のものが入っているケースがあるので，それは取り除く
    df_mesh_daytype_out = add_mesh_latlon_col(df_group, mesh_col_index=0)

    path = f"{output_dir}/trip_end_mesh{mesh_level}_day.csv"
    df_mesh_daytype_out.to_csv(path, index=False)
    print(f"export population per mesh-daytype:{path}")

