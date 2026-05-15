import pandas as pd
import os
import glob

def merge_and_aggregate(input_dir="CProfileResults", lib="scipy"):
    path = os.path.join(input_dir, f"{lib}_*.xlsx")
    files = glob.glob(path)
    
    if not files:
        return

    all_data = []
    for file in files:
        df = pd.read_excel(file)
        all_data.append(df)

    combined_df = pd.concat(all_data, ignore_index=True)
    index_cols = ["File Path", "Line", "Function"]
    value_cols = ["Total Calls", "Primitive Calls"]
    aggregated_df = combined_df.groupby(index_cols)[value_cols].sum().reset_index()
    aggregated_df = aggregated_df.sort_values(by="Total Calls", ascending=False)

    aggregated_df.to_excel(f"{input_dir}/{lib}_Agg.xlsx", index=False, engine='openpyxl')
    
if __name__ == "__main__":
    merge_and_aggregate(lib="scipy")