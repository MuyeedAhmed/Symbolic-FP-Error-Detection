import pandas as pd
import glob
import os
import argparse

def merge_results(pattern, output):
    files = glob.glob(pattern)
    print(files)
    if not files:
        return

    all_dfs = []
    for f in files:
        try:
            df = pd.read_excel(f)
            all_dfs.append(df)
        except Exception as e:
            print(f"Error loading {f}: {e}")

    if not all_dfs:
        return

    combined = pd.concat(all_dfs, ignore_index=True)

    agg_dict = {
        'Call Count': 'sum',
        'Total Executions': 'sum',
        'Function Size': 'max',
        'File Path': 'first'
    }
    
    merged = combined.groupby(['Method Path', 'Target']).agg(agg_dict).reset_index()
    merged = merged.sort_values(by='Call Count', ascending=False)
    
    Filtered = FilterResults(merged)

    Filtered.to_excel(output, index=False)

def FilterResults(df):
    return df[~df["Method Path"].str.contains(r"\.tests\.", na=False, regex=True)]

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--pattern", default="*_FanIn.xlsx")
    parser.add_argument("--output", default="Dynamic_FanIn_Merged.xlsx")
    args = parser.parse_args()
    merge_results(args.pattern, args.output)
