import pandas as pd
import glob
import os
import argparse

def merge_results(pattern, output):
    files = glob.glob(pattern)
    if not files:
        print(f"No files found matching {pattern}")
        return

    all_dfs = []
    for f in files:
        try:
            df = pd.read_excel(f)
            all_dfs.append(df)
            print(f"Loaded {f}")
        except Exception as e:
            print(f"Error loading {f}: {e}")

    if not all_dfs:
        return

    combined = pd.concat(all_dfs, ignore_index=True)
    
    # Group and aggregate
    # We sum Call Count and Total Executions
    # We take the max of Function Size (should be identical)
    agg_dict = {
        'Call Count': 'sum',
        'Total Executions': 'sum',
        'Function Size': 'max',
        'File Path': 'first'
    }
    
    merged = combined.groupby(['Method Path', 'Target']).agg(agg_dict).reset_index()
    merged = merged.sort_values(by='Call Count', ascending=False)
    
    merged.to_excel(output, index=False)
    print(f"Merged results saved to {output}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--pattern", default="*_FanIn.xlsx")
    parser.add_argument("--output", default="Dynamic_FanIn_Merged.xlsx")
    args = parser.parse_args()
    merge_results(args.pattern, args.output)
