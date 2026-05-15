import pandas as pd
import glob
import os
import argparse

def merge_comprehensive(pattern, output):
    files = glob.glob(pattern)
    if not files:
        print(f"No files found matching {pattern}")
        return

    edge_dfs = []
    flat_dfs = []
    
    for f in files:
        try:
            # Read both sheets
            edge_df = pd.read_excel(f, sheet_name='Edge Profile')
            flat_df = pd.read_excel(f, sheet_name='Flat Profile')
            edge_dfs.append(edge_df)
            flat_dfs.append(flat_df)
            print(f"Loaded {f}")
        except Exception as e:
            print(f"Error loading {f}: {e}")

    if not edge_dfs and not flat_dfs:
        return

    # --- Merge Edge Profiles ---
    if edge_dfs:
        combined_edge = pd.concat(edge_dfs, ignore_index=True)
        agg_edge = {
            'Call Count': 'sum',
            'Callee Total Executions': 'sum',
            'Callee Size': 'max',
            'Caller File': 'first',
            'Callee File': 'first'
        }
        merged_edge = combined_edge.groupby(['Caller Path', 'Callee Path']).agg(agg_edge).reset_index()
        merged_edge = merged_edge.sort_values(by='Call Count', ascending=False)
    else:
        merged_edge = pd.DataFrame()

    # --- Merge Flat Profiles ---
    if flat_dfs:
        combined_flat = pd.concat(flat_dfs, ignore_index=True)
        agg_flat = {
            'Execution Count': 'sum',
            'Function Size': 'max',
            'File Path': 'first',
            'Line No': 'first'
        }
        merged_flat = combined_flat.groupby(['Method Path']).agg(agg_flat).reset_index()
        merged_flat = merged_flat.sort_values(by='Execution Count', ascending=False)
    else:
        merged_flat = pd.DataFrame()

    # --- Save merged result ---
    try:
        with pd.ExcelWriter(output) as writer:
            if not merged_edge.empty:
                merged_edge.to_excel(writer, sheet_name='Edge Profile', index=False)
            if not merged_flat.empty:
                merged_flat.to_excel(writer, sheet_name='Flat Profile', index=False)
        print(f"Merged comprehensive results saved to {output}")
    except Exception as e:
        print(f"Error saving merged Excel: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--pattern", default="ComprehensiveResults/*_Profile.xlsx")
    parser.add_argument("--output", default="ComprehensiveResults/Comprehensive_Profile_Merged.xlsx")
    args = parser.parse_args()
    merge_comprehensive(args.pattern, args.output)
