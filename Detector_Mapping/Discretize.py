import json
import pandas as pd
import re
import argparse
import yaml
import os
import sys

def main():
    parser = argparse.ArgumentParser(
                        prog='Discretize',
                        description='Assigns discrete grid indices to physical coordinates based on a YAML config.')
    parser.add_argument('--in_file', required=True)
    parser.add_argument('--out_file')
    parser.add_argument('--in_dir', default="./")
    parser.add_argument('--out_dir')
    parser.add_argument('--type', required=True, choices=['pad', 'strip'])
    parser.add_argument('--config', default="./data/base-pad-discretiser.yaml")
    parser.add_argument('--snappy', action='store_true', 
                        help='If set, resolves overlaps by keeping the closest pin and moving others to the nearest empty virtual slot instead of throwing a fatal error.')

    args = parser.parse_args()

    # Path handling
    in_dir = args.in_dir
    out_dir = args.out_dir if args.out_dir is not None else in_dir
    in_file = args.in_file

    if args.out_file is None:
        base_name, _ = os.path.splitext(os.path.basename(in_file))
        out_file = f"{base_name}_indexed.json"
    else:
        out_file = args.out_file

    print(f"INFO: Processing file as a {args.type}-detector...")

    # Load YAML Config
    try:
        with open(args.config, 'r') as f:
            config = yaml.safe_load(f)
    except FileNotFoundError:
        sys.exit(f"ERROR: Configuration file {args.config} not found.")

    # Read and clean JSON
    input_path = os.path.join(in_dir, in_file)
    with open(input_path, 'r') as f:
        raw_data = f.read()

    clean_data = re.sub(r'(?i)\b(nan|none)\b', 'null', raw_data)
    data = json.loads(clean_data)

    # Explode lists into rows
    df = pd.DataFrame(data)
    explode_cols = ["Pin_id", "Pad_pos_x", "Pad_pos_y", "Pad_area", "Track_id"]
    df = df.explode([c for c in explode_cols if c in df.columns], ignore_index=True)

    # Convert numeric fields for math
    df['Pad_pos_x'] = pd.to_numeric(df['Pad_pos_x'], errors='coerce')
    df['Pad_pos_y'] = pd.to_numeric(df['Pad_pos_y'], errors='coerce')
    df['Track_id'] = pd.to_numeric(df['Track_id'], errors='coerce')

    processed_areas = []

    # Process each area defined in the YAML config
    for area_name, area_params in config.get('areas', {}).items():
        # Match area
        area_mask = (df["Pad_area"] == area_name) | (df["Pad_area"] == str(area_name)) | (df["Pad_area"] == int(area_name) if str(area_name).isdigit() else False)
        area_df = df[area_mask].copy()

        if area_df.empty:
            print(f"WARNING: Area '{area_name}' found in config but not in input data. Skipping.")
            continue

        valid_tracks = area_df['Track_id'] != -1

        # --- "CLOSEST POINT" GRID & DUPLICATE CHECKS ---
        if args.type == "pad":
            x0, y0 = area_params['x0'], area_params['y0']
            dx, dy = area_params['del_x'], area_params['del_y']

            area_df.loc[valid_tracks, 'Pad_id_x'] = ((area_df.loc[valid_tracks, 'Pad_pos_x'] - x0) / dx).round()
            area_df.loc[valid_tracks, 'Pad_id_y'] = ((area_df.loc[valid_tracks, 'Pad_pos_y'] - y0) / dy).round()
            
            # Check for duplicates among valid tracks
            dupes_mask = area_df[valid_tracks].duplicated(subset=['Pad_id_x', 'Pad_id_y'], keep=False)
            if dupes_mask.any():
                colliding_pins = area_df[valid_tracks][dupes_mask].sort_values(by=['Pad_id_x', 'Pad_id_y']).copy()
                
                # Calculate the virtual position based on the assigned grid index
                colliding_pins['Virt_pos_x'] = x0 + (colliding_pins['Pad_id_x'] * dx)
                colliding_pins['Virt_pos_y'] = y0 + (colliding_pins['Pad_id_y'] * dy)
                
                if not args.snappy:
                    print(f"\n[!] FATAL ERROR: Duplicate Pad indices detected in Area '{area_name}'.")
                    print("The following pins mapped to the same virtual index:")
                    print(colliding_pins[['Pin_id', 'Pad_pos_x', 'Pad_pos_y', 'Virt_pos_x', 'Virt_pos_y', 'Pad_id_x', 'Pad_id_y']].to_string(index=False))
                    sys.exit(f"\nPlease check your config grid (x0, y0, del_x, del_y) for Area '{area_name}' or run with --snappy to auto-resolve.")
                else:
                    print(f"\n[!] WARNING: Duplicate Pad indices detected in Area '{area_name}'. (--snappy is active)")
                    print("The following pins mapped to the same virtual index and will be dynamically re-routed:")
                    print(colliding_pins[['Pin_id', 'Pad_pos_x', 'Pad_pos_y', 'Virt_pos_x', 'Virt_pos_y', 'Pad_id_x', 'Pad_id_y']].to_string(index=False))
                    
                    # Snappy Resolution: Keep the closest, move the rest to nearest empty slot
                    occupied = set(zip(area_df.loc[valid_tracks, 'Pad_id_x'], area_df.loc[valid_tracks, 'Pad_id_y']))
                    duplicate_groups = colliding_pins.groupby(['Pad_id_x', 'Pad_id_y'])
                    
                    for (idx_x, idx_y), group in duplicate_groups:
                        virt_x = x0 + idx_x * dx
                        virt_y = y0 + idx_y * dy
                        
                        # Calculate euclidean distance to ideal center
                        group['dist'] = ((group['Pad_pos_x'] - virt_x)**2 + (group['Pad_pos_y'] - virt_y)**2)**0.5
                        group = group.sort_values('dist')
                        
                        # First one (closest) keeps its spot, others get moved
                        for i in range(1, len(group)):
                            row_idx = group.index[i]
                            r = 1
                            found = False
                            
                            # Expanding square search (spiral) for the nearest empty slot
                            while not found:
                                for dx_i in range(-r, r+1):
                                    for dy_i in range(-r, r+1):
                                        if abs(dx_i) == r or abs(dy_i) == r:
                                            test_x, test_y = idx_x + dx_i, idx_y + dy_i
                                            if (test_x, test_y) not in occupied:
                                                area_df.loc[row_idx, 'Pad_id_x'] = test_x
                                                area_df.loc[row_idx, 'Pad_id_y'] = test_y
                                                occupied.add((test_x, test_y))
                                                found = True
                                                break
                                    if found: break
                                r += 1

            area_df.loc[~valid_tracks, ['Pad_id_x', 'Pad_id_y']] = -1
            area_df['Pad_id_x'] = area_df['Pad_id_x'].fillna(-1).astype(int)
            area_df['Pad_id_y'] = area_df['Pad_id_y'].fillna(-1).astype(int)

        elif args.type == "strip":
            start = area_params['start']
            pitch = area_params['pitch']
            axis_val = area_params['axis'].upper()
            axis_col = 'Pad_pos_x' if axis_val == 'X' else 'Pad_pos_y'

            area_df.loc[valid_tracks, 'strip_id'] = ((area_df.loc[valid_tracks, axis_col] - start) / pitch).round()
            
            # Check for duplicates among valid tracks
            dupes_mask = area_df[valid_tracks].duplicated(subset=['strip_id'], keep=False)
            if dupes_mask.any():
                colliding_pins = area_df[valid_tracks][dupes_mask].sort_values(by=['strip_id']).copy()
                
                # Calculate the virtual position based on the assigned grid index
                virt_col_name = f"Virt_pos_{axis_val.lower()}"
                colliding_pins[virt_col_name] = start + (colliding_pins['strip_id'] * pitch)

                if not args.snappy:
                    print(f"\n[!] FATAL ERROR: Duplicate Strip indices detected in Area '{area_name}'.")
                    print(f"The following pins mapped to the same virtual strip_id:")
                    print(colliding_pins[['Pin_id', axis_col, virt_col_name, 'strip_id']].to_string(index=False))
                    sys.exit(f"\nPlease check your config grid (start, pitch) for Area '{area_name}' or run with --snappy to auto-resolve.")
                else:
                    print(f"\n[!] WARNING: Duplicate Strip indices detected in Area '{area_name}'. (--snappy is active)")
                    print(f"The following pins mapped to the same virtual strip_id and will be dynamically re-routed:")
                    print(colliding_pins[['Pin_id', axis_col, virt_col_name, 'strip_id']].to_string(index=False))
                    
                    # Snappy Resolution: Keep the closest, move the rest to nearest empty slot
                    occupied = set(area_df.loc[valid_tracks, 'strip_id'])
                    duplicate_groups = colliding_pins.groupby('strip_id')
                    
                    for s_id, group in duplicate_groups:
                        virt_pos = start + s_id * pitch
                        
                        # Calculate 1D distance to ideal center
                        group['dist'] = abs(group[axis_col] - virt_pos)
                        group = group.sort_values('dist')
                        
                        # First one (closest) keeps its spot, others get moved
                        for i in range(1, len(group)):
                            row_idx = group.index[i]
                            offset = 1
                            found = False
                            
                            # Linear +1/-1 ping-pong search for nearest empty strip
                            while not found:
                                for sign in [1, -1]:
                                    test_id = s_id + sign * offset
                                    if test_id not in occupied:
                                        area_df.loc[row_idx, 'strip_id'] = test_id
                                        occupied.add(test_id)
                                        found = True
                                        break
                                offset += 1
            
            area_df.loc[~valid_tracks, 'strip_id'] = -1
            area_df['strip_id'] = area_df['strip_id'].fillna(-1).astype(int)
            
            area_df['coordinate'] = axis_val

        processed_areas.append(area_df)

    if not processed_areas:
        sys.exit("ERROR: No configured areas matched the input data. Cannot proceed.")

    # --- SINGLE OUTPUT COMBINATION ---
    full_df = pd.concat(processed_areas, ignore_index=True)
    full_df = full_df.sort_values(by=["Connector_id", "Pin_id"])
    
    agg_dict = {col: lambda x: x.tolist() for col in full_df.columns if col != 'Connector_id'}
    grouped_df = full_df.groupby('Connector_id', sort=False).agg(agg_dict).reset_index()

    # --- CUSTOM JSON FORMATTING ---
    json_lines = ["[\n"]
    records = grouped_df.to_dict(orient='records')
    
    for i, rec in enumerate(records):
        json_lines.append("  {\n")
        keys = list(rec.keys())
        for j, k in enumerate(keys):
            val_str = json.dumps(rec[k]) 
            comma = "," if j < len(keys) - 1 else ""
            json_lines.append(f'    "{k}": {val_str}{comma}\n')
            
        comma2 = "," if i < len(records) - 1 else ""
        json_lines.append(f"  }}{comma2}\n")
    json_lines.append("]")
    
    json_str = "".join(json_lines)

    # Save the single output file
    out_path = os.path.join(out_dir, out_file)
    with open(out_path, "w") as f:
        f.write(json_str)
    print(f"SUCCESS: Saved compiled config to {out_path}")

if __name__ == "__main__":
    main()
