import json
import argparse
import pandas as pd
import re
import os

def main():
    parser = argparse.ArgumentParser(
                        prog='MapFormatter',
                        description='Formatting the mapping for Corryvreckan (strip detectors ONLY!)')
    parser.add_argument('--in_file', required=True, help='Input JSON file')
    parser.add_argument('--out_file', help='Output config file (optional)')

    args = parser.parse_args()

    in_file = args.in_file
    
    # 1. Handle Output File Naming (Replacing extension with .config)
    if args.out_file is None:
        base_name, _ = os.path.splitext(in_file)
        out_file = f"{base_name}.config"
    else:
        out_file = args.out_file

    # 2. Read and Clean Data
    with open(in_file, 'r') as f:
        raw_data = f.read()

    # Clean 'nan'/'none' text to valid JSON 'null'
    clean_data = re.sub(r'(?i)\b(nan|none)\b', 'null', raw_data)
    data = json.loads(clean_data)

    df = pd.DataFrame(data)
    
    # 3. Explode lists to rows
    # Include the expected "strip_id" and new "coordinate" variables. 
    explode_cols = [
        "Pin_id", "Pad_id_x", "Pad_id_y", "Pad_pos_x", "Pad_pos_y", 
        "Pad_area", "Track_id", "strip_id", "coordinate"
    ]
    
    # Only explode columns that actually exist in the dataframe to avoid KeyErrors
    explode_cols = [col for col in explode_cols if col in df.columns]
    
    if explode_cols:
        df = df.explode(explode_cols, ignore_index=True)

    # 4. Generate the Output
    with open(out_file, "w") as f:
        for _, row in df.iterrows():
            # Extract variables safely, providing defaults if missing
            apv = row.get("Connector_id", 0)
            channel = row.get("Pin_id", 0)
            strip = row.get("strip_id", 0) 
            
            # Fetch the axis directly from the input data instead of a dictionary map
            axis = row.get("coordinate", "U") 

            # Handle potential NaNs safely before formatting
            try:
                apv_val = int(float(apv)) if pd.notnull(apv) else 0
                ch_val = int(float(channel)) if pd.notnull(channel) else 0
                strp_val = int(float(strip)) if pd.notnull(strip) else 0
            except ValueError:
                # Fallback to string if they contain non-numeric characters
                apv_val, ch_val, strp_val = apv, channel, strip

            # 5. Apply the specific aesthetic formatting
            # {ch_val:>3} ensures the channel number is right-aligned to 3 spaces
            line = f"apv={apv_val}, axis={axis} , channel={ch_val:>3} ,  strip={strp_val}"
            f.write(line + "\n")
            
    print(f"Successfully wrote mapped config to: {out_file}")

if __name__ == "__main__":
    main()
