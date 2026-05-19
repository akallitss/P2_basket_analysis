import json
import argparse
import pandas as pd
import re

def fmt(x, width):
    """Convert value to string and pad to fixed width."""
    s = "" if pd.isna(x) else str(x)
    return s.ljust(width)  # left-aligned; use rjust for right alignment

parser = argparse.ArgumentParser(
                    prog='MapTransform',
                    description='Mapping transformations')
parser.add_argument('--in_file',required=True)
parser.add_argument('--out_file')

args=parser.parse_args()

in_file=args.in_file
out_file=""
out_dir="./"
if args.out_file is None:
	if in_file.rfind("/")<0:
		out_file=in_file
	else:
		out_file=in_file[in_file.rfind("/")+1:-1]
		out_dir=in_file[:in_file.rfind("/")]

	out_file=out_file[:out_file.index(".")]+"_TBreco.txt"
else:
	out_file=args.out_file

# Read data
with open(args.in_file) as f:
    raw_data=f.read()

# Clean 'nan'/'none' text to valid JSON 'null'
clean_data = re.sub(r'(?i)\b(nan|none)\b', 'null', raw_data)
data = json.loads(clean_data)

df = pd.DataFrame(data)
# Explode lists to rows
df = df.explode(["Pin_id","Pad_id_x","Pad_id_y","Pad_pos_x","Pad_pos_y","Pad_area","Track_id"], ignore_index=True)

# --- MODIFIED SECTION ---

# OLD LOGIC:
# df['n'] = df.groupby('Pin_id').cumcount()
# wide = df.pivot(index="Pin_id", columns="n", values=["Pad_id_x", "Pad_id_y"])

# NEW LOGIC:
# Use Connector_id directly as the pivot column.
# This assumes (Pin_id, Connector_id) pairs are unique.
wide = df.pivot(index="Pin_id", columns="Connector_id", values=["Pad_id_x", "Pad_id_y"])

# Flatten column names using the Connector_id number
# We cast 'i' to int to ensure we get suffixes like '_1' instead of '_1.0'
wide.columns = [f"{name}_{int(i)}" for (name, i) in wide.columns]

# ------------------------

# 4. Reorder columns as x_ID, y_ID, x_ID, y_ID...
cols = wide.columns

# Extract numeric suffixes (Connector IDs) correctly
n_values = sorted({
    int(c.split("_")[-1])
    for c in cols
    if c.startswith("Pad_id_")
})

reordered_cols = []
for i in n_values:
    if f"Pad_id_x_{i}" in cols:
        reordered_cols.append(f"Pad_id_x_{i}")
    if f"Pad_id_y_{i}" in cols:
        reordered_cols.append(f"Pad_id_y_{i}")

# Add Pin_id column first
reordered_cols = ["Pin_id"] + reordered_cols

# Apply ordering and reset index to make Pin_id a normal column
wide = wide.reset_index()[reordered_cols]

col_width = 12  # adjust as needed

with open(out_dir+"/"+out_file, "w") as f:
    # Write header
    header = "".join(fmt(col, col_width) for col in wide.columns)
    f.write(header + "\n")

    # Write rows
    for _, row in wide.iterrows():
        line = "".join(fmt(val, col_width) for val in row.values)
        f.write(line + "\n")
