import json
import argparse
import pandas as pd
import re

parser = argparse.ArgumentParser(
                    prog='MapTransform',
                    description='Mapping transformations')
parser.add_argument('--rotate_pads')
parser.add_argument('--flip_pads',action="store_true")
parser.add_argument('--rotate_connectors',action="store_true")
parser.add_argument('--in_file',required=True)
parser.add_argument('--out_file')
parser.add_argument('--out_dir')

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

	out_file=out_file[:out_file.index(".")]+"_Transformed.json"
else:
	out_file=args.out_file



def rotate_pin(pin):
	return 127-pin


with open(args.in_file) as f:
	raw_data=f.read()

# (?i)   -> Case-insensitive flag (matches NaN, nan, None, none)
# \b     -> Word boundary (start of word)
# (nan|none) -> Match group: either 'nan' OR 'none'
# \b     -> Word boundary (end of word)
clean_data = re.sub(r'(?i)\b(nan|none)\b', 'null', raw_data)
data = json.loads(clean_data)

#print(json.dumps(data, indent=4))

df=pd.DataFrame(data)
df=df.explode(["Pin_id","Pad_id_x","Pad_id_y","Pad_pos_x","Pad_pos_y","Pad_area","Track_id"],ignore_index=True)

if args.rotate_connectors:
	df["Pin_id"]=df["Pin_id"].apply(rotate_pin)
	
if args.rotate_pads is not None:
	for i in range(int(args.rotate_pads)):
		maximum=df["Pad_id_x"].max()
		for idx, row in df.iterrows():
			x=row["Pad_id_x"]
			y=row["Pad_id_y"]

			df.at[idx,"Pad_id_x"]=y
			df.at[idx,"Pad_id_y"]=maximum-x

if args.flip_pads:
	for idx, row in df.iterrows():
		x=row["Pad_id_x"]
		y=row["Pad_id_y"]

		df.at[idx,"Pad_id_x"]=y
		df.at[idx,"Pad_id_y"]=x

group_keys = ["Connector_id"]

df = pd.DataFrame(
    df.groupby(group_keys, sort=False)
      .apply(lambda g: {
          **{k: v for k, v in zip(group_keys, g.name if isinstance(g.name, tuple) else (g.name,))},
          **{col: g[col].tolist() for col in df.columns if col not in group_keys}
      })
      .tolist()
)

df['Pin_id'] = df['Pin_id'].astype(str)
df['Pad_id_x'] = df['Pad_id_x'].astype(str)
df['Pad_id_y'] = df['Pad_id_y'].astype(str)
df['Pad_pos_x'] = df['Pad_pos_x'].astype(str)
df['Pad_pos_y'] = df['Pad_pos_y'].astype(str)
df['Pad_area'] = df['Pad_area'].astype(str)
df['Track_id'] = df['Track_id'].astype(str)

json_str=df.to_json(orient="records",indent=2)
json_str=re.sub(r'\"\[',' [',json_str)
json_str=re.sub(r']\"',']',json_str)

with open(out_dir+"/"+out_file, "w") as f:
	f.write(json_str)
