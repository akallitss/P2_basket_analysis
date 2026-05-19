###########################################################
#Developed and maintained by Romano Orlandini (21/11/2025)
#romano.orlandini@cern.ch
###########################################################
import sys
sys.path.insert(0, "./includes")

import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import cdist
from scipy.spatial import KDTree
from matplotlib.widgets import RectangleSelector
import pandas as pd
import json, re
import tkinter as tk
from tkinter import ttk
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.text import Text
from matplotlib.path import Path
import math
import argparse
import os
import pathlib

import gerber
from gerber.primitives import Line, Circle, Rectangle
import traceback

#Custom includes
import gbr_reader
import GUI_utils

import yaml

parser = argparse.ArgumentParser(
                    prog='Romapping',
                    description='(Semi)Automatic mapping of gerber file for MPGD')
parser.add_argument('--track_files', nargs='+', required=True)
parser.add_argument('--pin_file', required=True)
parser.add_argument('--out_file')
parser.add_argument('--n_connectors', type=int, required=True)
parser.add_argument('--n_pad_areas', type=int, default=1)
parser.add_argument('--checks',action="store_true")
parser.add_argument('--in_dir')
parser.add_argument('--out_dir', default="./outputs")
parser.add_argument('--config', default="./data/base-config.yaml")

args=parser.parse_args()
#print(args.n_pad_areas)

in_dir="./"+args.in_dir if args.in_dir is not None else "./"
out_dir="./"+args.out_dir if args.out_dir is not None else "./"

out_file=""
if args.out_file is None:
	out_path=pathlib.Path(in_dir)
	out_file=out_path.name+"_match.json" if args.in_dir is not None else "out_match.json"
else:
	out_file=args.out_file

n_connectors=args.n_connectors
n_pad_areas=args.n_pad_areas

#Utilities
colours=['b','g','r','c','m','y','b','k']

#Variables
merged_tracks=[]

print("############################################################################")
print("Your Romapping experience is starting...")
print("############################################################################")

# --- MONKEY PATCH TO FIX PCB-TOOLS AMOutlinePrimitive CRASH ---
# Save the original, broken __init__ function
original_init = gerber.am_statements.AMOutlinePrimitive.__init__

# Define a safe version that catches the empty list
def safe_init(self, code, exposure, start_point, points, rotation):
    if not points:
        # If the CAD tool output a weird macro with no points, 
        # fake a closed 1-point shape so the library doesn't crash.
        points = [start_point] 
    
    # Pass the sanitized data back to the original function
    original_init(self, code, exposure, start_point, points, rotation)

# Inject our safe version into the library
gerber.am_statements.AMOutlinePrimitive.__init__ = safe_init
# --------------------------------------------------------------

# ---LOAD CONFIGURATION FILE ---
print("INFO: LOADING CONFIG FILE...")
def load_config(file_path):
  try:
    with open(file_path, 'r') as file:
      # Use safe_load to turn the YAML into a Python dictionary
      config = yaml.safe_load(file)
      return config
  except FileNotFoundError:
    print("ERROR: The configuration file was not found.")
    return None

config_data = load_config(args.config)
if config_data:
  MAX_TRACK_WIDTH_MM=config_data.get('MAX_TRACK_WIDTH_MM')
  MAX_MERGE_DIST_MM = config_data.get('MAX_MERGE_DIST_MM')

  MAX_PIN_L1=config_data.get("MAX_PIN_L1")
  MAX_PIN_L2=config_data.get("MAX_PIN_L2")

  print(f"INFO: Track Width: {MAX_TRACK_WIDTH_MM} mm")
  print(f"INFO: Merge Distance: {MAX_MERGE_DIST_MM} mm") 
  print(f"INFO: Max pin length: {MAX_PIN_L1} mm")
  print(f"INFO: Max pin height: {MAX_PIN_L2} mm")

  print()
# --------------------------------------------------------------

#READING TRACK FILES#################################################################
# Configuration parameter for your loose criteria

#MAX_TRACK_WIDTH_MM = 0.18
##MAX_MERGE_DIST_MM = 0.3
#MAX_MERGE_DIST_MM = 0.1

print("INFO: Reading track files...")
for i, track_file in enumerate(args.track_files, start=1):
	print(f"INFO: Reading file {i}/{len(args.track_files)} ({track_file})...")

	tracks = []
	current_segment = []
	filepath = os.path.join(in_dir, track_file)

	try:
		# 1. Open the file manually using standard 'r' mode
		with open(filepath, 'r') as f:
			gerber_data = f.read()

		# 2. Use loads() instead of read() to parse the raw string
		camfile = gerber.loads(gerber_data)
		camfile.to_metric() 

	except Exception as e:
		print(f"ERROR: Could not parse {track_file}. Reason: {e}")
		traceback.print_exc() # This will print the exact lines where pcb-tools failed
		continue

	# Iterate through all primitives parsed by the library
	for primitive in camfile.primitives:

		# We only care about drawn lines (equivalent to your D01 tracking)
		if isinstance(primitive, Line):

			# Apply your criteria: Circular aperture AND diameter < threshold
			if isinstance(primitive.aperture, Circle) and primitive.aperture.diameter < MAX_TRACK_WIDTH_MM:

				# Extract coordinates as lists to match your exact data structure
				start_pt = list(primitive.start)
				end_pt = list(primitive.end)

				# Reconstruct the continuous polyline segments
				if not current_segment:
					current_segment.extend([start_pt, end_pt])
				# Check if the start of this line touches the end of the previous one
				elif math.isclose(current_segment[-1][0], start_pt[0], abs_tol=1e-5) and math.isclose(current_segment[-1][1], start_pt[1], abs_tol=1e-5):
					current_segment.append(end_pt)
				else:
					# The trace was broken (equivalent to a D02/D03 move). Save and start new.
					tracks.append(current_segment)
					current_segment = [start_pt, end_pt]

	# Save the very last segment if it exists
	if current_segment:
		tracks.append(current_segment)

	print("INFO: Track segments collected!")
	print("INFO: Merging segments...")

	# Filter out tiny artifacts (keeping your existing logic exactly as-is)
	tracks_temp = []
	for m in tracks:
		dist = 0
		for j in range(len(m)-1):
			dist += math.dist(m[j], m[j+1])

		if dist < 0.05:
			continue
		tracks_temp.append(m)

	tracks = tracks_temp

	# Call to your custom KDTree merger
	#merged_tracks += gbr_reader.merge_segments_sampled_kdtree(tracks, sample_spacing=2, kdtree_threshold=1, merge_threshold=0.1)
	merged_tracks += gbr_reader.merge_segments_sampled_kdtree(tracks, sample_spacing=2, kdtree_threshold=1, merge_threshold=MAX_TRACK_WIDTH_MM)

	print("INFO: Segments merged!\n")


#READING PIN FILES#################################################################
print("INFO: Reading pin file...")

rect_segments = []
pin_tracks = []
current_segment = []

filepath = os.path.join(in_dir, args.pin_file)

try:
	# Open the file manually and parse to bypass the open() bug
	with open(filepath, 'r') as f:
		gerber_data = f.read()

	camfile = gerber.loads(gerber_data)
	camfile.to_metric()

except Exception as e:
	print(f"ERROR: Could not parse pin file {args.pin_file}.")
	traceback.print_exc()
	exit()

# Iterate through all primitives parsed by the library
for primitive in camfile.primitives:
	# GET RECTANGLES (i.e. pins)
	# Flashed D03 rectangular pads directly become 'Rectangle' primitives
	if isinstance(primitive, Rectangle) and ( (primitive.width<=MAX_PIN_L1 and primitive.height<=MAX_PIN_L2) or (primitive.width<=MAX_PIN_L2 and primitive.height<=MAX_PIN_L1) ) :
		# primitive.position is an (x, y) tuple, convert to list to match your format
		rect_segments.append(list(primitive.position))

	# GET CIRCLES (i.e. pin tracks)
	# Drawn D01 lines
	elif isinstance(primitive, Line):
		# Safe check to ensure the line has an aperture with a diameter below our threshold
		if hasattr(primitive, 'aperture') and hasattr(primitive.aperture, 'diameter'):
			if primitive.aperture.diameter < MAX_TRACK_WIDTH_MM:
				start_pt = list(primitive.start)
				end_pt = list(primitive.end)

				# Reconstruct the continuous polyline segments
				if not current_segment:
					current_segment.extend([start_pt, end_pt])
				elif math.isclose(current_segment[-1][0], start_pt[0], abs_tol=1e-5) and math.isclose(current_segment[-1][1], start_pt[1], abs_tol=1e-5):
					current_segment.append(end_pt)
				else:
					pin_tracks.append(current_segment)
					current_segment = [start_pt, end_pt]

# Save the very last track segment if it exists
if current_segment:
	pin_tracks.append(current_segment)

# Remove duplicates from rect_segments
rect_segments = list(map(list, {tuple(lst) for lst in rect_segments}))

print("INFO: Pin segments collected!")
print("INFO: Merging segments...")

# Call to your custom KDTree merger
#pin_merged_tracks = gbr_reader.merge_segments_sampled_kdtree(pin_tracks, sample_spacing=2, kdtree_threshold=1, merge_threshold=0.1)
pin_merged_tracks = gbr_reader.merge_segments_sampled_kdtree(pin_tracks, sample_spacing=2, kdtree_threshold=1, merge_threshold=MAX_MERGE_DIST_MM)
pin_merged_tracks = [gbr_reader.extend_track(track) for track in pin_merged_tracks]

print("INFO: Segments merged!\n")
print()

print("INFO: Linking pins to initial pin tracks...")
pin_track_links=gbr_reader.assign_points_to_tracks(rect_segments,pin_merged_tracks,0.4) #Index of the track corresponding to pin at list position
#print(pin_track_links)
print("INFO: Linked!")
print()

print("INFO: Reversing pin tracks to start at pin...")
for i, link in enumerate(pin_track_links):
	if link is not None:
	#if link >0:
		d0 = math.dist(rect_segments[i],pin_merged_tracks[link][0])
		d1 = math.dist(rect_segments[i],pin_merged_tracks[link][-1])
		
		if d1<d0:	
			pin_merged_tracks[link].reverse()
print("INFO: Reversed!")
print()

print("INFO: Merging pad tracks to pin tracks (this might take a while)...")
pin_track_ends=[t[-1] for t in pin_merged_tracks]
pin_track_to_pad_track_links=gbr_reader.assign_points_to_track_endpoints(pin_track_ends,merged_tracks,2) #Index of the track corresponding to pin at list position

full_merged_tracks=[]
for i, link in enumerate(pin_track_to_pad_track_links):
	if link is not None:	
	#if link >0:	
		d0 = math.dist(pin_merged_tracks[i][-1], merged_tracks[link][0])
		d1 = math.dist(pin_merged_tracks[i][-1], merged_tracks[link][-1])

		if d1<d0:	
			merged_tracks[link].reverse()

		full_merged_tracks.append(pin_merged_tracks[i]+merged_tracks[link])

print("INFO: Segments merged!")
print()

###################################################################################################
#GUI & VISUALIZATION (Connector selection)
###################################################################################################

#Visualise tracks (optional)########################################################################
if args.checks:

	fig_checks, ax_checks = plt.subplots(figsize=(10,7))
	#Pins
	x,y=[],[]
	for r in rect_segments:
		x.append(r[0])
		y.append(r[1])	
	ax_checks.plot(x,y,'s',markersize=3)

	#Full tracks
	i=0
	x,y=[],[]
	for r in full_merged_tracks:
		x=[i[0] for i in r]
		y=[i[1] for i in r]
		ax_checks.plot(x,y,color=colours[i%5],linestyle='-',linewidth=1)
		i+=1
#		x2.append(r[0][0])
#		y2.append(r[0][1])
#		x3.append(r[-1][0])
#		y3.append(r[-1][1])
	
#	ax_checks.plot(x2,y2,'v',markersize=5,color='orange')
#	ax_checks.plot(x3,y3,'^',markersize=5,color='lime')
	ax_checks.set_title("Tracks visualisation (just for checks...)")
	ax_checks.set_xlabel('X [mm]')
	ax_checks.set_ylabel('Y [mm]')

	#Pin tracks (only head and tail)
	i=0
	x2,y2=[],[]
	x3,y3=[],[]
	for r in pin_merged_tracks:
	#	x=[i[0] for i in r]
	#	y=[i[1] for i in r]
	#	ax_checks.plot(x,y,color=colours[i%5],linestyle='-',linewidth=1)
	#	i+=1
	
		x2.append(r[0][0])
		y2.append(r[0][1])
		x3.append(r[-1][0])
		y3.append(r[-1][1])
	
	ax_checks.plot(x2,y2,'v',markersize=5,color='orange')
	ax_checks.plot(x3,y3,'^',markersize=5,color='lime')
	
	#Pad tracks (only head and tail)
	i=0
	x2,y2=[],[]
	x3,y3=[],[]
	for r in merged_tracks:
	#	x=[i[0] for i in r]
	#	y=[i[1] for i in r]
	#	ax_checks.plot(x,y,color=colours[i%5],linestyle='-',linewidth=1)
	#	i+=1
	
		x2.append(r[0][0])
		y2.append(r[0][1])
		x3.append(r[-1][0])
		y3.append(r[-1][1])
	
	ax_checks.plot(x2,y2,'v',markersize=5,color='red')
	ax_checks.plot(x3,y3,'^',markersize=5,color='green')
	
	plt.show()

#Draw pins
fig, ax = plt.subplots(figsize=(10,7))
x,y=[],[]
for r in rect_segments:
	x.append(r[0])
	y.append(r[1])
	i+=1

ax.plot(x,y,'s',markersize=3)
ax.set_title("Select connector (1/"+str(n_connectors)+")")
ax.set_xlabel('X [mm]')
ax.set_ylabel('Y [mm]')

#START GUI########################################
root = tk.Tk()
root.title("Connectors selection GUI")

control_frame = ttk.Frame(root)
control_frame.pack(side=tk.TOP, fill=tk.X)

canvas = FigureCanvasTkAgg(fig, master=root)
canvas_widget = canvas.get_tk_widget()
canvas_widget.pack(fill=tk.BOTH, expand=True)

selector = GUI_utils.RectSelectorInteractive(ax)

def set_move(): selector.set_mode('move'); mode_var.set("Move")
def set_rotate(): selector.set_mode('rotate'); mode_var.set("Rotate")
def Flip(): selector.flip(); #flipped=True

ttk.Button(control_frame, text="Move", command=set_move).pack(side=tk.LEFT, padx=4)
ttk.Button(control_frame, text="Rotate", command=set_rotate).pack(side=tk.LEFT, padx=4)
ttk.Button(control_frame, text="Filp", command=Flip).pack(side=tk.LEFT, padx=4)
mode_var = tk.StringVar(value="Move")
ttk.Label(control_frame, textvariable=mode_var, width=8).pack(side=tk.LEFT, padx=6)

output = tk.Text(root, height=6)
output.pack(fill=tk.X)

areas=[]
def on_enter(event=None):
	verts, ang = selector.get_rectangle_parameters()
	if verts is None:
		output.insert(tk.END, "No rectangle.\n")
		return

	areas.append({
		"angle": ang,
		"flipped": selector.flipped,
		"vertices": verts.copy()
	})

	next_connector=len(areas)+1
	ax.set_title(f"Select connector ({next_connector}/{n_connectors})")
	canvas.draw()

	if len(areas) >= n_connectors:
		root.quit()
		root.destroy()

root.bind("<Return>", on_enter)
root.mainloop()

# --- Continue rest of script after plot window closes ---

#ENUMERATE PINS##################################################
mapping_df=pd.DataFrame(columns=["Track_id","Connector_id","Pin_id","Pad_pos_x","Pad_pos_y","Pad_area"])

#Loop over connectors
for i,pa in enumerate(areas,start=1):  #Loop over connectors
#	if i>=n_connectors break

	print("Connector "+str(i)+":")
	print("	Vertices: "+str(pa['vertices']))
	print("	Angle: "+str(pa["angle"]))
	print("	Flipped: "+str(pa["flipped"]))
	print()

	#Selection boundaries (lax)
	x_left = min(item[0] for item in pa['vertices'])
	x_right = max(item[0] for item in pa['vertices'])

	y_low = min(item[1] for item in pa['vertices'])
	y_high = max(item[1] for item in pa['vertices'])

	active_pins=[]
	active_pin_ids=[]
	for j, p in enumerate(rect_segments): #they are actually the pins
		x=p[0]
		y=p[1]

		if x>x_right or x<x_left or y>y_high or y<y_low: continue

		active_pins.append(p)
		active_pin_ids.append(j)


# Calculate angle in radians once
	rad_angle = pa['angle'] * math.pi / 180

	# Determine orientation: If cos^2 is 1, it's horizontal (0 or 180 deg)
	is_horizontal = round(math.cos(rad_angle) * math.cos(rad_angle)) == 0

	# Explicitly assign our working axes based on orientation
	if is_horizontal:
			split_axis = 1	# Y-axis separates the top row from bottom row
			sort_axis = 0	 # X-axis sorts the sequence of pins
	else:
			split_axis = 0	# X-axis separates the left row from right row
			sort_axis = 1	 # Y-axis sorts the sequence of pins

	# Get mid separation using the correct split axis
	coords_for_split = [t[split_axis] for t in active_pins]
#	mid_connector = sum(coords_for_split) / len(coords_for_split)
	min_split = min(coords_for_split)
	max_split = max(coords_for_split)
	mid_connector = (min_split + max_split) / 2.0
	print(f"INFO: Splitting connectors at {mid_connector}")

	# Sign logic (determining connector orientation)
	left_sign = round(np.sign(math.sin(rad_angle + math.pi/4))) #Determines if left side coordinates are smaller or grater than right side
	flip_sign = 1 if not pa['flipped'] else -1

	left_pins = []
	left_pin_ids = []
	right_pins = []
	right_pin_ids = []

	# Split pins into left/right (or top/bottom) using the split_axis
	for j, p in enumerate(active_pins):
			if flip_sign * left_sign * p[split_axis] < flip_sign * left_sign * mid_connector:
					left_pins.append(p)
					left_pin_ids.append(active_pin_ids[j])
					print(f"Pin {j} put on left. Split axis: {p[split_axis]}. Sort axis: {p[sort_axis]}")
			else:
					right_pins.append(p)
					right_pin_ids.append(active_pin_ids[j])
					print(f"Pin {j} put on right. Split axis: {p[split_axis]}. Sort axis: {p[sort_axis]}")

	# Safety check to prevent crashing if a side is empty
	if not left_pins or not right_pins:
			print(f"ERROR: Empty pin list on Connector {i}. Left: {len(left_pins)}, Right: {len(right_pins)}")
			continue 

	# Ascending logic (determines if pin index increase with increasing or decreasing pin coordinate)
	ascending_sign = round(np.sign(math.sin(rad_angle+3*math.pi/4)))
	is_reversed = (ascending_sign == -1)

	# Sort left side using the sort_axis
	comb = list(zip(left_pins, left_pin_ids))
	comb.sort(key=lambda x: x[0][sort_axis], reverse=is_reversed)
	left_pins, left_pin_ids = zip(*comb)

	# Sort right side using the sort_axis
	comb = list(zip(right_pins, right_pin_ids))
	comb.sort(key=lambda x: x[0][sort_axis], reverse=is_reversed)
	right_pins, right_pin_ids = zip(*comb)


	for j, pin_id in enumerate(left_pin_ids):		
#		print(j)
#		print(f"Left pins order: {left_pins[j]}. Index: {j*2}. Pin_id: {pin_id}. Track id: {pin_track_links[pin_id]}")
		if j==64: continue
		if j>64:
			print(f"ERROR: There are too many pins! {j}/64")

		if pin_track_links[pin_id] is None or pin_track_links[pin_id]<0:
			Track_id=-1
		else:
			Track_id=pin_track_to_pad_track_links[pin_track_links[pin_id]]

		if Track_id is None:
			Track_id=-1
			
		mapping_df.loc[len(mapping_df)] = {'Track_id': Track_id, 'Connector_id': i, 'Pin_id': j*2} #Left tracks are even pins???
		mapping_df.loc[mapping_df['Track_id'] == Track_id, 'Pad_pos_x'] = merged_tracks[Track_id][-1][0]
		mapping_df.loc[mapping_df['Track_id'] == Track_id, 'Pad_pos_y'] = merged_tracks[Track_id][-1][1]

	for j, pin_id in enumerate(right_pin_ids):
#		print(f"Right pins order: {right_pins[j]}. Index: {2*(j-1)+3}. Pin_id: {pin_id}. Track id: {pin_track_links[pin_id]}")
		if j==0: continue
		if j>64:
			print("ERROR: There are too many pins!")

		if pin_track_links[pin_id] is None or pin_track_links[pin_id]<0:
			Track_id=-1
		else:
			Track_id=pin_track_to_pad_track_links[pin_track_links[pin_id]]

		if Track_id is None:
			Track_id=-1

		mapping_df.loc[len(mapping_df)] = {'Track_id': Track_id, 'Connector_id': i, 'Pin_id': 2*(j-1) + 1} #Right tracks are odd pins???

		mapping_df.loc[mapping_df['Track_id'] == Track_id, 'Pad_pos_x'] = merged_tracks[Track_id][-1][0]
		mapping_df.loc[mapping_df['Track_id'] == Track_id, 'Pad_pos_y'] = merged_tracks[Track_id][-1][1]

#SECOND GUI########################################
#Reset plot
ax.clear()

x,y=[],[]
for r in full_merged_tracks:
	x.append(r[-1][0])
	y.append(r[-1][1])

ax.plot(x,y,'s',markersize=3)
ax.set_title(f"Select pad area (1/{n_pad_areas})")
ax.set_xlabel('X [mm]')
ax.set_ylabel('Y [mm]')

x,y=[],[]
for idx, row in mapping_df.iterrows():
	if row["Track_id"] is None: continue
	track_id=int(row["Track_id"])

	track=merged_tracks[track_id]
	x.append(track[-1][0])
	y.append(track[-1][1])

ax.plot(x,y,'v',markersize=3,color='orange')

#Start GUI
root = tk.Tk()
root.title("Pad area selection GUI")

control_frame = ttk.Frame(root)
control_frame.pack(side=tk.TOP, fill=tk.X)

canvas = FigureCanvasTkAgg(fig, master=root)
canvas_widget = canvas.get_tk_widget()
canvas_widget.pack(fill=tk.BOTH, expand=True)

selector = GUI_utils.RectSelectorPads(ax)

def set_move(): selector.set_mode('move'); mode_var.set("Move")
def set_rotate(): selector.set_mode('rotate'); mode_var.set("Rotate")
def Flip(): selector.flip(); #flipped=True

ttk.Button(control_frame, text="Move", command=set_move).pack(side=tk.LEFT, padx=4)
ttk.Button(control_frame, text="Rotate", command=set_rotate).pack(side=tk.LEFT, padx=4)
ttk.Button(control_frame, text="Filp", command=Flip).pack(side=tk.LEFT, padx=4)
mode_var = tk.StringVar(value="Move")
ttk.Label(control_frame, textvariable=mode_var, width=8).pack(side=tk.LEFT, padx=6)

output = tk.Text(root, height=6)
output.pack(fill=tk.X)

areas=[]
def on_enter(event=None):
	verts, ang = selector.get_rectangle_parameters()
	if verts is None:
		output.insert(tk.END, "No rectangle.\n")
		return

	areas.append({
		"angle": ang,
		"flipped": selector.flipped,
		"vertices": verts.copy()
	})

	next_area=len(areas)+1
	ax.set_title(f"Select pad area ({next_area}/{n_pad_areas})")
	canvas.draw()


	if len(areas) >= n_pad_areas:
		root.quit()
		root.destroy()

root.bind("<Return>", on_enter)
root.mainloop()


#Loop over pad areas
for i,pa in enumerate(areas,start=1):  #Loop over connectors

	print("Pad area "+str(i)+":")
	print("	Vertices: "+str(pa['vertices']))
	print("	Angle: "+str(pa["angle"]))
	print("	Flipped: "+str(pa["flipped"]))
	print()

	#Selection boundaries (lax)
	x_left = min(item[0] for item in pa['vertices'])
	x_right = max(item[0] for item in pa['vertices'])

	y_low = min(item[1] for item in pa['vertices'])
	y_high = max(item[1] for item in pa['vertices'])

	for idx, row in mapping_df.iterrows():
		track_id=row["Track_id"]
		if track_id is None:
			print("ERROR: Track is none!")
			exit()
		if not pd.isna(row["Pad_area"]): continue
		if track_id<0:
			mapping_df.at[idx, "Pad_pos_x"]=None
			mapping_df.at[idx, "Pad_pos_y"]=None
		else:
			x=row["Pad_pos_x"]
			y=row["Pad_pos_y"]

			if x>x_right or x<x_left: continue
			if y>y_high or y<y_low: continue
	
			mapping_df.loc[mapping_df["Connector_id"] == row["Connector_id"], "Pad_area"] = i

#	print(mapping_df)

#SAVE OUTPUT#########################################################################
out_df=mapping_df
#print(out_df)
out_df = out_df.sort_values(by=['Connector_id','Pin_id'])

group_keys = ["Connector_id"]

out_df = pd.DataFrame(
    out_df.groupby(group_keys, sort=False)
      .apply(lambda g: {
          **{k: v for k, v in zip(group_keys, g.name if isinstance(g.name, tuple) else (g.name,))},
          **{col: g[col].tolist() for col in out_df.columns if col not in group_keys}
      })
      .tolist()
)

out_df['Pin_id'] = out_df['Pin_id'].astype(str)
out_df['Pad_pos_x'] = out_df['Pad_pos_x'].astype(str)
out_df['Pad_pos_y'] = out_df['Pad_pos_y'].astype(str)
out_df['Pad_area'] = out_df['Pad_area'].astype(str)
out_df['Track_id'] = out_df['Track_id'].astype(str)

json_str=out_df.to_json(orient="records",indent=2)
json_str=re.sub(r'\"\[',' [',json_str)
json_str=re.sub(r']\"',']',json_str)

with open(out_dir+"/"+out_file, "w") as f:
	f.write(json_str)
