import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.text import Text
from matplotlib.path import Path
from scipy.spatial import KDTree
import math

def dot(a, b):
	return a[0]*b[0] + a[1]*b[1]

def subtract(a, b):
	return (a[0]-b[0], a[1]-b[1])

def length(v):
	return math.hypot(v[0], v[1])

def clamp(x, a, b):
	return max(a, min(b, x))

def sample_polyline(seg, spacing):
	"""Return [points] sampled every `spacing` distance along the polyline."""
	if len(seg) == 1:
		return [tuple(seg[0])]

	sampled = [tuple(seg[0])]
	acc = 0.0

	for i in range(len(seg) - 1):
		x1, y1 = seg[i]
		x2, y2 = seg[i+1]

		dx = x2 - x1
		dy = y2 - y1
		seg_len = math.hypot(dx, dy)

		if seg_len == 0:
			continue

		# Walk along this segment at fixed spacing
		t = spacing - acc
		while t < seg_len:
			px = x1 + dx * (t / seg_len)
			py = y1 + dy * (t / seg_len)
			sampled.append((px, py))
			t += spacing

		acc = seg_len - (t - spacing)

		# Add original vertex
		sampled.append((x2, y2))

	return sampled

def build_kdtree_sampled(segments, sample_spacing):
	points = []
	seg_map = []

	for idx, seg in enumerate(segments):
		for p in sample_polyline(seg, spacing=sample_spacing):
			points.append(p)
			seg_map.append(idx)

	tree = KDTree(points)
	return tree, seg_map

def inch_to_mm(x):
	return x*25.4

def point_segment_distance(p, a, b):
	"""Distance from point p to line segment ab."""
	ax, ay = a
	bx, by = b
	px, py = p

	AB = (bx - ax, by - ay)
	AP = (px - ax, py - ay)
	AB_len2 = AB[0]*AB[0] + AB[1]*AB[1]

	if AB_len2 == 0:
		return math.dist(p, a)

	t = (AP[0]*AB[0] + AP[1]*AB[1]) / AB_len2
	t = max(0, min(1, t))

	cx = ax + t * AB[0]
	cy = ay + t * AB[1]
	return math.dist(p, (cx, cy))


def polyline_distance(seg1, seg2):
	"""Minimum distance between two polylines."""
	min_d = float("inf")
	for i in range(len(seg1) - 1):
		a1, a2 = seg1[i], seg1[i + 1]
		for j in range(len(seg2) - 1):
			b1, b2 = seg2[j], seg2[j + 1]

			min_d = min(
				min_d,
				point_segment_distance(b1, a1, a2),
				point_segment_distance(b2, a1, a2),
				point_segment_distance(a1, b1, b2),
				point_segment_distance(a2, b1, b2)
			)
	return min_d

def merge_two(seg1, seg2):
	endpoint_pairs = [
		(seg1[0],  seg2[0],  's1_front', 's2_front'),
		(seg1[0],  seg2[-1], 's1_front', 's2_back'),
		(seg1[-1], seg2[0],  's1_back',  's2_front'),
		(seg1[-1], seg2[-1], 's1_back',  's2_back'),
	]

	p1, p2, t1, t2 = min(
		endpoint_pairs,
		key=lambda t: math.dist(t[0], t[1])
	)

	# Orient seg1 so that the matching endpoint is at the end
	if t1 == 's1_front':
		seg1 = seg1[::-1]
	# Orient seg2 so matching endpoint is at the start
	if t2 == 's2_back':
		seg2 = seg2[::-1]

	# Remove duplicate if needed
	if seg1[-1] == seg2[0]:
		return seg1 + seg2[1:]
	else:
		return seg1 + seg2

def merge_two_closest_endpoints(seg1, seg2):
	"""
	Merge two polyline segments so that the closest endpoints touch.
	Always returns: seg1_tail → seg2_head concatenation.
	"""

	# Extract endpoints
	a1, a2 = seg1[0], seg1[-1]
	b1, b2 = seg2[0], seg2[-1]

	# List all endpoint pairings
	candidates = [
		(a1, b1, 's1_front', 's2_front'),
		(a1, b2, 's1_front', 's2_back'),
		(a2, b1, 's1_back',  's2_front'),
		(a2, b2, 's1_back',  's2_back'),
	]

	# Pick the pair with minimum distance
	p1, p2, t1, t2 = min(candidates, key=lambda t: math.dist(t[0], t[1]))

	# --- Orient seg1 so that p1 is its *tail* ---
	if t1 == 's1_front':		  # p1 is seg1 head → need to flip
		seg1 = list(reversed(seg1))

	# --- Orient seg2 so that p2 is its *head* ---
	if t2 == 's2_back':		   # p2 is seg2 tail → need to flip
		seg2 = list(reversed(seg2))

	# Now seg1[-1] and seg2[0] are the closest endpoints.

	# Avoid duplicate vertex if they match
	if seg1[-1] == seg2[0]:
		return seg1 + seg2[1:]
	else:
		return seg1 + seg2

def build_kdtree(segments):
	"""Returns KD-tree and mapping from KD-point index → segment index."""
	points = []
	seg_index = []

	for i, seg in enumerate(segments):
		# Use both endpoints
		points.append(seg[0])
		seg_index.append(i)

		points.append(seg[-1])
		seg_index.append(i)

	tree = KDTree(points)
	return tree, seg_index

def merge_segments_kdtree(
	segments,
	kdtree_threshold,
	merge_threshold
):
	segments = [list(seg) for seg in segments]  # shallow copy
	changed = True

	while changed:
		changed = False

		# Build KD-tree of endpoints
		tree, index_map = build_kdtree(segments)

		used = set()
		new_segments = []

		for i in range(len(segments)):
			if i in used:
				continue

			merged = segments[i]

			# Endpoints of this segment
			endpoints = [merged[0], merged[-1]]

			# Candidate segments whose endpoints are within KD-tree threshold
			candidates = set()
			for p in endpoints:
				for idx in tree.query_ball_point(p, r=kdtree_threshold):
					seg_id = index_map[idx]
					if seg_id != i:
						candidates.add(seg_id)

			# Now perform *true* geometric distance test on candidates
			for j in sorted(candidates):
				if j in used:
					continue

				if polyline_distance(merged, segments[j]) < merge_threshold:
					merged = merge_two(merged, segments[j])
					used.add(j)
					changed = True

			new_segments.append(merged)

		segments = new_segments

	return segments

def merge_segments_kdtree_point(segments, threshold=0.2, verbose=False):
	"""
	Merge segments or polylines whose endpoints are closer than `threshold`,
	using a KD-tree for efficient neighbor search.

	Handles:
	  - Segments with 2 or more points
	  - Empty or invalid segments (silently skipped)
	  - All endpoint connection types (head↔tail, head↔head, etc.)

	Parameters
	----------
	segments : list of list of [x, y]
		Each element is a list of at least two coordinates.
	threshold : float
		Distance threshold for connecting endpoints.
	verbose : bool
		If True, prints info about filtered or merged segments.

	Returns
	-------
	merged_tracks : list of list of [x, y]
	"""

	# --- Step 1: Filter invalid segments ---
	valid_segments = []
	for i, s in enumerate(segments):
		if not isinstance(s, (list, tuple)) or len(s) < 2:
			if verbose:
				print(f"⚠️ Skipping invalid segment at index {i}: {s}")
			continue
		# Ensure each point is at least 2D numeric
		s_clean = [np.array(p[:2], dtype=float).tolist() for p in s if len(p) >= 2]
		if len(s_clean) >= 2:
			valid_segments.append(s_clean)
		elif verbose:
			print(f"⚠️ Skipping degenerate segment at index {i}: {s}")

	segments = valid_segments
	n = len(segments)
	if n == 0:
		if verbose:
			print("No valid segments to merge.")
		return []

	# --- Step 2: Extract heads and tails (first and last points) ---
	heads = np.array([s[0] for s in segments], dtype=float)
	tails = np.array([s[-1] for s in segments], dtype=float)

	# --- Step 3: Build KD-tree for endpoints ---
	endpoints = np.vstack([
		np.hstack([heads, np.arange(n)[:, None], np.zeros((n, 1))]),  # head=0
		np.hstack([tails, np.arange(n)[:, None], np.ones((n, 1))])	# tail=1
	])
	coords = endpoints[:, :2]
	seg_ids = endpoints[:, 2].astype(int)

	tree = KDTree(coords)

	# --- Step 4: Build adjacency graph (segments connected by close endpoints) ---
	adj = {i: set() for i in range(n)}
	for idx, (x, y, seg_id, _) in enumerate(endpoints):
		seg_id = int(seg_id)
		nearby = tree.query_ball_point([x, y], threshold)
		for j in nearby:
			other_seg = int(seg_ids[j])
			if other_seg != seg_id:
				adj[seg_id].add(other_seg)
				adj[other_seg].add(seg_id)

	# --- Step 5: Find connected components ---
	visited = set()
	groups = []
	for i in range(n):
		if i not in visited:
			stack = [i]
			component = []
			while stack:
				node = int(stack.pop())
				if node not in visited:
					visited.add(node)
					component.append(node)
					stack.extend(int(x) for x in adj[node] - visited)
			groups.append(component)

	# --- Step 6: Merge connected segments into ordered tracks ---
	merged_tracks = []
	for group in groups:
		used = set()
		pts = list(segments[group[0]])
		used.add(group[0])
		changed = True
		while changed:
			changed = False
			for i in group:
				i = int(i)
				if i in used:
					continue
				s = list(segments[i])

				# Compute endpoint distances
				d_tt = np.linalg.norm(np.array(pts[-1]) - np.array(s[0]))  # tail→head
				d_th = np.linalg.norm(np.array(pts[-1]) - np.array(s[-1])) # tail→tail
				d_ht = np.linalg.norm(np.array(pts[0]) - np.array(s[0]))   # head→head
				d_hh = np.linalg.norm(np.array(pts[0]) - np.array(s[-1]))  # head→tail

				if d_tt <= threshold:
					pts.extend(s[1:])
				elif d_th <= threshold:
					pts.extend(s[::-1][1:])
				elif d_ht <= threshold:
					pts = s[::-1][:-1] + pts
				elif d_hh <= threshold:
					pts = s[:-1] + pts
				else:
					continue
				used.add(i)
				changed = True
		merged_tracks.append(pts)

	if verbose:
		print(f"✅ Merged {n} segments into {len(merged_tracks)} tracks.")
	return merged_tracks

def merge_segments_sampled_kdtree(
	segments,
	sample_spacing,
	kdtree_threshold,
	merge_threshold
):
	segments = [list(seg) for seg in segments]
	changed = True

	while changed:
		changed = False

		# Build KD-tree of sampled points
		tree, seg_map = build_kdtree_sampled(segments, sample_spacing)

		used = set()
		new_segments = []

		for i in range(len(segments)):
			if i in used:
				continue

			merged = segments[i]

			# Get candidate segments by sampled points proximity
			samples = sample_polyline(merged, spacing=sample_spacing)
			candidate_ids = set()

			for p in samples:
				nearby = tree.query_ball_point(p, r=kdtree_threshold)
				for idx in nearby:
					sid = seg_map[idx]
					if sid != i:
						candidate_ids.add(sid)

			# Now check true minimum distance on candidates
			for j in sorted(candidate_ids):
				if j in used:
					continue
				if polyline_distance(merged, segments[j]) < merge_threshold:
					#merged = merge_two(merged, segments[j])
					merged = merge_two_closest_endpoints(merged, segments[j])
					used.add(j)
					changed = True

			new_segments.append(merged)

		segments = new_segments

	return segments

def build_track_vertex_kdtree(tracks):
	"""
	Build a KD-tree of all track vertices.
	Returns KDTree + parallel array mapping vertex index -> track index.
	"""
	all_pts = []
	pt_to_track = []

	for ti, track in enumerate(tracks):
		for v in track:
			all_pts.append((v[0], v[1]))
			pt_to_track.append(ti)

	tree = KDTree(all_pts)
	return tree, pt_to_track


def assign_points_to_tracks(points, tracks, threshold):
	"""
	For each point, find the nearest track vertex within threshold.
	Returns: list of track indices or None for each point.
	"""

	tree, pt_to_track = build_track_vertex_kdtree(tracks)
	result = []

	for p in points:
		# get candidate vertices within threshold
		candidates = tree.query_ball_point(p, threshold)

		if not candidates:
			#result.append(-1)
			result.append(None)
			continue

		# find closest among the candidates
		min_dist = float('inf')
		best_track = None

		for idx in candidates:
			vx, vy = tree.data[idx]
			d = math.hypot(vx - p[0], vy - p[1])
			if d < min_dist:
				min_dist = d
				best_track = pt_to_track[idx]

		result.append(best_track)

	return result

def build_track_endpoint_kdtree(tracks):
	pts = []		   # KD-tree coordinates
	endpoint_info = [] # (track_index, endpoint_type) for safety check

	for ti, track in enumerate(tracks):
		if len(track) == 0:
			continue

		# First endpoint
		pts.append(tuple(track[0]))
		endpoint_info.append((ti, 'start'))

		# Last endpoint
		if len(track) > 1:
			pts.append(tuple(track[-1]))
			endpoint_info.append((ti, 'end'))

	tree = KDTree(pts)
	return tree, endpoint_info


def assign_points_to_track_endpoints(points, tracks, threshold):
	tree, endpoint_info = build_track_endpoint_kdtree(tracks)
	results = []

	for p in points:
		nearby_idx = tree.query_ball_point(p, threshold)

		if not nearby_idx:
			#results.append(-1)
			results.append(None)
			continue

		# Find the closest valid endpoint
		best_track = None
		best_dist = float('inf')

		for idx in nearby_idx:
			# We *know* this is an endpoint
			track_id, endpoint_type = endpoint_info[idx]

			ex, ey = tree.data[idx]
			dist = math.hypot(ex - p[0], ey - p[1])

			if dist < best_dist:
				best_dist = dist
				best_track = track_id

		results.append(best_track)

	return results


def extend_track(points):
	"""
	Reorder a list of 2D points so that:
	- The endpoints are the two farthest-apart points in the track.
	- The interior points are ordered by greedy nearest-neighbor path.
	"""

	n = len(points)
	if n <= 2:
		return points[:]   # Already simplest possible track

	# --- 1. Find the farthest pair of points (track endpoints) ---
	max_dist = -1
	start_idx, end_idx = 0, 1

	for i in range(n):
		for j in range(i + 1, n):
			d = math.dist(points[i], points[j])
			if d > max_dist:
				max_dist = d
				start_idx, end_idx = i, j

	start = points[start_idx]
	end = points[end_idx]

	# --- 2. Build path: greedy nearest-neighbor from the start ---
	unvisited = set(range(n))
	unvisited.remove(start_idx)
	ordered = [start]
	current = start

	while unvisited:
		# pick the nearest unvisited point
		next_idx = min(
			unvisited,
			key=lambda k: math.dist(current, points[k])
		)
		ordered.append(points[next_idx])
		unvisited.remove(next_idx)
		current = points[next_idx]

	# --- 3. Ensure the final point in the greedy path is the far endpoint ---
	# If not, reverse the whole track.
	if ordered[-1] != end:
		# If the far endpoint is somewhere in the sequence,
		# reversing will make the track start at 'end' and end at 'start'.
		ordered.reverse()

	return ordered
