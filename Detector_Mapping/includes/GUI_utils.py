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

class RectSelectorInteractive:
		def __init__(self, ax):
				self.ax = ax
				self.canvas = ax.figure.canvas
				self.mode = 'move'
				self.center = None
				self.w = None
				self.h = None
				self.angle = 0.0
				self.patch = None
				self.corner_labels = []
				self.dragging = False
				self.creating = False
				self.press_event = None
				self.start_center = None
				self.start_angle = None
				self.start_vec = None
				self.cid_press = self.canvas.mpl_connect('button_press_event', self.on_press)
				self.cid_release = self.canvas.mpl_connect('button_release_event', self.on_release)
				self.cid_motion = self.canvas.mpl_connect('motion_notify_event', self.on_motion)
				self.flipped=False
				self.names = ["0","1","127","126"]


		def flip(self):
			#print("Flippin logic")
			if not self.flipped:
				self.names = ["1","0","126","127"]
				self.flipped=True
				#print("flip yes")
			else:
				self.names = ["0","1","127","126"]
				self.flipped=False
				#print("flip no")
			
			self.update_labels()
			self.canvas.draw_idle()


		def contains_point(self, event):
				from matplotlib.path import Path
				if self.patch is None:
						return False
				if event.xdata is None or event.ydata is None:
						return False
				verts = self._rect_vertices()
				return Path(verts).contains_point((event.xdata, event.ydata))

		def _rect_vertices(self, center=None, w=None, h=None, angle=None):
				center = self.center if center is None else center
				w = self.w if w is None else w
				h = self.h if h is None else h
				angle = self.angle if angle is None else angle
				dx = w/2
				dy = h/2
				corners = np.array([[-dx,-dy],[dx,-dy],[dx,dy],[-dx,dy]])
				c = np.cos(angle)
				s = np.sin(angle)
				R = np.array([[c,-s],[s,c]])
				return corners @ R.T + center

		# ---- Drawing & labels ----
		def clear_rect(self):
				if self.patch is not None:
						self.patch.remove()
				for lbl in self.corner_labels:
						lbl.remove()
				self.patch = None
				self.corner_labels = []
				self.canvas.draw_idle()

		def draw_rect(self):
				verts = self._rect_vertices()
				if self.patch is None:
						#self.patch = Polygon(verts, closed=True, fill=False, linewidth=2)
						self.patch = Polygon(verts, closed=True, fill=False, linewidth=2, zorder=999, color='black')
						self.ax.add_patch(self.patch)
				else:
						self.patch.set_xy(verts)
				self.update_labels()
				self.canvas.draw_idle()

		def update_labels(self):
		    verts = self._rect_vertices()
		    names = self.names
		
		    # Offset before rotation (relative to corner)
		    # You might need to tweak 0.02 depending on your axis scale
		    base_offset = np.array([0.02, 0.02]) 
		    c = np.cos(self.angle)
		    s = np.sin(self.angle)
		    R = np.array([[c, -s],[s, c]])
		
		    rotated_offsets = [R @ base_offset for _ in range(4)]
		
		    # Create labels if they don't exist
		    if len(self.corner_labels) == 0:
		        for (x, y), label, off in zip(verts, names, rotated_offsets):
		            t = self.ax.text(
		                x + off[0], y + off[1], label,
		                ha='center', va='center', fontsize=9,
		                fontweight='bold',       # Make text thicker
		                color='black',           # Text font color
		                zorder=1000,             # HIGHER than the rectangle (999) so it sits on top
		                bbox=dict(               # This creates the "PowerPoint-like" fill
		                    boxstyle="round,pad=0.3", # Rounded corners, padding around text
		                    facecolor="yellow",       # Background color (yellow, white, cyan, etc.)
		                    #edgecolor="black",        # Border color of the box
		                    linewidth=1,
		                    alpha=0.8                 # 0.0 is transparent, 1.0 is opaque
		                )
		            )
		            self.corner_labels.append(t)
		            
		    # Update positions if they already exist
		    else:
		        for lbl, (x, y), name, off in zip(self.corner_labels, verts, names, rotated_offsets):
		            lbl.set_position((x + off[0], y + off[1]))
		            lbl.set_text(name)
		            # Optional: Update rotation of the text box itself if you want it to spin with the rect
		            # lbl.set_rotation(np.rad2deg(self.angle))



		def on_press(self, event):
				if event.inaxes != self.ax or event.button != 1:
						return

				# If rectangle exists but click is OUTSIDE → remake new rectangle
				if self.patch is not None and not self.contains_point(event):
						# clear old
						self.clear_rect()
						# start creating new
						self.creating = True
						self.press_event = event
						return

				# If no rectangle exists — create
				if self.patch is None:
						self.creating = True
						self.press_event = event
						return

				# If clicking inside existing rect — move/rotate
				self.dragging = True
				self.press_event = event
				self.start_center = self.center.copy()
				self.start_angle = self.angle
				self.start_vec = np.array([event.xdata - self.center[0],
																	 event.ydata - self.center[1]])
				if np.linalg.norm(self.start_vec) < 1e-8:
						self.start_vec = np.array([1e-8,0])(self, event)
				if event.inaxes != self.ax or event.button != 1:
						return

				# If no rectangle exists — start creating one
				if self.patch is None:
						self.creating = True
						self.press_event = event
						return

				# If clicking inside existing rect — move/rotate
				if self.contains_point(event):
						self.dragging = True
						self.press_event = event
						self.start_center = self.center.copy()
						self.start_angle = self.angle
						self.start_vec = np.array([event.xdata - self.center[0],
																			 event.ydata - self.center[1]])
						if np.linalg.norm(self.start_vec) < 1e-8:
								self.start_vec = np.array([1e-8,0])


		def on_motion(self, event):
				if event.inaxes != self.ax:
						return

				# --- Creating new rectangle by drag ---
				if self.creating:
						x0,y0 = self.press_event.xdata, self.press_event.ydata
						x1,y1 = event.xdata, event.ydata
						cx = (x0+x1)/2
						cy = (y0+y1)/2
						w = abs(x1-x0)
						h = abs(y1-y0)
						self.center = np.array([cx,cy])
						self.w = w
						self.h = h
						self.angle = 0.0
						self.draw_rect()
						return

				# --- Move or rotate existing rectangle ---
				if not self.dragging:
						return
				if event.xdata is None:
						return
				dx = event.xdata - self.press_event.xdata
				dy = event.ydata - self.press_event.ydata

				if self.mode == 'move':
						self.center = self.start_center + np.array([dx,dy])
				else:
						cur_vec = np.array([event.xdata - self.center[0],
																event.ydata - self.center[1]])
						if np.linalg.norm(cur_vec) < 1e-8:
								return
						cross = np.cross(self.start_vec, cur_vec)
						dot = np.dot(self.start_vec, cur_vec)
						da = np.arctan2(cross, dot)
						self.angle = self.start_angle + da
				self.draw_rect()

		def on_release(self, event):
				if self.creating:
						self.creating = False
						# Finalize by erasing old rectangle if any
						# (clear_rect already handled before creation)
						return
				if self.dragging:
						self.dragging = False

		# ---- API ----
		def set_mode(self, mode):
				self.mode = mode

		def get_rectangle_parameters(self):
				if self.patch is None:
						return None, None
				verts = self._rect_vertices()
				return verts, (np.rad2deg(self.angle) % 360)

class RectSelectorPads:
		def __init__(self, ax):
				self.ax = ax
				self.canvas = ax.figure.canvas
				self.mode = 'move'
				self.center = None
				self.w = None
				self.h = None
				self.angle = 0.0
				self.patch = None
				self.corner_labels = []
				self.dragging = False
				self.creating = False
				self.press_event = None
				self.start_center = None
				self.start_angle = None
				self.start_vec = None
				self.cid_press = self.canvas.mpl_connect('button_press_event', self.on_press)
				self.cid_release = self.canvas.mpl_connect('button_release_event', self.on_release)
				self.cid_motion = self.canvas.mpl_connect('motion_notify_event', self.on_motion)
				self.flipped=False
				self.names = ["0","1","127","126"]

		def contains_point(self, event):
				from matplotlib.path import Path
				if self.patch is None:
						return False
				if event.xdata is None or event.ydata is None:
						return False
				verts = self._rect_vertices()
				return Path(verts).contains_point((event.xdata, event.ydata))

		def _rect_vertices(self, center=None, w=None, h=None, angle=None):
				center = self.center if center is None else center
				w = self.w if w is None else w
				h = self.h if h is None else h
				angle = self.angle if angle is None else angle
				dx = w/2
				dy = h/2
				corners = np.array([[-dx,-dy],[dx,-dy],[dx,dy],[-dx,dy]])
				c = np.cos(angle)
				s = np.sin(angle)
				R = np.array([[c,-s],[s,c]])
				return corners @ R.T + center

		# ---- Drawing & labels ----
		def clear_rect(self):
				if self.patch is not None:
						self.patch.remove()
				for lbl in self.corner_labels:
						lbl.remove()
				self.patch = None
				self.corner_labels = []
				self.canvas.draw_idle()

		def draw_rect(self):
				verts = self._rect_vertices()
				if self.patch is None:
						#self.patch = Polygon(verts, closed=True, fill=False, linewidth=2)
						self.patch = Polygon(verts, closed=True, fill=False, linewidth=2, zorder=999, color='red',linestyle='--')
						self.ax.add_patch(self.patch)
				else:
						self.patch.set_xy(verts)
				self.canvas.draw_idle()

		def on_press(self, event):
				if event.inaxes != self.ax or event.button != 1:
						return

				# If rectangle exists but click is OUTSIDE → remake new rectangle
				if self.patch is not None and not self.contains_point(event):
						# clear old
						self.clear_rect()
						# start creating new
						self.creating = True
						self.press_event = event
						return

				# If no rectangle exists — create
				if self.patch is None:
						self.creating = True
						self.press_event = event
						return

				# If clicking inside existing rect — move/rotate
				self.dragging = True
				self.press_event = event
				self.start_center = self.center.copy()
				self.start_angle = self.angle
				self.start_vec = np.array([event.xdata - self.center[0],
																	 event.ydata - self.center[1]])
				if np.linalg.norm(self.start_vec) < 1e-8:
						self.start_vec = np.array([1e-8,0])(self, event)
				if event.inaxes != self.ax or event.button != 1:
						return

				# If no rectangle exists — start creating one
				if self.patch is None:
						self.creating = True
						self.press_event = event
						return

				# If clicking inside existing rect — move/rotate
				if self.contains_point(event):
						self.dragging = True
						self.press_event = event
						self.start_center = self.center.copy()
						self.start_angle = self.angle
						self.start_vec = np.array([event.xdata - self.center[0],
																			 event.ydata - self.center[1]])
						if np.linalg.norm(self.start_vec) < 1e-8:
								self.start_vec = np.array([1e-8,0])


		def on_motion(self, event):
				if event.inaxes != self.ax:
						return

				# --- Creating new rectangle by drag ---
				if self.creating:
						x0,y0 = self.press_event.xdata, self.press_event.ydata
						x1,y1 = event.xdata, event.ydata
						cx = (x0+x1)/2
						cy = (y0+y1)/2
						w = abs(x1-x0)
						h = abs(y1-y0)
						self.center = np.array([cx,cy])
						self.w = w
						self.h = h
						self.angle = 0.0
						self.draw_rect()
						return

				# --- Move or rotate existing rectangle ---
				if not self.dragging:
						return
				if event.xdata is None:
						return
				dx = event.xdata - self.press_event.xdata
				dy = event.ydata - self.press_event.ydata

				if self.mode == 'move':
						self.center = self.start_center + np.array([dx,dy])
				else:
						cur_vec = np.array([event.xdata - self.center[0],
																event.ydata - self.center[1]])
						if np.linalg.norm(cur_vec) < 1e-8:
								return
						cross = np.cross(self.start_vec, cur_vec)
						dot = np.dot(self.start_vec, cur_vec)
						da = np.arctan2(cross, dot)
						self.angle = self.start_angle + da
				self.draw_rect()

		def on_release(self, event):
				if self.creating:
						self.creating = False
						# Finalize by erasing old rectangle if any
						# (clear_rect already handled before creation)
						return
				if self.dragging:
						self.dragging = False

		# ---- API ----
		def set_mode(self, mode):
				self.mode = mode

		def get_rectangle_parameters(self):
				if self.patch is None:
						return None, None
				verts = self._rect_vertices()
				return verts, (np.rad2deg(self.angle) % 360)

