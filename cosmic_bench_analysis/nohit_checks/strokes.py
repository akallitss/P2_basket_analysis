import re, numpy as np, sys
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
B='/local/home/ak271430/Documents/PostDocSaclay/Detector_Drawings/Version_Apr26/Insulation_masks/V1/'
def strokes(f):
    s=open(f).read().split('\n'); ap={}; cur=None; x=y=0; segs=[]
    for l in s:
        m=re.match(r'%ADD(\d+)C,([\d.]+)',l)
        if m: ap['D'+m.group(1)]=float(m.group(2)); continue
        m=re.match(r'^(D\d+)\*$',l)
        if m and m.group(1) in ap: cur=m.group(1); continue
        m=re.match(r'^X(-?\d+)Y(-?\d+)(?:I(-?\d+)J(-?\d+))?D0([123])\*',l)
        if m:
            nx,ny=int(m.group(1))*1e-6,int(m.group(2))*1e-6
            if m.group(5)=='1' and m.group(3) is None: segs.append((cur,ap.get(cur),x,y,nx,ny))
            x,y=nx,ny
    return segs
for f,ax_title in [(B+'P2_BASKET-Mask_M1_V1.gbr','M1'),(B+'P2_BASKET-Mask_M2_V1.gbr','M2')]:
    S=strokes(f); 
    import collections; c=collections.Counter((a,w) for a,w,*_ in S); print(ax_title, c)
    L=np.array([[s[2],s[3],s[4],s[5]] for s in S]); W=np.array([s[1] for s in S])
    ln=np.hypot(L[:,2]-L[:,0],L[:,3]-L[:,1]); print('  stroke length median', np.median(ln).round(3), 'total length m', ln.sum()/1000, 'x', L[:,0].min().round(1), L[:,0].max().round(1), 'y', L[:,1].min().round(1), L[:,1].max().round(1))
    fig,axs=plt.subplots(1,2,figsize=(14,7))
    for ax,lim in zip(axs,[None,(300,330,160,190)]):
        for w,col in zip(sorted(set(W)),['k','r','b','g']):
            k=W==w; ax.add_collection(LineCollection(L[k].reshape(-1,2,2),colors=col,linewidths=0.6,label=f'pen {w} mm ({k.sum()})'))
        if lim: ax.set_xlim(lim[0],lim[1]); ax.set_ylim(lim[2],lim[3])
        else: ax.autoscale()
        ax.set_aspect('equal'); ax.legend(fontsize=8); ax.set_title(ax_title+' strokes')
    fig.savefig(f'/local/home/ak271430/.claude-dylan/jobs/0b89599e/tmp/strokes_{ax_title}.png',dpi=90)
