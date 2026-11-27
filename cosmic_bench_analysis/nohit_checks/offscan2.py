import sys; sys.path.insert(0,'.')
import numpy as np, pandas as pd, p2_study_cache as sc, p2_study_common as cm, p2_mapping as pmap
from scipy.spatial import cKDTree
hits, rays, g, cfg = sc.load('det2_long_repos')
d34=cfg.out_dir('34_clustering_matching'); A=np.load(d34+'/alignment.npz'); T=cm.make_T(g,R=A['R'],s=A['s'],pad_mean=A['pad_mean'],m3_mean=A['m3_mean'])
ra=pd.read_pickle(d34+'/rays_aligned.pkl'); xy=ra[['x','y']].to_numpy(); th=ra.theta.to_numpy()
nohit=~ra.eventId.isin(hits.eventId).to_numpy()
pil=pmap.load_pillars(cfg.MASK_GBR_PATH); px,py=T.apply(pil.x.to_numpy(),pil.y.to_numpy())
sm=np.c_[px,py][~pil.big.to_numpy()]; big=np.c_[px,py][pil.big.to_numpy()]
db=cKDTree(big).query(xy)[0]; tr=cKDTree(sm)
sel0=(ra.sd>=4).to_numpy()&(db>10)&(th<8)
xm,ym=np.median(xy[sel0],0)
grid=np.arange(-2.0,2.01,0.1)
for qn,(sx,sy) in enumerate([(-1,-1),(-1,1),(1,-1),(1,1)]):
    sel=sel0&(np.sign(xy[:,0]-xm)==sx)&(np.sign(xy[:,1]-ym)==sy)
    P=xy[sel]; nh=nohit[sel]; base=nh.mean(); best=(-9,0,0,0,0)
    zs=[]
    for dx in grid:
        for dy in grid:
            d,_=tr.query(P-np.array([dx,dy])); k=d<0.5
            if k.sum()<20: continue
            z=(nh[k].mean()-base)/np.sqrt(base*(1-base)/k.sum()); zs.append(z)
            if z>best[0]: best=(z,dx,dy,nh[k].mean(),k.sum())
    print(f'quadrant {qn} (x{"+" if sx>0 else "-"} y{"+" if sy>0 else "-"}): n {sel.sum()} base {100*base:.2f}%  best z {best[0]:.1f} at ({best[1]:+.1f},{best[2]:+.1f}) rate {100*best[3]:.1f}% n {best[4]}; z sd {np.std(zs):.2f}')
