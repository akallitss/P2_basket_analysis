import sys; sys.path.insert(0,'.')
import numpy as np, pandas as pd, p2_study_cache as sc, p2_study_common as cm, p2_mapping as pmap
from scipy.spatial import cKDTree
hits, rays, g, cfg = sc.load('det2_long_repos')
d34=cfg.out_dir('34_clustering_matching'); A=np.load(d34+'/alignment.npz'); T=cm.make_T(g,R=A['R'],s=A['s'],pad_mean=A['pad_mean'],m3_mean=A['m3_mean'])
ra=pd.read_pickle(d34+'/rays_aligned.pkl'); xy=ra[['x','y']].to_numpy(); th=ra.theta.to_numpy()
nohit=~ra.eventId.isin(hits.eventId).to_numpy()
pil=pmap.load_pillars(cfg.MASK_GBR_PATH); px,py=T.apply(pil.x.to_numpy(),pil.y.to_numpy())
sm=np.c_[px,py][~pil.big.to_numpy()]; big=np.c_[px,py][pil.big.to_numpy()]
db=cKDTree(big).query(xy)[0]
res=[]
for lab,sel in [('vert<6',(ra.sd>=4).to_numpy()&(db>10)&(th<6)),('incl>=12',(ra.sd>=4).to_numpy()&(db>10)&(th>=12))]:
    P=xy[sel]; nh=nohit[sel]; base=nh.mean()
    grid=np.arange(-2.0,2.01,0.1); M=np.full((len(grid),len(grid)),np.nan); N=np.zeros_like(M)
    tr=cKDTree(sm)
    for i,dx in enumerate(grid):
        for j,dy in enumerate(grid):
            d,_=tr.query(P-np.array([dx,dy])); k=d<0.5
            M[i,j]=nh[k].mean(); N[i,j]=k.sum()
    z=(M-base)/np.sqrt(base*(1-base)/N)
    i,j=np.unravel_index(np.nanargmax(z),z.shape)
    print(f'{lab}: n {sel.sum()} base {100*base:.2f}%  best offset ({grid[i]:+.1f},{grid[j]:+.1f}) rate {100*M[i,j]:.1f}% (n {int(N[i,j])}) z={z[i,j]:.1f}; at (0,0): {100*M[20,20]:.1f}% z={z[20,20]:.1f}; z distribution: median {np.nanmedian(z):.2f} sd {np.nanstd(z):.2f}; frac z>3 {np.mean(z>3):.3f}')
    res.append((lab,grid,z,M))
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig,axs=plt.subplots(1,2,figsize=(12,5))
for ax,(lab,grid,z,M) in zip(axs,res):
    im=ax.imshow(z.T,origin='lower',extent=[grid[0],grid[-1],grid[0],grid[-1]],cmap='RdBu_r',vmin=-5,vmax=5); plt.colorbar(im,ax=ax,label='z of no-hit excess within 0.5 mm of shifted pillars')
    ax.set(title=f'{lab}: Gerber pillar pattern shifted by (dx,dy)',xlabel='dx [mm]',ylabel='dy [mm]')
fig.tight_layout(); fig.savefig(cfg.out_dir('39_nohit_waveforms')+'/pillar_offset_scan.png',dpi=120)
