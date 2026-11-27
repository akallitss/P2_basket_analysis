import sys; sys.path.insert(0,'.')
import numpy as np, pandas as pd, p2_study_cache as sc, p2_study_common as cm, p2_mapping as pmap
from scipy.spatial import cKDTree
hits, rays, g, cfg = sc.load('det2_long_repos')
d34=cfg.out_dir('34_clustering_matching'); A=np.load(d34+'/alignment.npz')
T=cm.make_T(g,R=A['R'],s=A['s'],pad_mean=A['pad_mean'],m3_mean=A['m3_mean'])
ra=pd.read_pickle(d34+'/rays_aligned.pkl'); xy=ra[['x','y']].to_numpy()
nohit=~ra.eventId.isin(hits.eventId).to_numpy(); th=ra.theta.to_numpy()
B='/local/home/ak271430/Documents/PostDocSaclay/Detector_Drawings/Version_Apr26/'
masks={'Saclay V1':B+'Insulation_masks/V1/P2_BASKET-Mask_M2_V1.gbr','Saclay V2':B+'Insulation_masks/V2/P2_BASKET-Mask_M2_V2.gbr','CERN':B+'Bulk_CERN/P2_Mask2.gbr'}
for nm,p in masks.items():
    pil=pmap.load_pillars(p); 
    if not len(pil): print(nm,'none'); continue
    px,py=T.apply(pil.x.to_numpy(),pil.y.to_numpy()); sm=np.c_[px,py][~pil.big.to_numpy()]; big=np.c_[px,py][pil.big.to_numpy()]
    dp,_=cKDTree(sm).query(xy); db,_=cKDTree(big).query(xy)
    sel=(ra.sd.to_numpy()>=4)&(db>10)
    pitch=np.median(cKDTree(sm).query(sm,k=2)[0][:,1])
    print(f'{nm}: {len(sm)} small pillars r={pil[~pil.big].r.median():.3f} pitch {pitch:.2f}; big at same place as V1? ')
    for lab,s in [('<6deg',sel&(th<6)),('6-12',sel&(th>=6)&(th<12)),('>=12',sel&(th>=12))]:
        e=[0,0.25,0.5,0.75,1.0,1.25,1.5,2.0,3.0]
        out=[]
        for lo,hi in zip(e[:-1],e[1:]):
            k=s&(dp>=lo)&(dp<hi)
            out.append(f'{lo:.2f}-{hi:.2f}: {100*nohit[k].mean():4.1f}% ({k.sum()})' if k.sum()>20 else f'{lo:.2f}-{hi:.2f}:   -  ')
        print('   ',lab,' | '.join(out))
