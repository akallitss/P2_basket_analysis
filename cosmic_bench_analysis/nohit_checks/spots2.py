import sys; sys.path.insert(0,'.')
import numpy as np, pandas as pd, p2_study_cache as sc, p2_study_common as cm
from scipy.spatial import cKDTree
from scipy.stats import binom
hits, rays, g, cfg = sc.load('det2_long_repos')
d34=cfg.out_dir('34_clustering_matching'); ra=pd.read_pickle(d34+'/rays_aligned.pkl')
A=np.load(d34+'/alignment.npz'); T=cm.make_T(g,R=A['R'],s=A['s'],pad_mean=A['pad_mean'],m3_mean=A['m3_mean'])
act=ra.sd.to_numpy()>=2; r=ra[act].reset_index(drop=True)
r['nohit']=~r.eventId.isin(hits.eventId)
xy=r[['x','y']].to_numpy(); tr=cKDTree(xy); p0=r.nohit.mean()
for R in (1.0,1.5,2.0):
    nb=tr.query_ball_point(xy,R)
    n=np.array([len(i)-1 for i in nb]); k=np.array([r.nohit.to_numpy()[i].sum() for i in nb])-r.nohit.to_numpy()
    r[f'n{R}']=n; r[f'k{R}']=k
    r[f'p{R}']=binom.sf(k-1,n,p0)   # P(>=k | n, p0)
R=1.5
spot=(r[f'p{R}']<1e-3)
print(f'global no-hit {100*p0:.2f}%  active {len(r)}')
print(f'tracks in spots (neighbour no-hit rate p<1e-3 within {R} mm): {spot.sum()} ({100*spot.mean():.2f}% of tracks), their no-hit rate {100*r.nohit[spot].mean():.1f}%')
print(f'  no-hits inside spots: {r.nohit[spot].sum()} of {r.nohit.sum()} ({100*r.nohit[spot].sum()/r.nohit.sum():.0f}%)')
print(f'  efficiency outside spots: no-hit {100*r.nohit[~spot].mean():.2f}%')
# group spot tracks into spots
S=r[spot]; G=cKDTree(S[['x','y']].to_numpy()).query_pairs(3.0,output_type='ndarray')
from scipy.sparse import coo_matrix; from scipy.sparse.csgraph import connected_components
n=len(S); M=coo_matrix((np.ones(len(G)),(G[:,0],G[:,1])),shape=(n,n)); nc,lab=connected_components(M,directed=False)
S=S.assign(lab=lab)
sp=S.groupby('lab').agg(x=('x','mean'),y=('y','mean'),n=('nohit','size'),nohit=('nohit','sum'),sx=('x','std'),sy=('y','std')).sort_values('nohit',ascending=False)
# pad under each spot
pads=cm.apply_T(T,g['pad_xy']); j=cKDTree(pads).query(sp[['x','y']].to_numpy())[1]
sp['pad']=g['pad_id'][j]; sp['d_padc']=np.hypot(*(pads[j]-sp[['x','y']].to_numpy()).T)
print(f'{nc} spots; top 25:'); print(sp.head(25).round(2).to_string())
# angle dependence outside spots
th=r.theta.to_numpy()
for lo,hi in [(0,4),(4,8),(8,12),(12,16),(16,90)]:
    k=(th>=lo)&(th<hi)
    print(f'theta {lo}-{hi}: no-hit all {100*r.nohit[k].mean():.2f}%  outside spots {100*r.nohit[k&~spot].mean():.2f}%  inside {100*r.nohit[k&spot].mean():.1f}% (n {int((k&spot).sum())})')
r.to_pickle(cfg.out_dir('39_nohit_waveforms')+'/tracks_spots.pkl'); sp.to_csv(cfg.out_dir('39_nohit_waveforms')+'/spots.csv')
