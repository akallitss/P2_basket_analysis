import sys; sys.path.insert(0,'.')
import numpy as np, pandas as pd, p2_study_cache as sc
from scipy.spatial import cKDTree
hits, rays, g, cfg = sc.load('det2_long_repos')
d34=cfg.out_dir('34_clustering_matching'); ra=pd.read_pickle(d34+'/rays_aligned.pkl')
F=pd.read_csv(cfg.out_dir('39_nohit_waveforms')+'/pad_waveform_features.csv')
best=F.loc[F.groupby(['kind','ray']).mf_sig.idxmax()]
empty_ev=set(best[(best.kind=='nohit')&(best.mf_sig<3.5)].eventId)
nohit=~ra.eventId.isin(hits.eventId).to_numpy()
act=(ra.sd.to_numpy()>=2)
xy=ra[['x','y']].to_numpy(); th=ra.theta.to_numpy()
rng=np.random.default_rng(1)
for lab,sel in [('all no-hit',act&nohit),('empty no-hit',act&ra.eventId.isin(empty_ev).to_numpy()),('empty, <8deg',act&ra.eventId.isin(empty_ev).to_numpy()&(th<8))]:
    P=xy[sel]; n=len(P); pool=xy[act]
    for r in (1.0,1.5,2.5):
        obs=len(cKDTree(P).query_pairs(r))
        exp=[len(cKDTree(pool[rng.choice(len(pool),n,replace=False)]).query_pairs(r)) for _ in range(200)]
        print(f'{lab:15s} n={n:5d} r<{r}: pairs {obs:4d}  random {np.mean(exp):6.1f} +- {np.std(exp):4.1f}  (p(>=obs)={np.mean(np.array(exp)>=obs):.3f})')
# local gain: hits amplitude (cluster charge) of HIT tracks within 3 mm of each no-hit track vs of random tracks
clu=pd.read_pickle(cfg.out_dir('34_study_cache')+'/clusters_cc20.pkl')
q=clu.groupby('eventId').q.max()
hit_tr=ra[act&~nohit].copy(); hit_tr['q']=hit_tr.eventId.map(q)
T=cKDTree(hit_tr[['x','y']].to_numpy())
def localq(P):
    idx=T.query_ball_point(P, r=3.0); return np.array([np.median(hit_tr.q.to_numpy()[i]) if len(i)>=3 else np.nan for i in idx])
for lab,sel in [('around empty no-hit', act&ra.eventId.isin(empty_ev).to_numpy()),('around random active', act)]:
    P=xy[sel]; P=P[rng.choice(len(P),min(len(P),1200),replace=False)]
    lq=localq(P); print(f'{lab:22s}: median of local median cluster charge {np.nanmedian(lq):.0f}  (n={np.isfinite(lq).sum()})  frac local<500: {np.nanmean(lq<500):.3f}')
