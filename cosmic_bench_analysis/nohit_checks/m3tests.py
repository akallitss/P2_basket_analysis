import sys, glob; sys.path.insert(0,'.')
import numpy as np, pandas as pd, uproot, awkward as ak
f39='/local/home/ak271430/Documents/PostDocSaclay/data/Cosmic_Bench/Analysis/det2/p2_det2_long_run_430_repositioned_11-25-26/long_run_det2_430_730/39_nohit_waveforms/pad_waveform_features.csv'
import p2_study_cache as sc, p2_study_common as cm
hits, rays, g, cfg = sc.load('det2_long_repos')
R = pd.read_csv(f39)
best = R.loc[R.groupby(['kind','ray']).mf_sig.idxmax()]
cls = best[best.kind!='noise'][['kind','eventId','mf_sig']].copy()
cls['cls'] = np.where(cls.kind=='hit','hit', np.where(cls.mf_sig>3.5,'nohit_small','nohit_empty'))
fs=sorted(glob.glob(cfg.m3_tracking_dir+'/*.root'))
br=['evn','evttime','Chi2X','Chi2Y','NClusX','NClusY']+[f'm{a}{i}' for a in 'XY' for i in range(4)]+[f'z{a}{i}' for a in 'XY' for i in range(4)]
parts=[]
for f in fs:
    a=uproot.open(f+':T').arrays(br, library='ak')
    one = ak.num(a['mX0'])==1
    d={'eventId':ak.to_numpy(a['evn'][one]).astype(np.int64),'t':ak.to_numpy(a['evttime'][one])*1e-8}
    for k in br[2:]: d[k]=ak.to_numpy(ak.firsts(a[k][one]))
    parts.append(pd.DataFrame(d))
M=pd.concat(parts).drop_duplicates('eventId')
# time since previous trigger (all M3 events, incl. those without tracks)
allt=[]
for f in fs:
    a=uproot.open(f+':T').arrays(['evn','evttime'],library='np'); allt.append(pd.DataFrame({'eventId':a['evn'].astype(np.int64),'t':a['evttime']*1e-8}))
A=pd.concat(allt).drop_duplicates('eventId').sort_values('t'); A['dt_prev']=A.t.diff(); A['dt_next']=-A.t.diff(-1)
M=M.merge(A[['eventId','dt_prev','dt_next']],on='eventId',how='left')
M=M.merge(cls[['eventId','cls']],on='eventId')
for a in 'XY':
    M[f'sl_top_{a}']=(M[f'm{a}0']-M[f'm{a}1'])/(M[f'z{a}0']-M[f'z{a}1'])
    M[f'sl_bot_{a}']=(M[f'm{a}2']-M[f'm{a}3'])/(M[f'z{a}2']-M[f'z{a}3'])
    M[f'sl_all_{a}']=(M[f'm{a}0']-M[f'm{a}3'])/(M[f'z{a}0']-M[f'z{a}3'])
    M[f'kink_{a}']=np.degrees(np.arctan(M[f'sl_top_{a}'])-np.arctan(M[f'sl_bot_{a}']))
M['theta']=np.degrees(np.arctan(np.hypot(M.sl_all_X,M.sl_all_Y)))
M['kink']=np.hypot(M.kink_X,M.kink_Y)
pd.set_option('display.width',200)
def summ(d):
    return pd.Series({'n':len(d),'kink_med':d.kink.median(),'kink>1deg%':100*(d.kink>1).mean(),'kink>3deg%':100*(d.kink>3).mean(),
                      'chi2X_med':d.Chi2X.median(),'chi2Y_med':d.Chi2Y.median(),'nclusX':d.NClusX.mean(),'dt_prev_med_ms':1e3*d.dt_prev.median(),
                      'dt_prev<1ms%':100*(d.dt_prev<1e-3).mean(),'dt_next<1ms%':100*(d.dt_next<1e-3).mean()})
M['thb']=pd.cut(M.theta,[0,6,12,90])
print(M.groupby(['thb','cls'],observed=True).apply(summ).round(3).to_string())
M.to_pickle(cfg.out_dir('39_nohit_waveforms')+'/m3_track_details.pkl')
