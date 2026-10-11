import time, numpy as np, xc
import task2b_ablation as B
nc = xc.NetCfg('ladder_6', 0, 'thin')
ref = B.setup_net(999,'thin','ladder_6',0)
print('corr equal', np.array_equal(nc.corr, ref[-1]))
ops = xc.Ops(nc.STA, nc.corr, xc.make_grid()); print('build', ops.t_build)
PW = xc.PKDictW(ops.BD, 4)
st = xc.window_stream(nc, 3001)
for w_i in range(4):
    w = next(st)
    t=time.time(); r1,_ = xc.bomp_pk(ops.PD, w['MFn'], w['truth']); t1=time.time()-t
    r2,_ = xc.bomp_pk(PW, w['MFn'], w['truth'])
    t=time.time(); r3,l3,_ = xc.tssfw(ops, w, w['truth'], f=0.01); t3=time.time()-t
    t=time.time(); r4 = xc.mtfw(ops, w['MFn'], w['truth']); t4=time.time()-t
    t=time.time(); o = xc.oracle_scores(ops, w); t5=time.time()-t
    print(len(w['truth']), r1==r2, 'bomp',np.round(r1,1).tolist(), f'{t1:.2f}s', 'ts',np.round(r3,1).tolist(), f'{t3:.2f}s', f'mtfw {t4:.2f}s', 'oracle', np.round(o,1), f'{t5:.2f}s')
