"""Bisect the ray family (phi = pi + d) for the boundary between rays that hit Sigma and rays
that turn back; the limiting ray is the heteroclinic to the joint saddle (3, 0.606), whose
action is W* = inf over Sigma of W_c.  Saves the path and prints its profile."""
import sys
EPSa,RTOLa,NSIGa=sys.argv[1:4]; VERB=sys.argv[4:]
sys.argv=["x",EPSa,RTOLa,NSIGa]+VERB
import numpy as np
exec(open('main_sigma.py').read().split("if __name__")[0])
ds=-np.geomspace(1e-10,1e-1,181); hits=[ray(d)[0] is not None for d in ds]
k=max(i for i in range(180) if hits[i] and not hits[i+1]); a,b=ds[k],ds[k+1]
best=None; HITS=[]
for it in range(70):
    m=0.5*(a+b); y,h=ray(m)
    if y is not None:
        a=m; best=(m,y,h); HITS.append((y[0],np.hypot(y[2],y[3]),y[4]))
        if VERB: print(f"   it{it} x1c={y[0]:.5f} x2c={y[1]:.6f} |p|={np.hypot(y[2],y[3]):.2e} S={y[4]:.10f}")
    else: b=m
    if b-a<1e-20*a+1e-300 and False: break
m,y,h=best
print(f'EPS={EPSa} RTOL={RTOLa} NSIG={NSIGa} d*={m:.15e} end x=({y[0]:.6f},{y[1]:.6f}) |p|={np.hypot(y[2],y[3]):.1e} W*={y[4]:.10f} H={h:.1e} Gamma={V1(S1)-y[4]:.6e}')
# path profile
p=np.array([-1.0,m]); dx=sh.Cm@p; dx*=np.sqrt(2*EPS/(dx@sh.Z@dx))
x0=sh.xs+dx; p0=sh.Z@dx
s=solve_ivp(sh.rhs,(0,3000),[x0[0],x0[1],p0[0],p0[1],EPS],method='DOP853',rtol=RTOL,atol=RTOL*1e-3,events=[ev_sigma,sh.ev_s1],dense_output=True)
t=np.linspace(0,s.t[-1],20001); Y=s.sol(t)
np.save(f'optpath{TAG}_eps{EPSa}_rtol{RTOLa}.npy',np.vstack([t,Y]))
k=np.argmin(Y[0]); print(f'  min x1 on path = {Y[0,k]:.6f} at x2={Y[1,k]:.4f}, S so far={Y[4,k]:.6f}')
for x2t in [2.5,2.0,1.4037,1.2,1.0,0.8,0.7,0.65]:
    j=np.argmax(Y[1]<x2t); print(f'  x2 crosses {x2t}: x1={Y[0,j]:.5f} S={Y[4,j]:.6f} p2={Y[3,j]:+.4f} t={t[j]:.3f}')
# where x2 crosses instantaneous mid root / ghost
g=np.array([Y[1,j]-mid_or_ghost(Y[0,j]) for j in range(len(t))])
j=np.argmax(g<0); print(f'  x2 crosses mid/ghost curve at x1={Y[0,j]:.5f} x2={Y[1,j]:.5f} S={Y[4,j]:.6f}')
# quadratic extrapolation of hit actions to |p| -> 0 (hits near the saddle only)
HALL=np.array(HITS)
print(f'  min S over all Sigma hits (upper bound on W*) = {HALL[:,2].min():.10f} at x1={HALL[np.argmin(HALL[:,2]),0]:.4f}; Gamma_lowerbound={V1(S1)-HALL[:,2].min():.6e}')
HA=np.array([h for h in HITS if h[0]>2.9]).reshape(-1,3)
if len(HA) < 4:
    print('  extrapolated: INSUFFICIENT near-saddle hits'); sys.exit()
c=np.polyfit(HA[:,1]**2,HA[:,2],1)
print(f'  extrapolated W*(|p|->0) = {c[1]:.10f} (slope {c[0]:.2f}, {len(HA)} hits, min|p|={HA[:,1].min():.1e}); Gamma_extrap={V1(S1)-c[1]:.6e}')
