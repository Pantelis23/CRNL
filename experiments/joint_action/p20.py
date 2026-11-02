"""Free-x2 (p2=0) ray: first crossing of Sigma and its Gamma bound V1(s1)-V1(x_dag). Usage: p20.py EPS RTOL NSIG"""
import sys
import numpy as np
exec(open('main_sigma.py').read().split("if __name__")[0])
y, h = ray(0.0)
print(f'p2=0 ray: x_dag=({y[0]:.6f},{y[1]:.6f}) S={y[4]:.10f} V1(x_dag)={V1(y[0]):.10f} Gamma_bound={V1(S1)-y[4]:.4e} H={h:.1e}')
