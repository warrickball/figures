#!/usr/bin/env python3

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-n', type=float, nargs='+',
                    default=[1,2,3,4,5],
                    help="polytropic indices for which to solve "
                    "(default=1--5)")
parser.add_argument('--kind', default='rho',
                    choices=['rho', 'theta', 'phi', 'uv'])
parser.add_argument('--ximax', type=float, default=10,
                    help="maximum value of xi for integrations, "
                    "which matters for n >= 5 (default=10)")
parser.add_argument('--max-step', type=float, default=1e-2,
                    help="maximum step for ODE solver (default=1e-2)")
parser.add_argument('--vectorize', action='store_true',
                    help="vectorize the function passed to the ODE solver")
parser.add_argument('--method', type=str, default='RK45',
                    help="solver choice for ODE solver (default='RK45')")
args = parser.parse_args()

import numpy as np
import matplotlib.pyplot as pl
from scipy.integrate import solve_ivp

if args.vectorize:
    def lee(t, y, n):
        return np.vstack([-y[1]/t**2, y[0]**n*t**2])
else:
    def lee(t, y, n):
        return (-y[1]/t**2, y[0]**n*t**2)

def surface(t, y, n):
    return y[0]

surface.terminal = True

xi0 = 1e-6
theta0 = 1-xi0**2/6

for n in args.n:
    sol = solve_ivp(lee, (xi0, args.ximax), (theta0, 0), vectorized=args.vectorize, args=(n,),
                    max_step=args.max_step, events=surface, method=args.method)
    x = sol.t
    
    if args.kind == 'rho':
        y = sol.y[0]**n
    elif args.kind == 'theta':
        y = sol.y[0]
    elif args.kind == 'phi':
        y = sol.y[1]
    elif args.kind == 'uv':
        x = sol.t**3*sol.y[0]**n/sol.y[1]
        y = (n+1)*sol.y[1]/sol.t/sol.y[0]
        
    pl.plot(x, y, label=n)

pl.xlabel(r'$\xi$')

if args.kind == 'rho':
    pl.ylabel(r'$\rho/\rho_c$')
elif args.kind == 'theta':
    pl.ylabel(r'$\theta$')
elif args.kind == 'phi':
    pl.ylabel(r'$\phi$')
elif args.kind == 'uv':
    pl.xlabel('U')
    pl.ylabel('V')
    pl.axis([0,3,0,6])
    
pl.legend()
pl.show()
