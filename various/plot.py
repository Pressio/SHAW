#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import sys, re
from argparse import ArgumentParser

# radius of the earth
earthRadius = 6371. #km

# thickness of the mantle
cmbDepth = 2891. # km

# the radius of the core-mantle boundary
cmbRadius = earthRadius - cmbDepth

def ak135f(dv):
  n = len(dv)
  rho, vs = np.zeros(n), np.zeros(n)
  for i in range(n):
    d = dv[i]
    dSq = d*d
    dCu = d*d*d

    if(d >= 0 and d <= 3.):
      rho[i] = 1.02
      vs[i]  = 1.

    elif(d > 3. and d <= 3.3):
      rho[i] = 2.
      vs[i]  = 1.

    elif(d > 3.3 and d <= 10.):
      rho[i] = 2.6
      vs[i]  = 3.2

    elif(d > 10. and d <= 18.):
      rho[i] = 2.92
      vs[i]  = 3.9

    elif(d > 18. and d <= 80.):
      rho[i] = 3.688908 - 0.002755944*d + 5.244987e-6*dSq
      vs[i]  = 4.479938 - 2.838134e-4*d - 3.537925e-6*dSq

    elif(d > 80. and d <= 120.):
      rho[i] = 3.6524 - 0.00188*d
      vs[i]  = 4.47 + 0.00025*d

    elif(d > 120. and d <= 210.):
      rho[i] = 3.561983 - 0.001138889*d
      vs[i]  = 4.4754   + 0.00020444*d

    elif(d > 210. and d <= 410.):
      rho[i] = 3.130252 + 0.0009128*d
      vs[i]  = 4.151532 + 0.0017548*d

    elif(d > 410. and d <= 660.):
      rho[i] = 3.948468 - 4.548571e-5*d
      vs[i]  = 4.211150 + 2.120343e-3*d

    elif(d > 660. and d <= 2740.):
      rho[i] = 3.789334 + 8.533642e-4*d - 9.455671e-8*dSq
      vs[i]  = 5.530732 + 9.439965e-4*d - 1.202153e-7*dSq

    elif(d > 2740. and d <= 2891.5):
      rho[i] = 4.268296 + 0.0005202*d
      vs[i]  = 6.648918 + 0.0002188*d

    elif(d > 2891.5 and d <= 5153.5):
      rho[i] = 3.523060 + 2.916881e-3*d - 2.421664e-7*dSq
      vs[i]  = 0.0

    elif(d > 5153.5 and d <= 6371):
      rho[i] = 4.565499 + 2.651449e-3*d - 2.080745e-7*dSq
      vs[i]  = -0.8068098 + 1.405390e-3*d - 1.103537e-7*dSq

  return [rho, vs]


def prem(rv):
  n = len(rv)
  rho, vs = np.zeros(n), np.zeros(n)
  for i in range(n):
    rKm = rv[i]
    x   = rKm/earthRadius;
    xSq = x*x;
    xCu = x*x*x;

    if(rKm >= 6356.0):
      rho[i] = 2.6
      vs[i] = 3.2

    elif(rKm >= 6346.6 and rKm < 6356.0):
      rho[i] = 2.9
      vs[i] = 3.9

    elif(rKm >= 6291.0 and rKm < 6346.6):
      rho[i] = 2.691  + 0.6924*x
      vs[i]  = 2.1519 + 2.3481*x

    elif(rKm >= 6151.0 and rKm < 6291.0):
      rho[i] = 2.691  + 0.6924*x
      vs[i]  = 2.1519 + 2.3481*x

    elif(rKm >= 5971.0 and rKm < 6151.0):
      rho[i] = 7.1089 - 3.8045*x
      vs[i]  = 8.9496 - 4.4597*x

    elif(rKm >= 5771.0 and rKm < 5971.0):
      rho[i] = 11.2494 - 8.0298*x
      vs[i]  = 22.3512 - 18.5856*x

    elif(rKm >= 5701.0 and rKm < 5771.0):
      rho[i] = 5.3197 - 1.4836*x
      vs[i]  = 9.9839 - 4.9324*x

    elif(rKm >= 5600.0 and rKm < 5701.0):
      rho[i] = 7.9565 - 6.4761*x  + 5.5283*xSq - 3.0807*xCu
      vs[i] = 22.3459 - 17.2473*x - 2.0834*xSq + 0.9783*xCu

    elif(rKm >= 3630.0 and rKm < 5600.0):
      rho[i] = 7.9565 - 6.4761*x  + 5.5283*xSq  - 3.0807*xCu
      vs[i] = 11.1671 - 13.7818*x + 17.4575*xSq - 9.2777*xCu

    elif(rKm >= 3480.0 and rKm < 3630.0):
      rho[i] = 7.9565 - 6.4761*x + 5.5283*xSq - 3.0807*xCu
      vs[i]  = 6.9254 + 1.4672*x - 2.0834*xSq + 0.9783*xCu

    elif(rKm >= 1221.5 and rKm < 3480.0):
      rho[i] = 12.5815 - 1.2638*x - 3.6426*xSq - 5.5281*xCu
      vs[i] = 0.0

    elif(rKm < 1221.5):
      rho[i] = 13.0885 - 8.8381*xSq
      vs[i]  = 3.6678  - 4.4475*xSq

  return [rho, vs]


def getPrem():
  rr = np.linspace(cmbRadius, earthRadius, 50000)
  [rho, vs] = prem(rr)
  return [rho, vs, rr]

def getAk135f():
  rr = np.linspace(cmbRadius, earthRadius, 50000)
  [rho, vs] = ak135f(rr)
  return [rho, vs, np.flipud(rr)]


###############################
if __name__== "__main__":
###############################
  parser = ArgumentParser()
  parser.add_argument("-model", "--model",
                      dest="model", default="empty",
                      help="Model: prem or ak135f. Must be set.")
  # parse args
  args = parser.parse_args()
  assert(args.model != "empty")

  if args.model == "prem":
    [rho, vs, rr] = getPrem()
  else:
    [rho, vs, rr] = getAk135f()

  fig = plt.figure(1)
  ax = fig.add_subplot(111)

  h1 = plt.plot(rho, rr, '-k', linewidth=1.5, label=r'$\rho$')
  h2 = plt.plot(vs, rr, '--k', linewidth=1.5, label=r'$v_s$')

  x = [0, 15]
  zcmb = [cmbRadius,cmbRadius]

  plt.text(7, 1600, "Core", fontsize=10)

  xfill = np.linspace(0, 15, 100)
  yfill1 = cmbRadius*np.ones(100)
  yfill2 = earthRadius*np.ones(100)
  ax.fill_between(xfill, yfill1, yfill2, facecolor='gray', alpha=0.4)
  #plt.text(7, cmbRadius-300, "Core-mantle boundary", fontsize=10)

  # plt.text(11.2, 500, "inner core", fontsize=10)
  # plt.text(11, 2200, "outer core", fontsize=10)
  if (args.model == 'prem'):
    zcore = [1221.5, 1221.5]
    plt.text(11, 4500, "Mantle", fontsize=10)
    plt.text(6.5, 5950, "Upper mantle and crust", fontsize=10)
    plt.plot(x, [5701, 5701], 'k--', lw=0.5)
  elif args.model == 'ak135f':
    zcore = [1217.5, 1217.5]
    plt.text(11, 4500, "Mantle", fontsize=10)
    plt.text(6.5, 5950, "Upper mantle and crust", fontsize=10)
    plt.plot(x, [5711, 5711], 'k--', lw=0.5)

  #plt.plot(x, zcore, 'k--', lw=1)
  plt.plot(x, zcmb, 'k--', lw=0.5)
  plt.plot(x, [earthRadius, earthRadius], 'k--', lw=0.5)

  ax.set_xlabel(r'$\rho~[g/cm^3]$, $v_s~[km/s]$', fontsize=15)
  ax.set_xlim(0,15)
  ax.set_xticks(np.linspace(0, 15, 16, endpoint=True))

  ax.set_ylabel(r'$r~[km]$', fontsize=15)
  ax.set_ylim(0,earthRadius+1000)
  ax.set_yticks(np.array([0, 500, 1000, 1500, 2000, 2500, 3000, cmbRadius,
                          4000, 4500, 5000, 5500, 6000, earthRadius]))

  ax.tick_params(axis='both', which='major', labelsize=11)
  plt.legend(loc="upper right", ncol=2)
  ax.set_aspect(0.0025)

  fig.savefig(args.model+".pdf", format="pdf", bbox_inches='tight', dpi=300)
  plt.show()
