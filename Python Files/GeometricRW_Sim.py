# -*- coding: utf-8 -*-
"""
Created on Sat Dec 14 16:46:24 2022

@author: leele2
"""
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import binom

def Simulate(prev, increment):
    return (prev*(mu*delta_t + (delta_t**0.5)*sigma*increment) + prev)

# Parameters
mu    = 0.10 #Drift (Percentage increase per unit of time)
sigma = 0.15 #Scale (Percentage increase/decrease per step)
S_0   = 10.0 #Initial point
p     = 0.50 #Probability of +1; P(Y = -1) = (1-p)
alpha = 1.00 #Step up
beta  = -1.0 #Step down

# Creating Random Walk Data
N       = 3    #Number of simulations to run
T       = 1    #Final end time for each simulation
T_N     = 200  #How many increments to split time range into
delta_t = (T)/T_N

sim = {}
for i in range(1, N + 1):
    sim[i] = []
    sim[i].append(S_0)
    for j, k in enumerate(binom.rvs(1, p, size = T_N)):
        if k == 1:
            sim[i].append(Simulate(sim[i][j], alpha))
        else:
            sim[i].append(Simulate(sim[i][j], beta))

# Plotting Simulations
plt.rc('font', family='serif')
fig, ax1 = plt.subplots()
plt.gcf().set_size_inches(6, 4.5)
for i in range(1, N + 1):
    ax1.plot([x * delta_t for x in range(0, T_N + 1)], sim[i], linewidth=1)
ax1.set_title(str(N) + ' Simulations of Geometric Brownian Motion $Y_i \sim \mathcal{B}$er$_m$($p=0.5$)')
ax1.set_xlabel('$t$')
ax1.set_ylabel('$S_n$')
plt.tight_layout()
plt.show()
# Saving plot
cwd = str(Path(__file__).parent.parent.absolute())
fig.savefig(cwd + "\Latex_Files\Main\Chapters\C1\plots\BM_Simulations.png", bbox_inches='tight')