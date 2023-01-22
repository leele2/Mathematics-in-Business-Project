# -*- coding: utf-8 -*-
"""
Created on Sat Dec 14 16:46:24 2022

@author: leele2
"""
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import norm

def Simulate(prev, increment):
    return (prev*(mu*delta_t + (delta_t**0.5)*sigma*increment) + prev)

# Parameters
mu    = 0.10 #Drift (Percentage increase per unit of time)
sigma = 0.15 #Scale (Percentage increase/decrease per step)
S_0   = 10.0 #Initial point

# Creating Random Walk Data
N       = 3    #Number of simulations to run
T       = 1    #Final end time for each simulation
T_N     = 200  #How many increments to split time range into
delta_t = (T)/T_N

sim = {}
for i in range(1, N + 1):
    sim[i] = []
    sim[i].append(S_0)
    for j, k in enumerate(norm.rvs(0, 1, size = T_N)):
        sim[i].append(Simulate(sim[i][j], k))

# Plotting Simulations
plt.rc('font', family='serif')
fig, ax1 = plt.subplots()
plt.gcf().set_size_inches(6, 4.5)
for i in range(1, N + 1):
    ax1.plot([x * delta_t for x in range(0, T_N + 1)], sim[i], linewidth=1)
ax1.set_title(str(N) + ' Simulations of Geometric Brownian Motion $Y_i \sim \mathcal{N}(0,1)$')
ax1.set_xlabel('$t$')
ax1.set_ylabel('$S_n$')
plt.tight_layout()
plt.show()
# Saving plot
cwd = str(Path(__file__).parent.parent.absolute())
fig.savefig(cwd + "\Latex_Files\Main\Chapters\C1\plots\BM_Norm_Simulations.png", bbox_inches='tight')