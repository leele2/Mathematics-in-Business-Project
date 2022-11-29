# -*- coding: utf-8 -*-
"""
Created on Sat Nov 19 17:11:24 2022

@author: leele2
"""
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import binom

# Parameters
mu    = 0.0 #Drift
sigma = 0.0 #Scale
S_0   = 0.0 #Initial point
p     = 0.5 #Probability of +1; P(Y = -1) = (1-p)

# Creating Random Walk Data
N   = 5   #Number of simulations to run
T   = 100 #How many forward steps to add to each run

sim = {}
for i in range(1, N + 1):
    sim[i] = []
    sim[i].append(S_0)
    for j, k in enumerate(binom.rvs(1, p, size = T)):
        if k == 1:
            sim[i].append(sim[i][j] + 1)
        else:
            sim[i].append(sim[i][j] - 1)

# Plotting Simulations
plt.rc('font', family='serif')
fig, ax1 = plt.subplots()
plt.gcf().set_size_inches(6, 4.5)
for i in range(1, N + 1):
    ax1.plot(sim[i], linewidth=2)
ax1.set_title('Five Symetric Random Walk Simulations')
ax1.set_xlabel('$n$')
ax1.set_ylabel('$S_n$')
plt.tight_layout()
plt.show()
# Saving plot
cwd = str(Path(__file__).parent.parent.absolute())
fig.savefig(cwd + "\Latex_Files\Main\Chapters\C1\plots\RW_Simulations.png", bbox_inches='tight')