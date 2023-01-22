# -*- coding: utf-8 -*-
"""
Created on Sat Nov 19 17:11:24 2022

@author: leele2
"""
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import binom

def Simulate(prev, increment):
    return (prev + mu + sigma * increment)

# Parameters
mu    = +0.0 #Drift
sigma = +1.0 #Scale
S_0   = +10. #Initial point
p     = +0.5 #Probability of +1; P(Y = -1) = (1-p)
alpha = +1.0 #Step up
beta  = -1.0 #Step down

# Creating Random Walk Data
N   = 3   #Number of simulations to run
T   = 100 #How many forward steps to add to each run

sim = {}
for i in range(1, N + 1):
    sim[i] = []
    sim[i].append(S_0)
    for j, k in enumerate(binom.rvs(1, p, size = T)):
        if k == 1:
            sim[i].append(Simulate(sim[i][j], alpha))
        else:
            sim[i].append(Simulate(sim[i][j], beta))

# Plotting Simulations
plt.rc('font', family='serif')
fig, ax1 = plt.subplots()
plt.gcf().set_size_inches(6, 4.5)
for i in range(1, N + 1):
    ax1.plot(sim[i], linewidth=1.5)
ax1.set_title(str(N) + ' Simulations of Symetric Random Walk')
ax1.set_xlabel('$n$')
ax1.set_ylabel('$S_n$')
plt.tight_layout()
plt.show()
# Saving plot
cwd = str(Path(__file__).parent.parent.absolute())
fig.savefig(cwd + "\Latex_Files\Main\Chapters\C1\plots\RW_Simulations.png", bbox_inches='tight')