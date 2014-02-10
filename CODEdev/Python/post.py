#!/usr/bin/env python
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

File3 = '../bin/Pr03/RMSlog.dat'
File4 = '../bin/Pr04/RMSlog.dat'
File5 = '../bin/Pr05/RMSlog.dat'

data3 = np.loadtxt(File3)
data4 = np.loadtxt(File4)
data5 = np.loadtxt(File5)
nIter3 = data3[:,0]
RMS3 = data3[:,1]
nIter4 = data4[:,0]
RMS4 = data4[:,1]
nIter5 = data5[:,0]
RMS5 = data5[:,1]

MinX = 0
MaxX = 250
MinY = 1.0e-6
MaxY = max(RMS3)

p = plt.plot(nIter3,RMS3, 'k-', label='Grid #3')
p = plt.plot(nIter4,RMS4, 'b--', label='Grid #4')
p = plt.plot(nIter5,RMS5, 'r-', label='Grid #5')
plt.setp(p, linewidth='1.0')
plt.axis([MinX,MaxX, MinY, MaxY])
plt.xscale('linear')
plt.yscale('log')
plt.xlabel('Number of iteration', fontsize=22)
plt.ylabel('RMS residual', fontsize=22)
plt.grid(True)
ax = plt.gca()
xlabels = plt.getp(ax, 'xticklabels')
ylabels = plt.getp(ax, 'yticklabels')
plt.setp(xlabels, fontsize=18)
plt.setp(ylabels, fontsize=18)
plt.legend(
          loc='upper right',
          borderpad=0.25,
          handletextpad=0.25,
          borderaxespad=0.25,
          labelspacing=0.0,
          handlelength=2.0,
          numpoints=1)
legendText = plt.gca().get_legend().get_texts()
plt.setp(legendText, fontsize=18)
legend = plt.gca().get_legend()
legend.draw_frame(False)

pltFile = 'RMSlog.png'
fig = plt.gcf()
fig.set_size_inches(8,5)
plt.tight_layout()
plt.savefig(pltFile, format='png')
plt.close()

print "RMSlog.png DONE!!"

