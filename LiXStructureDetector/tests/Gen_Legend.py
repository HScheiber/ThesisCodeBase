# -*- coding: utf-8 -*-
"""
Created on Fri May 13 19:17:32 2022

@author: Hayden
"""

#%% Generate a series of line plots: partitions vs time
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
labels = ["Liquid","Rocksalt","Wurtzite","5-5","NiAs","Sphalerite","$\\beta$-BeO","AntiNiAs","CsCl"]
n_classes = len(labels)
# Generate the histogram along one dimension

y = np.zeros(n_classes)

# colors
cols = sns.color_palette("muted",n_classes)

# Set up the matplotlib figure
plt.rcParams['text.latex.preamble'] = [
       r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
       r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
       r'\usepackage{amsfonts}',    # set the normal font here
       r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
       r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
]

font = {'family' : 'sans-serif',
        'size'   : 18}

matplotlib.rc('font', **font)

# Set up the matplotlib figure
for idx,y_dat in enumerate(y[:]):
    if idx > -1:
        fig = plt.scatter(0, 0, s=30, c=cols[idx], label=labels[idx])
        #plt.plot(0.5, 0.5, label=labels[idx], color=cols[idx], alpha=0.8, linewidth=3, linestyle='solid')
        
plt.axis('off')
plt.legend(loc=10,bbox_to_anchor=(0.,0.,0.5,1.),ncol=1,markerscale=2,fancybox=False,
            shadow=False,framealpha=1,edgecolor='0.75',mode="expand",
            columnspacing=1.0,handletextpad=0)
# plt.legend(ncol=5,markerscale=1,fancybox=True,
#            shadow=False,framealpha=1,edgecolor='0.75',mode="expand")
plt.margins(x=0)

plt.savefig('Legend.svg', format='svg',bbox_inches='tight')
