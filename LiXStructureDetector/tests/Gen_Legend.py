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

y = np.zeros(9)

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
        'size'   : 21}

matplotlib.rc('font', **font)

# Set up the matplotlib figure
for idx,y_dat in enumerate(y[:]):
    if idx > -1:
        fig = plt.scatter(0.5, 0.5, s=30, c=cols[idx], label=labels[idx])
plt.axis('off')
plt.legend(loc=10,bbox_to_anchor=(0.,0.,1.,1.),ncol=2,markerscale=2,fancybox=False,
           shadow=False,framealpha=1,edgecolor='0.',mode="expand",
           columnspacing=1.0,handletextpad=0.2)
plt.margins(x=0)

plt.savefig('Legend.svg', format='svg',bbox_inches='tight')
