import matplotlib.pyplot as plt
import numpy as np
import sys
import copy
import matplotlib
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
from scipy import interpolate
import tools

Base = '../' 
binfolderList=['/binned-new/']
#simname=str(sys.argv[0])

#list_of_sims = ['LCDM', 'EXP001']
#snapshots = [92,56,41,35]
#list_of_reds = ['_z_000','_z_055','_z_100','_z_161']
list_of_sims = ['LCDM']
snapshots = [92,56]
list_of_reds = ['z=0','z=0.55']
list_particles = ['_']
#fill the snapshot numbers with zeros, to match filenames
list_snaps = [str(snap).zfill(3) for snap in snapshots ]

pklabels=['','$P(k)$', '$\Delta(k)$']
pkindex=1
#pkindex=1 -> dimensionful P(k),  pkindex=2 -> dimensionless Delta(k)

simu_set = len(list_of_sims)


#choose which file to plot
#simulation index
simindx=0
#snapshot (redshift) index
snaindx=1

#choose data and labels
list_of_files_and_labels = (Base+list_of_sims[simindx]+binfolderList[0]+list_of_sims[simindx]+list_particles[0]+str(list_snaps[snaindx])+'.txt', 
        'CoDECS '+pklabels[pkindex]+' '+list_of_sims[simindx]+' '+list_of_reds[snaindx])

datalist = [np.loadtxt(list_of_files_and_labels[0])[:,(0,pkindex)], list_of_files_and_labels[1]]



G = gridspec.GridSpec(4,4)

fig=plt.figure(1, figsize=(20,12), dpi=80,facecolor='w')

collist=['b','g','r','c']
collist2=['Indigo','Olive','OrangeRed','SkyBlue']

axes1 = fig.add_subplot(G[:,:])

minorLocator   = ticker.LogLocator(base=10,subs=np.linspace(1.0,10,10))

axes1.loglog(datalist[0][:,0],datalist[0][:,1], '-p', label=datalist[1], color=collist[1], lw=3, ms=3, markevery=20, alpha=1.0 )

xmini,xmaxi,ymini,ymaxi = axes1.axis('tight')
axes1.axis(ymax=ymaxi*1.2, xmin=xmini*0.9, xmax=xmaxi*1.5, ymin=ymini)
axes1.legend(loc='best',markerscale=1.5,prop={'size':12},numpoints=3,handlelength=3)
axes1.grid(True,which="major",ls=":")
axes1.tick_params(which='both',length=6, width=1, labeltop=True)
axes1.set_ylabel(pklabels[pkindex],size='x-large')
axes1.xaxis.set_minor_locator(minorLocator)
axes1.text(0.02,0.1,'Non-Linear power spectrum '+list_of_sims[simindx]+', '+list_of_reds[snaindx]
        ,verticalalignment='top',horizontalalignment='left',transform=axes1.transAxes,fontsize=14,style='italic') 
axes1.set_xlabel("k in $h/Mpc$",size='x-large')
minorLocator   = ticker.LogLocator(0.1,[0.01,0.02])

#axes2 = fig.add_subplot(G[2:,:])



fig.savefig('nonlinear-Pk_'+list_of_sims[simindx]+'_'+list_of_reds[snaindx]+'.png',dpi=200)

plt.show()
