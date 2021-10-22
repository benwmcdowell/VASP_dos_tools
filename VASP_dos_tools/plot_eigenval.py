import numpy as np
import matplotlib.pyplot as plt
import sys

def parse_eigenval(ifile):
    with open(ifile) as file:
        for i in range(6):
            line=file.readline().split()
        kpts=int(line[1])
        nstates=int(line[2])
        eigenval=np.zeros((nstates*kpts))
        for j in range(kpts):
            for i in range(2):
                file.readline()
            for i in range(nstates):
                line=file.readline().split()
                if len(line)==5:
                    #if calc is spin polarized
                   eigenval[i+j*nstates]+=(float(line[1])+float(line[2]))/2
                elif len(line)==3:
                    eigenval[i+j*nstates]+=float(line[1])
        
    return eigenval,nstates

def plot_eigenval(ifile,**args):
    if 'nbins' in args:
        nbins=args['nbins']
    else:
        nbins=1000
        
    eigenval,nstates=parse_eigenval(ifile)
        
    if 'doscar' in args:
        ef=parse_doscar(args['doscar'])
        eigenval-=ef
    
    fig,axs=plt.subplots(1,1,tight_layout=True)
    axs.hist(eigenval,bins=nbins)
    fig.show()
    
#reads DOSCAR
def parse_doscar(filepath):
    with open(filepath,'r') as file:
        line=file.readline().split()
        for i in range(5):
            line=file.readline().split()
        ef=float(line[3])
    return ef

