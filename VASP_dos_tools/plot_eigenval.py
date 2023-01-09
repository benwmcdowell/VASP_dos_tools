import numpy as np
import matplotlib.pyplot as plt
import sys

def parse_eigenval(ifile):
    with open(ifile) as file:
        for i in range(6):
            line=file.readline().split()
        nkpts=int(line[1])
        nstates=int(line[2])
        eigenval=np.zeros((nstates*nkpts))
        kpts=np.zeros((nstates*nkpts,3))
        for j in range(nkpts):
            for i in range(2):
                line=file.readline()
            kpts[i+j*nstates]=np.array([float(k) for k in line.split[:3]])
            for i in range(nstates):
                line=file.readline().split()
                if len(line)==5:
                    #if calc is spin polarized
                   eigenval[i+j*nstates]+=(float(line[1])+float(line[2]))/2
                elif len(line)==3:
                    eigenval[i+j*nstates]+=float(line[1])
        
    return eigenval,nstates,kpts

def plot_eigenval(ifile,**args):
    if 'nbins' in args:
        nbins=args['nbins']
    else:
        nbins=1000
        
    eigenval,nstates=parse_eigenval(ifile)[:2]
        
    if 'doscar' in args:
        ef=parse_doscar(args['doscar'])
        eigenval-=ef
    
    fig,axs=plt.subplots(1,1,tight_layout=True)
    axs.hist(eigenval,bins=nbins)
    fig.show()
    
def plot_dispersion(ifile,kvec):
    eigenval,nstates,kpts=parse_eigenval(ifile)
    e=[]
    k=[]
    for i,j in zip(kpts,eigenval):
        if np.dot(i,kvec)/np.linalg.norm(i)/np.linalg.norm(kvec)==1:
            e.append(j)
            k.append(np.dot(i,kvec)/np.linalg.norm(kvec))
    plt.figure()
    plt.scatter(k,e)
    plt.show()
    
#reads DOSCAR
def parse_doscar(filepath):
    with open(filepath,'r') as file:
        line=file.readline().split()
        for i in range(5):
            line=file.readline().split()
        ef=float(line[3])
    return ef

