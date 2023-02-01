import numpy as np
import matplotlib.pyplot as plt
import os

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
            current_k=np.array([float(k) for k in line.split()[:3]])
            for i in range(nstates):
                kpts[i+j*nstates]=current_k
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
    
def plot_dispersion(fpath,kvec):
    os.chdir(fpath)
    eigenval,nstates,kpts=parse_eigenval('./EIGENVAL')
    lv=parse_poscar('./POSCAR')[0]
    ef=parse_doscar('./DOSCAR')
    
    kvec=np.dot(kvec,lv)
    eigenval-=ef
    
    e=np.zeros((nstates,int(len(kpts)/nstates)))
    k=np.zeros((np.shape(e)[1],3))
    for i in range(nstates):
        for j in range(np.shape(e)[1]):
            if i==0:
                k[j]=np.dot(kpts[j*nstates+i],kvec)/np.linalg.norm(kvec)
            e[i,j]=eigenval[j*nstates+i]
            
    plt.figure()
    for i in range(nstates):
        plt.plot(k,e[i])
    plt.show()
    
    return k,e
    
#reads DOSCAR
def parse_doscar(filepath):
    with open(filepath,'r') as file:
        line=file.readline().split()
        for i in range(5):
            line=file.readline().split()
        ef=float(line[3])
    return ef

#reads POSCAR
def parse_poscar(ifile):
    with open(ifile, 'r') as file:
        lines=file.readlines()
        sf=float(lines[1])
        latticevectors=[float(lines[i].split()[j])*sf for i in range(2,5) for j in range(3)]
        latticevectors=np.array(latticevectors).reshape(3,3)
        atomtypes=lines[5].split()
        atomnums=[int(i) for i in lines[6].split()]
        if 'Direct' in lines[7] or 'Cartesian' in lines[7]:
            start=8
            mode=lines[7].split()[0]
        else:
            mode=lines[8].split()[0]
            start=9
            seldyn=[''.join(lines[i].split()[-3:]) for i in range(start,sum(atomnums)+start)]
        coord=np.array([[float(lines[i].split()[j]) for j in range(3)] for i in range(start,sum(atomnums)+start)])
        if mode!='Cartesian':
            for i in range(sum(atomnums)):
                for j in range(3):
                    while coord[i][j]>1.0 or coord[i][j]<0.0:
                        if coord[i][j]>1.0:
                            coord[i][j]-=1.0
                        elif coord[i][j]<0.0:
                            coord[i][j]+=1.0
                coord[i]=np.dot(coord[i],latticevectors)
            
    #latticevectors formatted as a 3x3 array
    #coord holds the atomic coordinates with shape ()
    try:
        return latticevectors, coord, atomtypes, atomnums, seldyn
    except NameError:
        return latticevectors, coord, atomtypes, atomnums