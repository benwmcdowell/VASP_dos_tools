import numpy as np
import matplotlib.pyplot as plt
import sys

def parse_procar(ifile,**args):
    if 'kpoints' in args:
        kpts=args['kpts']
    else:
        kpts=0
    with open(ifile) as file:
        for i in range(2):
            line=file.readline().split()
        nkpts=int(line[3])
        if kpts==0:
            kpts=[i for i in range(nkpts)]
        nbands=int(line[7])
        nions=int(line[11])
        eigenval=np.zeros((nkpts*nbands))
        total_weight=np.zeros((nkpts*nbands,nions))
        k_weights=np.zeros((nkpts*nbands))
        for i in range(nkpts):
            for j in range(2):
                line=file.readline()
            if i in kpts:
                for j in range(nbands):
                    k_weights[i*nbands+j]=float(line.split()[-1])
            for j in range(nbands):
                for k in range(2):
                    line=file.readline()
                eigenval[i*nbands+j]=float(line.split()[4])
                for k in range(2):
                    line=file.readline()
                for k in range(nions):
                    line=file.readline().split()
                    total_weight[i*nbands+j,k]+=float(line[-1])
                line=file.readline().split()
                total_weight[i*nbands+j,:]/=float(line[-1])
            file.readline()
                
    return eigenval,total_weight,k_weights

#cutoff is the minimum total weight for which specific eigenstates are selected
def plot_atom_eigenval(atoms_to_plot,cutoff,procar,poscar,doscar):
    atomtypes,atomnums=parse_poscar(poscar)[2:4]
    atoms=[]
    for i in range(len(atomnums)):
        for j in range(atomnums[i]):
            atoms.append(atomtypes[i])
    atom_filter=np.zeros((len(atoms)))
    for i in range(len(atoms)):
        if atoms[i] in atoms_to_plot:
            atom_filter[i]+=1.0
    ef=parse_doscar(doscar)
    total_eigenval,weights,total_kweights=parse_procar(procar)
    total_eigenval-=ef
    eigenval=[]
    kweights=[]
    for i in range(len(total_eigenval)):
        if sum(atom_filter*weights[i])>cutoff:
            eigenval.append(total_eigenval[i])
            kweights.append(total_kweights[i])
            
    eigenval=np.array(eigenval)
    kweights=np.array(kweights)
    kweights/=sum(kweights)
    
    return eigenval,kweights
    
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