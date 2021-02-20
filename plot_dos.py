from numpy import array,dot
import sys
import matplotlib.pyplot as plt
import getopt
from os.path import exists

def plot_dos(doscar,poscar,**args):
    dos, energies, ef = parse_doscar(doscar)
    atomtypes, atomnums = parse_poscar(poscar)[2:4]
    
    if 'irange' in args and len(args['irange'])>0:
        integrate_dos=True
        irange=args['irange']
        erange=[]
        for i in range(len(energies)):
            if energies[i]>min(irange) and energies[i]<max(irange):
                erange.append(energies[i])
            elif energies[i]<min(irange):
                emin=i
            elif energies[i]>max(irange):
                emax=i
                break
    else:
        integrate_dos=False
    
    if 'nums' in args:
        nums=args['nums']
    else:
        nums=[]
    if 'types' in args:
        types=args['types']
    else:
        types=[]
    selected_atoms=[]
    if len(types)>0 and len(nums)>0:
        counter=1
        for i in atomtypes:
            if i in types:
                for j in range(atomnums[atomtypes.index(i)]):
                    if counter in nums:
                        selected_atoms.append(counter)
                    counter+=1
            else:
                counter+=atomnums[atomtypes.index(i)]
    elif len(nums)>0:
        selected_atoms=nums
    elif len(types)>0:
        counter=0
        for i in atomtypes:
            if i in types:
                for j in range(atomnums[atomtypes.index(i)]):
                    selected_atoms.append(counter)
                    counter+=1
            else:
                counter+=atomnums[atomtypes.index(i)]
    else:
        selected_atoms=[i for i in range(sum(atomnums))]
    
    plt.figure()
    for i in selected_atoms:
        for j in range(len(atomnums)):
            if i<sum(atomnums[:j+1]):
                atomlabel=atomtypes[j]
                break
        if not integrate_dos:
            tempy=array([0.0 for j in range(len(energies))])
            for j in range(len(dos[i+1])):
                tempy+=dos[i+1][j]
            plt.plot(energies,tempy,label='{} #{}'.format(atomlabel,i-sum(atomnums[:atomtypes.index(atomlabel)])))
        else:
            tempy=array([0.0 for j in range(len(erange))])
            for j in range(len(dos[i+1])):
                tempy+=dos[i+1][j][emin+1:emax]
            for j in range(1,len(tempy)):
                tempy[j]+=tempy[j-1]
            plt.plot(erange,tempy,label='{} #{}'.format(atomlabel,i-sum(atomnums[:atomtypes.index(atomlabel)])))
    plt.xlabel('energy - $E_f$ / eV')
    if not integrate_dos:
        plt.ylabel('DOS / states $eV^{-1}$')
    else:
        plt.ylabel('integrated DOS / # states')
    plt.legend()
    plt.show()
    
def parse_doscar(filepath):
    with open(filepath,'r') as file:
        line=file.readline().split()
        atomnum=int(line[0])
        for i in range(5):
            line=file.readline().split()
        nedos=int(line[2])
        ef=float(line[3])
        dos=[]
        energies=[]
        for i in range(atomnum+1):
            if i!=0:
                line=file.readline()
            for j in range(nedos):
                line=file.readline().split()
                if i==0:
                    energies.append(float(line[0]))
                if j==0:
                    temp_dos=[[] for k in range(len(line)-1)]
                for k in range(len(line)-1):
                    try:
                        temp_dos[k].append(float(line[k+1]))
                    except:
                        temp_dos[k].append(float('e-'.join(line[k+1].split('-'))))
            dos.append(temp_dos)
    energies=array(energies)-ef
        
    #dos is formatted as [[total dos],[atomic_projected_dos for i in range(atomnum)]]
    #total dos has a shape of (4,nedos): [[spin up],[spin down],[integrated, spin up],[integrated spin down]]
    #atomic ldos have shapes of (6,nedos): [[i,j] for j in [spin up, spin down] for i in [s,p,d]]
    #energies has shape (1,nedos) and contains the energies that each dos should be plotted against
    return dos, energies, ef

def parse_poscar(ifile):
    with open(ifile, 'r') as file:
        lines=file.readlines()
        sf=float(lines[1])
        latticevectors=[float(lines[i].split()[j])*sf for i in range(2,5) for j in range(3)]
        latticevectors=array(latticevectors).reshape(3,3)
        atomtypes=lines[5].split()
        atomnums=[int(i) for i in lines[6].split()]
        if 'Direct' in lines[7] or 'Cartesian' in lines[7]:
            start=8
            mode=lines[7].split()[0]
        else:
            mode=lines[8].split()[0]
            start=9
            seldyn=[''.join(lines[i].split()[-3:]) for i in range(start,sum(atomnums)+start)]
        coord=array([[float(lines[i].split()[j]) for j in range(3)] for i in range(start,sum(atomnums)+start)])
        if mode!='Cartesian':
            for i in range(sum(atomnums)):
                for j in range(3):
                    while coord[i][j]>1.0 or coord[i][j]<0.0:
                        if coord[i][j]>1.0:
                            coord[i][j]-=1.0
                        elif coord[i][j]<0.0:
                            coord[i][j]+=1.0
                coord[i]=dot(coord[i],latticevectors)
            
    #latticevectors formatted as a 3x3 array
    #coord holds the atomic coordinates with shape ()
    try:
        return latticevectors, coord, atomtypes, atomnums, seldyn
    except NameError:
        return latticevectors, coord, atomtypes, atomnums
    
if __name__=='__main__':
    irange=[]
    doscar='./DOSCAR'
    if exists('./CONTCAR'):
        poscar='./CONTCAR'
    else:
        poscar='./POSCAR'
    atomnums=[]
    atomtypes=[]
    try:
        opts,args=getopt.getopt(sys.argv[1:],'a:t:i:h',['atomnums=','types=','integrated=','help'])
    except getopt.GetoptError:
        print('error in command line syntax')
        sys.exit(2)
    for i,j in opts:
        if i in ['-a','--atomnums']:
            atomnums=[int(k) for k in j.split(',')]
        if i in ['-t','--types']:
            atomtypes=[str(k) for k in j.split(',')]
        if i in ['-i','--integrated']:
            irange=[float(k) for k in j.split(',')]
        if i in ['-h','--help']:
            print('''
plotting options:
-a, --atomnums          specify site projected DOS to plot by the index of atoms: 1,2,3,etc...
-t, --types             specify which site projected DOS to plot by atom type: Au,C,etc...
-i, --integrated        integrate the DOS between the range specified. ie -i 0,3 will plot the integrated DOS from 0 to 3 eV above the Fermi level

help options:
-h, --help               display this help message
                  ''')
            sys.exit()
    plot_dos(doscar,poscar,nums=atomnums,types=atomtypes,irange=irange)
