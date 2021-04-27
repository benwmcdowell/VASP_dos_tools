from numpy import array,dot,shape
import sys
import matplotlib.pyplot as plt
import getopt
from os.path import exists

def plot_dos(doscar,poscar,**args):
    dos, energies, ef, orbitals = parse_doscar(doscar)
    atomtypes, atomnums = parse_poscar(poscar)[2:4]
    
    if 'full' in args:
        full_dos_only=args['full']
    else:
        full_dos_only=False
    
    if 'orbitals' in args and len(args['orbitals'])>0:
        orbitals_to_plot=args['orbitals']
    else:
        orbitals_to_plot=orbitals
    
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
                    if j+1 in nums:
                        selected_atoms.append(counter)
                    counter+=1
            else:
                counter+=atomnums[atomtypes.index(i)]
    elif len(nums)>0:
        selected_atoms=nums
    elif len(types)>0:
        counter=1
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
        if full_dos_only:
            plt.plot(energies,array(dos[0][0])+array(dos[0][1]),label='total dos')
            break
        for j in range(len(atomnums)):
            if i-1<sum(atomnums[:j+1]):
                atomlabel=atomtypes[j]
                break
        if not integrate_dos:
            tempy=array([0.0 for j in range(len(energies))])
            for j in range(len(dos[i])):
                for k in orbitals_to_plot:
                    if k in orbitals[j]:
                        tempy+=dos[i][j]
                        break
            plt.plot(energies,tempy,label='{} #{}'.format(atomlabel,i-sum(atomnums[:atomtypes.index(atomlabel)])))
        else:
            tempy=array([0.0 for j in range(len(erange))])
            for j in range(len(dos[i])):
                for k in orbitals_to_plot:
                    if k in orbitals[j]:
                        tempy+=dos[i][j][emin+1:emax]
                        break
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
    
#reads DOSCAR
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
                    temp_dos[k].append(float(line[k+1]))
            dos.append(temp_dos)
    energies=array(energies)-ef
    
    #orbitals contains the type of orbital found in each array of the site projected dos
    num_columns=shape(dos[1:])[1]
    if num_columns==3:
        orbitals=['s','p','d']
    elif num_columns==6:
        orbitals=['s_up','s_down','p_up','p_down','d_up','d_down']
    elif num_columns==9:
        orbitals=['s','p_y','p_z','p_x','d_xy','d_yz','d_z2','d_xz','d_x2-y2']
    elif num_columns==18:
        orbitals=['s_up','s_down','p_y_up','p_y_down','p_z_up','p_z_down','p_x_up','p_x_down','d_xy_up','d_xy_down','d_yz_up','d_yz_down','d_z2_up','d_z2_down','d_xz_up','d_xz_down','d_x2-y2_up','d_x2-y2_down']
        
    #dos is formatted as [[total dos],[atomic_projected_dos for i in range(atomnum)]]
    #total dos has a shape of (4,nedos): [[spin up],[spin down],[integrated, spin up],[integrated spin down]]
    #atomic ldos have shapes of (6,nedos): [[i,j] for j in [spin up, spin down] for i in [s,p,d]]
    #energies has shape (1,nedos) and contains the energies that each dos should be plotted against
    return dos, energies, ef, orbitals

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
    full=False
    orbitals_to_plot=[]
    try:
        opts,args=getopt.getopt(sys.argv[1:],'a:t:i:hfo:',['atomnums=','types=','integrated=','help','full','orbitals='])
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
        if i in ['-f','--full']:
            full=True
        if i in ['-o','--orbitals']:
            orbitals_to_plot=j.split(',')
        if i in ['-h','--help']:
            print('''
plotting options:
-a, --atomnums          specify site projected DOS to plot by the index of atoms: 1,2,3,etc...
-t, --types             specify which site projected DOS to plot by atom type: Au,C,etc...
-i, --integrated        integrate the DOS between the range specified. ie -i 0,3 will plot the integrated DOS from 0 to 3 eV above the Fermi level
-f, --full              plot the total DOS instead of site projected DOS
-o, --orbitals          plot only the contributions from specific orbital projections
                        for example, to plot only s up and px down use: -o s_up,p_x_down

help options:
-h, --help               display this help message
                  ''')
            sys.exit()
    if exists(doscar):
        plot_dos(doscar,poscar,nums=atomnums,types=atomtypes,irange=irange,full=full,orbitals=orbitals_to_plot)
