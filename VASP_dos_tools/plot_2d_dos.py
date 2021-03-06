from numpy import array,dot,shape,zeros,outer
import sys
import matplotlib.pyplot as plt
import getopt
from os.path import exists

def plot_2d_dos(doscar1,doscar2,poscar1,poscar2,**args):
    dos=[]
    energies=[]
    ef=[]
    orbitals=[]
    atomtypes=[]
    atomnums=[]
    for i,j in zip([doscar1,doscar2],[poscar1,poscar2]):
        tempvar=parse_doscar(i)
        dos.append(tempvar[0])
        energies.append(tempvar[1])
        ef.append(tempvar[2])
        orbitals.append(tempvar[3])
        tempvar=parse_poscar(j)[2:4]
        atomtypes.append(tempvar[0])
        atomnums.append(tempvar[1])
    
    if 'orbitals' in args and len(args['orbitals'])>0:
        orbitals_to_plot=args['orbitals']
    else:
        orbitals_to_plot=orbitals[0]
    contributing_orbitals=[]
    for k in orbitals[0]:
        for l in orbitals_to_plot:
            if l in k:
                contributing_orbitals.append(k)
                break
    if len(contributing_orbitals)==len(orbitals[0]):
        contributing_orbitals=['all']
    
    if 'nums' in args:
        nums=args['nums']
    else:
        nums=[]
        
    if 'types' in args:
        types=args['types']
    else:
        types=[]
        
    if 'average' in args:
        average_dos=True
    else:
        average_dos=False
        
    if 'cmap' in args:
        cmap=args['cmap']
    else:
        cmap='jet'
        
    start=[0,0]
    end=[len(i) for i in energies]
    if 'energy_range' in args:
        for i in range(2):
            if min(energies[i])>args['energy_range'][0] or max(energies[i])<args['energy_range'][1]:
                print('min or max of energy range exceeds bounds of doscar{}'.format(i))
                sys.exit()
            for j in range(len(energies[i])):
                if energies[i][j]<args['energy_range'][0]:
                    start[i]=j
                if energies[i][j]>args['energy_range'][1]:
                    end[i]=j
                    break
        
    selected_atoms=[[],[]]
    for k in range(2):
        if len(types)>0 and len(nums)>0:
            counter=1
            for i in atomtypes[k]:
                if i in types:
                    for j in range(atomnums[k][atomtypes[k].index(i)]):
                        if j+1 in nums:
                            selected_atoms[k].append(counter)
                        counter+=1
                else:
                    counter+=atomnums[k][atomtypes[k].index(i)]
        elif len(nums)>0:
            selected_atoms[k]=nums
        elif len(types)>0:
            counter=1
            for i in atomtypes[k]:
                if i in types:
                    for j in range(atomnums[k][atomtypes[k].index(i)]):
                        selected_atoms[k].append(counter)
                        counter+=1
                else:
                    counter+=atomnums[k][atomtypes[k].index(i)]
        else:
            selected_atoms[k]=[i for i in range(sum(atomnums))]
    
    #plots the 2d dos spectra as individual figures
    tempx=array([energies[0][start[0]:end[0]] for i in range(end[1]-start[1])])
    tempy=array([energies[1][start[1]:end[1]] for i in range(end[0]-start[0])]).transpose()
    if average_dos:
            fig,axs=plt.subplots(2,2,sharex='col',sharey='row',gridspec_kw={'width_ratios': [4,1], 'height_ratios': [4,1]})
            fig.delaxes(axs[1,1])
            proj_dos=zeros((end[1]-start[1],end[0]-start[0]))
            temp_dosx=zeros((3,end[0]-start[0]))
            temp_dosy=zeros((3,end[1]-start[1]))
    for i,j in zip(selected_atoms[0],selected_atoms[1]):
        if not average_dos:
            fig,axs=plt.subplots(2,2,sharex='col',sharey='row',gridspec_kw={'width_ratios': [4,1], 'height_ratios': [4,1]})
            fig.delaxes(axs[1,1])
            proj_dos=zeros((end[1]-start[1],end[0]-start[0]))
            temp_dosx=zeros((3,end[0]-start[0]))
            temp_dosy=zeros((3,end[1]-start[1]))
        for k in range(len(atomnums[0])):
            if i-1<sum(atomnums[0][:k+1]):
                atomlabel=atomtypes[0][k]
                break
        for k in range(len(orbitals[0])):
            for l in orbitals_to_plot:
                if l in orbitals[0][k]:
                    proj_dos+=outer(dos[1][j][k][start[1]:end[1]],dos[0][i][k][start[0]:end[0]])
                    if k<2:
                        m=0
                    if k>1 and k<8:
                        m=1
                    if k>7:
                        m=2
                    temp_dosx[m]+=dos[0][i][k][start[0]:end[0]]
                    temp_dosy[m]+=dos[1][j][k][start[1]:end[1]]
                    break
                
        if not average_dos:
            fig.suptitle('{} | contrbuting orbitals: {}'.format('{} #{}'.format(atomlabel,i-sum(atomnums[0][:atomtypes[0].index(atomlabel)])),', '.join(contributing_orbitals)))
            dos2d=axs[0,0].pcolormesh(tempx,tempy,proj_dos,shading='nearest',cmap=cmap)
            axs[0,0].plot([energies[0][start[0]:end[0]][k] for k in [0,-1]],[energies[1][start[1]:end[1]][k] for k in [0,-1]],color='red',linestyle='dashed')
            axs[0,1].set_ylabel('energy - $E_f$ / eV')
            axs[1,0].set_xlabel('energy - $E_f$ / eV')
            for k,l in zip(range(3),['s','p','d']):
                axs[0,1].plot(temp_dosy[k],energies[1][start[1]:end[1]],label=l)
                axs[1,0].plot(energies[0][start[0]:end[0]],temp_dosx[k],label=l)
            axs[0,1].set_xlim(axs[0,1].get_xlim()[1],axs[0,1].get_xlim()[0])
            axs[1,0].set_yticklabels([])
            axs[0,1].set_xticklabels([])
            axs[0,1].yaxis.set_label_position('right')
            axs[0,0].yaxis.tick_right()
            axs[0,1].yaxis.tick_right()
            axs[0,1].tick_params(axis='x', which='both', top=False, bottom=False)
            axs[1,0].tick_params(axis='y', which='both', top=False, bottom=False)
            axs[0,1].tick_params(labelright=True)
            axs[0,0].tick_params(labelright=False)
            handles, labels = axs[1,0].get_legend_handles_labels()
            fig.legend(handles, labels, bbox_to_anchor=(0.8,0.15), loc='lower right')
            fig.tight_layout()
            fig.canvas.draw()
            
    if average_dos:
        fig.suptitle('averaged over: {}\n contrbuting orbitals: {}'.format(', '.join(['{} #{}'.format(i,j) for i,j in zip(types,nums)]),', '.join(contributing_orbitals)))
        dos2d=axs[0,0].pcolormesh(tempx,tempy,proj_dos,shading='nearest',cmap=cmap)
        axs[0,0].plot([energies[0][start[0]:end[0]][k] for k in [0,-1]],[energies[1][start[1]:end[1]][k] for k in [0,-1]],color='red',linestyle='dashed')
        axs[0,1].set_ylabel('energy - $E_f$ / eV')
        axs[1,0].set_xlabel('energy - $E_f$ / eV')
        for k,l in zip(range(3),['s','p','d']):
            axs[0,1].plot(temp_dosy[k],energies[1][start[1]:end[1]],label=l)
            axs[1,0].plot(energies[0][start[0]:end[0]],temp_dosx[k],label=l)
        axs[0,1].set_xlim(axs[0,1].get_xlim()[1],axs[0,1].get_xlim()[0])
        axs[1,0].set_yticklabels([])
        axs[0,1].set_xticklabels([])
        axs[0,1].yaxis.set_label_position('right')
        axs[0,0].yaxis.tick_right()
        axs[0,1].yaxis.tick_right()
        axs[0,1].tick_params(axis='x', which='both', top=False, bottom=False)
        axs[1,0].tick_params(axis='y', which='both', top=False, bottom=False)
        axs[0,1].tick_params(labelright=True)
        axs[0,0].tick_params(labelright=False)
        handles, labels = axs[1,0].get_legend_handles_labels()
        fig.legend(handles, labels, loc='lower right')
        fig.tight_layout()
        fig.canvas.draw()
    
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
        if i in ['-o','--orbitals']:
            orbitals_to_plot=j.split(',')
        if i in ['-h','--help']:
            print('''
plotting options:
-a, --atomnums          specify site projected DOS to plot by the index of atoms: 1,2,3,etc...
-t, --types             specify which site projected DOS to plot by atom type: Au,C,etc...
-o, --orbitals          plot only the contributions from specific orbital projections
                        for example, to plot only s up and px down use: -o s_up,p_x_down

help options:
-h, --help               display this help message
                  ''')
            sys.exit()
    if exists(doscar):
        plot_dos(doscar,poscar,nums=atomnums,types=atomtypes,irange=irange,full=full,orbitals=orbitals_to_plot)
