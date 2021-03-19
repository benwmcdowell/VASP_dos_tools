from numpy import array,dot,zeros
import sys
from getopt import getopt

def calc_energy_shift(doscar,doscar_ref):
    dos, energies, ef = parse_doscar(doscar)
    dos_ref, energies_ref, ef_ref = parse_doscar(doscar_ref)
    
    efindex=[0,0]
    mindiff=[energies[-1]-energies[0] for i in range(2)]
    for i,j in zip(range(2),[energies,energies_ref]):
        for k in range(len(j)):
            if abs(energies[k])<mindiff[i]:
                efindex[i]=k
                mindiff[i]=abs(energies[k])
                
    total_dos=[zeros(len(energies[:efindex[0]])),zeros(len(energies_ref[:efindex[1]]))]
    for i,j in zip([dos,dos_ref],range(2)):
        for k in range(1,len(i)):
            for l in i[k]:
                total_dos[j]+=array(l[:efindex[j]])
                
    max_overlap=0.0
    for i in range(-len(energies[:efindex[0]]),len(energies[:efindex[0]])+1):
        if i<0:
            overlap=dot(total_dos[0][:i],total_dos[1][-i:])
        elif i>0:
            overlap=dot(total_dos[0][i:],total_dos[1][:-i])
        elif i==0:
            overlap=dot(total_dos[0],total_dos[1])
        if overlap>max_overlap:
            eshift=i
            max_overlap=overlap
    
    eshift=eshift*(energies[-1]-energies[0])/len(energies)
    
    return eshift    

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

if __name__=='__main__':
    short_opts='h'
    long_opts=['help']
    try:
        doscar=sys.argv[1]
        doscar_ref=sys.argv[2]
    except IndexError:
        print('missing required arguments. exiting...')
        sys.exit()
    try:
        opts,args=getopt(sys.argv[3:],short_opts,long_opts)
    except IndexError:
        print('error specifying optional arguments')
        sys.exit()
    for i,j in opts:
        if i in ['-h','--help']:
            print('''
help options:
    -h, --help                    the energy shift is calculated for the first argument, relative to the second
                                  to calculate the shift in DOS energies for a charged cell, use the charged cell as the first argument and the neutral cell as the second
''')
            sys.exit()
    try:
        eshift=calc_energy_shift(doscar,doscar_ref)
        print('shift in energies is {} eV for {} relative to {}'.format(eshift,doscar,doscar_ref))
    except NameError:
        print('incorrect specification of files. exiting...')
        sys.exit()
