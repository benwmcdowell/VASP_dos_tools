import numpy as np
import matplotlib.pyplot as plt

def parse_eigenval(ifile):
    with open(ifile) as file:
        for i in range(6):
            line=file.readline().split()
        nstates=int(line[0])
        for i in range(2):
            file.readline()
        eigenval=np.zeros((nstates))
        for i in range(nstates):
            line=file.readline().split()
            eigenval[i]+=(float(line[1])+float(line[2]))/2
        
    return eigenval,nstates
