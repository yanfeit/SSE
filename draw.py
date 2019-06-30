import os
from itertools import izip
import matplotlib.pyplot as plt


def num(s):
    """This function is used to convert string to int or float."""
    try:
        return int(s)
    except ValueError:
        return float(s)

def CollectData(filename = "lx20t00.06t10.01mu0.02beta20.dat"):
    """This function is used to collect informatin from files."""
    entry = {}
    with open(filename, 'r') as ifile:
        for line in ifile:
            args = line.split()
            key = args[0]
            if (len(args) == 2):
                value = num(args[1])
            else:
                value = [num(args[1]), num(args[2])]
            entry[key] = value
    return entry
    
def CreateDataSet():
    """This function is used to collect all files' data. Return a list"""
    data = []
    files = [f for f in os.listdir('.') if os.path.isfile(f) and f[-3:] == 'dat']
    for file in files:
        data.append(CollectData(file))
        
    return data
    
def CreateLine(data, key = "t_parallel"):
    axis = []
    for i in range(len(data)):
        axis.append(data[i][key])
    return axis

    
if __name__ == "__main__":
    
    data = CreateDataSet()
    mu = CreateLine(data, "mu")
    density = CreateLine(data, "density")
    rhos = CreateLine(data, "rhos")
    rhosTT = CreateLine(data, "rhosTT")
    rhosTB = CreateLine(data, "rhosTB")
    compressibility = CreateLine(data, "compressibility")
    
    sorted_lists = sorted(izip(density, rhos, rhosTT, rhosTB, compressibility, mu),\
    key = lambda x: x[5])
    density, rhos, rhosTT, rhosTB, compressibility, mu = [[x[i] for x in sorted_lists\
    ] for i in range(6)]
    
    fig, ax1 = plt.subplots()
    ax1.set_xlabel(r'$\mu/U$', fontsize = 'x-large')
    ax1.set_ylabel(r'$\rho$  $\kappa$  $\tilde{\rho}^{TT}_{s}$  $\tilde{\rho}^{TB}_{s}$',\
    fontsize = 'x-large')
    ax2 = ax1.twinx()
    ax2.set_ylabel(r'$\rho_{s}$', fontsize = 'x-large')
    
    
    l1 = ax1.errorbar(mu, [i[0] for i in density], yerr = [i[1] for i in density],\
    color = 'green', label = r'$\rho$')
    l3 = ax1.errorbar(mu, [i[0] for i in compressibility], yerr = [i[1] for i in compressibility],\
    color = 'blue', label = r'$\kappa$')
    l4 = ax1.errorbar(mu, [i[0] for i in rhosTT], yerr = [i[1] for i in rhosTT],\
    color = 'burlywood', label = r'$\tilde{\rho}^{TT}_{s}$')
    l5 = ax1.errorbar(mu, [i[0] for i in rhosTB], yerr = [i[1] for i in rhosTB],\
    color = 'purple', label = r'$\tilde{\rho}^{TB}_{s}$')        
    l2 = ax2.errorbar(mu, [i[0] for i in rhos], yerr = [i[1] for i in rhos],\
    color = 'red', label = r'$\rho_{s}$')
    
    lines = [l1, l2, l3, l4, l5]
    ax1.legend(lines, [l.get_label() for l in lines])
   
    plt.title(r'$t_{\parallel}/U = 0.02$, $t_{\perp}/U = 0.01$')
    plt.show()
    
    
    
    
    
    
    




