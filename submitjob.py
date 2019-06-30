import random
import subprocess
import numpy as np

#seed = random.randrange(1, 1000000000, 2)
#parameter = {'lx':8, 'ly':8,
#            'beta': 10, 'mpps': 4, 'warmsteps': 10000, 'mcsteps': 100000,
#           't_parallel': 0.06, 't_perp': 0.01, 'interaction': 0.00, 'chempot': 0.5,
#          'seed': seed
#         }

def WriteInputFile(parameter = {'lx':8, 'ly':8,
                 'beta': 10.0, 'mpps': 4, 'warmstpes': 10000, 'mcsteps': 100000,
                 't_paralle': 0.06, 't_perp': 0.01, 'interaction': 0.00, 'chempot': 0.5,
                 'seed': 434535
                 }):
    """This function is used to generate input file for bilayer program."""
        
    context = "{0}\t{1}\n{2}\t{3}\t{4}\t{5}\n{6}\t{7}\t{8}\t{9}\n{10}\n\n\nlx ly\nbeta mpps warmsteps mcsteps\nt_paralle t_perp interaction chempot\n(random number seed)".format(\
    parameter['lx'], parameter['ly'], parameter['beta'], parameter['mpps'],\
    parameter['warmsteps'], parameter['mcsteps'], parameter['t_parallel'], parameter['t_perp'],\
    parameter['interaction'], parameter['chempot'], parameter['seed'])
    
    filename = 'lx' + str(parameter['lx']) + 't0' + str(parameter['t_parallel']) +\
    't1' + str(parameter['t_perp']) + 'mu' + str(parameter['chempot']) + 'beta' +\
    str(parameter['beta']) + '.in'
    
    ifile = open(filename, 'wb')
    ifile.write(context)
    
    return filename
    
def WriteQsubFile(parameter = {'lx':8, 'ly':8,
                 'beta': 10.0, 'mpps': 4, 'warmstpes': 10000, 'mcsteps': 100000,
                 't_paralle': 0.06, 't_perp': 0.01, 'interaction': 0.00, 'chempot': 0.5,
                 'seed': 434535
                 }, walltime = '10:00:00'):    
    """This function is used to generate shell script"""
    
    filename = 'lx' + str(parameter['lx']) + 't0' + str(parameter['t_parallel']) +\
    't1' + str(parameter['t_perp']) + 'mu' + str(parameter['chempot']) + 'beta' +\
    str(parameter['beta']) + '.sh'
    
    context = "#PBS -lwalltime={0}\n#PBS -lnodes=1\n#PBS -j oe\ndate\ncd $PBS_O_WORKDIR\n~/bin/bilayer \t{1}\ndate\n\
    ".format(walltime, WriteInputFile(parameter))
    
    ifile = open(filename,'wb')
    ifile.write(context)
    
    return filename
    
def SubmitFile(filename):
    subprocess.call(['qsub', filename])
    
    


if __name__ == '__main__':
    
    for t_perp in np.linspace(0, 0.02, 101):
        seed = random.randrange(1, 10000000000, 2)
        parameter = {'lx':20, 'ly':20,
                     'beta': 20, 'mpps': 4, 'warmsteps': 10000, 'mcsteps': 100000,
                     't_parallel': 0.02, 't_perp': t_perp, 'interaction': 0.00, 'chempot': 0.9,
                     'seed': seed
                     }
        SubmitFile(WriteQsubFile(parameter))
         
    
    
                
    
    
    


