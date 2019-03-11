import numpy as np               # http://www.numpy.org/
import matplotlib.pyplot as plt    # https://matplotlib.org/
from pathlib import Path         # https://docs.python.org/3.5/library/pathlib.html#module-pathlib
import molsim_utilities as mu
import pandas as pd

##fixing numbers, opposite charge : neg/pos


strtype = ['m40ncen', 'm20m20', 'm10m20m10', 'm10m10m10m10', 'm5m10m5', 'm5m20m5', 'm10m10m10', 'm20m10', 'm5m30m5', 'm15m10m15', 'm30m10', 'm10m5m10', 'm15m5m15' ] 
chargeratio = ['0', '1', '1', '1', '1', '2', '2', '2', '3', '3', '3' , '4', '6' ]
monomernumber = ['40', '40', '40', '40', '20', '30', '30', '30', '40', '40', '40', '25', '35' ]
polymer = ['0', '2', '3', '4', '3', '3', '3', '3', '2', '3', '3', '2', '3', '3' ] 
#strtype = ['m30m10', 'm5m30m5','m15m10m15']



pathstring  = []


ph = []
#ree = pd.DataFrame(data = [0,0,0,0,0,0,0], columns = ['Monomer index', 'Prob', 'Err', 'pH', 'Structure', 'Charge ratio', 'R-type'])
rg = [] 
#linecolor = ['wheat', 'gold', 'orange', 'darkorange', 'chocolate', 'sienna', 'saddlebrown']
#linecolor2 = ['powderblue', 'lightblue', 'lightskyblue', 'royalblue', 'blue', 'mediumblue', 'darkblue']

for i in range(len(strtype)):
    pathstring.append('alsi/polymer/quenchedpoly/{:s}'.format(str(chargeratio[i])) + '/{:s}/'.format(str(strtype[i])))

reesum = pd.DataFrame
startcount = 1
for j in range(len(pathstring)):
    for p in sorted(Path(str(pathstring[j])).glob('zero/ph*/cylinder_shell.list')):
        print(p)  
        reesingle = mu.getDistribution(p, 'ree pa')  
        rgsingle = mu.getDistribution(p, 'rg pa') 
        #f = mu.getDistribution(p, 'rg pa')

        labelstring = 'ph = {:d}'.format(int(p.parts[-2][2])) 
        print(p.parts)
        print(p.parts[-2])
        print(p.parts[-2][2])
        ph.append(p.parts[-2][2])
                
        ree = pd.DataFrame(reesingle, columns = ['Monomer index', 'Prob', 'Err']) 
        ree['pH'] = int(p.parts[-2][2])
        ree['Structure'] = str(p.parts[-4])
        ree['Polymer type'] = polymer[j]
        ree['Charge ratio'] = str(p.parts[-5])
        ree['Monomer number'] = monomernumber[j]
        ree['R-type'] = 'Ree'  
            
        rg = pd.DataFrame(rgsingle, columns = ['Monomer index', 'Prob', 'Err'])
        rg['pH'] = int(p.parts[-2][2])
        rg['Structure'] = str(p.parts[-4])
        rg['Polymer type'] = polymer[j]
        rg['Charge ratio'] = str(p.parts[-5])
        rg['Monomer number'] = monomernumber[j]
        rg['R-type'] = 'Rg' 

        local = str(p.parts[-4])
        if local.find('cen') == True: ##이 부분이 안됨!!!!! 
            ree['Location'] = 'Inside'
            rg['Location'] = 'Inside'
        else:
            ree['Location'] = 'Outside'
            rg['Location'] = 'Outside'
    
    
        if startcount == 1:

            #reesum = ree
            ree.to_csv("ree.csv", mode="w")
            rg.to_csv("ree.csv", mode="a", header=False)
        else:    
            ree.to_csv("ree.csv", mode="a", header=False)
            rg.to_csv("ree.csv", mode="a", header=False)

            #pd.merge(reesum, ree, on = 'Monomer index')

        startcount += 1

            
 
'''         
        plt.figure(1)
        plt.plot(d[:,0], d[:,1], label = labelstring, color = '{:s}'.format(linecolor[int(p.parts[-2][2])-2])) 
        plt.title('{:s}'.format(str(strtype[j])) + 'Ree')
        
        plt.figure(2)
        plt.plot(e[:,0], e[:,1], label = labelstring, color = '{:s}'.format(linecolor[int(p.parts[-2][2])-2])) 
        plt.title('{:s}'.format(str(strtype[j])) + 'Rg')

    plt.figure(1)
    plt.ylim(0.000, 1.000)
    plt.ylabel('Prob(Ree)')
    plt.xlabel('Ree [Å]')
    plt.legend()
    
    plt.figure(2)
    plt.ylim(0.000, 1.000)
    plt.ylabel('Prob(Rg)')
    plt.xlabel('Rg [Å]')
    plt.legend()
            
    plt.show()                                                              

'''