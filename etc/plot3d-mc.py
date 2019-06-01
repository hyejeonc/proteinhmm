# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 14:12:25 2019

@author: HYEJEONG
>2P4FA
$ sn:6945 prot_size:286 total_dr:93 num_dr:2 size_dr:81,12 #1-81 #152-163
EDRYKEKLLQKAKAEGVESIEELKKRLADQIEEKKKELNKIDPLRELEQHLNAGSRIHTNKEHKTTKMSNKSNEKSGNVLPKDKPYKTLDDYLKLDKIKDLSKQEVEFLWRAKWSNRDDSLVAVVPYVKTFQGMYKYAVKNPLFVLPLPRENAADGNKADKDSVPVELQYVQWQFAGPNTVHCLITSLAEYKLHQDFAKPHTTIQFHLDLANDKDMVLMNGQVESDSNVSLQDAQLLLLNVQRFYGAMGSETSIAKERIQLLEDFNKGSQNFDINKLIQLAQSMEN
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++----------------------------------------------------------------------++++++++++++---------------------------------------------------------------------------------------------------------------------------

"""
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

seq = "EDRYKEKLLQKAKAEGVESIEELKKRLADQIEEKKKELNKIDPLRELEQHLNAGSRIHTNKEHKTTKMSNKSNEKSGNVLPKDKPYKTLDDYLKLDKIKDLSKQEVEFLWRAKWSNRDDSLVAVVPYVKTFQGMYKYAVKNPLFVLPLPRENAADGNKADKDSVPVELQYVQWQFAGPNTVHCLITSLAEYKLHQDFAKPHTTIQFHLDLANDKDMVLMNGQVESDSNVSLQDAQLLLLNVQRFYGAMGSETSIAKERIQLLEDFNKGSQNFDINKLIQLAQSMEN" 
l = len(seq)
print(seq[0])
print(seq[1])




fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

Axes3D.scatter(xs, ys, zs=0, zdir='z', s=20, c=None, depthshade=True, *args, **kwargs)