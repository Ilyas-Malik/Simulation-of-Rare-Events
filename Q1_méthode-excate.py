#######################################################################################
###                                                                                 ###
###          """ Projet SNA : Dynamique de population, modele de Galton-Watson"""   ###
###          """ Ismail Atourki && Ilyas Malik """                                  ###
###                                                                                 ###
#######################################################################################

## Méthode exact de calcul de probabilité de non extinction

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import scipy.stats as sps
import numpy.random as npr
plt.close("all")

n=int(40)
moy=0.7
p=1/(moy+1)
r=int(5E2)
initmin=20

################# Méthode itérative

def Probaextinction():
    a=[0]
    for k in range (n):
        a.append(1/(1+moy*(1-a[-1])))
    a=1-np.array(a)**initmin    
    return a[-1]

print("Méthode exacte 1",Probaextinction())

############## Méthode avec la matrice de transition
def P1(x,y):
    if x==y==0:
        return 1
    return p**x*(1-p)**y*sp.special.binom(x+y-1,x-1)

Y1=np.ones((r,r))

for i in range(r):
    for j in range(r):
        Y1[i,j]=P1(i,j)

Y1=np.linalg.matrix_power(Y1,n)
print(" ")
print("Méthode exacte 2: {} ".format(1-Y1[initmin,0]))
print("\nDifférence:",1-Y1[initmin,0]-Probaextinction())





