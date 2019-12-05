                #######################################################################################
                ###                                                                                 ###  
                ###          """ Projet SNA : Dynamique de population, modele de Galton-Watson"""   ###  
                ###          """ Ismail Atourki && Ilyas Malik """                                  ###  
                ###                                                                                 ###  
                #######################################################################################
                
"""Cas sous critique """         

       
""" Q1: la probabilite de non-extinction"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sps
import scipy as sp
import numpy.random as npr
from timeit import default_timer as timer
plt.close("all")
###############################
""" Données du problème """
###############################

n=int(30)
m=int(1E3)
moyb=.1       # Le nombre moyen d'enfant par individu
moyh=1.5
pb=.5
moyprim=.97
alpha=0.05
q=sps.norm.ppf(.975)
initmin=20    # La loi de la population initiale suit une loi 
initmax=21    # uniforme dans [initmin,initmax]


#########################################################
""" Méthode de Monte-Carlo naif"""
#########################################################
def moyenne(pb):
    return npr.choice(a=[moyb,moyh],p=[pb,1-pb])


def mat(pb):
    Y=[]
    for j in range(m):
#        if (int(4*j/m)==4*j/m):
        print("tourne : ",int(100*j/m),"%")
            
        X=[npr.randint(initmin,initmax)]
        for i in range (n):
            y=np.random.geometric(1/(moyenne(pb)+1),size=X[-1])-1
            X.append(np.sum(y))
        Y.append(X)
    Y=np.asarray(Y)
    return Y

Y=mat(.4)

Z = (Y[:,-1]>0)
v=np.mean(Z)
s=np.std(Z)
d=s/np.sqrt(m)*q
print("Probabilité de non extinction:{}".format(v))
print("l'intervalle de confiance à 95% est [{:.6f},{:.6f}]".format(v-d,v+d))



nb=6

for i in range (m):
    if nb==0:
        break
    if Y[i,n]>0:
        nb-=1
        plt.plot(Y[i],label="mb={},mh={},pb={}")
plt.xlabel("Générations")


















