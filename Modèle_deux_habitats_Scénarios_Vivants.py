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
m=int(1E2)
alpha=0.05
q=sps.norm.ppf(.975)
initmin=20    # La loi de la population initiale suit une loi 
initmax=21    # uniforme dans [initmin,initmax]
moyb=0.00001
moyh=1.5
moys=0.9
pb=0.1
f=moys*( (moyh-moys) -pb*(moyh-moyb) )/( (moys-moyb)*(moyh-moys) )


#########################################################
""" Méthode de Monte-Carlo naif"""
#########################################################
def moyenne(pb,moyb,moyh):
    return npr.choice(a=[moyb,moyh],p=[pb,1-pb])


def mat(pb,moyb,moyh,moys,f):          #f est la proportion de la population dans l'habitat 1
    Y=[]
    for j in range(m):
#        if (int(4*j/m)==4*j/m):
#        print("tourne : ",int(100*j/m),"%")
            
        X=[npr.randint(initmin,initmax)]
        for i in range (n):
            f1=int(f*X[-1])
            y1=np.sum(np.random.geometric(1/(moyenne(pb,moyb,moyh)+1),size=f1)-1)
            y2=np.sum(np.random.geometric(1/(moys+1),size=X[-1]-f1)-1)
            X.append(y1+y2)
        Y.append(X)
    Y=np.asarray(Y)
    return Y

x=[i for i in range(n+1)]
Y=mat(pb,moyb,moyh,moys,f)
nb=6
for i in range(m):
    if nb==0:
        break
    if Y[i,-1]!=0:
        nb-=1
        plt.plot(x,Y[i,:])
plt.title("scénarios non éteints pour f_optimale")
plt.xlabel("Générations n")
plt.ylabel("Xn")



















