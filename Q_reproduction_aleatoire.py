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

n=int(20)
m=int(1E2)
moyb=.4       # Le nombre moyen d'enfant par individu
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
#            print("tourne : ",int(100*j/m),"%")
            
        X=[npr.randint(initmin,initmax)]
        for i in range (n):
            y=np.random.geometric(1/(moyenne(pb)+1),size=X[-1])-1
            X.append(np.sum(y))
        Y.append(X)
    Y=np.asarray(Y)
    Y = (Y[:,-1]>0)
    v=np.mean(Y)
    s=np.std(Y)
    d=s/np.sqrt(m)*q
    return v,d

x=[.1*i for i in range(11)]
Y=np.zeros(11)
D=np.zeros(11)

moyb=.8
moyh=1.2

for i in range(len(x)):
    print(i)
    Y[i],D[i]=mat(x[i])
plt.subplot(2,2,1)
plt.plot(x,Y,label="mb={}, mh={}".format(moyb,moyh))
plt.plot(x,Y+D,"r")
plt.plot(x,Y-D,"r")
plt.ylabel("P(Xn>0)")
plt.legend(loc="best")

moyb=.6
moyh=1.3

for i in range(len(x)):
    print(i)
    Y[i],D[i]=mat(x[i])
plt.subplot(2,2,2)
plt.plot(x,Y,label="mb={}, mh={}".format(moyb,moyh))
plt.plot(x,Y+D,"r")
plt.plot(x,Y-D,"r")
plt.ylabel("P(Xn>0)")
plt.legend(loc="best")

moyb=.4
moyh=1.4

for i in range(len(x)):
    print(i)
    Y[i],D[i]=mat(x[i])
plt.subplot(2,2,3)
plt.plot(x,Y,label="mb={}, mh={}".format(moyb,moyh))
plt.plot(x,Y+D,"r")
plt.plot(x,Y-D,"r")
plt.ylabel("P(Xn>0)")
plt.legend(loc="best")


moyb=0.2
moyh=1.8

for i in range(len(x)):
    print(i)
    Y[i],D[i]=mat(x[i])
plt.subplot(2,2,4)
plt.plot(x,Y,label="mb={}, mh={}".format(moyb,moyh))
plt.plot(x,Y+D,"r")
plt.plot(x,Y-D,"r")
plt.ylabel("P(Xn>0)")
plt.legend(loc="best")















