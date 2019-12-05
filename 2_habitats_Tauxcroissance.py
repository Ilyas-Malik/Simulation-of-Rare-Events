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





def tauxcroissance(pb,moyb,moyh,moys,f):
    return ( f*moyb+(1-f)*moys )**pb*( f*moyh+(1-f)*moys )**(1-pb)


x=np.linspace(0,1,100)

y=np.zeros(100)
for i in range(100):
    y[i]=tauxcroissance(pb,moyb,moyh,moys,x[i])

plt.plot(x,y,label="Taux de croissance")
plt.vlines(f,.5,1.2,"r",label="valeur optimale de f (75%)")
plt.hlines(1,0,1,"y",label="1")
plt.legend(loc="best")
plt.xlabel("fraction f")



















