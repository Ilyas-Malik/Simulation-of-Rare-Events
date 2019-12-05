#######################################################################################
###                                                                                 ###
###          """ Projet SNA : Dynamique de population, modele de Galton-Watson"""   ###
###          """ Ismail Atourki && Ilyas Malik """                                  ###
###                                                                                 ###
#######################################################################################

## Methode exact de calcul de probabilite de non extinction

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sps
import numpy.random as npr
plt.close("all")

n=int(40)

moy=0.8
initmin=20


l=[1,2,5,10,20,40,100]
m=[.1,.3,.5,.7,.9,.99]
res=[]

############## Influence de X0
for x in l:
    loc=[]
    for i in range(n):
        loc.append(Probaextinction(.9,i,x))
    res.append(loc)
res=np.asarray(res)

for i in range (len(l)):
    plt.plot(res[i],label="X0={}".format(l[i]))
    plt.legend(loc="best")
plt.xlabel("Instant n")
plt.ylabel("Probabilité d'extinction P(Xn)")


########### Influence de la moyenne m
plt.figure()
res=[]
for x in m:
    loc=[]
    for i in range(n):
        loc.append(Probaextinction(x,i,20))
    res.append(loc)
res=np.asarray(res)

for i in range (len(m)):
    plt.plot(res[i],label="m={}".format(m[i]))
    plt.legend(loc="best")
plt.xlabel("Instant n")
plt.ylabel("Probabilité d'extinction P(Xn)")


