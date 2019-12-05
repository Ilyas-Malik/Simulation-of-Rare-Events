                #######################################################################################
                ###                                                                                 ###  
                ###          """ Projet SNA : Dynamique de population, modele de Galton-Watson"""   ###  
                ###          """ Ismail Atourki && Ilyas Malik """                                  ###  
                ###                                                                                 ###  
                #######################################################################################
                
"""Cas sous critique """         

       


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
m=int(1e4)
moy=.7  # Le nombre moyen d'enfant par individu
p=1/(moy+1)
moyprim=.97
alpha=0.05
q=sps.norm.ppf(.975)
initmin=20    # La loi de la population initiale suit une loi 
initmax=21    # uniforme dans [initmin,initmax]
x=150


""" Q4.1 Estimation de P(somme de X_n>x)"""

#########################################################
""" Méthode de Monte-Carlo naif"""
#########################################################
#

print("Estimation de P(somme de X_n>x)")
def mat(mo):
    Y=[]
    for j in range(m):
#        if (int(4*j/m)==4*j/m):
#            print("tourne : ",int(100*j/m),"%")
            
        X=[npr.randint(initmin,initmax)]
        for i in range (n):
            y=np.random.geometric(1/(mo+1),size=X[-1])-1
            X.append(np.sum(y))
        Y.append(X)
    Y=np.asarray(Y)
    return Y

Y=mat(moy)

## pour tracer les scenario des resources consommees
plt.figure()
Z = (np.sum(Y,axis=1)>x)
A=np.cumsum(Y,axis=1)
l=np.arange(0,len(A[0]),1)
for k in range (40):
    
    plt.plot(l,A[k])
    
plt.ylabel("Somme des X_n")
plt.xlabel("generations")
plt.legend(loc="best")

plt.hlines(y=x,xmin=0,xmax=50,color='black',label="x=".format(x))
v=np.mean(Z)
s=np.std(Z)
d=s/np.sqrt(m)*q

# Calcul de la probabilite et des iuntervalle de confiance
print("\nMonte-Carlo naif: ",m,"simulations\n\nValeur: {}".format(v))

print(" ")
#intervalle de confiance
print("\nIntervalle de confiance à 95% : [{:.8f},{:.8f}]".format(v-d ,v+d ))

    
    
    
#print(" ")
###########################################################
"""Methode Echantillonnage preferentiel"""
###########################################################

def P(x,y):
    temp=np.logical_and((x==y),(x==0))
    r= p**(x+temp)*(1-p)**y*sp.special.binom(x+temp+y-1,x+temp-1)
    r*=1-temp
    r+=temp
    return r


def h(X):
    return np.sum(X,axis=1)>x

def Q(x,y,m):
    temp=np.logical_and((x==y),(x==0))
    p=1/(m+1)
    r= p**(x+temp)*(1-p)**y*sp.special.binom(x+temp+y-1,x+temp-1)
    r*=1-temp
    r+=temp
    return r


def Lp(X,m):
    c=np.ones(np.size(X,axis=0))
    for k in range(1,np.size(X,axis=1)):
        c*=P(X[:,k-1],X[:,k])/Q(X[:,k-1],X[:,k],m)    
    return c

def mat(mo):
    Y=[]
    for j in range(m):
#        if (int(4*j/m)==4*j/m):
#            print("tourne : ",int(100*j/m),"%")
            
        X=[npr.randint(initmin,initmax)]
        for i in range (n):
            y=np.random.geometric(1/(mo+1),size=X[-1])-1
            X.append(np.sum(y))
        Y.append(X)
    Y=np.asarray(Y)
    return Y

Y=mat(moyprim)
pm=Lp(Y,moyprim)*h(Y)
v=np.mean(pm)
s=np.std(pm)
d=s/np.sqrt(m)*q

print("\n\nEchantillionnage preferentiel: ",m,"simulations \n\nValeur: {}".format(v))

print(" ")
#intervalle de confiance
print("\nIntervalle de confiance à 95% : [{:.8f},{:.8f}]".format(v-d ,v+d ))

print(" ")



       
""" Q4.2 Estimation de P( max X_i>x)"""

x=22
##########################################################
#""" Méthode de Monte-Carlo naif"""
##########################################################
print("Estimation de P(max de X_n>x)")
def mat(mo):
    Y=[]
    for j in range(m):
#        if (int(4*j/m)==4*j/m):
#            print("tourne : ",int(100*j/m),"%")
            
        X=[npr.randint(initmin,initmax)]
        for i in range (n):
            y=np.random.geometric(1/(mo+1),size=X[-1])-1
            X.append(np.sum(y))
        Y.append(X)
    Y=np.asarray(Y)
    return Y

Y=mat(moy)
Z = (np.max(Y,axis=1)>x)

v=np.mean(Z)
s=np.std(Z)
d=s/np.sqrt(m)*q


## scenario pour maxX_i
plt.figure()
A=[]

    
A=np.maximum.accumulate(Y,axis=1)
l=np.arange(0,len(A[0]),1)
bb=0
max=3
for k in range (10):
    plt.plot(l,A[k])
        
for k in range (m):
    if A[k,-1]>x and bb<max:
        plt.plot(l,A[k])
        bb+=1
    
plt.ylabel("Max des X_n")
plt.xlabel("generations")
plt.hlines(y=x,xmin=0,xmax=50,color='black')
plt.legend(loc="best")

print(" ")
##probabilite et intervalle de confiance
print("\nMonte-Carlo naif: ",m,"simulations\n\nValeur: {}".format(v))

print(" ")
#intervalle de confiance
print("\nIntervalle de confiance à 95% : [{:.8f},{:.8f}]".format(v-d ,v+d ))

    
    
    
##print(" ")
###########################################################
"""Methode Echantillonnage preferentiel"""
###########################################################
def P(x,y):
    temp=np.logical_and((x==y),(x==0))
    r= p**(x+temp)*(1-p)**y*sp.special.binom(x+temp+y-1,x+temp-1)
    r*=1-temp
    r+=temp
    return r


def h(X):
    return np.max(X,axis=1)>x

def Q(x,y,m):
    temp=np.logical_and((x==y),(x==0))
    p=1/(m+1)
    r= p**(x+temp)*(1-p)**y*sp.special.binom(x+temp+y-1,x+temp-1)
    r*=1-temp
    r+=temp
    return r


def Lp(X,m):
    c=np.ones(np.size(X,axis=0))
    for k in range(1,np.size(X,axis=1)):
        c*=P(X[:,k-1],X[:,k])/Q(X[:,k-1],X[:,k],m)    
    return c

def mat(mo):
    Y=[]
    for j in range(m):
#        if (int(4*j/m)==4*j/m):
#            print("tourne : ",int(100*j/m),"%")
            
        X=[npr.randint(initmin,initmax)]
        for i in range (n):
            y=np.random.geometric(1/(mo+1),size=X[-1])-1
            X.append(np.sum(y))
        Y.append(X)
    Y=np.asarray(Y)
    return Y

Y=mat(moyprim)
pm=Lp(Y,moyprim)*h(Y)
v=np.mean(pm)
s=np.std(pm)
d=s/np.sqrt(m)*q


print(" ")
print("\n\nEchantillionnage preferentiel: ",m,"simulations \n\nValeur: {}".format(v))

print(" ")
#intervalle de confiance
print("\nIntervalle de confiance à 95% : [{:.8f},{:.8f}]".format(v-d ,v+d ))
















#
