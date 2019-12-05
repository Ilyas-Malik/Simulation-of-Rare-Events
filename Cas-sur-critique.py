

                #######################################################################################
                ###                                                                                 ###  
                ###          """ Projet SNA : Dynamique de population, modele de Galton-Watson"""   ###  
                ###          """ Ismail Atourki && Ilyas Malik """                                  ###  
                ###                                                                                 ###  
                #######################################################################################

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sps
import scipy as sp
from timeit import default_timer as timer

plt.close("all")

x=100
alpha=.05
n=int(40)
m=int(1e4)
moy=1.2 # Le nombre moyen d'enfant par individu
p=1/(moy+1)
moyprim=.97
alpha=0.05
q=sps.norm.ppf(.975)
x0_min=1
x0_max=2
c=1.08
""" Methode de Monte-Carlo naïf P(0<X<c^n)"""
#
# cette fonction calculer les valeur de la probabilité pour diffirents valeurs de m 
#ainsi lque les bornes des intervalles de confiance  afin de tracer les courbe de convergence vers la valeur exact
def monte_carlo(m):
        
    X=np.random.randint(x0_min,x0_max,size=m)
    Z=np.zeros(shape=(n,m))
    Z[0]=X
    for k in range(1,n):
        for i in range(m):
            y=np.random.geometric(1/(moy+1),size=int(Z[k-1,i]))-1
            Z[k,i]=np.sum(y)
    xz=np.logical_and(Z[-1,:]>0,Z[-1,:]<=c**n)
        
    s=np.std(xz)
    q_alpha=sps.norm.ppf(1-alpha/2)
    b=s*q_alpha/np.sqrt(m)
    p=np.mean(xz)
    return p,p-b,p+b
     

l_monte=np.arange(1000,10000,1000)
monte_p=[]
inf_b=[]
sup_b=[]
for j in l_monte:
    #print(j)
    monte_p.append(monte_carlo(j)[0])
    inf_b.append(monte_carlo(j)[1])
    sup_b.append(monte_carlo(j)[2])
    
plt.plot(l_monte,monte_p,color="b",label="Probabilite numerique ")
plt.plot(l_monte,sup_b,color="r",label="Borne_sup de IC")
plt.plot(l_monte,inf_b,color="r",label='Borne_inf de IC')


#### Calcul d'un seul valeur exacte
#
#X=np.random.randint(x0_min,x0_max,size=m)
#Z=np.zeros(shape=(n,m))
#Z[0]=X
#for k in range(1,n):
#    for i in range(m):
#            y=np.random.geometric(1/(moy+1),size=int(Z[k-1,i]))-1
#            Z[k,i]=np.sum(y)
#    xz=np.logical_and(Z[-1,:]>0,Z[-1,:]<=c**n)
#    p=np.mean(xz)  
#xz=np.logical_and(Z[-1,:]>0,Z[-1,:]<=c**n)
#p=np.mean(xz)
#print("probabilite P(0<X<c^n) = {}".format(p))
#
#
### Scenario pour ou 0<X_n<c^n
#bb=0
#max=50
#l=np.arange(0,n,1)
#for k in range(m):
#    if  bb<=max:
#        plt.plot(l,Z[:,k])
#        bb+=1
#exp=c**l
#
#plt.plot(l,exp,color="black",label="Fonction de croissance pour c={} et m={}".format(c,moy))
#plt.legend(loc='best')
#plt.ylim([0,20])
#
#
##intervalle de confiance
#
#s=np.std(xz)
#q_alpha=sps.norm.ppf(1-alpha/2)
#b=s*q_alpha/np.sqrt(m)
#
#print("Intervalle de confiance pour la probabilite p au niveau {} : [{:.4f}, "
#      "{:.4f}]".format(1 - alpha, p - b, p + b))
#    
#    
###    
#
###########################################################
#"""Methode Echantillonnage preferentiel  P(0<X<c^n)"""
###########################################################
#
#
#def P(x,y):
#    temp=np.logical_and((x==y),(x==0))
#    r= p**(x+temp)*(1-p)**y*sp.special.binom(x+temp+y-1,x+temp-1)
#    r*=1-temp
#    r+=temp
#    return r
#
#
#def h(X):
#    return np.logical_and(X[:,-1]>0,X[:,-1]<=c**n)
#
#def Q(x,y,m):
#    temp=np.logical_and((x==y),(x==0))
#    p=1/(m+1)
#    r= p**(x+temp)*(1-p)**y*sp.special.binom(x+temp+y-1,x+temp-1)
#    r*=1-temp
#    r+=temp
#    return r
#
#
#def Lp(X,m):
#    c=np.ones(np.size(X,axis=0))
#    for k in range(1,np.size(X,axis=1)):
#        c*=P(X[:,k-1],X[:,k])/Q(X[:,k-1],X[:,k],m)    
#    return c
#
#def mat(mo):
#    Y=[]
#    for j in range(m):
##        if (int(4*j/m)==4*j/m):
##            print("tourne : ",int(100*j/m),"%")
#            
#        X=[np.random.randint(x0_min,x0_max)]
#        for i in range (n):
#            y=np.random.geometric(1/(mo+1),size=X[-1])-1
#            X.append(np.sum(y))
#        Y.append(X)
#    Y=np.asarray(Y)
#    return Y
#
#Y=mat(moyprim)
#pm=Lp(Y,moyprim)*h(Y)
#v=np.mean(pm)
#s=np.std(pm)
#d=s/np.sqrt(m)*q
#
#print("\n\nEchantillionnage preferentiel: ",m,"simulations \n\nValeur: {}".format(v))
#
#
##intervalle de confiance
#print("\nIntervalle de confiance à 95% : [{:.8f},{:.8f}]".format(v-d ,v+d ))
#
#print(" ")






    
    
