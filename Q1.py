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

n=int(40)
m=int(1E5)
moy=.7       # Le nombre moyen d'enfant par individu
p=1/(moy+1)
moyprim=.97
alpha=0.05
q=sps.norm.ppf(.975)
initmin=20    # La loi de la population initiale suit une loi 
initmax=21    # uniforme dans [initmin,initmax]


#########################################################
""" Methode exacte """### : Matrice de transition de la chaine de Markov associe
#########################################################

def Probaextinction():
    a=[0]
    for k in range (n):
        a.append(1/(1+moy*(1-a[-1])))
    a=1-np.array(a)**initmin    
    return a[-1]

ve=Probaextinction()
print("\nValeur exacte:{:.8f}".format(ve))


#########################################################
""" Méthode de Monte-Carlo naif"""
#########################################################

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

Y = (Y[:,-1]>0)
v=np.mean(Y)
s=np.std(Y)
d=s/np.sqrt(m)*q

print("\nMonte-Carlo naif: ",m,"simulations\n\nValeur: {}".format(v))
print("Erreur relative: {:.2f} %".format(np.abs(v-ve)/ve*100))

#intervalle de confiance
print("\nIntervalle de confiance à 95% : [{:.8f},{:.8f}]".format(v-d ,v+d ))
if (v!=0):
    print("Largeur relative: {:.2f} %".format(d/v*100))
print(" ")


#########################################################
"""Methode de selection-mutation """
#########################################################

##fonction G et f:
#fonction G et f:
lamda= .01
M=int(1E4)
L=20
def G(x):
    return np.exp(lamda*x)

def G1(x,i):
    return G(x[i])


def f(x):
    return (x>0)

  
#r=np.arange(1,60)
#
#s=np.random.choice(a=r,size=100000,p=G(r)/np.sum(G(r)))
#
#c=np.random.choice(a=s,size=100000,p=G(s)/np.sum(G(s)))
#
#
#plt.hist(c,bins=np.size(r),normed=True)
def SM(M):
    cte=1
    X0=npr.randint(initmin,initmax,size=M)
    X=np.zeros((n+2,M),dtype=int)
    X[1]=X0
    for i in range(1,n+1):
    #    print("tourne : ",int(100*i/n),"%")
        cte*=np.mean(G1(X,i))
        I=np.random.choice(a=np.arange(M),size=M,p=G1(X,i)/np.sum(G1(X,i)))
        
        X=X[:,I]      # Selection
        
        z=np.random.geometric(1/(moy+1),size=(int(np.max(X[i])),M))-1
        
        for k in range(M):                              # Mutation
            X[i+1][k]=np.sum(z[:int(X[i][k]),k])        
    return X,cte

def selection_mutation(L):   
    pm=[]
    for l in range(L):
        X,cte=SM(M)
        esp1=np.mean(f(X[n+1])/   np.product(  G(X[:-1,:]) ,0 )    )*cte
        pm.append(esp1)        
    return pm
    
#start = timer()

pm=selection_mutation(L)    
#end = timer()
#print(end - start)

#print(pm)
#intervalle de confiance 
mean_pro=np.mean(pm)

print("\nSelection mutation: ",L,"fois avec",M,"simulations\n\nValeur:  {:.8f}".format(mean_pro))
s_l=sps.t.ppf(0.975,df=L-1)
b=s_l*np.std(pm)/np.sqrt(L)
print("Erreur relative: {:.2f} %".format(np.abs(mean_pro-ve)/ve*100))

print("\nIntervalle de confiance à 95%: [{:.8f}, "
     "{:.8f}]".format(mean_pro - b, mean_pro + b))
print("Largeur relative: {:.2f} %".format(b/mean_pro*100))
        

#########################################################
"""Methode Echantillonnage preferentiel"""
#########################################################

def P(x,y):
    temp=np.logical_and((x==y),(x==0))
    r= p**(x+temp)*(1-p)**y*sp.special.binom(x+temp+y-1,x+temp-1)
    r*=1-temp
    r+=temp
    return r


def h(X):
    return (X[:,-1]>0)

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

print("\n\nEchantillonnage preferentiel: ",m,"simulations \n\nValeur: {}".format(v))
print("Erreur relative: {:.2f} %".format(np.abs(v-ve)/ve*100))

#intervalle de confiance
print("\nIntervalle de confiance à 95% : [{:.8f},{:.8f}]".format(v-d ,v+d ))
print("Largeur relative: {:.2f} %".format(d/v*100))
print(" ")














