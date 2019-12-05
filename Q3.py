
                #######################################################################################
                ###                                                                                 ###  
                ###          """ Projet SNA : Dynamique de population, modele de Galton-Watson"""   ###  
                ###          """ Ismail Atourki && Ilyas Malik """                                  ###  
                ###                                                                                 ###  
                #######################################################################################
                
"""Cas sous critique """         

       
""" Q3: Scenario typique ou la population n'est pas eteinte """

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sps
from timeit import default_timer as timer

plt.close("all")

alpha=.05
n=int(40)
m=int(1e4)
moy=.8  # Le nombre moyen d'enfant par individu
initmin=20    # La loi de la population initiale suit une loi 
initmax=21    # uniforme dans [initmin,initmax]


#########################################################
""" Méthode de Monte-Carlo naif"""
#########################################################

def mat(mo):
    Y=[]
    for j in range(m):
#        if (int(4*j/m)==4*j/m):
#            print("tourne : ",int(100*j/m),"%")
            
        X=[np.random.randint(initmin,initmax)]
        for i in range (n):
            y=np.random.geometric(1/(mo+1),size=X[-1])-1
            X.append(np.sum(y))
        Y.append(X)
    Y=np.asarray(Y)
    return Y

Y=mat(moy)
        
    
"""Scenario typique ou la population n'est pas éteinte """  
l=np.arange(0,n+1,1)
max=10
bb=0
for k in range(m):
    if Y[k,-1]>0 and bb<max :
        plt.plot(l,Y[k,:])
        bb+=1


plt.title("Scenario ou la popultation n'est pas eteinte")


##########################################################
#"""Methode de selection-mutation """
##########################################################
#
###fonction G et f:
##fonction G et f:
#lamda=0.01
#M=int(1E3)
#L=20
#def G(X):
#    return np.exp(lamda*X)
#
#def G1(X,i):
#    return G(X[i])
#
#
#def f(X):
#    return (np.sum(X,axis=0)>x)
#
#  
##r=np.arange(1,60)
##
##s=np.random.choice(a=r,size=100000,p=G(r)/np.sum(G(r)))
##
##c=np.random.choice(a=s,size=100000,p=G(s)/np.sum(G(s)))
##
##
##plt.hist(c,bins=np.size(r),normed=True)
#def SM(M):
#    cte=1
#    X0=npr.randint(initmin,initmax,size=M)
#    X=np.zeros((n+2,M),dtype=int)
#    X[1]=X0
#    for i in range(1,n+1):
#    #    print("tourne : ",int(100*i/n),"%")
#        cte*=np.mean(G1(X,i))
#        I=np.random.choice(a=np.arange(M),size=M,p=G1(X,i)/np.sum(G1(X,i)))
#        
#        X=X[:,I]      # Selection
#        
#        z=np.random.geometric(1/(moy+1),size=(int(np.max(X[i])),M))-1
#        
#        for k in range(M):                              # Mutation
#            X[i+1][k]=np.sum(z[:int(X[i][k]),k])        
#    return X,cte
#
