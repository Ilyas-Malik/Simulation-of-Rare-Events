
                
"""Cas sous critique """         

       
""" Q2: la loi de X_n conditionnellement a la survie A_n"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sps
import scipy as sp
import numpy.random as npr
plt.close("all")

n=int(20)
m=int(1e4)
moy=.8        # Le nombre moyen d'enfant par individu
p=1/(moy+1)
q=sps.norm.ppf(.975)
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
            
        X=[npr.randint(initmin,initmax)]
        for i in range (n):
            y=np.random.geometric(1/(mo+1),size=X[-1])-1
            X.append(np.sum(y))
        Y.append(X)
    Y=np.asarray(Y)
    return Y

Y=mat(moy)
Z=Y[:,-1]
Z=Z[Z!=0]

plt.figure()
plt.hist(Z,bins=int(m**(1./3)),normed= True,label="monte carlo naive sachant non extinction")

loc,scale=sps.expon.fit(Z,floc=0)

x=np.linspace(0,6*scale,100)
expo=sps.expon.pdf(x,scale=scale)
plt.plot(x,expo,"r",label="loi expo de moyenne={}".format(round(scale,2)))

plt.legend(loc="best")



#########################################################
"""Methode de selection-mutation """
#########################################################

##fonction G et f:
#fonction G et f:
lamda= 0.01
#lamda=0.1
def G(x):
    return np.exp(lamda*x)

def G1(x,i):
    return G(x[i])


def f(x):
    return (x>0)

## histogramme de la fonction poids G
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


plt.figure()
X=SM(m)[0]
Z=X[-1,:]
Z=Z[Z!=0]

plt.hist(Z,bins=int(m**(1./3)),normed= True,label="selection mutation sachant non extinction")

loc,scale=sps.expon.fit(Z,floc=0)


x=np.linspace(0,6*scale,100)
expo=sps.expon.pdf(x,scale=scale)

lgauss,sgauss=sps.norm.fit_loc_scale(Z)
plt.plot(x,sps.norm.pdf(x,loc=lgauss,scale=sgauss),"g",label="Gaussienne de m={} sigma={}".format(round(lgauss,2),round(sgauss,2)))


plt.plot(x,expo,"r",label="loi expo de moyenne={}".format(round(scale,2)))
plt.legend(loc="best")
    


















