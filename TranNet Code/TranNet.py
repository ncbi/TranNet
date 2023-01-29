import math
import numpy as np
import pandas as pd
from sklearn.linear_model import Lasso
def Fun(X, y, w): # Compute the objective function value and return the fitness error as a positive scalar
    return np.linalg.norm(X.dot(w)-y)**2

def Grad(X, y, w): # Compute the gradient of the objective function and return a vector of the descent direction
    return 2*X.T.dot(X.dot(w)-y)

def Alpha(X, y, w, g): # Bisection search for choosing optimal size of step and return a positive scalar value
    a=0; b=1; gr = (math.sqrt(5) + 1) / 2
    c=b-(b-a)/gr; d =a+(b-a)/gr
    while abs(b - a) > 1e-15:
        if Fun(X, y, w-c*g) < Fun(X, y, w-d*g):
            b=d
        else:
            a=c
        c=b-(b-a)/gr; d = a + (b - a) / gr
    return (b + a) / 2

def Project(v, up): # Project a vector (current solution) onto l_1 norm space, and return a unit-norm vector
    n, = v.shape
    s = np.abs(v)
    if s.sum() <= up:
        return v
    else:
        u = np.sort(s)[::-1]
        cSum = np.cumsum(u)
        rho = np.nonzero(u * np.arange(1, n+1) > (cSum - up))[0][-1]
        theta = float(cSum[rho] - up) / (rho+1)
        x = (s - theta).clip(min=0)
        x *= np.sign(v)
    return x
def InitialSol(X, y): # Initialize starting point with LASSO regression and return a initial point close to optimal
    dense_lasso = Lasso(alpha=0.1, fit_intercept=False, max_iter=100, tol=1e-5)
    dense_lasso.fit(X, y)
    return dense_lasso.coef_

# Gradient projection method for find optimal transition matrix over l_1 norm space and return transition matrix between given two matrix
def ProjectedGradient(RESPON, EXPLAN):
    W=np.zeros((len(RESPON.columns), len(EXPLAN.columns)))
    X=np.array(EXPLAN); Y=np.array(RESPON); Error=0
    for i in range(len(RESPON.columns)):
        W[i,:]=InitialSol(X, Y[:,i]) # Use LASSO to obtain an initial solution
        k=0; g=Grad(X, Y[:,i], W[i,:])
        FunOld=10;
        W[i,:]=Project(W[i,:]-Alpha(X, Y[:,i], W[i,:], g)*g, 1)
        FunNew=Fun(X, Y[:,i], W[i,:])
        while abs(FunOld-FunNew) > 1e-10 and k<100:
            FunOld=FunNew
            g=Grad(X, Y[:,i], W[i,:])
            xk=W[i,:]-Alpha(X, Y[:,i], W[i,:], g)*g
            W[i,:]=Project(xk, 1)
            FunNew=Fun(X, Y[:,i], W[i,:])
            k=k+1
        Error=Error+FunNew
    return pd.DataFrame(W,columns=EXPLAN.columns, index=RESPON.columns)
nEXP=pd.read_csv('NormalData.csv', index_col=0) #Load Normal expression data
tEXP=pd.read_csv('TumorData.csv', index_col=0)  #Load Tumor expression data
M=ProjectedGradient(tEXP, nEXP)                 #Find optimal transition matrix
(M.T).to_csv('Transition_Matrix_M.csv')
RP=(abs(M)).sum(axis=0).sort_values(ascending=False)  #Compute scores of regulatory potentials for genes
RP.to_frame(name='RP Score').to_csv('Regulatory_Potential.csv')
