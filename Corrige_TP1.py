import numpy as np
import numpy.linalg as alg
import matplotlib.pyplot as plt

#Exercice 1
#Question 1
def Vandermonde(x) :
    N = np.size(x)
    Vander = np.zeros((N,N))
    for j in range (N) :
        for i in range (N) :
            Vander[i,j]= (x[i])**j
    return Vander

#On peut aussi le faire sans boucle pour économiser de la mémoire et temps de calcul :

def matrice_vandermonde_mieux(x):
    X = np.array(x)
    N = len(X)
    V = X[:, np.newaxis] ** np.arange(N)
    return V

#Question 2
def solvand(x,y) :
    N = np.size(x)
    NN = np.size(y)
        if N != NN :
       print("Probleme de taille entre les vecteurs")
    Vander = Vandermonde(x)
    a = alg.solve(Vander, y)
    return a

#Question 3
x=np.array([0, 1, 2, 3])
y=np.array([1, 2, 9, 28])

rep=solvand(x,y)
print(rep)

#Question 4
