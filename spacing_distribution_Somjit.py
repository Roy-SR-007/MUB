"""
@author: Somjit Roy

Simulating the Distribution of the Spacings/Distances of the Args of Eigen Values
in randomly generated unitary matrices distributed with Haar Measure.

haar_measure(n):= this routine takes into account the random generation of 
                  unitary matrices with dimension 'n', which are distributed 
                  haar measure.

density_spacing(s):= returns the value of the conjectured distribution of the
                     spacings of the arguments at a specified value, say 's'.
                     
eigen_args(n):= this routines makes a call to the function haar_measure() and 
                determines the 'n' arguments of the rigen values of randomly 
                generated unitary matrix over haar measure.

spacing_dist(dim,r):= this routine basically takes input for the choice of 
                      dimensions of the matrices, along with the # simulations.
                      Thereby, simulating the distribution of the 
                      spacings/consecutive difference of the arguments (theta) 
                      of the eigen values.
"""

import numpy as np
import math
import matplotlib.pyplot as plt

"""
Generating a random unitary matrix distributed with Haar Measure.
"""
def haar_measure(n):
# n:= dimension of the matrices.
    z1 = np.random.randn(n,n) + 1j*np.random.randn(n,n)
    z = (1/math.sqrt(2))*z1
    Q,R = np.linalg.qr(z)
    Lambda = np.diag(np.diag(R)/abs(np.diag(R)))
    return np.dot(Lambda,Q)
    
"""
Returning the value of the density (pdf) at a specified value of the distance, say s.
"""
def density_spacing(s):
# s:= the specified value at which the distributional value of the spacing is to 
#     be computed. s can be an array of inputs/values.
    res = []
    for i in range(0,len(s)):
        res.append((32/((math.pi)*(math.pi)))*(s[i]*s[i])*(math.exp(-4*s[i]*s[i]*(1/math.pi))))
    
    return res

"""
Determining the n arguments of the eigen values of a complex matrix.
"""
def eigen_args(n):
# n:= dimension of the matrices.
    A = haar_measure(n)
    w,v = np.linalg.eig(A)
    ev = w
    real = np.real(ev)
    imag = np.imag(ev)
    theta = []
    
    for i in range(real.size):
        theta.append(math.atan2(imag[i],real[i]))
    
    theta = np.array(theta)
    return theta

"""
Simulating the Spacing Distribution of the Args.
"""
def spacing_dist(dim,r):
# dim:= an array which takes in the choice of the dimensions of the matrix, we want to simulate.
# r:= an integer, representing the # simulations, we want to run.
    s_arr = []
    s = []
    theta_sorted = []
    k = 331
    
    for d in range(dim[0],(dim[len(dim)-1]+1)):
        for sim in range(1,r+1):
            theta_sorted = []
            theta_sorted.append(sorted(eigen_args(d)))
            
            for i in range(len(theta_sorted[0])-1):
                x = len(theta_sorted[0])/(2*3.142)
                x = x*(theta_sorted[0][i+1] - theta_sorted[0][i])
                s.append(x)
            
            s_arr.append(s)
            s = []
            
        s_arr = np.array(s_arr)
        plt.subplot(k)
        plt.tight_layout()
        plt.title('d = '+str(d), fontsize=8, color='black', loc='left', style='italic')
        plt.hist(s_arr.flatten(),edgecolor="k",density=True,bins=20,alpha=0.4,color="red")
        x = np.linspace(0,3,1000)
        y = density_spacing(x)
        plt.plot(x,y,'-',color='blue')
        plt.xlabel('Spacing/Distance')
        plt.ylabel('Density')
        k = k + 1
        s_arr = []

    plt.suptitle('Spacing (Normalized) Distribution of Arguments', y=1.02)
    

"""
A test case with dimensions varying in the set {3,4,5,6,...,11} 
and number of simulations r = 1000.
"""

dim = range(3,12)
r = 10000
spacing_dist(dim,r)


