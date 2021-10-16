import numpy as np
import matplotlib.pyplot as plt

p = int(input('Enter the Prime Number'))

m = int(input('Enter the Power of the Prime Number'))

q = p**m

t = np.arange(q)
t2 = t**2

mt = np.arange(q+1)

Phs = np.zeros((q,q,q), dtype=np.complex64)
MUBs = np.zeros((q,q,q), dtype=np.complex64)


for i in t:
    for k in t:
        Phs[i,k,:]=( 2*np.pi*1j/q )*(i*t2 + k*t)

MUBs = np.exp(Phs)
#ABS_MUBs = np.abs(MUBs)

#print(MUBs)
#print(ABS_MUBs)

MUB1 = MUBs[3,:,:]
#TMUB1 = np.transpose(MUB1)
#TMUB = np.conjugate(TMUB1)

r = 4
MUB2 = MUBs[r,:,:]
TMUB2 = np.transpose(MUB2)
TMUB = np.conjugate(TMUB2)

print(np.abs(MUB1@TMUB))
#print(np.transpose(MUBs[1,:,:]))
#print(np.exp(MUBs))
#Tp = MUBs[1,:,:]@np.transpose(MUBs[1,:,:]
#print(Tp)


