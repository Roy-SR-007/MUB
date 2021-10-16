import numpy as np
import matplotlib.pyplot as plt

q = int(input('Enter the Order of Affine Plane'))

#m = int(input('Enter the Power of the Prime Number'))

#q = p**m

t = np.arange(q)
t2 = t**2

#mt = np.arange(q+1)

Y = np.zeros((q,q,q), dtype=np.int64)
Blk = np.zeros((q+1,q,q), dtype=np.int64)
MUBs = np.zeros((q+1,q**2,q**2), dtype=np.complex128)


for i in t:
    for k in t:
        Y[i,k,:]= (i*t +k)%q


#print(Y)                  
        

for i in t:
    for k in t:
        Blk[i,k,:]= Y[i,k,:]*q + t

Blk[q,:,:] = np.transpose(Blk[0,:,:])

print(Blk)

dftmtx = np.fft.fft(np.eye(q))/q

for i in t:
    for k in t:
        MUBs[i,Blk[i,k,:],q*k:q*(k+1)]= dftmtx[:,:]

for k in t:
    MUBs[q,Blk[q,k,:],q*k:q*(k+1)]= dftmtx[:,:]

#print(MUBs)


MUB1 = MUBs[1,:,:]
print(MUB1)
TMUB1 = np.transpose(MUB1)
TMUB = np.conjugate(TMUB1)

r = 3
MUB2 = MUBs[r,:,:]
print(MUB2)
#TMUB2 = np.transpose(MUB2)
#TMUB = np.conjugate(TMUB2)

print(np.abs(TMUB @ MUB2))

