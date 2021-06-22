
import numpy as np
import matplotlib.pyplot as plt


L = int(8)
J = 1
s = np.zeros((L,L,L))

def total_E(s,J,L):
    E=0
    for i in range(L):
        for j in range(L):
            for k in range(L):
                E = E + (s[i,j,k] * (s[(i+1)%L,j,k] + s[i,(j+1)%L,k] + s[i,j,(k+1)%L])) 
    return  -J*E 

def metropolis(s,J,L):
    monteCarlo_Steps = 10
    
    for m in range(monteCarlo_Steps):

        for i in range(L):
            for j in range(L):
                for k in range(L):
                    ii = np.random.randint(0,L)
                    jj = np.random.randint(0,L)
                    kk = np.random.randint(0,L)
                    
                    p = np.random.rand()
                    
                    e_o=s[ii,jj,kk]*(s[(ii+1)%L,jj,kk]+s[(ii-1)%L,jj,kk]+
                    s[ii,(jj-1)%L,kk]+s[ii,(jj+1)%L,kk]+
                    s[ii,jj,(kk-1)%L]+s[ii,jj,(kk+1)%L])
                    delE = 2*J*e_o
                    if(delE <= 0):
                        s[ii][jj][kk] = -s[ii][jj][kk]
                    elif(p<np.exp(-beta*delE)):
                        s[ii][jj][kk] = -s[ii][jj][kk]
                    
        
    return s

def magnetization(s,T,L):
    metropolis(s, J, L)
    
    m = 0
    m2 = 0
    E = total_E(s, J, L)
    Ebar = 0
    E2bar = 0
    Ebar += E
    E2bar += E**2
    cnt = 0
    for n in range(100):
        
        m += np.abs(s.mean())
        m2 += (s.mean())**2
        cnt +=1
        
        for i in range(L):
            for K in range(L):
                for j in range(L):
                    ii = np.random.randint(0,L)
                    jj = np.random.randint(0,L)
                    kk = np.random.randint(0,L)
                    p = np.random.rand()
    
                    e_o=s[ii,jj,kk]*(s[(ii+1)%L,jj,kk]+s[(ii-1)%L,jj,kk]+
                    s[ii,(jj-1)%L,kk]+s[ii,(jj+1)%L,kk]+
                    s[ii,jj,(kk-1)%L]+s[ii,jj,(kk+1)%L])
                    delE = 2*J*e_o                
                    
                    if(delE <= 0):
                        s[ii][jj][kk] = -s[ii][jj][kk]
                        E += delE
                    elif(p<np.exp(-beta*delE)):
                        s[ii][jj][kk] = -s[ii][jj][kk]
                        E+= delE
        Ebar += E
        E2bar += E**2
    return m/cnt,Ebar/(cnt*L*L),(m2/cnt - (m/cnt)**2)/(T)

Temp = np.arange(0.1,5,0.5)
Ebar = []
m = []
chi = []
cv = []
for T in Temp:
    beta = 1.0/T
    for i in range(L):
        for k in range(L):
            for j in range(L):
                s[i][j][k] = np.random.choice([1])
                
    m.append(magnetization(s, T, L)[0])
    Ebar.append(magnetization(s, T, L)[1])
    chi.append(magnetization(s, T, L)[2])

plt.plot(Temp,chi)
plt.xlabel('Temp')
plt.ylabel('chi')
