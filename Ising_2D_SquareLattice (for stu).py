
# Ising_model_in_two_dimensions

import numpy as np
import matplotlib.pyplot as plt
Ls=[4,8,16,32]
for pp in Ls:
    
    L = pp
    J = 1
    s = np.zeros((L,L))
    
    def total_E(s,J,L):
        E = 0
        for i in range (L):
            for j in range (L):
                E -=s[i][j]*(s[(i+1)%L][j]+s[i-1][j]+s[i][(j+1)%L]+s[i][j-1])
        return 0.5*J*E
    
    def metropolis(s,J,L):
        monteCarlo_Steps = 100
        
        for m in range(monteCarlo_Steps):
    
            for i in range(L):
                for j in range(L):
                    ii = np.random.randint(0,L)
                    jj = np.random.randint(0,L)
                    p = np.random.rand()
                    
                    delE = 2*J*s[ii][jj]*(s[(ii+1)%L][jj]+s[ii-1][jj]+s[ii][(jj+1)%L]+s[ii][jj-1])
                    if(delE <= 0):
                        s[ii][jj] = -s[ii][jj]
                    elif(p<np.exp(-beta*delE)):
                        s[ii][jj] = -s[ii][jj]
                    
        
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
                for j in range(L):
                    ii = np.random.randint(0,L)
                    jj = np.random.randint(0,L)
                    p = np.random.rand()
                    
                    delE = 2*J*s[ii][jj]*(s[(ii+1)%L][jj]+s[ii-1][jj]+s[ii][(jj+1)%L]+s[ii][jj-1])
                    if(delE <= 0):
                        s[ii][jj] = -s[ii][jj]
                        E += delE
                    elif(p<np.exp(-beta*delE)):
                        s[ii][jj] = -s[ii][jj]
                        E+= delE
            Ebar += E
            E2bar += E**2
        return m/cnt,Ebar/(cnt*L*L),(m2/cnt - (m/cnt)**2)/(T)
    
    Temp = np.arange(0.1,7,0.1)
    Ebar = []
    m = []
    chi = []
    cv = []
    for T in Temp:
        beta = 1.0/T
        for i in range(L):
            for j in range(L):
                s[i][j] = np.random.choice([1])
                
        m.append(magnetization(s, T, L)[0])
        #Ebar.append(magnetization(s, T, L)[1])
        #chi.append(magnetization(s, T, L)[2])
    
    #plt.scatter(Temp,m,label=str(pp)) #scatter plot
    plt.plot(Temp,m,label='L='+str(pp))
plt.xlabel('Temperature')
plt.ylabel('Magnetization')
plt.legend()

