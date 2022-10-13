import numpy as np
import scipy.constants as sc
import sys

#Function for reading parameters from input
def ReadParams(filename):
    f = open(filename,"r")
    i=0
    m=[]
    line=f.readline()
    while line:
        sep = line.split()
        appo=float(sep[1])
        m.append(appo)
        i+1
        line=f.readline()
    f.close()
    return m

#Function for reading the cores
def ReadInput(filename):
    f = open(filename,"r")
    a=[]
    b=[]
    data=f.read().splitlines()
    for i in range(1,np.size(data)):
        line=data[i]
        sep = line.split()
        a.append(float(sep[1]))
        b.append(float(sep[2]))
    f.close()
    a1=np.array(a)
    b1=np.array(b)
    return a1,b1

#Function for reading flux from input
def ReadFlux(filename, Nbin):
    TextFile = filename
    f = open(TextFile,"r")
    i=0
    v=[]
    line=f.readline()
    while line:
        v.append(float(line))
        i+1
        line=f.readline()
    f.close()
    if(np.size(v)!=Nbin):
        print("The dimension of flux input file differs from the numbers of bins")
        sys.exit(1)
    return np.array(v)

def Gauss(e,sigma): 
    g=np.array(1/(np.sqrt(2*sc.pi)*sigma)*np.exp(-0.5*(e)**2/(sigma**2)))
    return g    

#Oscillation probability NO
def Osc_prob(x,d, O, T13, T12, M21, M3l):
    # L=d*10**3
    # y=x
    hbar=sc.hbar/sc.e
    L=d*10**3/(hbar*sc.c)
    y=x*10**6    
    if(O==True):
        M31=M3l
        M32=M31-M21
    else:
        M32=M3l
        M31=M32+M21

    # D21=1.27*M21*L/y
    # D31=1.27*M31*L/y
    # D32=1.27*M32*L/y
    D21=M21*L/y/4
    D31=M31*L/y/4
    D32=M32*L/y/4

    p=1 - 4*(1-T13)**2*T12*(1-T12)*np.sin(D21)**2 - 4*(1-T12)*T13*(1-T13)*np.sin(D31)**2 - 4*T12*T13*(1-T13)*np.sin(D32)**2
    return p

# Asimov flux without background
def Flux_teo(x, flux, ord, T13, T12, M21, M3l, dist, w, N, A, B,C): #x in MeV
    bin=x[2]-x[1]
    #Theorical flux of anti-nu
    F=flux
    #IBD cross section
    Ee=x-1.293
    m_e=0.511
    pe=np.sqrt(Ee**2-m_e**2)
    s=10**(-43)*pe*Ee*(x**(-0.07056+0.02018*np.log(x)-0.001953*np.log(x)**3))
    #Oscillation probability
    Y=np.zeros(np.size(F))
    if(np.size(dist)==1):
        Y=s*F*Osc_prob(x,dist, ord, T13, T12, M21, M3l)
    else:
        for i in range(0,np.size(dist)):
            FO=s*F*Osc_prob(x,dist[i], ord, T13, T12, M21, M3l)
            Y+=FO*w[i]  
    #Normalized flux
    tot=np.sum(Y)
    Y/=(tot*np.size(x)*bin)
    #Energy resolution
    out=[]
    Evis=np.array(x-0.8)
    deltaE=np.array(np.sqrt((A**2/Evis)+B**2+(C/Evis)**2)*Evis)
    for i in range(0,np.size(x)):
        E0=Evis[i]
        gi=Gauss(Evis-E0,deltaE)
        conv=np.sum(gi*Y)*np.size(x)*bin
        out.append(conv)
    yconv=np.array(out)
    yconv*=N/np.sum(yconv)
    return yconv

#Asimov flux with background
def Flux_teo_bg(x, flux, bkg, ord, T13, T12, M21, M3l, dist, w, N, A, B,C): #x in MeV
    return Flux_teo(x, flux, ord, T13, T12, M21, M3l, dist, w, N, A, B,C)+bkg

#Flux with statistical fluctuations
def Flux_Real(x, flux, bkg, ord, T13, T12, M21, M3l, dist, w, N, A, B, C):
    bin=x[2]-x[1]
    #Theorical flux of anti-nu
    F=flux
    #IBD cross section
    Ee=x-1.293
    m_e=0.511
    pe=np.sqrt(Ee**2-m_e**2)
    s=10**(-43)*pe*Ee*(x**(-0.07056+0.02018*np.log(x)-0.001953*np.log(x)**3))
    #Oscillation probability
    Y=np.zeros(np.size(F))
    if(np.size(dist)==1):
        Y=s*F*Osc_prob(x,dist, ord, T13, T12, M21, M3l)
    else:
        for i in range(0,np.size(dist)):
            FO=s*F*Osc_prob(x,dist[i], ord, T13, T12, M21, M3l)
            Y+=FO*w[i]  
    #Normalized flux
    Y*=N/np.sum(Y)
    #Poisson fluctuations
    new=np.array(np.random.poisson(Y),dtype='float')
    new/=(np.sum(new)*np.size(x)*bin)
    #Energy resolution
    out=[]
    Evis=np.array(x-0.8)
    deltaE=np.array(np.sqrt((A**2/Evis)+B**2+(C/Evis)**2)*Evis)
    for i in range(0,np.size(x)):
        E0=Evis[i]
        gi=Gauss(Evis-E0,deltaE)
        conv=np.sum(gi*new)*np.size(x)*bin
        out.append(conv)
    yconv=np.array(out)
    yconv*=N/np.sum(yconv)
    return yconv+bkg
