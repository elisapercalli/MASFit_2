#Main code that takes from input both the anti-nu spectrum and all the backgrounds
#Input background order: geo-nu, accidentals, alpha, fast neutrons, Li9
#run with "source runMASFit.sh"

#Checking computational time
import time as tm
start_tm=tm.time()

from matplotlib import pyplot as plt
import numpy as np
from iminuit import Minuit
import sys
import MASFit_func_2 as fn
import seaborn as sn

#Reading parameters from file
InputTextFile = sys.argv[1]
InputFluxFile = sys.argv[2]
InputCoresFile= sys.argv[3]
InputBkgFile1 = sys.argv[4]
InputBkgFile2 = sys.argv[5]
InputBkgFile3 = sys.argv[6]
InputBkgFile4 = sys.argv[7]
InputBkgFile5 = sys.argv[8]
InputBkgFile6 = sys.argv[9]

#Energy and other parameters
v=fn.ReadParams(InputTextFile)
Dist,Pw=fn.ReadInput(InputCoresFile)
Weights=Pw/Dist**2

Fluc=v[0]
if (Fluc==True):
    M=int(v[1])
else:
    M=1
Nbin=int(v[2])
Emin=v[3]
Emax=v[4]
time=1800
Ntot=int(v[5])
#Anti-nu oscillations parameters
Theta12=v[6]
Theta13_NO=v[7]
Theta13_IO=v[8]
DeltaM21=v[9]
DeltaM31_NO=v[10]
DeltaM32_IO=v[11]

#Energy resolution and systematic uncertainties
a=v[12]/100
b=v[13]/100
c=v[14]/100
sigma_a=v[15]/100
sigma_b=v[16]/100
sigma_c=v[17]/100
sigma_alphaC=v[18]/100
sigma_alphaD=v[19]/100
sigma_b2b=v[20]/100
sigma_alphaR=v[21]/100
sigma_B_geo=v[22]/100
sigma_sh_geo=v[23]/100
sigma_B_acc=v[24]/100
sigma_sh_acc=v[25]/100
sigma_B_alpha=v[26]/100
sigma_sh_alpha=v[27]/100
sigma_B_fn=v[28]/100
sigma_sh_fn=v[29]/100
sigma_B_Li9=v[30]/100
sigma_sh_Li9=v[31]/100
sigma_B_WR=v[32]/100
sigma_sh_WR=v[33]/100
Scan=bool(v[34])
Corr=bool(v[35])
#Fit parameters
Fix_M21=bool(v[36])
Fix_M3l=bool(v[37])
Fix_T13=bool(v[38])
Fix_T12=bool(v[39])
Fix_N=bool(v[40])
Fix_a=bool(v[41])
Fix_b=bool(v[42])
Fix_c=bool(v[43])

#Rate for normalization
IBD_r=v[44]
Norm_geo=v[45]/IBD_r*Ntot
Norm_acc=v[46]/IBD_r*Ntot
Norm_Li9=v[47]/IBD_r*Ntot
Norm_fn=v[48]/IBD_r*Ntot
Norm_alpha=v[49]/IBD_r*Ntot
Norm_WR=v[50]/IBD_r*Ntot

#Check for real or ideal baseline
if(v[51]!=0):
    Dist=v[51]

#Reading anti-nu flux from file
F=fn.ReadFlux(InputFluxFile, Nbin)
Bg_geo=fn.ReadFlux(InputBkgFile1, Nbin)
Bg_acc=fn.ReadFlux(InputBkgFile2, Nbin)
Bg_alpha=fn.ReadFlux(InputBkgFile3, Nbin)
Bg_fn=fn.ReadFlux(InputBkgFile4, Nbin)
Bg_Li9=fn.ReadFlux(InputBkgFile5, Nbin)
Bg_WR=fn.ReadFlux(InputBkgFile6, Nbin)

#Normalizaion of backgrounds
Bg_geo*=(Norm_geo/np.sum(Bg_geo))
Bg_acc*=(Norm_acc/np.sum(Bg_acc))
Bg_alpha*=(Norm_alpha/np.sum(Bg_alpha))
Bg_fn*=(Norm_fn/np.sum(Bg_fn))
Bg_Li9*=(Norm_Li9/np.sum(Bg_Li9))
Bg_WR*=(Norm_WR/np.sum(Bg_WR))
Bg_tot=Bg_acc+Bg_alpha+Bg_fn+Bg_geo+Bg_Li9+Bg_WR

#Creation of anti-nu NO flux
bin=(Emax-Emin)/Nbin
E2=[]
for i in range(0,Nbin):
    E2.append(Emin+bin*(i+0.5))
e=np.array(E2)

xe=e-0.8
D_CHI=[]
bins=int(1+3*np.log(M))

for i in range(0,M):
    if(i%10==0 and i!=0):
        print(i)

    #without statistical fluctuations
    if(Fluc==False):
        yMis=fn.Flux_teo_bg(e, F, Bg_tot, 1, Theta13_NO, Theta12, DeltaM21, DeltaM31_NO, Dist, Weights, Ntot, a, b, c)  #NO
        #yMis=fn.Flux_teo_bg(e, F, Bg_tot, 0, Theta13_IO, Theta12, DeltaM21, DeltaM32_IO, Dist, Weights, Ntot, a, b, c) #IO

    #with statistical fluctuations
    else:
        yMis=fn.Flux_Real2(e, F, Bg_tot, 1, Theta13_NO, Theta12, DeltaM21, DeltaM31_NO, Dist, Weights, Ntot, a, b, c)    #NO      
        #yMis=fn.Flux_Real2(e, F, Bg_tot, 0, Theta13_IO, Theta12, DeltaM21, DeltaM32_IO, Dist, Weights, Ntot, a, b, c)    #IO

    #Computing of errors
    err=np.sqrt(yMis)
    for k in range(0,Nbin):
        if err[k]==0:
            err[k]=1

    #Definition of the Chi squared
    def least_squares_syst(T13, T12, M21, M3l, ord, N, A, B, C,alpha_C,alpha_D,a_R1,a_R2,a_R3,a_R4,a_R5,a_R6,a_R7,a_R8,a_R9,\
         alpha_B1,alpha_B2,alpha_B3,alpha_B4,alpha_B5,alpha_B6):
        yT = fn.Flux_teo(e, F, ord, T13, T12, M21, M3l, Dist, Weights, N, A, B, C)    #without background
        weig_sum=Weights[0]*a_R1+Weights[1]*a_R2+Weights[2]*a_R3+Weights[3]*a_R4+Weights[4]*a_R5+Weights[5]*a_R6+Weights[6]*a_R7+Weights[7]*a_R8+Weights[8]*a_R9
        z_num=(yMis-yT*(1+alpha_C+alpha_D+weig_sum)-Bg_geo*(1+alpha_B1)-Bg_acc*(1+alpha_B2)-Bg_alpha*(1+alpha_B3)-Bg_fn*(1+alpha_B4)-Bg_Li9*(1+alpha_B5)-Bg_WR*(1+alpha_B6))**2
        z_den=yMis+(sigma_b2b*yT)**2+(sigma_sh_geo*Bg_geo)**2+(sigma_sh_acc*Bg_acc)**2+(sigma_sh_alpha*Bg_alpha)**2+(sigma_sh_fn*Bg_fn)**2+(sigma_sh_Li9*Bg_Li9)**2+(sigma_sh_WR*Bg_WR)**2
        CHI=np.sum(z_num/z_den)
        pull_a=((A-a)/sigma_a)**2
        pull_b=((B-b)/sigma_b)**2
        pull_c=((C-c)/sigma_c)**2
        pull_syst=(alpha_C/sigma_alphaC)**2+(alpha_D/sigma_alphaD)**2
        pull_syst_cores=(a_R1**2+a_R2**2+a_R3**2+a_R4**2+a_R5**2+a_R6**2+a_R7**2+a_R8**2+6*a_R9**2)/sigma_alphaR**2
        pull_bkg=(alpha_B1/sigma_B_geo)**2+(alpha_B2/sigma_B_acc)**2+(alpha_B3/sigma_B_alpha)**2+(alpha_B4/sigma_B_fn)**2+(alpha_B5/sigma_B_Li9)**2+(alpha_B6/sigma_B_WR)**2
        return CHI+pull_a+pull_b+pull_c+pull_syst+pull_syst_cores+pull_bkg
    least_squares_syst.errordef = Minuit.LEAST_SQUARES

    #Fit NO on NO
    m = Minuit(least_squares_syst, T13=Theta13_NO,T12=Theta12,M21=DeltaM21,M3l=DeltaM31_NO, ord=1, N=Ntot,A=a,B=b,C=c,alpha_C=0,alpha_D=0,\
        a_R1=0, a_R2=0, a_R3=0, a_R4=0, a_R5=0, a_R6=0, a_R7=0, a_R8=0, a_R9=0, alpha_B1=0,alpha_B2=0,alpha_B3=0,alpha_B4=0,alpha_B5=0,alpha_B6=0)

    #Fit parameters
    m.fixed["ord"]=True
    m.fixed["M21"]=Fix_M21
    m.fixed["M3l"]=Fix_M3l
    m.fixed["T13"]=Fix_T13
    m.fixed["T12"]=Fix_T12
    m.fixed["N"]=Fix_N
    m.fixed["A"]=Fix_a
    m.fixed["B"]=Fix_b
    m.fixed["C"]=Fix_c
    m.migrad()
    m.hesse()

    #Fit IO su NO
    n = Minuit(least_squares_syst, T13=Theta13_IO,T12=Theta12,M21=DeltaM21,M3l=DeltaM32_IO, ord=0, N=Ntot,A=a,B=b,C=c,alpha_C=0,alpha_D=0, \
        a_R1=0, a_R2=0, a_R3=0, a_R4=0, a_R5=0, a_R6=0, a_R7=0, a_R8=0, a_R9=0, alpha_B1=0,alpha_B2=0,alpha_B3=0,alpha_B4=0,alpha_B5=0,alpha_B6=0)

    #Fit parameters
    n.fixed["ord"]=True
    n.fixed["M21"]=Fix_M21
    n.fixed["M3l"]=Fix_M3l
    n.fixed["T13"]=Fix_T13
    n.fixed["T12"]=Fix_T12
    n.fixed["N"]=Fix_N
    n.fixed["A"]=Fix_a
    n.fixed["B"]=Fix_b
    n.fixed["C"]=Fix_c
    n.migrad()
    n.hesse()

    D_CHI.append(n.fval-m.fval)

if(Fluc==False):
    #Print on file of fit reconstructed parameters
    file_par = open("MASFit_parameters.txt", "w")
    file_par.write("Free parameters in fit NO")
    file_par.write("\nParameter\tInj_value\tRec_value\tError\tBias(%)")
    if(Fix_M21==False):
        file_par.write(f"\n\nDeltaM21\t{DeltaM21:.3e}\t{m.values[2]:.3e}\t{m.errors[2]:.3e}\t{(m.values[2]-DeltaM21)/DeltaM21*100:.3f}")
    if(Fix_M3l==False):
        file_par.write(f"\nDeltaM31\t{DeltaM31_NO:.3e}\t{m.values[3]:.3e}\t{m.errors[3]:.3e}\t{(m.values[3]-DeltaM31_NO)/DeltaM31_NO*100:.3f}")
    if(Fix_T13==False):
        file_par.write(f"\nSin^2(T_13)\t{Theta13_NO}\t{m.values[0]:.5f}\t{m.errors[0]:.5f}\t{(m.values[0]-Theta13_NO)/Theta13_NO*100:.3f}")
    if(Fix_T12==False):
        file_par.write(f"\nSin^2(T_12)\t{Theta12}\t{m.values[1]:.3f}\t{m.errors[1]:.3f}\t{(m.values[1]-Theta12)/Theta12*100:.3f}")
    if(Fix_N==False):
        file_par.write(f"\nN\t{Ntot}\t{m.values[5]}\t{m.errors[5]:.0f}\t{(m.values[5]-Ntot)/Ntot*100:.3f}")
    if(Fix_a==False):
        file_par.write(f"\na\t{a:.4f}\t{m.values[6]:.4f}\t{m.errors[6]:.4f}\t{(m.values[6]-a)/a*100:.3f}")
    if(Fix_b==False):
        file_par.write(f"\nb\t{b:.4f}\t{m.values[7]:.4f}\t{m.errors[7]:.4f}\t{(m.values[7]-b)/b*100:.3f}")
    if(Fix_c==False):
        file_par.write(f"\nc\t{c:.4f}\t{m.values[8]:.4f}\t{m.errors[8]:.4f}\t{(m.values[8]-c)/c*100:.3f}")
    file_par.write(f"\nChi^2 = {m.fval:.2f}")

    file_par.write("\n\nFree parameters in fit IO")
    file_par.write("\nParameter\tInj_value\tRec_value\tError\tBias(%)")
    if(Fix_M21==False):
        file_par.write(f"\n\nDeltaM21\t{DeltaM21:.3e}\t{n.values[2]:.3e}\t{n.errors[2]:.3e}\t{(n.values[2]-DeltaM21)/DeltaM21*100:.3f}")
    if(Fix_M3l==False):
        file_par.write(f"\nDeltaM32\t{DeltaM32_IO:.3e}\t{n.values[3]:.3e}\t{n.errors[3]:.3e}\t{(n.values[3]-DeltaM32_IO)/DeltaM32_IO*100:.3f}")
    if(Fix_T13==False):
        file_par.write(f"\nSin^2(T_13)\t{Theta13_IO}\t{n.values[0]:.5f}\t{n.errors[0]:.5f}\t{(n.values[0]-Theta13_IO)/Theta13_IO*100:.3f}")
    if(Fix_T12==False):
        file_par.write(f"\nSin^2(T_12)\t{Theta12}\t{n.values[1]:.3f}\t{n.errors[1]:.3f}\t{(n.values[1]-Theta12)/Theta12*100:.3f}")
    if(Fix_N==False):
        file_par.write(f"\nN\t{Ntot}\t{n.values[5]}\t{n.errors[5]:.0f}\t{(n.values[5]-Ntot)/Ntot*100:.3f}")
    if(Fix_a==False):
        file_par.write(f"\na\t{a:.4f}\t{n.values[6]:.4f}\t{n.errors[6]:.4f}\t{(n.values[6]-a)/a*100:.3f}")
    if(Fix_b==False):
        file_par.write(f"\nb\t{b:.4f}\t{n.values[7]:.4f}\t{n.errors[7]:.4f}\t{(n.values[7]-b)/b*100:.3f}")
    if(Fix_c==False):
        file_par.write(f"\nc\t{c:.4f}\t{n.values[8]:.4f}\t{n.errors[8]:.4f}\t{(n.values[8]-c)/c*100:.3f}")
    file_par.write(f"\nChi^2 = {n.fval:.2f}")
    file_par.write(f"\n\nDeltaChi^2 = {n.fval-m.fval:.2f}")
    file_par.close()

    #Plot parameters scan
    if(Scan==True):
        fig_Z3,ax_Z3 = plt.subplots(1,1)
        xm,ym,rm = m.mnprofile('M3l',bound=5, subtract_min=False)
        xn,yn,rn = n.mnprofile('M3l',bound=5, subtract_min=False)
        ax_Z3.plot(xm, ym, label="NO", markersize=2, color='blue')
        ax_Z3.plot(np.abs(xn), yn, label="IO", markersize=2, color='red')
        leg_Z3 = ax_Z3.legend();
        ax_Z3.grid()
        plt.title(f"Scan $\\chi^2$")
        plt.tight_layout()
        plt.xlabel("$|\\Delta m^2_3l| [10^{-3} eV^2]$")
        plt.ylabel("$\\chi^2$")
        fig_Z3.set_size_inches(10, 8)
        plt.savefig("MASFit_scan_2.pdf")

    #Correlation matrix
    if (Corr==True):
        plt.figure(figsize=(14,10))
        lab=["T12","M21","M3l","A","B","C","a_C","a_D","a_R1","a_R2","a_R3","a_R4","a_R5","a_R6","a_R7","a_R8","a_R9","a_B1","a_B2","a_B3","a_B4","a_B5","a_B6"]
        corr1 = m.covariance.correlation()
        corr2=np.delete(corr1,(0,4,5),axis=0)
        corr3=np.delete(corr2,(0,4,5),axis=1)
        sn.heatmap(corr3, xticklabels=lab, yticklabels=lab,fmt=".2f", cmap="BrBG",linewidths=.5,annot=True)
        plt.savefig("corr_matrix.png")

    #Plot
    plt.figure(figsize=(12, 8))
    plt.title("Fit assuming NO",fontsize=16)
    plt.errorbar(xe, yMis, yerr=err, markersize=3, fmt="o", label="data", alpha=0.4)
    plt.plot(xe, fn.Flux_teo_bg(e,F,Bg_tot, 1,m.values[0],m.values[1],m.values[2],m.values[3], Dist, Weights,\
        m.values[5],m.values[6],m.values[7],m.values[8]), label="fit NO", color="seagreen")
    plt.plot(xe, fn.Flux_teo_bg(e,F,Bg_tot, 0, n.values[0],n.values[1],n.values[2],n.values[3], Dist, Weights,\
        n.values[5],n.values[6],n.values[7],n.values[8]), label="fit IO", color="coral")
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)

    # Legend and fit info
    CHI_NO=m.fval/(len(e)-m.nfit)
    CHI_IO=n.fval/(len(e)-n.nfit)
    DeltaCHI=n.fval-m.fval
    plt.legend(title=f"$\\Delta$ $\\chi^2$ = {DeltaCHI:.3f}", title_fontsize=14,fontsize=14)
    plt.rc('legend', fontsize=16)
    plt.grid()
    plt.xlabel("E_vis [MeV]",fontsize=16)
    plt.ylabel("dN/dE [#/ 20 keV]",fontsize=16)
    plt.savefig("MASFit_plot.pdf")

else:
    #checking for outliers and deleting them
    mean_CHI=np.mean(D_CHI)
    med_CHI=stat.median(D_CHI)
    std_CHI=np.std(D_CHI)

    i=0
    j=0
    while (i<(M-j)):
        if (D_CHI[i]>med_CHI+5*std_CHI or D_CHI[i]<med_CHI-5*std_CHI):
            print("out",D_CHI[i])
            D_CHI.pop(i)
            med_CHI=stat.median(D_CHI)
            std_CHI=np.std(D_CHI)
            i=-1
            j+=1
        i+=1

    print(stat.median(D_CHI), np.mean(D_CHI),np.std(D_CHI))
    file_CHI = open("CHI_NO_allN.txt", "w")
    #file_CHI = open("CHI_IO_allN.txt", "w")
    for i in range(0,np.size(D_CHI)):
        file_CHI.write(f"\n{D_CHI[i]}")
    file_CHI.close()
    
    #Plot Chi squared
    plt.figure(figsize=(12, 8))
    plt.title("Fit su NO")
    plt.hist(D_CHI,bins)
    plt.grid()
    plt.xlabel(f"$\\Delta$ $\\chi^2$")
    plt.ylabel("N events")
    plt.savefig("MASFit_CHI.pdf")

print("Computational time: ",f"{(tm.time()-start_tm):.2f}"," s")
if(v[52]==True):
    plt.show()
