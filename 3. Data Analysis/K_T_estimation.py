import functions as Funk
import strain_gauge as DMSM

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import scipy.integrate as integrate
import scipy.odr as odr



############# Constant definition
c = 299792458.0 #m/s
Lambda = 1310.0  #1550.0# * 10**-9 #m  1310

T0 = 293

#Luna Parameter
K_S_C = 0.78 
K_S_M = 0.81
K_S_M_std = K_S_M*0.01
K_T_C =   6.45*10**(-6) #1/C
k_S_C = -6.668 #-Lambda/K_S_C/c *10**6  # -6.67#muStrain/GHz # used are -6.668
k_T_C = Lambda/K_T_C/c # -0.801 #C/GHz


c = 299792458.0 #m/s
Lambda = 1310.0
k_S_C = -6.668 #-Lambda/K_S_C/c *10**6  # -6.67#muStrain/GHz # used are -6.668
St_to_dvv = Lambda/c/k_S_C  # [10**9]  ## Umrechnung von Strain in norm. Shift




K_e_26 = [0.8308, 0.8301 ]
K_e_26_err = [0.0038,0.0037]
K_e_LN2 = [0.7816, 0.7893 ]
K_e_LN2_err = [0.0009,0.0032]


KE_G_std = 0.0008048738306906666
KE_B_std =  0.002304131801250109

DT = 30

def K_e(T,GoB,  Terr = 0,IfM = False, IfKerr = False):  
    
    if GoB ==  0:
        m =  0.00017223390701612612
        m_std =  8.210704853122784e-06
        m_sys =  1.0261889194373178e-06
        c =  0.7774684204746407
        c_std =  0.0017515057764383268
        c_sys =  0.004683672698285901
        cov =  [[ 1.17029148e-06 ,-4.68630234e-09], [-4.68630234e-09 , 2.57176793e-11]]
        Chi =  2.6213747153650107
        
        K = m * T + c
        if IfKerr == True:
            Kerr_stat = m *Terr        
            Kerr_sys = np.sqrt( np.abs(Chi*(T**2 *cov[1][1]**2 + cov[0][0]**2 + T*(cov[0][1] + cov[1][0])))+ (T*m_sys)**2 + c_sys**2  )
        
    if GoB ==  1:
        m =  0.00016388412854575812
        m_std =  1.0667120580183699e-05
        m_sys =  9.891838347601725e-07
        c =  0.7788708491034005
        c_std =  0.002612237670551766
        c_sys =  0.004688057943942114
        cov =  [[ 1.17029148e-06, -4.68630234e-09],  [-4.68630234e-09,  2.57176793e-11]]
        Chi =  3.1360342372149033
                
        K = m * T + c
        if IfKerr == True:
            Kerr_stat = m *Terr
            Kerr_sys = np.sqrt( np.abs(Chi*(T**2 *cov[1][1]**2 + cov[0][0]**2 + T*(cov[0][1] + cov[1][0]))) + (T*m_sys)**2 + c_sys**2  )
    if IfM == True:
        if IfKerr == False:
            return m
        if IfKerr == True:
            return np.sqrt(m_std **2 + m_sys**2) 
        
    
    if IfKerr == True:
        return Kerr_stat, Kerr_sys
    else: 
        return K
       
def PreanalysedData(Name):
    Data = np.loadtxt(Name)
    dv_G_M = Data[0] 
    dv_G_M_std = Data[1]
    dv_B_M = Data[2] 
    dv_B_M_std = Data[3]
    T_M = Data[4]
    T_M_std = Data[5]
    DMS_M = Data[6]
    DMS_M_std = Data[7]
    DMS_MK_sys = Data[8]
    T_MK_sys = np.zeros(len(T_M))+0.1
    
    delete_list=[]
    for i in delete_list:
        T_M = np.delete(T_M, i)
        T_M_std = np.delete(T_M_std, i)
        dv_G_M = np.delete(dv_G_M, i)
        dv_G_M_std = np.delete(dv_G_M_std, i)
        dv_B_M = np.delete(dv_B_M, i)
        dv_B_M_std = np.delete(dv_B_M_std, i)
        DMS_M = np.delete(DMS_M, i)
        DMS_M_std = np.delete(DMS_M_std, i)
        DMS_MK_sys =np.delete(DMS_MK_sys, i)
        T_MK_sys = np.delete(T_MK_sys, i)
        
        if len(Data) >10:
            DMS_MT = np.delete(DMS_MT, i)
            DMS_MT_std = np.delete(DMS_MT_std, i)
    
    return dv_G_M, dv_G_M_std, dv_B_M, dv_B_M_std, T_M, T_M_std, DMS_M, DMS_M_std, DMS_MK_sys, T_MK_sys


dv_G_MK, dv_G_MK_std, dv_B_MK, dv_B_MK_std, T_MK, T_MK_std, DMS_MK, DMS_MK_std ,DMS_MK_sys, T_MK_sys = PreanalysedData("step_test_time_dv_temp_strain_with_error_bor_and_ger.out")

Zusatz_Name = "GB_IQR_und_BStufe"

T_MK_Ref = T_MK[0]
T_MK_Ref_std = T_MK_std[0]
T_MK_Ref_sys = T_MK_sys[0]


def Sig_Pol4(B,Bs,x):
    b1 = B[0] + (B[1] + Bs[1]) * x + B[2]*x**2 + B[3]*x**3  + B[4]*x**4
    b2 = B[0] + (B[1] - Bs[1]) * x + B[2]*x**2 + B[3]*x**3  + B[4]*x**4
    
    c1 = B[0] + B[1] * x + (B[2] + Bs[2])*x**2 + B[3]*x**3  + B[4]*x**4
    c2 = B[0] + B[1] * x + (B[2] - Bs[2])*x**2 + B[3]*x**3  + B[4]*x**4
    
    d1 = B[0] + B[1] * x + B[2]*x**2 + (B[3] + Bs[3])*x**3  + B[4]*x**4
    d2 = B[0] + B[1] * x + B[2]*x**2 + (B[3] - Bs[3])*x**3  + B[4]*x**4
    
    e1 = B[0] + B[1] * x + B[2]*x**2 + B[3]*x**3  + (B[4] + Bs[4])*x**4
    e2 = B[0] + B[1] * x + B[2]*x**2 + B[3]*x**3  + (B[4] - Bs[4])*x**4
    
    return np.sqrt(  Bs[0]**2 + (np.abs(b1)-np.abs(b2))**2 /4  + (np.abs(c1)-np.abs(c2))**2 /4  + (np.abs(d1) - np.abs(d2))**2 /4  + (np.abs(e1) - np.abs(e2))**2 /4 )

def K_T(T,FitL,Fit_cor = None, Fit_syserr = None):
    B = FitL[1]
    C = FitL[2]
    D = [0,0]
    if len(FitL) >= 4:
        D = FitL[3]
    E = [0,0]
    if len(FitL) == 5:
        E = FitL[4]
    K = B[0]  +  2*C[0]*T + 3*D[0]*T**2 + 4*E[0]*T**3
    
    if Fit_cor != None:
        Kerr = KT_Error(T,Fit_cor) 
        if Fit_syserr != None:
            Kerr = KT_Error(T,Fit_cor,Fitsyserr = Fit_syserr ) 
            
    else:
        Kerr = Sig_Pol4([B[0],C[0],D[0],E[0],0],[B[1],C[1],D[1],E[1],0],T)
    return [K, Kerr]

def KT_Error(x, fitobj, Fitsyserr = None):
    confprob = 68.2
    p = fitobj.beta
    covscale = fitobj.sum_square/(len(x) - len(p))
    dof = len(p)
    
    if dof == 3:       
        dfdp = [np.zeros(len(x)), np.zeros(len(x)) + 1 , x*2]
        covscale = fitobj.sum_square/(len(fitobj.delta) - len(p) - 1)
        dof = len(p) - 1
        
    if dof == 4:       
        dfdp = [np.zeros(len(x)), np.zeros(len(x)) + 1 , x*2, 3*x**2]
        covscale = fitobj.sum_square/(len(fitobj.delta) - len(p) - 1)
        dof = len(p) - 1
        
    if dof == 5:       
        dfdp = [np.zeros(len(x)), np.zeros(len(x)) + 1 , x*2, 3*x**2, 4*x**3]
        covscale = fitobj.sum_square/(len(fitobj.delta) - len(p) - 1)
        dof = len(p) - 1

    abswei=False
    from scipy.stats import t
   # Given the confidence probability confprob = 100(1-alpha)
   # we derive for alpha: alpha = 1 - confprob/100 
    alpha = 1 - confprob/100.0
    prb = 1.0 - alpha/2
    tval = t.ppf(prb, dof)
   
    C = fitobj.cov_beta
    n = len(fitobj.beta) 
             # Number of parameters from covariance matrix
    
    N = len(x)
    if abswei:
       covscale = 1.0
    else:
       covscale = covscale
      
      
    df2 = np.zeros(N)
    for j in range(n):
       for k in range(n):
          df2 += dfdp[j]*dfdp[k]*C[j,k]
    df = np.sqrt(covscale*df2)
    delta = tval * df   
   
    if Fitsyserr != None:
        delta = np.sqrt(delta**2 + (Fitsyserr[0] + Fitsyserr[1]*x + Fitsyserr[2]*x**2)**2)
        if len(Fitsyserr) >3:
            delta = np.sqrt(delta**2 + (Fitsyserr[3] + Fitsyserr[4]*x + Fitsyserr[5]*x**2)**2)
  
    
    return delta

def K_T_2F(T,FitL):
    B = FitL[1]
    C = FitL[2]
    TC = FitL[-2]
    m = FitL[-1]   
    D = [0,0]
    if len(FitL) == 6:
        D = FitL[3]
        
    K = []
    Kerr = []
    for i in range(len(T)):
        if T[i] < TC[0]:
            
            K.append(0)
            Kerr.append(m[1])
            
        if T[i] >= TC[0]:
            K.append( B[0]  +  2*C[0]*T[i] + 3*D[0]*T[i]**2 )
            Kerr.append(Sig_Pol4([B[0],C[0],D[0],0,0],[B[1],C[1],D[1],0,0],T[i]))
            
    
    return [np.array(K), np.array(Kerr)]



def K_TN(T,FitL):
    A = FitL[0]
    B = FitL[1]
    C = FitL[2]

    return [A[0]*B[0]*(T - C[0])**(B[0]-1) , np.sqrt( (A[1]*B[0]*(T - C[0])**(B[0]-1))**2 )]

def F_T(T,T_0, FitL,Fit_cor = None,GoB = None,CombiKT = None, LinFit = None):
    B = FitL[1] 
    C = FitL[2] 
    D = [0,0]
    if len(FitL) >= 4:
        D = FitL[3]
    E = [0,0]
    if len(FitL) == 5:
        E = FitL[4]
    if Fit_cor != None:
        Ferr = Funk.Spesific_func_confidence_band(T,68.2,Fit_cor,8 + len(Fit_cor.beta))
        print(8 + len(Fit_cor.beta))
    else:
        Ferr = np.sqrt( (B[1]*(T-T_0) )**2 + ((T**2 - T_0**2)*C[1])**2 + ((T**3 - T_0**3)*D[1])**2 + ((T**4 - T_0**4)*E[1])**2) 
        
    F = ( B[0] * (T - T_0) + C[0] * (T**2-T_0**2) + D[0] * (T**3-T_0**3) + E[0]* (T**4 - T_0**4))
    
    
    Tc = 250
    if LinFit != None:
        BL = LinFit[1]
        CL = LinFit[2]
        DL = [0,0]
        if len(FitL) == 4:
            DL = LinFit[3]
        EL = [0,0]
        if len(FitL) == 5:
            EL = LinFit[4]
        
        if T_0 < Tc:    
            for i in range(len(T)):
                if T[i] > Tc:
                    F[i] = ( B[0] * (Tc - T_0) + C[0] * (Tc**2-T_0**2) + D[0] * (Tc**3-T_0**3) + E[0]* (Tc**4 - T_0**4)) + ( BL[0] * (T[i] - Tc) + CL[0] * (T[i]**2-Tc**2) + DL[0] * (T[i]**3-Tc**3) + EL[0]* (T[i]**4 - Tc**4))
        if T_0 >=Tc:
            for i in range(len(T)):
                if T[i] < Tc:
                    F[i] =  ( B[0] * (Tc - T_0) + C[0] * (Tc**2-T_0**2) + D[0] * (Tc**3-T_0**3) + E[0]* (Tc**4 - T_0**4)) + ( BL[0] * (T[i] - Tc) + CL[0] * (T[i]**2-Tc**2) + DL[0] * (T[i]**3-Tc**3) + EL[0]* (T[i]**4 - Tc**4))
                if T[i] >=Tc:
                    F[i] = ( BL[0] * (T[i] - T_0) + CL[0] * (T[i]**2-T_0**2) + DL[0] * (T[i]**3-T_0**3) + EL[0]* (T[i]**4 - T_0**4))
    
    
    Tc = 130
    if GoB != None:
        if T_0 < Tc:    
            for i in range(len(T)):
                if T[i] < Tc:
                    F[i] = 0
                if T[i] > Tc:
                    F[i] =  ( B[0] * (T[i] - Tc) + C[0] * (T[i]**2-Tc**2) + D[0] * (T[i]**3-Tc**3) + E[0]* (T[i]**4 - Tc**4))
        if T_0 >=Tc:
            for i in range(len(T)):
                if T[i] < Tc:
                    F[i] =  ( B[0] * (Tc - T_0) + C[0] * (Tc**2-T_0**2) + D[0] * (Tc**3-T_0**3) + E[0]* (Tc**4 - T_0**4))
            
    if CombiKT != None:
        B = FitL[0] 
        C = FitL[1] 
        D = [0,0]
        if len(FitL) >= 3:
            D = FitL[2]
        E = [0,0]
        if len(FitL) == 4:
            E = FitL[3]
        F = ( B[0] * (T - T_0) + C[0]/2 * (T**2-T_0**2) + D[0]/3 * (T**3-T_0**3) + E[0]/4* (T**4 - T_0**4))
        Ferr = np.zeros(len(T)) + 10*10**-6
        Tc = 130
        if GoB != None:
            if T_0 < Tc:    
                for i in range(len(T)):
                    if T[i] < Tc:
                        F[i] = 0
                    if T[i] > Tc:
                        F[i] =  ( B[0] * (T[i] - Tc) + C[0]/2 * (T[i]**2-Tc**2) + D[0]/3 * (T[i]**3-Tc**3) + E[0]/4* (T[i]**4 - Tc**4))
            if T_0 >=Tc:
                for i in range(len(T)):
                    if T[i] < Tc:
                        F[i] =  ( B[0] * (Tc - T_0) + C[0]/2 * (Tc**2-T_0**2) + D[0]/3 * (Tc**3-T_0**3) + E[0]/4* (Tc**4 - T_0**4))
                
            
        

    return [F,  Ferr]



def SysErrF(dv,dverr,T,Terr,Terrsys,T_R,T_R_err,T_R_sys,DMS,DMSerrsys,GoB,NFit):
    #Terrsys = 0.15+0.002*np.abs(T-273.15)
    N = 1000
    #4543
    
    m = K_e(0,GoB,IfM = True)
    IntDMST, DMST = DMSM.IntStrain_Al(T)
    IntDMST0, DMST0 = DMSM.IntStrain_Al(T_R)
    
    dvr = -dv - K_e(T, GoB)*DMS - DMST0*m*(T - T_R) + m*(IntDMST - IntDMST0)
    dvr_std = np.sqrt( dverr**2 + (K_e(T,GoB)*DMSM.CTE_Al(T) * Terr)**2 + (K_e(T_R,GoB)*DMSM.CTE_Al(T_R) * T_R_err)**2)
    
    Dat = [T,dvr ,Terr,dvr_std]
    Fit_P = Funk.Fit_Poli(Dat, NFit,SonderWertfK=[T_R,T_R_err])
    KP = 0
    
    if len(Fit_P.beta) == 3:
        KP = Fit_P.beta[1]  +  2*Fit_P.beta[2]*T  
    if len(Fit_P.beta) == 4:
        KP = Fit_P.beta[1]  +  2*Fit_P.beta[2]*T + 3*Fit_P.beta[3]*T**2
    if len(Fit_P.beta) == 5:
        KP = Fit_P.beta[1]  +  2*Fit_P.beta[2]*T + 3*Fit_P.beta[3]*T**2 + 4*Fit_P.beta[4]*T**3
            
    Prameter_Sys = [[],[],[],[],[]]
    Ksys = []

    if GoB ==  0:
        KE_sys = np.sqrt(0.0008**2 + 0.004**2)
        KE_2_sys = np.sqrt( 0.0012**2 +  0.0047**2)
        KE_22_sys = np.sqrt( 0.00081**2 +  0.0049**2)
        DatKE = [np.array([ 76.99078932,  77.00021746, 299.1261369 ]),np.array([0.7917514291367586, 0.7891880055332932, 0.8289872093707367]),np.array([0.00934783, 0.01503821, 0.03912767]),np.array([ 0.0012280583530470549, 0.000999730073825758, 0.0008171649643477754])]
        
    else:
        KE_sys = np.sqrt( 0.0016004631547986187**2 +  0.0047**2)
        KE_2_sys = np.sqrt(0.0014537841364968623**2 +  0.0047**2)
        KE_22_sys = np.sqrt(0.0007951016906917333**2 +  0.0049**2)
        DatKE = [np.array([ 76.99078932,  77.00021746, 299.1261369 ]),np.array([0.7935870117981854 , 0.789758467355037, 0.8278920440705645]),np.array([0.00934783, 0.01503821, 0.03912767]),np.array([ 0.0016004631547986187,0.0014537841364968623 , 0.0007951016906917333])]    
        
    
    for i in range(N):
        DatKE[1][0] += np.random.uniform(-1,1) * KE_sys
        DatKE[1][1] += np.random.uniform(-1,1) * KE_2_sys
        DatKE[1][2] += np.random.uniform(-1,1) * KE_22_sys
        
        DatKE[0][0] += np.random.uniform(-1,1) * 0.6#Terrsys
        DatKE[0][1] += np.random.uniform(-1,1) * 0.6#Terrsys
        DatKE[0][2] += np.random.uniform(-1,1) * 0.15#Terrsys
        
        Fit_KE= Funk.Fit_Poli(DatKE,1)
        c = Fit_KE.beta[0]
        c_stat = Fit_KE.sd_beta[0]
        m = Fit_KE.beta[1]
        m_stat = Fit_KE.sd_beta[1]   
        
        Ke_sysi = np.random.uniform(-1,1) * KE_sys #np.array(Ke_sysi)
        
        Kestat,Kesys = K_e(T,GoB,IfKerr = True)
        Ke_sysi =np.random.uniform(-1,1)*Kesys
        
        DMS_sysi = np.random.uniform(-1,1) * DMSerrsys
        Terrsysi = np.random.uniform(-1,1)*Terrsys
        T_R_errsysi = np.random.uniform(-1,1)*T_R_sys
        
        IntDMST, DMST = DMSM.IntStrain_Al(T + Terrsysi)
        IntDMST0, DMST0 = DMSM.IntStrain_Al(T_R + T_R_errsysi)
         
        
        dvr_std = np.sqrt( dverr**2 + ( (m*T + c)*DMSM.CTE_Al(T) * Terr)**2 + ( (m*T_R + c)*DMSM.CTE_Al(T_R) * T_R_err)**2)
        
        dvr = -dv - (K_e((T+Terrsysi + T_R + T_R_errsysi)/2,GoB) +  Ke_sysi   )*(DMS + DMS_sysi)
             
        
        Dat = [T + Terrsysi,dvr ,Terr,dvr_std]
        Fit_sys = Funk.Fit_Poli(Dat, NFit,SonderWertfK=[T_R + T_R_errsysi,T_R_err])
        
        
        for k in range(1,len(Fit_P.beta)):
            Prameter_Sys[k].append(Fit_sys.beta[k] - Fit_P.beta[k]) 
        

        if len(Fit_P.beta) == 3:
            Ksys.append( Fit_sys.beta[1]  +  2*Fit_sys.beta[2]*T  - KP )
        if len(Fit_P.beta) == 4:
            Ksys.append( Fit_sys.beta[1]  +  2*Fit_sys.beta[2]*T + 3*Fit_sys.beta[3]*T**2  - KP )
        if len(Fit_P.beta) == 5:
            Ksys.append(Fit_sys.beta[1]  +  2*Fit_sys.beta[2]*T + 3*Fit_sys.beta[3]*T**2 + 4*Fit_sys.beta[4]*T**3  - KP)
        
    
    Ksys = np.array(Ksys)
    N2 = np.linspace(0,len(T)-1,20)
    SysKT = []
    SysT = []

    for i in N2:

        SysKT.append( np.std(Ksys[:,int(i)],ddof = 1)*3) # Abschätzung nach oben da sonst zu kleiner Feher  
        SysT.append(T[int(i)])
        
    DatSys = [np.array(SysT),np.array(SysKT),np.zeros(len(SysT))+0.1,np.zeros(len(SysT))+0.001*10**-6]
    Fit_sys_Para = Funk.Fit_w_Plot(DatSys,[],2,"","","","Bilder/gmn-js-ramp--Sys_FT","Sys vs T",AT)
    
    Fit_P_sys =[ Fit_sys_Para.beta[0] ,Fit_sys_Para.beta[1],Fit_sys_Para.beta[2] ]#np.mean(SysKT)

    print(" SYS Parameter Error K: ", Fit_P_sys)
    
    return Fit_P_sys  


Breite = 5
DTS = 10
DTB = 20

Start = 0
Ende = len(T_MK) 

IntDMST, DMST = DMSM.IntStrain_Al(T_MK[Start:Ende])
IntDMST0, DMST0 = DMSM.IntStrain_Al(T_MK_Ref)

dvrK = -dv_G_MK[Start:Ende] - K_e((T_MK[Start:Ende] + T_MK_Ref)/2, 0)*DMS_MK[Start:Ende] 
dvrK_std = np.sqrt( dv_G_MK_std[Start:Ende]**2 + (K_e(T_MK[Start:Ende],0)*DMSM.CTE_Al(T_MK[Start:Ende]) * T_MK_std)**2 + (K_e(T_MK_Ref,0)*DMSM.CTE_Al(T_MK_Ref) * T_MK_Ref_std)**2)
Dat = [T_MK[Start:Ende],dvrK ,0.1,0.1]

Text1 = "Germanium Fiber"# Ger.-Data"  
Text2 = r"Fit quadratic $K_T$" 
Text3 = r"$Fit(T)\,=\,A\cdot (T_M-T_R) \,+\,\frac{B}{2}\cdot(T_M^2-T_R^2)$"
Name_PNG = "Bilder/gmn-js-ramp--ger-fit" + Zusatz_Name

AT = ['T [K]', r"$-\frac{\Delta \nu}{\hat{\nu}}_{Thermal} \times 10^{-6}$" , r'$(Data-Fit(T) \times 10^{-6}$']

Titel = ""
Fit_G_Klima = Funk.Fit_w_Plot(Dat,[],11,Text1,Text2,Text3,Name_PNG,Titel,AT,SonderWertfK=[T_MK_Ref,T_MK_Ref_std])

Fit_G_KlimaM = Fit_G_Klima
FP = Fit_G_Klima.beta
FPs = Fit_G_Klima.sd_beta
Fit_G_Klima = [[0,0], [FP[1],FPs[1]],[FP[2],FPs[2]],[0,0],[0,0]]

Start = 0
Ende = len(T_MK) 

IntDMST, DMST = DMSM.IntStrain_Al(T_MK[Start:Ende])
IntDMST0, DMST0 = DMSM.IntStrain_Al(T_MK_Ref)
dvrK = -dv_B_MK[Start:Ende] - K_e((T_MK[Start:Ende] + T_MK_Ref)/2, 1)*DMS_MK[Start:Ende] 
dvrK_std = np.sqrt( dv_B_MK_std[Start:Ende]**2 + (K_e(T_MK[Start:Ende],1)*DMSM.CTE_Al(T_MK[Start:Ende]) * T_MK_std)**2 + (K_e(T_MK_Ref,1)*DMSM.CTE_Al(T_MK_Ref) * T_MK_Ref_std)**2)


Dat = [T_MK[Start:Ende],dvrK ,0.1,0.1]

Text1 = "Boron Fiber"  
Text2 = r"Fit quadratic $K_T$" 
Text3 = r"$Fit(T)\,=\,A\cdot (T_M-T_R) \,+\,\frac{B}{2}\cdot(T_M^2-T_R^2)$" + "\n" + r"$ \,+\,\frac{C}{3}\cdot(T_M^3-T_R^3) $"

Name_PNG = "Bilder/gmn-js-ramp--boron-fit" + Zusatz_Name
AT = ['T [K]', r"$-\frac{\Delta \nu}{\hat{\nu}}_{Thermal} \times  10^{-6}$" , r"$(Data-Fit(T) \times 10^{-6}$"]
Titel = ""#r"Bor-Fiber $F(T_0,T)$ Pol-" + str(12-9) + " Fit, Climate Chamber"
Fit_B_Klima = Funk.Fit_w_Plot(Dat,[],12,Text1,Text2,Text3,Name_PNG,Titel,AT,SonderWertfK=[T_MK_Ref,T_MK_Ref_std])
Fit_B_KlimaM = Fit_B_Klima
FP = Fit_B_Klima.beta
FPs = Fit_B_Klima.sd_beta
Fit_B_Klima = [[0,0], [FP[1],FPs[1]],[FP[2],FPs[2]],[FP[3],FPs[3]],[0,0]]

Tx = np.linspace(77,290,300)
Tx2 = np.linspace(100,290,300)
Tx2B = np.linspace(77,100,300)
TP = np.linspace(50,77)

TxK = np.linspace(295,325,300)

Fit_B_K_sys = SysErrF(dv_B_MK[Start:Ende],dv_B_MK_std[Start:Ende],T_MK[Start:Ende],T_MK_std[Start:Ende],T_MK_sys[Start:Ende],T_MK_Ref,T_MK_Ref_std,T_MK_Ref_sys,DMS_MK[Start:Ende],DMS_MK_sys[Start:Ende],1,12)
Fit_G_K_sys = SysErrF(dv_G_MK[Start:Ende],dv_G_MK_std[Start:Ende],T_MK[Start:Ende],T_MK_std[Start:Ende],T_MK_sys[Start:Ende],T_MK_Ref,T_MK_Ref_std,T_MK_Ref_sys,DMS_MK[Start:Ende],DMS_MK_sys[Start:Ende],0,11)

################################################### Paper Plot KT

#IMPORTANT
plt.figure(194,dpi = 180) 

Tx = np.linspace(77,290,400)

plt.plot(TxK,K_T(TxK,Fit_B_Klima,Fit_B_KlimaM)[0]*10**6,label = r"$K_T$ - Boron Fiber",color ="tab:orange", linestyle='--')
y1 = (K_T(TxK,Fit_B_Klima,Fit_B_KlimaM)[0]-K_T(TxK,Fit_B_Klima,Fit_B_KlimaM,Fit_B_K_sys)[1])*10**6
y2 = (K_T(TxK,Fit_B_Klima,Fit_B_KlimaM)[0]+K_T(TxK,Fit_B_Klima,Fit_B_KlimaM,Fit_B_K_sys)[1])*10**6
plt.fill_between(TxK, y1, y2 , alpha = 0.5, color ="tab:orange")


plt.plot(TxK,K_T(TxK,Fit_G_Klima,Fit_G_KlimaM)[0]*10**6,label = r"$K_T$ - Germanium Fiber", color ="tab:blue", linestyle='--')

y1 = (K_T(TxK,Fit_G_Klima,Fit_G_KlimaM)[0]-K_T(TxK,Fit_G_Klima,Fit_G_KlimaM,Fit_G_K_sys)[1])*10**6
y2 = (K_T(TxK,Fit_G_Klima,Fit_G_KlimaM)[0]+K_T(TxK,Fit_G_Klima,Fit_G_KlimaM,Fit_G_K_sys)[1])*10**6
#
plt.fill_between(TxK, y1, y2 , alpha = 0.5, color ="tab:blue")



plt.xlabel("T [K]")
plt.ylabel(r"$K_T \times 10^{-6}$")
plt.legend(fontsize=7)
plt.xlim(295,325)
plt.grid()
plt.savefig("Bilder/gmn-js-ramp--KT_Optical_Fiber_Bor_Ger" + Zusatz_Name + ".eps")
plt.savefig("Bilder/gmn-js-ramp--KT_Optical_Fiber_Bor_Ger" + Zusatz_Name + ".pdf")
plt.show()


