import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import scipy.integrate as integrate
import scipy.odr as odr
from numpy import heaviside

def Add_Mean(Dat,Mean):
    for i in range(len(Dat)):
        Dat[i] += Mean
    return Dat
def Mean_from_Position(X,Y,p1,p2):
    m = []
    mx = []
    for i in range(len(X)):
        if X[i]>= p1 and X[i]<= p2:
            m.append(Y[i])
            mx.append(X[i])
    return [np.mean(mx), np.mean(m), np.std(m,ddof = 1)]

def G_Mean_from_Position(X,Y,YS,p1,p2):
    m = []
    ms = []
    mx = []
    for i in range(len(X)):
        if X[i]>= p1 and X[i]<= p2:
            m.append(Y[i])
            ms.append(YS[i])
            mx.append(X[i])
    return [np.mean(mx), Gewichtetes_Mittel(np.array(m),np.array(ms))[0], Gewichtetes_Mittel(np.array(m),np.array(ms))[1]]

def Std_from_Position(X,Y,p1,p2):
    m = []
    mx = []
    for i in range(len(X)):
        if X[i]>= p1 and X[i]<= p2:
            m.append(Y[i])
            mx.append(X[i])
    return [np.mean(mx), np.std(m,ddof = 1)] 

def ListMfP(X,Y,S,St,En):
    mmx = []
    mm = []
    P = St
    while P + S <= En:
        D = Mean_from_Position(X,Y,P,P+S)
        mmx.append(D[0])
        mm.append(D[1])
        P += S
    return [mmx,mm]
        
def Gewichtetes_Mittel(X,XS):
    Weights = 1/XS**2
    Mean = np.average(X,weights = Weights)
    STD = np.sqrt(1/sum(Weights))
    return [Mean,STD]
    

def Pol4(B,x):
    return B[0] + B[1] * x + B[2]*x**2 + B[3]*x**3 + B[4]* x**4 

def Sig_Pol4(B,x):
#    return np.sqrt(  (B[1])**2 + (B[2]*x*2)**2+ (3*B[3]*x**2)**2+ (4*B[4]*x**3)**2 )
    return np.sqrt(  B[0]**2+ (B[1]*x)**2 + (B[2]*x*2)**2+ (B[3]*x**3)**2+ (B[4]*x**4)**2 )

def Pol4_p_lin(B,x):
    
    miny = 1000000
    mini = 0
    for i in range(20,70):
        yi = Pol4(B,i)
        if yi <0:
            if yi <= miny:
                miny = yi
                mini = i
    mini =   -400
    re = np.array([])
    for i in range(len(x)):
        if x[i]<mini:
            re = np.append(re,miny)
        if x[i]>=mini:
            re = np.append(re,B[0] + B[1] * x[i] + B[2]*x[i]**2 + B[3]*x[i]**3 + B[4]* x[i]**4)
    return re
def Sig_Pol4_p_lin(B,x):
    miny = 1000000
    minys = 0
    mini = 0
    for i in range(20,70):
        yi = Pol4(B,i)
        if yi<=miny:
            miny = yi
            minys = Sig_Pol4(B,i)
            mini = i
    mini = -400
    re = np.array([])
    for i in range(len(x)):
        if x[i]<mini:
            re = np.append(re,minys)
        if x[i]>=mini:
            re = np.append(re,np.sqrt(  (B[1])**2 + (B[2]*x[i]*2)**2+ (3*B[3]*x[i]**2)**2+ (4*B[4]*x[i]**3)**2 ))
    return re

def Pol5(B,x):
    return B[0] + B[1] * x + B[2]*x**2 + B[3]*x**3 + B[4]* x**4 + B[5]* x**5 

def Pol3(B,x):
    return B[0] + B[1] * x + B[2]*x**2 + B[3]*x**3 
def Pol2(B,x):
    return B[0] + B[1] * x + B[2]*x**2 
def Pol1(B,x):
    return B[0] + B[1] * x 
def Pol0(B,x):
    return B[0]
# objective function
def objective5(x,a,b,c,d,e,f):
    return a + b * x + c*x**2 + d*x**3 + e*x**4 + f*x**5

def objective4(x,a,b,c,d,e):
    return a + b * x + c*x**2 + d*x**3 + e*x**4

def objective3(x, a, b, c, d ):
	return a  + b*x + c*x**2 + d*x**3 

def objective2(x, a, b, c):
	return a  + b*x + c*x**2 

def objective1(x, a, b):
	return a  + b*x 
def objective0(x, a):
	return a  




def LinFast(X,Y):
    m = (Y[1]-Y[0])/(X[1]-X[0])
    b = Y[0] - m*X[0]
    return [m,b]
def Lin_Min_Max(a,b,c):
    A = np.sqrt(b**2 - 4*a*c)/2/a
    return [A - b/2/a , -A - b/2/a]

def Dv_to_T_Lin(KGT,KGE,KBT,KBE,VG,VB,DT,TM):
    Alpha = KBE/KGE # G to B
    
    LG = LinFast(DT,KGT)
    LB = LinFast(DT,KBT)
    
    LG = [LG[0]*Alpha,LG[1]*Alpha]
    VG = VG*Alpha
    #print(DT)
    #print((-LG[0]*196+LG[1])*10**6)
    a = LG[0]-LB[0]
    b = LG[1]-LB[1]
    c = (-VG+VB)-a*TM**2-b*TM
    #print(a,b,c)
    return Lin_Min_Max(a,b,c)
    
def RTE_T(T,M):
    par = [25,0,0,0,0]
    T = T + 273.15
    if M == 0:
        par = [-419.67,  0.1313,  2.924e-3,  1.0796e-5, -1.9105e-8] 
    if M == 1:
        par = [ -1.1503e+3	 ,     1.0535	,    2.3016e-2	,  -6.8525e-5	, 7.7033e-8]
    return (par[0] + par[1]*T + par[2]*T*T + par[3]*T*T*T + par[4]*T*T*T*T)*10**-1

def CTE_Peek(T):
    par = [50 , -1.1503e+3	 ,     1.0535	,    2.3016e-2	,0]#,  -6.8525e-5	, 7.7033e-8]
    return (par[0] + par[1]*T + par[2]*T*T + par[3]*T*T*T + par[4]*T*T*T*T)*10**-1

def RTE_Peek(T,T0):
    S= np.array([])
    for i in T:
        S = np.append(S , integrate.quad(CTE_Peek, T0, i)[0])
    return S



#T = np.linspace(77,300,1000)
#plt.plot(T, RTE_Peek(T,293))
#plt.show()     

def CTE_T(T,M):
    par = [25,0,0,0,0]
    if M == 0:
        par = [-419.67,  0.1313,  2.924e-3,  1.0796e-5, -1.9105e-8] 
    if M == 1:
        par = [ -1.1503e+3	 ,     1.0535	,    2.3016e-2	,  -6.8525e-5	, 7.7033e-8]
    return  (par[1] + 2*par[2]*T + 3*par[3]*T*T + 4*par[4]*T*T*T)*10
#    return  (RTE_T(T+10,M)-RTE_T(T-10,M))/20

def Strain_Theo(T,TM,M):
    S = []
    for i in range(len(T)):
        S.append(RTE_T(T[i]+TM+273.15,M) - RTE_T(TM+273.15,M))
       # S.append(CTE_T(1,0)*T[i])
    return S
        
def T_Theo_Luna(Dv,par):
    if len(par) == 4:
        par = np.append(par,0)
    if len(par) == 3:
        par = np.append(par,[0,0])
    return par[0] + par[1]*Dv + par[2]*Dv**2 + par[3]*Dv**3 + par[4] * Dv**4




def T_Theo(Dv,par):
    if len(par) == 4:
        par = np.append(par,0)
    if len(par) == 3:
        par = np.append(par,[0,0])
    if len(par) == 2:
        par = np.append(par,[0,0,0])
    T = par[0] + par[1]*Dv + par[2]*Dv**2 + par[3]*Dv**3 + par[4] * Dv**4
#    Ts = np.sqrt( pars[0]**2 + (pars[1]*Dv)**2 + (pars[2]*Dv**2)**2+ (pars[3]*Dv**3)**2+ (pars[4]*Dv**4)**2 )
    return T

def T_Theo2(Dv,par):
    if len(par) == 4:
        par = np.append(par,0)
    if len(par) == 3:
        par = np.append(par,[0,0])
    if len(par) == 2:
        par = np.append(par,[0,0,0])

    Ts = np.sqrt( par[0]**2 + (par[1]*Dv)**2 + (par[2]*Dv**2)**2+ (par[3]*Dv**3)**2+ (par[4]*Dv**4)**2 )
    return Ts

def K_Theo2(T,par,pars):
    if len(par) == 4:
        par = np.append(par,0)
    if len(par) == 3:
        par = np.append(par,[0,0])
    dv = np.linspace(-400,1000,2000)
    dt = T_Theo_Luna(dv,par) 
    Ks = -7.32E-1 #linea approx
    for i in range(2000-1):
        if T<=dt[i] and T>=dt[i+1]:
            Ks = np.sqrt(  (pars[1])**2 + (pars[2]*dv[i]*2)**2+ (3*pars[3]*dv[i]**2)**2+ (4*pars[4]*dv[i]**3)**2 )
            break
    return Ks

def k_t(dv,T,S,par):
    kmin = []
    for i in range(len(T)):
        kmin.append(dv[i]/T[i] + S[i]/K_Theo(T[i],par)/T[i]  )
    K = 1/np.array(kmin)
    return K
    
    
    
def K_Theo_Luna(T,par):
    dv = np.linspace(290,-170,600)
    dt = T_Theo_Luna(dv,par) + 273.15
    K = -7.32E-1 #linea approx
    for i in range(600-1):
        if T<=dt[i] and T>=dt[i+1]:
            K = par[1] + par[2]*dv[i]*2 + 3*par[3]*dv[i]**2 + 4*par[4] * dv[i]**3
            break
    return K

def K_Theo(T,par):
    if len(par) == 4:
        par = np.append(par,0)
    if len(par) == 3:
        par = np.append(par,[0,0])
    dv = np.linspace(-400,1000,10000)
    dt = T_Theo_Luna(dv,par) 
  #  print(dt)
    K = -7.32E-1 #linea approx
    for i in range(10000-1):
        if T<=dt[i] and T>=dt[i+1]:
            K = par[1] + par[2]*dv[i]*2 + 3*par[3]*dv[i]**2 + 4*par[4] * dv[i]**3
         #   print("AA")
            break
    return K

def n(T):
    C = 0
    T0 = 293
    B = [  (0.691663 + 0.1107001 *C) , (0.0684043 + 0.000568306*C)**2 , (0.4079426 + 0.31021588*C) , (0.1162414 + 0.03772465*C)**2 , (0.8974749 - 0.043311091*C) , (9.896161 + 1.94577*C)**2]
    Lam = 1310
    n2 = 1+ B[0]*Lam**2/(Lam**2 - (B[1]*(T/T0)**2)**2) + B[2]*Lam**2/(Lam**2 - (B[3]*(T/T0)**2)**2) + B[4]*Lam**2/(Lam**2 - B[5]**2)
    
    return np.sqrt(n2)

def n_C_7980(T):
    Lam = 1310
    S1 = [1.10127E+00, -4.94251E-05 ,  5.27414E-07, -1.59700E-09 ,  1.75949E-12]
    S2 = [1.78752E-05 , 4.76391E-05, -4.49019E-07, 1.44546E-09, -1.57223E-12]
    S3 = [7.93552E-01, -1.27815E-03, 1.84595E-05, -9.20275E-08, 1.48829E-10]
    
    L1 = [-8.90600E-02,  9.08730E-06,   -6.53638E-08, 7.77072E-11  , 6.84605E-14 ]
    L2 = [2.97562E-01, -8.59578E-04 , 6.59069E-06 , -1.09482E-08,  7.85145E-13]
    L3 = [9.34454E+00, -7.09788E-03,  1.01968E-04, -5.07660E-07, 8.21348E-10]
    
    n2 = 1 + Pol4(S1,T)*Lam**2/(Lam**2 - Pol4(L1,T)**2) + Pol4(S2,T)*Lam**2/(Lam**2 - Pol4(L2,T)**2) + Pol4(S3,T)*Lam**2/(Lam**2 - Pol4(L3,T)**2)
    return np.sqrt(n2)

def dn(T):
    d = []
    Ti =  np.append(T,np.max(T)*1.01 )
    nn = n_C_7980(Ti)
    for i in range(len(T)):
        d.append(  (nn[i]-nn[i+1])/(Ti[i] - Ti[i+1])  )
    return np.array(d)
 
def dEE(T):
    A = 4.7
    B = 636
    E0 = 10.4
    return -1/(E0 - A*T**2/(T+B)) * A*T*(T+2*B)/(T+B)**2    

def dn2(T):
    Alpha = 0.62
    nU = 1.44
    Lam =1310
    Lam0 = 0.0684043
    K = (nU**2 - 1)
    G = -3*CTE_T(T,0)*K**2   
    R = Lam**2/(Lam**2 - Lam0**2)
    
    Eg = 10.4
    dEg = -2.5
   # H = -1/Eg*dEg*K**2
    H = -dEE(T)*K**2
    return (G*R + H*R**2)/2/n(T)

def T_K_T(T):
    Alpha = 0.62
    nU = 1.44
    Lam =1310
    Lam0 = 0.0684043
    K = (nU**2 - 1)
    G = -3*CTE_T(T,0)*K**2   
    R = Lam**2/(Lam**2 - Lam0**2)
    
    Eg = 10.4
    dEg = -2.5
    #H = -1/Eg*dEg*K**2
    H = -dEE(T)*K**2
    return  H*R**2/2/n(T)**2
    
def Alpha_Fiber(T):
    a = -4.22
    b = 35.5
    c = 0.335
    d = 1.253
    e = 535      
    return  (a*(b/T)**c * np.exp(b/T) * (np.exp(b/T) + 1)**(-2) + d*(e/T)**2 * np.exp(e/T)*(np.exp(e/T) - 1)**(-2))*10**-6

def Strain_Fiber(T,T0):
    S= np.array([])
    for i in T:
        S = np.append(S , integrate.quad(Alpha_Fiber, T0, i)[0])
    return S

  
def chisq(obs, exp, error):
    return np.sum((obs - exp) ** 2 / (error ** 2))

from scipy.stats import t


def prediction_interval( dfdp, x, y, yerr, signif, popt, pcov):
    """
    * Calculate Preduction Intervals
    *  adapted from Extraterrestrial Physics pdf
    *
    *  func    : function that we are using
    *  dfdp    : derivatives of that function (calculated using sympy.diff)
    *  x       : the x values (calculated using numpy.linspace)
    *  y       : the y values (calculated by passing the ODR parameters and x to func)
    *  y_err   : the maximum residual on the y axis
    *  signif  : the significance value (68., 95. and 99.7 for 1, 2 and 3 sigma respectively)
    *  popt    : the ODR parameters for the fit, calculated using scipy.odr.run()
    *  pcov    : the covariance matrix for the fit, calculated using scipy.odr.run()
    """
    
    # get number of fit parameters and data points
    npp = len(popt)
    n = len(x)

    # convert provided value to a significance level (e.g. 95. -> 0.05), then calculate alpha
    alpha = 1. - (1 - signif / 100.0) / 2

    # student’s t test
    tval = t.ppf(alpha, n - npp)
    print(tval)

    # process covarianvce matrix and derivatives
    
    d = np.zeros(n)
    for j in range(npp):
        for k in range(npp):
            d += dfdp[j] * dfdp[k] * pcov[j,k]
    print(d)
    print(yerr)
    # return prediction band offset for a new measurement
    return tval * np.sqrt(yerr**2 + d)


def confidence_band(x, dfdp, confprob, fitobj, covscale,dof, abswei=False):
   #----------------------------------------------------------
   # Given a value for x, calculate the error df in y = model(p,x)
   # This function returns for each x in a NumPy array, the
   # upper and lower value of the confidence interval. 
   # The arrays with limits are returned and can be used to
   # plot confidence bands.  
   # 
   #
   # Input:
   #
   # x        NumPy array with values for which you want
   #          the confidence interval.
   #
   # dfdp     A list with derivatives. There are as many entries in
   #          this list as there are parameters in your model.
   #
   # confprob Confidence probability in percent (e.g. 90% or 95%).
   #          From this number we derive the confidence level 
   #          (e.g. 0.05). The Confidence Band
   #          is a 100*(1-alpha)% band. This implies
   #          that for a given value of x the probability that
   #          the 'true' value of f falls within these limits is
   #          100*(1-alpha)%.
   # 
   # fitobj   The Fitter object from a fit with kmpfit
   #
   # f        A function that returns a value y = f(p,x)
   #          p are the best-fit parameters and x is a NumPy array
   #          with values of x for which you want the confidence interval.
   #
   # abswei   Are the weights absolute? For absolute weights we take
   #          unscaled covariance matrix elements in our calculations.
   #          For unit weighting (i.e. unweighted) and relative 
   #          weighting, we scale the covariance matrix elements with 
   #          the value of the reduced chi squared.
   #
   # Returns:
   #
   # y          The model values at x: y = f(p,x)
   # upperband  The upper confidence limits
   # lowerband  The lower confidence limits   
   #
   # Note:
   #
   # If parameters were fixed in the fit, the corresponding 
   # error is 0 and there is no contribution to the condidence
   # interval.
   #----------------------------------------------------------   
   from scipy.stats import t
   # Given the confidence probability confprob = 100(1-alpha)
   # we derive for alpha: alpha = 1 - confprob/100 
   alpha = 1 - confprob/100.0
   prb = 1.0 - alpha/2
   tval = t.ppf(prb, dof)
   
   C = fitobj.cov_beta
   n = len(fitobj.beta)              # Number of parameters from covariance matrix
#   p = fitobj.params
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
 #  y = f(p, x)
   delta = tval * df   
   #upperband = y + delta
   #lowerband = y - delta 
   return delta # y, upperband, lowerband


    
def Spesific_func_confidence_band(x, confprob, fitobj, Funk_Number, Yerr = None, abswei=False , Prediction_I = False):
    p = fitobj.beta
    covscale = fitobj.sum_square/(len(fitobj.delta) - len(p))
    dof = len(p)
    
    if Funk_Number == 11:       
        dfdp = [np.zeros(len(x)) + -p[1] - p[2]*2*p[0], x - p[0], x**2 -p[0]**2]
        covscale = fitobj.sum_square/(len(fitobj.delta) - len(p) - 1)
        dof = len(p) - 1
        
    if Funk_Number == 12:       
        dfdp = [np.zeros(len(x)) + -p[1] - p[2]*2*p[0] - p[3]*3*p[0]**2, x - p[0], x**2 -p[0]**2, x**3 - p[0]**3]
        covscale = fitobj.sum_square/(len(fitobj.delta) - len(p) - 1)
        dof = len(p) - 1
        
    if Funk_Number == 13:       
        dfdp = [np.zeros(len(x)) + -p[1] - p[2]*2*p[0] - p[3]*3*p[0]**2 - p[4]*4*p[0]**3, x - p[0], x**2 -p[0]**2, x**3 - p[0]**3, x**4 - p[0]**4]
        covscale = fitobj.sum_square/(len(fitobj.delta) - len(p) - 1)
        dof = len(p) - 1
        
    if Funk_Number == 112:  
        dfdp = [np.zeros(len(x)) + -p[1] - p[2]*2*p[0], x - p[0], x**2 -p[0]**2, KPol3H(fitobj.beta,x)*1/(1+np.exp(-10*(x-p[3]))) * (1-1/(1+np.exp(-10*(x-p[3])))), 1-1/(1+np.exp(-10*(x-p[3])))]
        covscale = fitobj.sum_square/(len(fitobj.delta) - len(p) - 1)
        dof = len(p) - 1
        
    if Funk_Number == 122:  
        dfdp = [np.zeros(len(x)) + -p[1] - p[2]*2*p[0] - p[3]*3*p[0]**2, x - p[0], x**2 -p[0]**2, x**3 - p[0]**3, KPol4H(fitobj.beta,x)*1/(1+np.exp(-10*(x-p[4]))) * (1-1/(1+np.exp(-10*(x-p[4])))), 1-1/(1+np.exp(-10*(x-p[4])))]
        covscale = fitobj.sum_square/(len(fitobj.delta) - len(p) - 1)
        dof = len(p) - 1
     
    if Funk_Number == 0:
        delta = 0
        return delta
    
    if Funk_Number < 10:
        print("Fit Number: " + str(Funk_Number))
        dfdp = [np.zeros(len(x)) + 1]
        for i in range(1,len(p)):
            dfdp.append(x**i)
            print("Number of parameter: " + str(len(p)) + " diff. : x^" + str(i) )
            
    if Prediction_I == True:
        return np.sqrt(confidence_band(x, dfdp, confprob, fitobj,covscale,dof)**2 + Yerr**2)
    
    return confidence_band(x, dfdp, confprob, fitobj,covscale,dof)


def format_zahl(number):
    if number == 0:
        return 0, 0
    exponent = 0
    if abs(number) < 1:
        while abs(number) < 1:
            number *= 10
            exponent -= 1
    else:
        while abs(number) >= 10:
            number /= 10
            exponent += 1
    return round(number, 2), exponent


def format_string(number, errr, nachkommerstellen=1):
    gerundete_zahl, exponent = format_zahl(number)
    gerundete_zahlerr, exponenterr = format_zahl(errr)
    gerundete_zahlerr = gerundete_zahlerr*10**(exponenterr-exponent)
    if exponent==0:
        return f"{gerundete_zahl:.{nachkommerstellen}f}" +r"$\, \pm \,$" + f"{gerundete_zahlerr:.{nachkommerstellen}f}"
    else:
       # nachkommerstellen -=1
        return f"({gerundete_zahl:.{nachkommerstellen}f}" +r"$\, \pm \,$" + f"{gerundete_zahlerr:.{nachkommerstellen}f})"  + r"$\times 10^{"+str(exponent)+"}$"


def Text_legend(myoutput,P):
    Text = []
    Num = []
    if P > 0:
        for i in range(len(myoutput.beta)):
            Ne = np.floor(np.log10(np.abs(myoutput.beta[i]))).astype(int) - np.floor(np.log10(np.abs(myoutput.sd_beta[i]))).astype(int)
            if P > 10 and i>0 :
                Ne = np.floor(np.log10(np.abs(myoutput.beta[i]/(i)))).astype(int) - np.floor(np.log10(np.abs(myoutput.sd_beta[i]/i))).astype(int)
            if Ne < 0:
                Ne = 0
            Text.append("{:."+str(Ne+1)+"e}")
            Num.append(Ne+1)
        
    if P == 0:
        Ne = np.floor(np.log10(np.abs(myoutput[0]))).astype(int) - np.floor(np.log10(np.abs(myoutput[1]))).astype(int)
        #NT = "{:."+str(Ne+1)+"e}" 
        #print(Ne)
        if Ne > 1:
            TextP = "\n" + "WM = " + format_string(myoutput[0],myoutput[1],Ne  )
            
        if Ne == 1:
            gerundete_zahl, exponent = format_zahl(myoutput[0])
            gerundete_zahl *= 10**exponent
            gerundete_zahlerr, exponenterr = format_zahl(myoutput[1])
            TextP = "\n" + "WM = " +  f"{gerundete_zahl:.{1}f}" +r"$\, \pm \,$" + f"{gerundete_zahlerr:.{1}f}"
        if Ne == 0:
            TextP = "\n" + "WM = " + format_string(myoutput[0],myoutput[1],Ne  )
            #print(TextP,Ne)
            #print(1/0)
            
        
    if P == 1:      
        #TextP = "\n" + "a = " + Text[0].format(myoutput.beta[0]) +r"$\, \pm \,$" +"{:.1e}".format( myoutput.sd_beta[0]) + "\n"+ "b = " + Text[1].format(myoutput.beta[1]) +r"$\, \pm \,$" +"{:.1e}".format( myoutput.sd_beta[1])
        TextP = "\n" + "a = " + format_string(myoutput.beta[0],myoutput.sd_beta[0],Num[0]) +  "\n"+ "b = " + format_string(myoutput.beta[1],myoutput.sd_beta[1],Num[1])
        
        
    if P == 2:
        TextP = "\n" + "a = " + Text[0].format(myoutput.beta[0]) +r"$\, \pm \,$" +"{:.1e}".format( myoutput.sd_beta[0]) + "\n"+ "b = " + Text[1].format(myoutput.beta[1]) +r"$\, \pm \,$" +"{:.1e}".format( myoutput.sd_beta[1])+ "\n"+ "c = " + Text[2].format(myoutput.beta[2]) +r"$\, \pm \,$" +"{:.1e}".format( myoutput.sd_beta[2]) 
    if P == 3:
        TextP = "\n" + "a = " + Text[0].format(myoutput.beta[0]) +r"$\, \pm \,$" +"{:.1e}".format( myoutput.sd_beta[0]) + "\n"+ "b = " + Text[1].format(myoutput.beta[1]) +r"$\, \pm \,$" +"{:.1e}".format( myoutput.sd_beta[1])+ "\n"+ "c = " + Text[2].format(myoutput.beta[2]) +r"$\, \pm \,$" +"{:.1e}".format( myoutput.sd_beta[2]) + "\n"+ "d = " + Text[3].format(myoutput.beta[3]) +r"$\, \pm \,$" +"{:.1e}".format( myoutput.sd_beta[3])
    if P == 4:
        TextP = "\n" + "a = " + Text[0].format(myoutput.beta[0]) +r"$\, \pm \,$" +"{:.1e}".format( myoutput.sd_beta[0]) + "\n"+ "b = " + Text[1].format(myoutput.beta[1]) +r"$\, \pm \,$" +"{:.1e}".format( myoutput.sd_beta[1])+ "\n"+ "c = " + Text[2].format(myoutput.beta[2]) +r"$\, \pm \,$" +"{:.1e}".format( myoutput.sd_beta[2]) + "\n"+ "d = " + Text[3].format(myoutput.beta[3]) +r"$\, \pm \,$" +"{:.1e}".format( myoutput.sd_beta[3])+ "\n"+ "e = " + Text[4].format(myoutput.beta[4]) +r"$\, \pm \,$" +"{:.1e}".format( myoutput.sd_beta[4])
    if P == 5:
        TextP = "\n" + "a = " + Text[0].format(myoutput.beta[0]) +r"$\, \pm \,$" +"{:.1e}".format( myoutput.sd_beta[0]) + "\n"+ "b = " + Text[1].format(myoutput.beta[1]) +r"$\, \pm \,$" +"{:.1e}".format( myoutput.sd_beta[1])+ "\n"+ "c = " + Text[2].format(myoutput.beta[2]) +r"$\, \pm \,$" +"{:.1e}".format( myoutput.sd_beta[2]) + "\n"+ "d = " + Text[3].format(myoutput.beta[3]) +r"$\, \pm \,$" +"{:.1e}".format( myoutput.sd_beta[3])+ "\n"+ "e = " + Text[4].format(myoutput.beta[4]) +r"$\, \pm \,$" +"{:.1e}".format( myoutput.sd_beta[4]) + "\n"+ "f = " + Text[5].format(myoutput.beta[5]) +r"$\, \pm \,$" +"{:.1e}".format( myoutput.sd_beta[5])
        
    if P == 11:
        TextP = "\n" + "A = " + Text[1].format(myoutput.beta[1]) +r"$\, \pm \,$" +"{:.1e}".format( myoutput.sd_beta[1]) + "\n"+ "B = " + Text[2].format(myoutput.beta[2]*2) +r"$\, \pm \,$" +"{:.1e}".format( myoutput.sd_beta[2]*2)
        TextP = "\n" + "A = " + format_string(myoutput.beta[1],myoutput.sd_beta[1],Num[1]) + "\n"+ "B = " + format_string(myoutput.beta[2]*2,myoutput.sd_beta[2]*2,Num[2])
        
    if P == 12: 
        TextP = "\n" + "A = " + Text[1].format(myoutput.beta[1]) +r"$\, \pm \,$" +"{:.1e}".format( myoutput.sd_beta[1]) + "\n"+ "B = " + Text[2].format(myoutput.beta[2]/2) +r"$\, \pm \,$" +"{:.1e}".format( myoutput.sd_beta[2]/2)+ "\n"+ "C = " + Text[3].format(myoutput.beta[3]/3) +r"$\, \pm \,$" +"{:.1e}".format( myoutput.sd_beta[3]/3)

        TextP = "\n" + "A = " + format_string(myoutput.beta[1],myoutput.sd_beta[1],Num[1]) + "\n"+ "B = " + format_string(myoutput.beta[2]*2,myoutput.sd_beta[2]*2,Num[2]) + "\n"+ "C = " + format_string(myoutput.beta[3]*3,myoutput.sd_beta[3]*3,Num[3])

        
    if P == 13:
        TextP = "\n" + "A = " + Text[1].format(myoutput.beta[1]) +r"$\, \pm \,$" +"{:.1e}".format( myoutput.sd_beta[1]) + "\n"+ "B = " + Text[2].format(myoutput.beta[2]*2) +r"$\, \pm \,$" +"{:.1e}".format( myoutput.sd_beta[2]*2)+ "\n"+ "C = " + Text[3].format(myoutput.beta[3]*3) +r"$\, \pm \,$" +"{:.1e}".format( myoutput.sd_beta[3]*3) + "\n"+ "D = " + Text[4].format(myoutput.beta[4]*4) +r"$\, \pm \,$" +"{:.1e}".format( myoutput.sd_beta[4]*4)
        TextP = "\n" + "A = " + format_string(myoutput.beta[1],myoutput.sd_beta[1],Num[1]) + "\n"+ "B = " + format_string(myoutput.beta[2]*2,myoutput.sd_beta[2]*2,Num[2]) + "\n"+ "C = " + format_string(myoutput.beta[3]*3,myoutput.sd_beta[3]*3,Num[3])+ "\n"+ "D = " + format_string(myoutput.beta[4]*4,myoutput.sd_beta[4]*4,Num[4])
    return TextP

def Nice_Plot(Dat1,Dat2,myoutput,Chi2,Pol,Text1,Text2,Text3,Name,Titel,AT,Text_Fit = None,XLim = False, YLim = False):
    
    NLegendSize = 10
    Fak = 1
    if np.abs(np.mean(Dat1[1])) < 10**-2:
        Fak = 10**6
        
    Fakx = 1
    if np.mean(Dat1[0]) < 10**-2:
        Fakx = 10**6
     
    TextP = ""
    if np.min(Dat1[0]) > 0:
        x_line = np.linspace(np.min(Dat1[0])*0.98 , np.max(Dat1[0])*1.02 , 200)
    if np.min(Dat1[0]) <= 0:
        x_line = np.linspace(np.min(Dat1[0])*1.1 , np.max(Dat1[0])*1.02 , 200)
    
    if Pol >0:    
        pl1 = Spesific_func_confidence_band(x_line,95.4,myoutput,Pol) #68.2 95.4
        print("Max Err: ", np.max(pl1), " mean Error: ", np.mean(pl1))
    #pl2 = Spesific_func_confidence_band(Dat1[0],68.2,myoutput,Pol, Prediction_I = True, Yerr = Dat1[3] ) 
    
    Text_CI = r"CI 0.95"
    
    if Pol ==5:
        p1 = Pol5(myoutput.beta,Dat1[0])
        if len(Dat2)>0:
            p2 = Pol5(myoutput.beta,Dat2[0])
        p3 = Pol5(myoutput.beta,x_line)
        p7 = Pol5(myoutput.beta,Dat1[0]+Dat1[2])
        p8 = Pol5(myoutput.beta,Dat1[0]-Dat1[2])
        #TextP = "\n" + "a = " + "{:.2e}".format(myoutput.beta[0]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[0]) + "\n"+ "b = " + "{:.2e}".format(myoutput.beta[1]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[1])+ "\n"+ "c = " + "{:.2e}".format(myoutput.beta[2]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[2]) + "\n"+ "d = " + "{:.2e}".format(myoutput.beta[3]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[3])+ "\n"+ "e = " + "{:.2e}".format(myoutput.beta[4]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[4]) + "\n"+ "f = " + "{:.2e}".format(myoutput.beta[5]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[5]) 
        TextP = Text_legend(myoutput,Pol)
        TextP = TextP + "\n" + r"$\frac{\chi^2}{Ndf} \, = \, $" + "{:.1f}".format(Chi2)
        
    if Pol ==4:
       
        p1 = Pol4(myoutput.beta,Dat1[0])
        if len(Dat2)>0:
            p2 = Pol4(myoutput.beta,Dat2[0])
        p3 = Pol4(myoutput.beta,x_line)
        p7 = Pol4(myoutput.beta,Dat1[0]+Dat1[2])
        p8 = Pol4(myoutput.beta,Dat1[0]-Dat1[2])
        #TextP = "\n" + "a = " + "{:.2e}".format(myoutput.beta[0]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[0]) + "\n"+ "b = " + "{:.2e}".format(myoutput.beta[1]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[1])+ "\n"+ "c = " + "{:.2e}".format(myoutput.beta[2]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[2]) + "\n"+ "d = " + "{:.2e}".format(myoutput.beta[3]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[3])+ "\n"+ "e = " + "{:.2e}".format(myoutput.beta[4]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[4]) 
        TextP = Text_legend(myoutput,Pol)
        TextP = TextP + "\n" + r"$\frac{\chi^2}{Ndf} \, = \, $" + "{:.1f}".format(Chi2)
        
    if Pol ==2:
        Text3 = r"$Fit(x)\,=\,a\,+\,b\cdot x\,+\, c \cdot x^2$"
       
        p1 = Pol2(myoutput.beta,Dat1[0])
        if len(Dat2)>0:
            p2 = Pol2(myoutput.beta,Dat2[0])
        p3 = Pol2(myoutput.beta,x_line)    
        p7 = Pol2(myoutput.beta,Dat1[0]+Dat1[2])
        p8 = Pol2(myoutput.beta,Dat1[0]-Dat1[2])
        #TextP = "\n" + "a = " + "{:.2e}".format(myoutput.beta[0]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[0]) + "\n"+ "b = " + "{:.2e}".format(myoutput.beta[1]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[1])+ "\n"+ "c = " + "{:.2e}".format(myoutput.beta[2]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[2]) 
        TextP = Text_legend(myoutput,Pol)
        TextP = TextP + "\n" + r"$\frac{\chi^2}{Ndf} \, = \, $" + "{:.1f}".format(Chi2)
    if Pol ==1:
        Text3 = r"$Fit(x)\,=\,a\,+\,b\cdot x$"
       
        p1 = Pol1(myoutput.beta,Dat1[0])
        if len(Dat2)>0:
            p2 = Pol1(myoutput.beta,Dat2[0])
        p3 = Pol1(myoutput.beta,x_line)    
        p7 = Pol1(myoutput.beta,Dat1[0]+Dat1[2])
        p8 = Pol1(myoutput.beta,Dat1[0]-Dat1[2])
        
        #TextP = "\n" + "a = " + "{:.2e}".format(myoutput.beta[0]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[0]) + "\n"+ "b = " + "{:.2e}".format(myoutput.beta[1]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[1])
        TextP = Text_legend(myoutput,Pol)
        TextP = TextP + "\n" + r"$\frac{\chi^2}{Ndf} \, = \, $" + "{:.1f}".format(Chi2)
            
        '''
        exponent1 = -math.floor(math.log10(np.abs(myoutput.beta[0]))) 
        exponent2 = -math.floor(math.log10(np.abs(myoutput.beta[0])))
        
        TextP = "\n" +"a = " + str(np.around(myoutput.beta[0]*10**exponent1,4)) + "e" + str(exponent1) + " $\pm$ " +  str(np.around(myoutput.sd_beta[0]*10**exponent1,4)) + "e" + str(exponent1)

        if exponent1 == 1 or exponent1 == -1 :
            TextP = "\n" +"a = " + str(np.around(myoutput.beta[0],4)) + " $\pm$ " +  str(np.around(myoutput.sd_beta[0],4)) 
        TextP2 = "\n"  +"b = " + str(np.around(myoutput.beta[1]*10**exponent2,4)) + "e" + str(exponent2) + " $\pm$ " +  str(np.around(myoutput.sd_beta[1]*10**exponent2,4)) + "e" + str(exponent2)
        if exponent2 == 1 or exponent2 == -1 :
            TextP2 = "\n"  +"b = " + str(np.around(myoutput.beta[1],4)) + " $\pm$ " +  str(np.around(myoutput.sd_beta[1],4)) 
        TextP += TextP2
        '''
    if Pol ==3:
        Text3 = r"$Fit(x)\,=\,a\,+\,b\cdot x \, + \, c\cdot x^2 \, + \, d\cdot x^3$"
        
        p1 = Pol3(myoutput.beta,Dat1[0])
        if len(Dat2)>0:
            p2 = Pol3(myoutput.beta,Dat2[0])
        p3 = Pol3(myoutput.beta,x_line)    
        p7 = Pol3(myoutput.beta,Dat1[0]+Dat1[2])
        p8 = Pol3(myoutput.beta,Dat1[0]-Dat1[2])
        #TextP = "\n" + "a = " + "{:.2e}".format(myoutput.beta[0]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[0]) + "\n"+ "b = " + "{:.2e}".format(myoutput.beta[1]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[1])+ "\n"+ "c = " + "{:.2e}".format(myoutput.beta[2]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[2]) + "\n"+ "d = " + "{:.2e}".format(myoutput.beta[3]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[3]) 
        TextP = Text_legend(myoutput,Pol)
        TextP = TextP + "\n" + r"$\frac{\chi^2}{Ndf} \, = \, $" + "{:.1f}".format(Chi2)
        
    if Pol ==0:
        pl1 = myoutput[1]#  np.std(Dat1[1])
        p1 = np.zeros(len(Dat1[0])) + myoutput[0]
        if len(Dat2)>0:
            p2 = np.zeros(len(Dat2[0])) + np.mean(Dat2[1])
        p3 = np.zeros(len(x_line)) + myoutput[0]
        p7 = 0#Pol1(myoutput.beta,Dat1[0]+Dat1[2])
        p8 =0# Pol1(myoutput.beta,Dat1[0]-Dat1[2])
        #Ne = np.floor(np.log10(np.abs(myoutput[0]))).astype(int) - np.floor(np.log10(np.abs(myoutput[1]))).astype(int)
        #NT = "{:."+str(Ne+1)+"e}"
        TextP =  Text_legend(myoutput,Pol)# "\n" +"Weighted Mean = " + NT.format(myoutput[0]) + " $\pm$ " +  "{:.2e}".format(myoutput[1]) 
        Text_CI = r"$1\sigma $"
        
        
    if Pol ==11:
        
        p1 = KPol3(myoutput.beta,Dat1[0])
        if len(Dat2)>0:
            p2 = KPol3(myoutput.beta,Dat2[0])
        p3 = KPol3(myoutput.beta,x_line)
        p7 = KPol3(myoutput.beta,Dat1[0]+Dat1[2])
        p8 = KPol3(myoutput.beta,Dat1[0]-Dat1[2])
        
        #TextP = "\n" + "A = " + "{:.2e}".format(myoutput.beta[1]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[1]) + "\n"+ "B = " + "{:.2e}".format(myoutput.beta[2]/2) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[2]/2)
        TextP = Text_legend(myoutput,Pol)
        TextP = TextP + "\n" + r"$\frac{\chi^2}{Ndf} \, = \, $" + "{:.1f}".format(Chi2)
        

    if Pol ==12:
        p1 = KPol4(myoutput.beta,Dat1[0])
        if len(Dat2)>0:
            p2 = KPol4(myoutput.beta,Dat2[0])
        p3 = KPol4(myoutput.beta,x_line)
        p7 = KPol4(myoutput.beta,Dat1[0]+Dat1[2])
        p8 = KPol4(myoutput.beta,Dat1[0]-Dat1[2])
        #TextP = "\n" + "A = " + "{:.2e}".format(myoutput.beta[1]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[1]) + "\n"+ "B = " + "{:.2e}".format(myoutput.beta[2]/2) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[2]/2)+ "\n"+ "C = " + "{:.2e}".format(myoutput.beta[3]/3) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[3]/3)
        TextP = Text_legend(myoutput,Pol)
        TextP = TextP + "\n" + r"$\frac{\chi^2}{Ndf} \, = \, $" + "{:.1f}".format(Chi2)

        
    if Pol ==13:
        p1 = KPol5(myoutput.beta,Dat1[0])
        if len(Dat2)>0:
            p2 = KPol5(myoutput.beta,Dat2[0])
        p3 = KPol5(myoutput.beta,x_line)
        p7 = KPol5(myoutput.beta,Dat1[0]+Dat1[2])
        p8 = KPol5(myoutput.beta,Dat1[0]-Dat1[2])
        #TextP = "\n" + "A = " + "{:.2e}".format(myoutput.beta[1]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[1]) + "\n"+ "B = " + "{:.2e}".format(myoutput.beta[2]/2) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[2]/2)+ "\n"+ "C = " + "{:.2e}".format(myoutput.beta[3]/3) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[3]/3) + "\n"+ "D = " + "{:.2e}".format(myoutput.beta[4]/4) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[4]/4)
        TextP = Text_legend(myoutput,Pol)
        TextP = TextP + "\n" + r"$\frac{\chi^2}{Ndf} \, = \, $" + "{:.1f}".format(Chi2)
        
        
    if Pol ==14:
        
        p1 = KPolN(myoutput.beta,Dat1[0])
        if len(Dat2)>0:
            p2 = KPolN(myoutput.beta,Dat2[0])
        p3 = KPolN(myoutput.beta,x_line)
        p7 = KPolN(myoutput.beta,Dat1[0]+Dat1[2])
        p8 = KPolN(myoutput.beta,Dat1[0]-Dat1[2])
        
    if Pol ==112:
        
        p1 = KPol3H(myoutput.beta,Dat1[0])
        if len(Dat2)>0:
            p2 = KPol3H(myoutput.beta,Dat2[0])
        p3 = KPol3H(myoutput.beta,x_line)
        p7 = KPol3H(myoutput.beta,Dat1[0]+Dat1[2])
        p8 = KPol3H(myoutput.beta,Dat1[0]-Dat1[2])
        
        TextP = "\n" + "a = " + "{:.2e}".format(myoutput.beta[1]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[1]) + "\n"+ "b = " + "{:.2e}".format(myoutput.beta[2]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[2])+ "\n"+ r"$T_C$ = " + "{:.2e}".format(myoutput.beta[3]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[3]) + "\n"+ "m = " + "{:.2e}".format(myoutput.beta[4]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[4]) 
        TextP = TextP + "\n" + r"$\frac{\chi^2}{Ndf} \, = \, $" + "{:.3f}".format(Chi2)

    if Pol ==122:
        
        p1 = KPol4H(myoutput.beta,Dat1[0])
        if len(Dat2)>0:
            p2 = KPol4H(myoutput.beta,Dat2[0])
        p3 = KPol4H(myoutput.beta,x_line)
        p7 = KPol4H(myoutput.beta,Dat1[0]+Dat1[2])
        p8 = KPol4H(myoutput.beta,Dat1[0]-Dat1[2])
        
        TextP = "\n" + "a = " + "{:.2e}".format(myoutput.beta[1]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[1]) + "\n"+ "b = " + "{:.2e}".format(myoutput.beta[2]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[2])+ "\n"+ "c = " + "{:.2e}".format(myoutput.beta[3]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[3]) + "\n"+ r"$T_C$ = " + "{:.2e}".format(myoutput.beta[4]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[4]) + "\n"+ "m = " + "{:.2e}".format(myoutput.beta[5]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[5])
        TextP = TextP + "\n" + r"$\frac{\chi^2}{Ndf} \, = \, $" + "{:.3f}".format(Chi2)
    
    if Text_Fit != None:
        Text3 = Text_Fit

    fig, axarray = plt.subplots(2, 1, figsize=(7,5), sharex=True, gridspec_kw={'height_ratios': [5, 2]}, dpi = 150)
    axarray[0].set_title(Titel)
    # Grafische Darstellung der Rohdaten
    axarray[0].errorbar(Dat1[0]*Fakx, Dat1[1]*Fak, xerr=Dat1[2], yerr=Dat1[3]*Fak, fmt='.', marker='.', label = Text1)
    
    if len(Dat2)>0:
        axarray[0].errorbar(Dat2[0]*Fakx, Dat2[1]*Fak, xerr=Dat2[2], yerr=Dat2[3]*Fak, color = "orange", fmt='.', marker='o',  markeredgecolor ="orange", label = Text2)
    
   # axarray[0].set_xlabel('$\Delta v$ / 1/GHz')
    #axarray[0].set_ylabel('T extern / °C')
    #axarray[0].set_xlabel(AT[0], fontsize = NLegendSize)
    axarray[0].set_ylabel(AT[1], fontsize = NLegendSize)
 #   axarray[0].set_ylabel('Strain DMS $10^{-6}$')
    axarray[0].grid()
    
    
    #sigR2  = ( np.abs(Dat1[1] - p7 ) + np.abs( Dat1[1] - p8 ) )/2
    sigR2  = np.abs( np.abs(p7 ) - np.abs( p8 ) )/2
    sigmaRes = np.sqrt(Dat1[3]**2 + sigR2**2)*Fak
    
    #axarray[0].plot(Y_A_V, Pol4(myoutput.beta,Y_A_V), color='green')
    axarray[0].plot(x_line*Fakx, p3*Fak, '--', color='red', label = Text3 + TextP)
    
    axarray[0].plot(x_line*Fakx, (p3 + pl1)*Fak, ':', color='red', alpha = 0.5,label = Text_CI) #r"$1 \sigma$")
    axarray[0].plot(x_line*Fakx, (p3 - pl1)*Fak, ':', color='red', alpha = 0.5)
    
    #axarray[0].plot(Dat1[0], (p1 + pl2)*Fak, ':', color='green', alpha = 0.5,label = r"PI 0.68") #r"$1 \sigma$")
    #axarray[0].plot(Dat1[0], (p1 - pl2)*Fak, ':', color='green', alpha = 0.5)
    
    #axarray[0].plot( x_line,  y_line4, '--', color='blue', label = "Pol. 4 Fit")
    axarray[0].legend(prop={'size':NLegendSize})   
    axarray[0].legend(fontsize = NLegendSize)  
    # Zunächst plotten wir eine gestrichelte Nulllinie, dann den eigentlichen Residuenplot:
    axarray[1].axhline(y=0., color='black', linestyle='--')
    axarray[1].errorbar(Dat1[0]*Fakx,( Dat1[1] - p1)*Fak, yerr=sigmaRes, fmt='.', marker='.')
    
    #axarray[1].plot(x_line, (pl1)*Fak, ':', color='red', alpha = 0.5)  # 1 sigma markierung im residuen plot
    #axarray[1].plot(x_line, (-pl1)*Fak, ':', color='red', alpha = 0.5)
    
    if len(Dat2)>0:
        axarray[1].errorbar(Dat2[0]*Fakx, (Dat2[1] - p2)*Fak , color = "orange" , yerr=sigmaRes[0:len(Dat2[0])], fmt='.', marker='.',  markeredgecolor ="orange")
    #axarray[1].errorbar(Y_A_V,Y_A_T - Pol4(par_F,np.sort(Y_A_V)), yerr=sigmaRes, color='red', fmt='.', marker='o', markeredgecolor='red')
    
    #axarray[1].set_xlabel('$\Delta v$ / 1/GHz')
    axarray[1].set_xlabel(AT[0],fontsize = NLegendSize)
   # axarray[1].set_ylabel('$(T-Fit(\Delta v)$ / °C')
    axarray[1].set_ylabel(AT[2],fontsize = NLegendSize)
    
    #axarray[1].set_ylabel('$(Strain-Fit(\Delta v)$ $10^{-6}$')
    
    # Wir sorgen dafür, dass die y-Achse beim Residuenplot symmetrisch um die Nulllinie ist:
    ymax = max([abs(x) for x in axarray[1].get_ylim()])
    axarray[1].set_ylim(-ymax, ymax)
    #axarray[1].set_ylim(-0.0002,0.0002)
    if XLim != False:
        axarray[0].set_xlim(XLim[0],XLim[1])
        axarray[1].set_xlim(XLim[0],XLim[1])
    if YLim != False:
        axarray[0].set_ylim(YLim[0],YLim[1])
        
    axarray[1].grid()
    
    plt.tight_layout()
    fig.subplots_adjust(hspace=0.050)
   # plt.tight_layout() 
    plt.savefig(Name + ".eps",dpi = 1200)
    plt.savefig(Name + ".pdf",dpi = 1200)
    plt.show()
    
    
def Nice_Plot_2YA(Dat1,Dat2,myoutput,Pol,Text1,Text2,Text3,Name,Titel,AT, YFaktor):
    Fak = 1# 10**6
    if np.min(Dat1[0]) > 0:
        x_line = np.linspace(np.min(Dat1[0])*0.1 , np.max(Dat1[0])*1.05 , 200)
    if np.min(Dat1[0]) <= 0:
        x_line = np.linspace(np.min(Dat1[0])*1.1 , np.max(Dat1[0])*1.05 , 200)
    if Pol ==4:
        p1 = Pol4(myoutput.beta,Dat1[0])
        if len(Dat2)>0:
            p2 = Pol4(myoutput.beta,Dat2[0])
        p3 = Pol4(myoutput.beta,x_line)
        p7 = Pol4(myoutput.beta,Dat1[0]+Dat1[2])
        p8 = Pol4(myoutput.beta,Dat1[0]-Dat1[2])
    if Pol ==2:
        p1 = Pol2(myoutput.beta,Dat1[0])
        if len(Dat2)>0:
            p2 = Pol2(myoutput.beta,Dat2[0])
        p3 = Pol2(myoutput.beta,x_line)    
        p7 = Pol2(myoutput.beta,Dat1[0]+Dat1[2])
        p8 = Pol2(myoutput.beta,Dat1[0]-Dat1[2])
    if Pol ==1:
        p1 = Pol1(myoutput.beta,Dat1[0])
        if len(Dat2)>0:
            p2 = Pol1(myoutput.beta,Dat2[0])
        p3 = Pol1(myoutput.beta,x_line)    
        p7 = Pol1(myoutput.beta,Dat1[0]+Dat1[2])
        p8 = Pol1(myoutput.beta,Dat1[0]-Dat1[2])
    if Pol ==3:
        p1 = Pol3(myoutput.beta,Dat1[0])
        if len(Dat2)>0:
            p2 = Pol3(myoutput.beta,Dat2[0])
        p3 = Pol3(myoutput.beta,x_line)    
        p7 = Pol3(myoutput.beta,Dat1[0]+Dat1[2])
        p8 = Pol3(myoutput.beta,Dat1[0]-Dat1[2])
    

    fig, axarray = plt.subplots(2, 1, figsize=(7,5), sharex=True, gridspec_kw={'height_ratios': [5, 2]}, dpi = 150)
    axarray[0].set_title(Titel)
    # Grafische Darstellung der Rohdaten
    axarray[0].errorbar(Dat1[0], Dat1[1]*Fak, xerr=Dat1[2], yerr=Dat1[3]*Fak, fmt='.', marker='o', label = Text1)
    
    if len(Dat2)>0:
        axarray[0].errorbar(Dat2[0], Dat2[1]*Fak, xerr=Dat2[2], yerr=Dat2[3]*Fak, color = "orange", fmt='.', marker='o',  markeredgecolor ="orange", label = Text2)
    
   # axarray[0].set_xlabel('$\Delta v$ / 1/GHz')
    #axarray[0].set_ylabel('T extern / °C')
    axarray[0].set_xlabel(AT[0])
    axarray[0].set_ylabel(AT[1])
 #   axarray[0].set_ylabel('Strain DMS $10^{-6}$')
    axarray[0].grid()
    
    
    #sigR2  = ( np.abs(Dat1[1] - p7 ) + np.abs( Dat1[1] - p8 ) )/2
    sigR2  = np.abs( np.abs(p7 ) - np.abs( p8 ) )/2
    sigmaRes = np.sqrt(Dat1[3]**2 + sigR2**2)*Fak
    #axarray[0].plot(Y_A_V, Pol4(myoutput.beta,Y_A_V), color='green')
    axarray[0].plot(x_line, p3*Fak, '--', color='red', label = Text3)
    #axarray[0].plot( x_line,  y_line4, '--', color='blue', label = "Pol. 4 Fit")
    axarray[0].legend()
    
    if len(AT) > 3:
        ymax = max([abs(x) for x in axarray[0].get_ylim()])
        ymin = min([abs(x) for x in axarray[0].get_ylim()])
        ax21 = axarray[0].twinx()    
        ax21.set_ylabel(AT[3])
        ax21.set_ylim(ymin*YFaktor,ymax*YFaktor)
        #ax21.plot(x_line, p3*Fak*YFaktor, '--', color='red')
        #x = np.array(np.linspace(int(ymin),int(ymax),int(ymax)+1))
        #print(x)
        #values = x * YFaktor
        #ax21.set_yticks(x,values)
        #ax21.grid()
        
    # Zunächst plotten wir eine gestrichelte Nulllinie, dann den eigentlichen Residuenplot:
    axarray[1].axhline(y=0., color='black', linestyle='--')
    axarray[1].errorbar(Dat1[0],( Dat1[1] - p1)*Fak, yerr=sigmaRes, fmt='.', marker='o')
    if len(Dat2)>0:
        axarray[1].errorbar(Dat2[0], (Dat2[1] - p2)*Fak , color = "orange" , yerr=sigmaRes[0:len(Dat2[0])], fmt='.', marker='o',  markeredgecolor ="orange")
    #axarray[1].errorbar(Y_A_V,Y_A_T - Pol4(par_F,np.sort(Y_A_V)), yerr=sigmaRes, color='red', fmt='.', marker='o', markeredgecolor='red')
    
    #axarray[1].set_xlabel('$\Delta v$ / 1/GHz')
    axarray[1].set_xlabel(AT[0])
   # axarray[1].set_ylabel('$(T-Fit(\Delta v)$ / °C')
    axarray[1].set_ylabel(AT[2])
    
    #axarray[1].set_ylabel('$(Strain-Fit(\Delta v)$ $10^{-6}$')
    
    # Wir sorgen dafür, dass die y-Achse beim Residuenplot symmetrisch um die Nulllinie ist:
    ymax = max([abs(x) for x in axarray[1].get_ylim()])
    axarray[1].set_ylim(-ymax, ymax)
    #axarray[1].set_ylim(-0.0002,0.0002)
    axarray[1].grid()
    
    plt.tight_layout()
    fig.subplots_adjust(hspace=0.0)
   # plt.tight_layout() 
    plt.savefig(Name + ".png")
    plt.show()   
    

        

def Nice_Plot_w_Lit(Dat1,Dat2,Dat3, Dat3_P,myoutput,Chi2,Pol,Text1,Text2,Text3,Text4, Text5,Name,Titel,AT, ThirdSigPlot = False):

    Fak = 1
    if np.mean(Dat1[1]) < 10**-2:
        Fak = 10**6
        
    NLegendSize = 10
    if np.min(Dat1[0]) > 0:
        x_line = np.linspace(np.min(Dat1[0])*0.98 , np.max(Dat1[0])*1.02 , 200)
    if np.min(Dat1[0]) <= 0:
        x_line = np.linspace(np.min(Dat1[0])*1.1 , np.max(Dat1[0])*1.02 , 200)
        
    pl1 = Spesific_func_confidence_band(x_line,95.4,myoutput,Pol) #68.2 95.4
    
    if Pol ==5:
        p1 = Pol5(myoutput.beta,Dat1[0])
        if len(Dat2)>0:
            p2 = Pol5(myoutput.beta,Dat2[0])
        p3 = Pol5(myoutput.beta,x_line)
        p7 = Pol5(myoutput.beta,Dat1[0]+Dat1[2])
        p8 = Pol5(myoutput.beta,Dat1[0]-Dat1[2])
        #TextP = "\n" + "a = " + "{:.2e}".format(myoutput.beta[0]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[0]) + "\n"+ "b = " + "{:.2e}".format(myoutput.beta[1]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[1])+ "\n"+ "c = " + "{:.2e}".format(myoutput.beta[2]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[2]) + "\n"+ "d = " + "{:.2e}".format(myoutput.beta[3]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[3])+ "\n"+ "e = " + "{:.2e}".format(myoutput.beta[4]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[4]) + "\n"+ "f = " + "{:.2e}".format(myoutput.beta[5]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[5])
        TextP = Text_legend(myoutput,Pol)
        TextP = TextP + "\n" + r"$\frac{\chi^2}{Ndf} \, = \, $" + "{:.1f}".format(Chi2)
        
        
    if Pol ==4:
        p1 = Pol4(myoutput.beta,Dat1[0])
        if len(Dat2)>0:
            p2 = Pol4(myoutput.beta,Dat2[0])
        p3 = Pol4(myoutput.beta,x_line)
        p7 = Pol4(myoutput.beta,Dat1[0]+Dat1[2])
        p8 = Pol4(myoutput.beta,Dat1[0]-Dat1[2])
        #TextP = "\n" + "a = " + "{:.2e}".format(myoutput.beta[0]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[0]) + "\n"+ "b = " + "{:.2e}".format(myoutput.beta[1]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[1])+ "\n"+ "c = " + "{:.2e}".format(myoutput.beta[2]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[2]) + "\n"+ "d = " + "{:.2e}".format(myoutput.beta[3]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[3])+ "\n"+ "e = " + "{:.2e}".format(myoutput.beta[4]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[4])
        TextP = Text_legend(myoutput,Pol)
        TextP = TextP + "\n" + r"$\frac{\chi^2}{Ndf} \, = \, $" + "{:.1f}".format(Chi2)
        
    if Pol ==2:
        p1 = Pol2(myoutput.beta,Dat1[0])
        if len(Dat2)>0:
            p2 = Pol2(myoutput.beta,Dat2[0])
        p3 = Pol2(myoutput.beta,x_line)    
        p7 = Pol2(myoutput.beta,Dat1[0]+Dat1[2])
        p8 = Pol2(myoutput.beta,Dat1[0]-Dat1[2])
        #TextP = "\n" + "a = " + "{:.2e}".format(myoutput.beta[0]) +r"$\, \pm \,$" +"{:.1e}".format( myoutput.sd_beta[0]) + "\n"+ "b = " + "{:.2e}".format(myoutput.beta[1]) +r"$\, \pm \,$" +"{:.1e}".format( myoutput.sd_beta[1])+ "\n"+ "c = " + "{:.2e}".format(myoutput.beta[2]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[2]) 
        TextP = Text_legend(myoutput,Pol)
        TextP = TextP + "\n" + r"$\frac{\chi^2}{Ndf} \, = \, $" + "{:.1f}".format(Chi2)
        
    if Pol ==1:
        p1 = Pol1(myoutput.beta,Dat1[0])
        if len(Dat2)>0:
            p2 = Pol1(myoutput.beta,Dat2[0])
        p3 = Pol1(myoutput.beta,x_line)    
        p7 = Pol1(myoutput.beta,Dat1[0]+Dat1[2])
        p8 = Pol1(myoutput.beta,Dat1[0]-Dat1[2])
        #TextP = "\n" + "a = " + "{:.2e}".format(myoutput.beta[0]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[0]) + "\n"+ "b = " + "{:.2e}".format(myoutput.beta[1]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[1])
        TextP = Text_legend(myoutput,Pol)
        TextP = TextP + "\n" + r"$\frac{\chi^2}{Ndf} \, = \, $" + "{:.1f}".format(Chi2)
        
    if Pol ==3:
        p1 = Pol3(myoutput.beta,Dat1[0])
        if len(Dat2)>0:
            p2 = Pol3(myoutput.beta,Dat2[0])
        p3 = Pol3(myoutput.beta,x_line)    
        p7 = Pol3(myoutput.beta,Dat1[0]+Dat1[2])
        p8 = Pol3(myoutput.beta,Dat1[0]-Dat1[2])
        
        #TextP = "\n" + "a = " + "{:.2e}".format(myoutput.beta[0]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[0]) + "\n"+ "b = " + "{:.2e}".format(myoutput.beta[1]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[1])+ "\n"+ "c = " + "{:.2e}".format(myoutput.beta[2]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[2]) + "\n"+ "d = " + "{:.2e}".format(myoutput.beta[3]) +r"$\, \pm \,$" +"{:.2e}".format( myoutput.sd_beta[3])
        TextP = Text_legend(myoutput,Pol)
        TextP = TextP + "\n" + r"$\frac{\chi^2}{Ndf} \, = \, $" + "{:.1f}".format(Chi2)
     
    NPlot = 2  
    NPlotSice = [5,2]
    if ThirdSigPlot == True and len(Dat3) == 4 :
        NPlot = 3
        NPlotSice = [5,2,2]

    fig, axarray = plt.subplots(NPlot, 1, figsize=(7,5), sharex=True, gridspec_kw={'height_ratios': NPlotSice}, dpi = 150)
   
    axarray[0].set_title(Titel)
    # Grafische Darstellung der Rohdaten
    axarray[0].errorbar(Dat1[0], Dat1[1]*Fak, xerr=Dat1[2], yerr=Dat1[3]*Fak, fmt='.', marker='.', label = Text1)
    
    if len(Dat2)>0:
        axarray[0].errorbar(Dat2[0], Dat2[1]*Fak, xerr=Dat2[2], yerr=Dat2[3]*Fak, color = "orange", fmt='.', marker='.',  markeredgecolor ="orange", label = Text2)
    if len(Dat3_P)>0:
        axarray[0].errorbar(Dat1[0][Dat3_P[0]:Dat3_P[1]] , Dat1[1][Dat3_P[0]:Dat3_P[1]]*Fak, xerr=Dat1[2][Dat3_P[0]:Dat3_P[1]], yerr=Dat1[3][Dat3_P[0]:Dat3_P[1]]*Fak, color = "black", fmt='.', marker='.',  markeredgecolor ="black", label = Text5)
    
   # axarray[0].set_xlabel('$\Delta v$ / 1/GHz')
    #axarray[0].set_ylabel('T extern / °C')
    #axarray[0].set_xlabel(AT[0])
    axarray[0].set_ylabel(AT[1],fontsize = NLegendSize)
 #   axarray[0].set_ylabel('Strain DMS $10^{-6}$')
    axarray[0].grid()
    
    
    #sigR2  = ( np.abs(Dat1[1] - p7 ) + np.abs( Dat1[1] - p8 ) )/2
    sigR2  = np.abs( np.abs(p7 ) - np.abs( p8 ) )/2
    sigmaRes = np.sqrt(Dat1[3]**2 + sigR2**2)*Fak
    
    
    #axarray[0].plot(Y_A_V, Pol4(myoutput.beta,Y_A_V), color='green')
    axarray[0].plot(x_line, p3*Fak, '--', color='red', label = Text3 + TextP)
    axarray[0].plot(x_line, (p3 + pl1)*Fak, ':', color='red', alpha = 0.5,label = r"CI 0.95") #r"$1 \sigma$")
    axarray[0].plot(x_line, (p3 - pl1)*Fak, ':', color='red', alpha = 0.5)
    
    if len(Dat3)>0:
        axarray[0].plot(Dat3[0], Dat3[1], color='green', label = Text4)
        if len(Dat3)>=3:
            axarray[0].fill_between(Dat3[0], (1-Dat3[2])*Dat3[1],(1+Dat3[2])*Dat3[1] , alpha = 0.5, color = "g")

    #axarray[0].plot( x_line,  y_line4, '--', color='blue', label = "Pol. 4 Fit")
    axarray[0].legend(fontsize = NLegendSize)
    # Zunächst plotten wir eine gestrichelte Nulllinie, dann den eigentlichen Residuenplot:
    axarray[1].axhline(y=0., color='black', linestyle='--')
    axarray[1].errorbar(Dat1[0],( Dat1[1] - p1)*Fak, yerr=sigmaRes, fmt='.', marker='.')
    if len(Dat2)>0:
        axarray[1].errorbar(Dat2[0], (Dat2[1] - p2)*Fak , color = "orange" , yerr=sigmaRes[0:len(Dat2[0])], fmt='.', marker='.',  markeredgecolor ="orange")
        
    if len(Dat3_P)>0:
        axarray[1].errorbar(Dat1[0][Dat3_P[0]:Dat3_P[1]], (Dat1[1][Dat3_P[0]:Dat3_P[1]] - p1[Dat3_P[0]:Dat3_P[1]] )*Fak , color = "black" , yerr=sigmaRes[Dat3_P[0]:Dat3_P[1]], fmt='.', marker='.',  markeredgecolor ="black")
 
    #axarray[1].errorbar(Y_A_V,Y_A_T - Pol4(par_F,np.sort(Y_A_V)), yerr=sigmaRes, color='red', fmt='.', marker='o', markeredgecolor='red')
    
    if ThirdSigPlot == False:
        axarray[1].set_xlabel(AT[0],fontsize = NLegendSize)
    axarray[1].set_ylabel(AT[2],fontsize = NLegendSize)
    
    #axarray[1].set_ylabel('$(Strain-Fit(\Delta v)$ $10^{-6}$')
    
    # Wir sorgen dafür, dass die y-Achse beim Residuenplot symmetrisch um die Nulllinie ist:
    ymax = max([abs(x) for x in axarray[1].get_ylim()])
    axarray[1].set_ylim(-ymax, ymax)
    #axarray[1].set_ylim(-0.0002,0.0002)
    axarray[1].grid()
    
    
    if ThirdSigPlot == True  and len(Dat3) == 4:
        
        axarray[2].axhline(y=0., color='black', linestyle='--')
        axarray[2].set_ylabel(r"$t$-Test" + "\n" +  "to Lit.",fontsize = NLegendSize)    #r"$\sigma\,-\,$diviation" + "\n" + "to Lit.")
        axarray[2].set_xlabel(AT[0],fontsize = NLegendSize)
        axarray[2].grid()
                
        Sigdiv = (Dat3[3] - Dat1[1])/ np.sqrt( (Dat3[2]*Dat3[3])**2 + (sigmaRes/Fak)**2)
        axarray[2].scatter(Dat1[0],Sigdiv,marker='.')
        ymax = max([abs(x) for x in axarray[2].get_ylim()]) + 1
        axarray[2].set_ylim(-ymax, ymax)
        
    plt.rc('font', size=12) 
    plt.rc('axes', labelsize=12)
    plt.tight_layout()
    fig.subplots_adjust(hspace=0.050)
   # plt.tight_layout() 
    plt.savefig(Name + ".png")
    plt.show()
    
    Errest = Spesific_func_confidence_band(x_line,68.2,myoutput,Pol)
    print("Mean Fit Error: ",np.mean(Errest), " maximum Fit Error: ", np.max(Errest), " Mean Fit error in %: ", np.mean(Errest/p3), " max % err: ",np.max(Errest/p3) )
 
#Calculating the Gaussian PDF values given Gaussian parameters and random variable X
def gaus(X,C,X_mean,sigma):
    return C*np.exp(-(X-X_mean)**2/(2*sigma**2))  
    
def Hist_Plot(x_data,NameFig = "Hist_All", NBin = False,xlabel = "",ylabel = "Probability"):
    NBins = int(max(x_data))
    if NBin != False:
        NBins = NBin
    hist, bin_edges = np.histogram(x_data,bins=NBins)
    hist=hist/sum(hist)
    
    n = len(hist)
    x_hist=np.zeros((n),dtype=float) 
    for ii in range(n):
        x_hist[ii]=(bin_edges[ii+1]+bin_edges[ii])/2
        
    y_hist=hist
           
    
    mean = sum(x_hist*y_hist)/sum(y_hist)                  
    sigma = sum(y_hist*(x_hist-mean)**2)/sum(y_hist) 
    
    #Gaussian least-square fitting process
    param_optimised,param_covariance_matrix = curve_fit(gaus,x_hist,y_hist,p0=[max(y_hist),mean,sigma],maxfev=5000)
    
    #print fit Gaussian parameters
    print("Fit parameters: ")
    print("=====================================================")
    print("C = ", param_optimised[0], "+-",np.sqrt(param_covariance_matrix[0,0]))
    print("X_mean =", param_optimised[1], "+-",np.sqrt(param_covariance_matrix[1,1]))
    print("sigma = ", param_optimised[2], "+-",np.sqrt(param_covariance_matrix[2,2]))
    print("\n")
    
    Text_Label = 'Gaussian fit' + "\n" + r"$Fit(x)\,=\,C\cdot e^{(-(x - \mu)^2/(2 \sigma^2))}$" 
    Text_Label  = Text_Label + "\n" + "c = " + "{:.3f}".format(param_optimised[0]) +r"$\, \pm \,$" +"{:.3f}".format(np.sqrt(param_covariance_matrix[0,0])) + "\n"+ r"$\mu$ = " + "{:.2f}".format(param_optimised[1]) +r"$\, \pm \,$" +"{:.2f}".format( np.sqrt(param_covariance_matrix[1,1]))+ "\n"+ r"$ \sigma$ = " + "{:.2f}".format(param_optimised[2]) +r"$\, \pm \,$" +"{:.2f}".format( np.sqrt(param_covariance_matrix[2,2])) 
    #Text_Label = ""
    
    #STEP 4: PLOTTING THE GAUSSIAN CURVE -----------------------------------------
    plt.figure(123232, dpi = 180)
    x_hist_2=np.linspace(np.min(x_hist),np.max(x_hist),500)
    plt.plot(x_hist_2,gaus(x_hist_2,*param_optimised),'r',label=Text_Label,linewidth = 5)
    plt.legend()
    
    #Normalise the histogram values
    weights = np.ones_like(x_data) / len(x_data)
    
    plt.hist(x_data, weights=weights,bins=NBins)
    
    
    #setting the label,title and grid of the plot
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid("on")
    plt.savefig(NameFig + ".png")
    plt.show()
    
def Fit_w_Plot(Dat1,Dat2,Nump,Text1,Text2,Text3,Name_PNG,Titel,AT,SonderWertfK=[0,0],Text_Fit=None,XLim = False,YLim = False):
    degrees_freedom = len(Dat1[0]) - Nump
    

    if Nump >0 and Nump < 11:
        if Nump ==5 :
            Pol_4 = odr.Model(Pol5)
            popt, cov = curve_fit(objective5, Dat1[0],Dat1[1]) 
            
        if Nump ==4 :
            Pol_4 = odr.Model(Pol4)
            popt, cov = curve_fit(objective4, Dat1[0],Dat1[1])
          
        if Nump ==3 :
            Pol_4 = odr.Model(Pol3)
            popt, cov = curve_fit(objective3, Dat1[0],Dat1[1])
        
        if Nump ==2 :
            Pol_4 = odr.Model(Pol2)
            popt, cov = curve_fit(objective2, Dat1[0],Dat1[1])
         
        if Nump ==1 :
            Pol_4 = odr.Model(Pol1)
            popt, cov = curve_fit(objective1, Dat1[0],Dat1[1])
                                   
        mydata = odr.RealData(Dat1[0],Dat1[1],sx = Dat1[2], sy = Dat1[3])
        myodr = odr.ODR(mydata,Pol_4,beta0 = popt)
        myoutput = myodr.run()
        #myoutput.pprint()
        Chi2 = myoutput.sum_square/(degrees_freedom - 1)
        #print("Chi2/ndf: ",Chi2)

    if Nump >10:
        if Nump ==11 : #Sonderfit
            Pol_4 = odr.Model(KPol3)
            popt, cov = curve_fit(objectiveK3, Dat1[0],Dat1[1],bounds=([SonderWertfK[0]-SonderWertfK[1],-np.inf,-np.inf], [SonderWertfK[0]+SonderWertfK[1],np.inf,np.inf]))
        #    print(popt)
            FixtPara = [0,1,1]
            degrees_freedom = len(Dat1[0]) -2
        
        if Nump == 12:
            Pol_4 = odr.Model(KPol4)
            popt, cov = curve_fit(objectiveK4, Dat1[0],Dat1[1],bounds=([SonderWertfK[0]-SonderWertfK[1],-np.inf,-np.inf,-np.inf], [SonderWertfK[0]+SonderWertfK[1],np.inf,np.inf,np.inf]))
        #    print(popt)
            FixtPara = [0,1,1,1]
            degrees_freedom = len(Dat1[0]) -3
            
        if Nump == 13:
            Pol_4 = odr.Model(KPol5)
            popt, cov = curve_fit(objectiveK5, Dat1[0],Dat1[1],bounds=([SonderWertfK[0]-SonderWertfK[1],-np.inf,-np.inf,-np.inf, -np.inf], [SonderWertfK[0]+SonderWertfK[1],np.inf,np.inf,np.inf, np.inf]))
            #print(popt)
            FixtPara = [0,1,1,1,1]
            degrees_freedom = len(Dat1[0]) -4
            
        if Nump == 14:
            Pol_4 = odr.Model(KPolN)
            popt, cov = curve_fit(objectiveKN, Dat1[0],Dat1[1],bounds=([SonderWertfK[0]-SonderWertfK[1],-np.inf,3,77], [SonderWertfK[0]+SonderWertfK[1],np.inf,4,200]))
            print(popt)
            FixtPara = [0,1,1,1]
            degrees_freedom = len(Dat1[0]) -3
            
        if Nump ==112 : #Sonderfit
            Pol_4 = odr.Model(KPol3H)
            popt, cov = curve_fit(objectiveK3H, Dat1[0],Dat1[1],bounds=([SonderWertfK[0]-SonderWertfK[1],-np.inf,-np.inf,SonderWertfK[2],-np.inf], [SonderWertfK[0]+SonderWertfK[1],np.inf,np.inf,SonderWertfK[3],np.inf]))
            print(popt)
            FixtPara = [0,1,1,0,1]
            
        if Nump ==122 : #Sonderfit
            Pol_4 = odr.Model(KPol4H)
            popt, cov = curve_fit(objectiveK4H, Dat1[0],Dat1[1],bounds=([SonderWertfK[0]-SonderWertfK[1],-np.inf,-np.inf,-np.inf,SonderWertfK[2],-np.inf], [SonderWertfK[0]+SonderWertfK[1],np.inf,np.inf,np.inf,SonderWertfK[3],np.inf]))
            print(popt)
            FixtPara = [0,1,1,1,0,1]
            
        
        mydata = odr.RealData(Dat1[0],Dat1[1],sx = Dat1[2], sy = Dat1[3])
        myodr = odr.ODR(mydata,Pol_4,beta0 = popt,ifixb=FixtPara)
        myoutput = myodr.run()
        #myoutput.pprint()
        Chi2 = myoutput.sum_square/degrees_freedom

       # print("Chi2/ndf: ",Chi2)

    if Nump == 0:
        Weights = 1/Dat1[3]**2
        Mean = np.average(Dat1[1],weights = Weights)
        STD = np.sqrt(1/sum(Weights)) #/np.sqrt(len(Dat1[1])
        myoutput = [Mean,STD]
        Chi2 = 1
 #       print("Mean = " + str(np.around(Mean,4)) + " $\pm$ " +  str(np.around(STD,4)) )
    

    Nice_Plot(Dat1,Dat2,myoutput,Chi2,Nump,Text1,Text2,Text3,Name_PNG,Titel,AT,Text_Fit = Text_Fit,XLim= XLim,YLim = YLim) #arguments: 'Dat1', 'Dat2', 'myoutput', 'Pol', 'Text1', 'Text2', 'Text3', and 'Name', Achsen
    return myoutput



def Fit_w_Plot_w_2YA(Dat1,Dat2,Nump,Text1,Text2,Text3,Name_PNG,Titel,AT,YFaktor):
    degrees_freedom = len(Dat1[0]) - Nump
    

    if Nump ==4 :
        Pol_4 = odr.Model(Pol4)
        popt, cov = curve_fit(objective4, Dat1[0],Dat1[1])
    if Nump ==3 :
        Pol_4 = odr.Model(Pol3)
        popt, cov = curve_fit(objective3, Dat1[0],Dat1[1])
    if Nump ==2 :
        Pol_4 = odr.Model(Pol2)
        popt, cov = curve_fit(objective2, Dat1[0],Dat1[1])
    if Nump ==1 :
        Pol_4 = odr.Model(Pol1)
        popt, cov = curve_fit(objective1, Dat1[0],Dat1[1])
        
    mydata = odr.RealData(Dat1[0],Dat1[1],sx = Dat1[2], sy = Dat1[3])
    myodr = odr.ODR(mydata,Pol_4,beta0 = popt)
    myoutput = myodr.run()
    myoutput.pprint()
    Chi2 = myoutput.sum_square/(degrees_freedom -1)
    print("Chi2/ndf: ",Chi2)
    
    Nice_Plot_2YA(Dat1,Dat2,myoutput,Nump,Text1,Text2,Text3,Name_PNG,Titel,AT,YFaktor) #arguments: 'Dat1', 'Dat2', 'myoutput', 'Pol', 'Text1', 'Text2', 'Text3', and 'Name', Achsen
    return myoutput

def Fit_Poli(Dat1,Nump,SonderWertfK = [0,0]):
    degrees_freedom = len(Dat1[0]) - Nump

    if Nump >0 and Nump < 11:
        if Nump ==4 :
            Pol_4 = odr.Model(Pol4)
            popt, cov = curve_fit(objective4, Dat1[0],Dat1[1])

        if Nump ==3 :
            Pol_4 = odr.Model(Pol3)
            popt, cov = curve_fit(objective3, Dat1[0],Dat1[1])

        if Nump ==2 :
            Pol_4 = odr.Model(Pol2)
            popt, cov = curve_fit(objective2, Dat1[0],Dat1[1])

        if Nump ==1 :
            Pol_4 = odr.Model(Pol1)
            popt, cov = curve_fit(objective1, Dat1[0],Dat1[1])
           
                
        mydata = odr.RealData(Dat1[0],Dat1[1],sx = Dat1[2], sy = Dat1[3])
        myodr = odr.ODR(mydata,Pol_4,beta0 = popt)
        myoutput = myodr.run()
        #myoutput.pprint()
        Chi2 = myoutput.sum_square/(degrees_freedom -1)
        print("Chi2/ndf: ",Chi2)
 
    if Nump >10:
        if Nump ==11 : #Sonderfit
            Pol_4 = odr.Model(KPol3)
            popt, cov = curve_fit(objectiveK3, Dat1[0],Dat1[1],bounds=([SonderWertfK[0]-SonderWertfK[1],-np.inf,-np.inf], [SonderWertfK[0]+SonderWertfK[1],np.inf,np.inf]))
            print(popt)
            FixtPara = [0,1,1]
        
        if Nump == 12:
            Pol_4 = odr.Model(KPol4)
            popt, cov = curve_fit(objectiveK4, Dat1[0],Dat1[1],bounds=([SonderWertfK[0]-SonderWertfK[1],-np.inf,-np.inf,-np.inf], [SonderWertfK[0]+SonderWertfK[1],np.inf,np.inf,np.inf]))
            print(popt)
            FixtPara = [0,1,1,1]
        
        
        if Nump == 13:
            Pol_4 = odr.Model(KPol5)
            popt, cov = curve_fit(objectiveK5, Dat1[0],Dat1[1],bounds=([SonderWertfK[0]-SonderWertfK[1],-np.inf,-np.inf,-np.inf, -np.inf], [SonderWertfK[0]+SonderWertfK[1],np.inf,np.inf,np.inf, np.inf]))
            print(popt)
            FixtPara = [0,1,1,1,1]
            
        if Nump == 14:
            Pol_4 = odr.Model(KPolN)
            popt, cov = curve_fit(objectiveKN, Dat1[0],Dat1[1],bounds=([SonderWertfK[0]-SonderWertfK[1],-np.inf,1,77], [SonderWertfK[0]+SonderWertfK[1],np.inf,5,300]))
            print(popt)
            FixtPara = [0,1,1,1]
            
        if Nump ==112 : #Sonderfit
            Pol_4 = odr.Model(KPol3H)
            popt, cov = curve_fit(objectiveK3H, Dat1[0],Dat1[1],bounds=([SonderWertfK[0]-SonderWertfK[1],-np.inf,-np.inf,0,-np.inf], [SonderWertfK[0]+SonderWertfK[1],np.inf,np.inf,300,np.inf]))
            print(popt)
            FixtPara = [0,1,1,0,1]
            
        if Nump ==122 : #Sonderfit
            Pol_4 = odr.Model(KPol4H)
            popt, cov = curve_fit(objectiveK4H, Dat1[0],Dat1[1],bounds=([SonderWertfK[0]-SonderWertfK[1],-np.inf,-np.inf,-np.inf,SonderWertfK[2],-np.inf], [SonderWertfK[0]+SonderWertfK[1],np.inf,np.inf,np.inf,SonderWertfK[3],np.inf]))
            print(popt)
            FixtPara = [0,1,1,1,0,1]
            
        mydata = odr.RealData(Dat1[0],Dat1[1],sx = Dat1[2], sy = Dat1[3])
        myodr = odr.ODR(mydata,Pol_4,beta0 = popt,ifixb=FixtPara)
        myoutput = myodr.run()
        #myoutput.pprint()
        Chi2 = myoutput.sum_square/degrees_freedom
        print("Chi2/ndf: ",Chi2)
        
    if Nump == 0:
        Weights = 1/Dat1[3]**2
        Mean = np.average(Dat1[1],weights = Weights)
        STD = np.sqrt(1/sum(Weights))
        myoutput = [Mean,STD]
        print("Mean = " + str(np.around(Mean,4)) + " $\pm$ " +  str(np.around(STD,4)) )
    
    
    return myoutput

def Plot_Diff_Data(Dat1,Dat2,Titel,Name,AT,GoB,Text1=None,Text2=None,TextBand = "Prediction range", Data2_konti = None, WSig = False,LitLin = None,LitLinLabel = "Standard Model",LitLin2 = None,LitLinLabel2 = "Standard Model"):
    color ="tab:blue"
    if GoB ==1:
        color ="tab:orange"
    Fak2 = 1
    if np.mean(Dat1[1]) < 10**-2:
        Fak2 = 10**6
        
    Fak = 100
    Diff = (Dat1[1] - Dat2[1])*Fak2
    Diff_err = np.sqrt(Dat1[3]**2 + Dat2[3]**2)*Fak2
    
    DiffP = (Dat1[1] - Dat2[1])/Dat1[1]*Fak/ Diff_err
    DiffP_err = np.sqrt( ((Dat1[3] /Dat1[1] ) + (Dat1[1] - Dat2[1])*Dat1[3]/Dat1[1]**2)**2 + (Dat2[3]/Dat1[1])**2)*Fak /Diff_err
    
    DiffP = Dat1[1]/Dat2[1] -1
    DiffP_err = np.sqrt( (Dat1[3]/Dat2[1])**2 + (Dat1[1]*Dat2[3]/Dat2[1]**2)**2 )
    
    DiffP = (Dat1[1] - Dat2[1])/np.sqrt(Dat1[3]**2 + Dat2[3]**2)
    DiffP_err = np.sqrt( ((Dat1[3] /Dat1[1] ) + (Dat1[1] - Dat2[1])*Dat1[3]/Dat1[1]**2)**2 + (Dat2[3]/Dat1[1])**2)*Fak /Diff_err
    
    
    #fig, axarray = plt.subplots(2, 1, figsize=(7,5), sharex=True, gridspec_kw={'height_ratios': [5, 2]}, dpi = 150)
    #axarray[0].set_title(Titel)
    
    # Grafische Darstellung der Rohdaten
    NP = 1
    if Text1 != None and WSig == True:
        fig, axarray = plt.subplots(3, 1, figsize=(7,5), sharex=True, gridspec_kw={'height_ratios': [5, 2, 2]}, dpi = 150)
        axarray[0].set_title(Titel)
        axarray[0].errorbar(Dat1[0], Dat1[1]*Fak2, xerr=Dat1[2], yerr=Dat1[3]*Fak2, fmt='.', marker='.', label = Text1)
        
        Dat21 = Dat2
        if Data2_konti !=None:
            Dat21 = Data2_konti
            
            
        axarray[0].plot(Dat21[0], Dat21[1]*Fak2, color = "black", label = Text2)
        axarray[0].fill_between(Dat21[0], (Dat21[1] + Dat21[3])*Fak2,(Dat21[1] - Dat21[3])*Fak2 , alpha = 0.25, color = "black",label = TextBand)
        axarray[1].axhline(y=0., color='black', linestyle='--')
        axarray[1].fill_between(Dat21[0], Dat21[3]*Fak2,- Dat21[3]*Fak2 , alpha = 0.25, color = "black")


        axarray[1].errorbar(Dat1[0], Diff, xerr=Dat1[2], yerr=Dat1[3]*Fak2, fmt='.', marker='.',color = "tab:blue",elinewidth =1)#, label = Text1)
        #Diff_err,
        

        axarray[0].set_ylabel(AT[1],fontsize = 12)
        axarray[1].set_ylabel(r"Diff.$\times 10^{-6}$",fontsize = 12)
        axarray[1].grid()
        NP = 2
        
    if WSig == False:
        fig, axarray = plt.subplots(2, 1, figsize=(7,5), sharex=True, gridspec_kw={'height_ratios': [5, 2]}, dpi = 150)
        axarray[0].set_title(Titel)
        axarray[0].errorbar(Dat1[0], Dat1[1]*Fak2, xerr=Dat1[2], yerr=Dat1[3]*Fak2, fmt='.', marker='.', label = Text1)
        
        Dat21 = Dat2
        if Data2_konti !=None:
            Dat21 = Data2_konti
            
            
        axarray[0].plot(Dat21[0], Dat21[1]*Fak2, color = "black", label = Text2)
        axarray[1].fill_between(Dat21[0], Dat21[3]*Fak2,- Dat21[3]*Fak2 , alpha = 0.25, color = "black")
        axarray[0].fill_between(Dat21[0], (Dat21[1] + Dat21[3])*Fak2,(Dat21[1] - Dat21[3])*Fak2 , alpha = 0.35, color = "black",label = TextBand)
        axarray[1].axhline(y=0., color='black', linestyle='--')


        axarray[1].errorbar(Dat1[0], Diff, xerr=Dat1[2], yerr=Dat1[3]*Fak2, fmt='.', marker='.',color = "tab:blue",elinewidth =1)#, label = Text1)
        #Diff_err,
        

        axarray[0].set_ylabel(AT[1],fontsize = 12)
        axarray[1].set_xlabel(AT[0],fontsize = 12)
        axarray[1].set_ylabel(r"Diff.$\times 10^{-6}$",fontsize = 12)
        axarray[1].grid()
        NP = 2
    
    else:
        fig, axarray = plt.subplots(2, 1, figsize=(7,5), sharex=True, gridspec_kw={'height_ratios': [5, 2]}, dpi = 150)
        axarray[0].set_title(Titel)
        axarray[0].errorbar(Dat1[0], Dat1[1]*Fak2, xerr=Dat1[2], yerr=Dat1[3]*Fak2, fmt='.', marker='.', label = Text1)

        #axarray[0].errorbar(Dat2[0], Dat2[1]*Fak2, xerr=Dat2[2], yerr=Dat2[3]*Fak2, fmt='.', marker='.',color = "black", label = Text2)
        axarray[0].plot(Dat2[0], Dat2[1]*Fak2, color = "black", label = Text2)
        axarray[0].fill_between(Dat2[0], (Dat2[1] + Dat2[3])*Fak2,(Dat2[1] - Dat2[3])*Fak2 , alpha = 0.5, color = "black")
        axarray[0].set_ylabel(AT[1], size=15)
        
        
    if LitLin != None:
        axarray[0].plot(LitLin[0],LitLin[1],color = "r", label =LitLinLabel )
        
    if LitLin2 != None:
        axarray[0].plot(LitLin2[0],LitLin2[1],color = "g", label =LitLinLabel2 )
        axarray[1].plot(LitLin2[0],LitLin2[1] -  Dat2[1]*Fak2 ,color = "g", label =LitLinLabel2 )

    axarray[0].grid()
    
    if Text1 != None:
        axarray[0].legend(fontsize = 9)     
    
    if WSig == True:
        # Zunächst plotten wir eine gestrichelte Nulllinie, dann den eigentlichen Residuenplot:
        axarray[NP].axhline(y=0., color='black', linestyle='--')
        
        #axarray[1].errorbar(Dat1[0],DiffP,xerr = Dat1[2], yerr=DiffP_err, fmt='.', marker='.')
        axarray[NP].errorbar(Dat1[0],DiffP,xerr = Dat1[2], fmt='.', marker='.', color = "tab:green")
        axarray[NP].set_xlabel(AT[0],fontsize = 15)
        axarray[NP].set_ylabel(r"$\sigma$ deviation",fontsize = 15)# AT[2]) #("t-Test \n to Pre.")#
     
        # Wir sorgen dafür, dass die y-Achse beim Residuenplot symmetrisch um die Nulllinie ist:
        ymax = max([abs(x) for x in axarray[NP].get_ylim()])
        axarray[NP].set_ylim(-ymax, ymax)
        axarray[NP].grid()
    
    plt.tight_layout()
    fig.subplots_adjust(hspace=0.050)
    plt.rc('font', size=12) 
    plt.rc('axes', labelsize=12)
   
    plt.savefig(Name + ".pdf")
    plt.show()
'''    
x = np.linspace(77,300,900)
xerr = np.zeros(900) + 0.1
y = np.zeros(900) + 100 + np.linspace(0,20,900)
yerr = np.zeros(900) + 1
y2 = np.zeros(900) + 10 + np.linspace(0,30,900)

Dat1 = [x,y,xerr,yerr]
Dat2 = [x,y2,xerr,yerr]
Plot_Diff_Data(Dat1,Dat2,"Test","TestDiff",["1","2","3"])
'''



def Sigmoid(x,b):
    return 1/(1+np.exp(-10*(x-b)))
    
def KPol4(B,x):
    return  B[1] * (x - B[0]) + B[2]*(x**2 - B[0]**2) + B[3]*(x**3 - B[0]**3)

def DifKPol4(B,x):
    return  B[1] + B[2]*x*2 + B[3]*x**2 * 3

def KPol5(B,x):
    return  B[1] * (x - B[0]) + B[2]*(x**2 - B[0]**2) + B[3]*(x**3 - B[0]**3) + B[4]*(x**4 - B[0]**4)

def DifKPol5(B,x):
    return  B[1] + B[2]*x*2  + B[3]*x**2 * 3 + B[4]*x**3 * 4
    

def KPolN(B,x):
    return B[1] *(x - B[3])**B[2] - B[1] *(B[0] - B[3])**B[2] 

def KPol3(B,x):
    return  B[1] * (x - B[0]) + B[2]*(x**2 - B[0]**2) 

def DifKPol3(B,x):
    return  B[1] + B[2]*x* 2

def KPol3H(B,x):
    return (B[1] * (x - B[0]) + B[2]*(x**2 - B[0]**2)) *  1/(1+np.exp(-10*(x-B[3]))) + B[4] *  (1-1/(1+np.exp(-10*(x-B[3]))))  

def KPol4H(B,x):
    return (B[1] * (x - B[0]) + B[2]*(x**2 - B[0]**2) + B[3]*(x**3 - B[0]**3)) *  1/(1+np.exp(-10*(x-B[4]))) + B[5] *  (1-1/(1+np.exp(-10*(x-B[4])))) 


def objectiveK3H(x,a,b,c,d,e):
    return (b * (x - a) + c*(x**2 - a**2)) *  1/(1+np.exp(-10*(x-d))) + e *  (1-1/(1+np.exp(-10*(x-d))))


def objectiveK4H(x,a,b,c,d,e,f):
    return (b * (x - a) + c*(x**2 - a**2) + d*(x**3 - a**3)) *  1/(1+np.exp(-10*(x-e))) + f *  (1-1/(1+np.exp(-10*(x-e))))


def objectiveK4(x, a, b, c, d ):
	return  b*(x - a) + c*(x**2 - a**2) + d*(x**3 - a**3)

def objectiveK5(x, a, b, c, d, e ):
	return  b*(x - a) + c*(x**2 - a**2) + d*(x**3 - a**3) + e*(x**4 - a**4)

def objectiveKN(x,T0,a,b,c):
    return a *(x - c)**2 - a *(T0 - c)**2 

def objectiveK3(x, a, b, c ):
	return  b*(x - a) + c*(x**2 - a**2) 

def Fit_w_Plot_w_Lit(Dat1,Dat2,Dat3,Dat3_P,Nump,Text1,Text2,Text3,Text4,Text5,Name_PNG,Titel,AT, ThirdSigPlot = False):
    
    degrees_freedom = len(Dat1[0]) - Nump
    if Nump ==5 :
        Pol_4 = odr.Model(Pol5)
        popt, cov = curve_fit(objective5, Dat1[0],Dat1[1]) 
    if Nump ==4 :
        Pol_4 = odr.Model(Pol4)
        popt, cov = curve_fit(objective4, Dat1[0],Dat1[1])
    if Nump ==3 :
        Pol_4 = odr.Model(Pol3)
        popt, cov = curve_fit(objective3, Dat1[0],Dat1[1])
    if Nump ==2 :
        Pol_4 = odr.Model(Pol2)
        popt, cov = curve_fit(objective2, Dat1[0],Dat1[1])
    if Nump ==1 :
        Pol_4 = odr.Model(Pol1)
        popt, cov = curve_fit(objective1, Dat1[0],Dat1[1])
        
    mydata = odr.RealData(Dat1[0],Dat1[1],sx = Dat1[2], sy = Dat1[3])
    myodr = odr.ODR(mydata,Pol_4,beta0 = popt)
    myoutput = myodr.run()
    myoutput.pprint()
    Chi2 = myoutput.sum_square/degrees_freedom
    print("Chi2/ndf: ",Chi2)
    
    Nice_Plot_w_Lit(Dat1,Dat2,Dat3,Dat3_P,myoutput,Chi2,Nump,Text1,Text2,Text3,Text4,Text5,Name_PNG,Titel,AT,ThirdSigPlot=ThirdSigPlot) #arguments: 'Dat1', 'Dat2', 'myoutput', 'Pol', 'Text1', 'Text2', 'Text3', and 'Name', Achsen
    return myoutput

def R_RTD(B,x): 
    return B[0] * ( 1.0 + B[1]*x + B[2]*x*x + heaviside(-x,0.0)*B[3]*x*x*x*(x-100.0) )

def R_RTD_O(x,a,b,c,d): 
    return a * ( 1.0 + b*x + c*x*x + heaviside(-x,0.0)*d*x*x*x*(x-100.0) )

def Fit_Cali_PT1000(Dat1,Dat2,Nump,Text1,Text2,Text3,Name_PNG,Titel,AT):
    degrees_freedom = len(Dat1[0]) - Nump
    
    VDu = odr.Model(R_RTD)
    popt, cov = curve_fit(R_RTD_O, Dat1[0],Dat1[1])
    
    mydata = odr.RealData(Dat1[0],Dat1[1],sx = Dat1[2], sy = Dat1[3])
    myodr = odr.ODR(mydata,VDu,beta0 = popt)
    myoutput = myodr.run()
    myoutput.pprint()
    Chi2 = myoutput.sum_square/degrees_freedom
    print("Chi2/ndf: ",Chi2)
    
    #Funk.Nice_Plot(Dat1,Dat2,myoutput,Nump,Text1,Text2,Text3,Name_PNG,Titel,AT) #arguments: 'Dat1', 'Dat2', 'myoutput', 'Pol', 'Text1', 'Text2', 'Text3', and 'Name', Achsen
    return myoutput
    
def Nice_Plot_lin(Dat1,par,Pol,Text1,Text2,Name):
    x_line = np.linspace(np.min(Dat1[0])*1.1 , np.max(Dat1[0])*1.1 , 200)

    if Pol ==1:
        p1 = Pol1(par,Dat1[0])
        p3 = Pol1(par,x_line)    
        p7 = Pol1(par,Dat1[0]+Dat1[2])
        p8 = Pol1(par,Dat1[0]-Dat1[2])

    

    fig, axarray = plt.subplots(2, 1, figsize=(7,5), sharex=True, gridspec_kw={'height_ratios': [5, 2]}, dpi = 150)
    # Grafische Darstellung der Rohdaten
    axarray[0].errorbar(Dat1[0], Dat1[1], xerr=Dat1[2], yerr=Dat1[3], fmt='.', marker='o', label = Text1)
    
    
    axarray[0].set_xlabel('$\Delta v$ / 1/GHz')
    axarray[0].set_ylabel('$\Delta T$ / °C')
    #axarray[0].set_ylabel('Strain DMS $10^{-6}$')
    axarray[0].plot(x_line, p3, '--', color='red', label = Text2)
    axarray[0].grid()
    
   # sigR2  = ( np.abs(Dat1[1] - p7 ) + np.abs( Dat1[1] - p8 ) )/2
    sigR2  = np.abs( np.abs(p7 ) - np.abs( p8 ) )/2
    sigmaRes = np.sqrt(Dat1[3]**2 + sigR2**2)
    #axarray[0].plot(Y_A_V, Pol4(myoutput.beta,Y_A_V), color='green')
    #axarray[0].plot( x_line,  y_line4, '--', color='blue', label = "Pol. 4 Fit")
    axarray[0].legend()
    # Zunächst plotten wir eine gestrichelte Nulllinie, dann den eigentlichen Residuenplot:
    axarray[1].axhline(y=0., color='black', linestyle='--')
    axarray[1].errorbar(Dat1[0], Dat1[1] - p1, yerr=sigmaRes, fmt='.', marker='o')
    
    axarray[1].set_xlabel('$\Delta v$ / 1/GHz')
    axarray[1].set_ylabel('$(\Delta T-Fit(\Delta v)$ / °C')
    
  #  axarray[1].set_ylabel('$(Strain-Fit(\Delta v)$ $10^{-6}$')
    
    # Wir sorgen dafür, dass die y-Achse beim Residuenplot symmetrisch um die Nulllinie ist:
    ymax = max([abs(x) for x in axarray[1].get_ylim()])
    axarray[1].set_ylim(-ymax, ymax)
    axarray[1].grid()
    
    plt.tight_layout()
    fig.subplots_adjust(hspace=0.0)
    #plt.tight_layout() 
    plt.savefig(Name + ".png")
    plt.show()
    
    
    
    
    
class StrainGauge: # Strain-Gauge Parametrisation

    Tref = 20.0  # reference temperature 20 oC
    Kfit =  2.0  # K for compensation-curve
    alFit = [-3.9811e-2+ 2* 9.2683e-4 + 3* -2.0261e-6 + 4* 1.7127e-9]
    
    def __init__(self,PAR=[ 350.0, 2.21, -285e-6, -1.7e-2, 10.8e-6, -9.56, 0.89, -1.82e-2, 4.35e-5 ]): # default parameters for HBM 1-LC11-6/350 (steel-compensated)
        self.R0 = PAR[0] # nominal resistance in Ohm
        self.K  = PAR[1] # K-factor for strain:           eps  = (R-R0)/R0 / K  
        self.KT = PAR[2] # temperature-coefficient of K:  K(T) = K + KT*(T-T0)  with  T in oC and T0 = 20 oC
        self.TR = PAR[3] # transverse sensitivity
        self.al = PAR[4] # temperature compensated strain
        self.A  = PAR[5] # deviation from compensation: eps_s(T) = 1e-6 * (A + B*T + C*T^2 + D*T^3) with T in oC
        self.B  = PAR[6]
        self.C  = PAR[7]
        self.D  = PAR[8]

    def eps_s(self,T): return 1e-6*(self.A + self.B*T + self.C*T*T + self.D*T*T*T)
    
    def eps_gauge(T): 
        alFit = [-3.9811e-2+ 2* 9.2683e-4 + 3* -2.0261e-6 + 4* 1.7127e-9]
        Tref = 20 + 273.15
        return alFit[0]*(T - Tref) + 1/2 * alFit[1]*(T**2 - Tref**2) + 1/3*alFit[2]*(T**3 - Tref**3) + 1/4*alFit[3]*(T**4 - Tref**4)
        
    def strain(self,R,T):
        K = self.K*(1+self.KT*(T-StrainGauge.Tref))  # temperature corrected K               with Tref = 20.0 oC
        eps   = (R-self.R0)/self.R0 / K              # measured strain
        eps_s = self.eps_s(T) * StrainGauge.Kfit/K   # compensation                          with Kfit = 2.0
        eps_f = eps - eps_s  
       # print(T)                        # eps_f = eps_substrate - eps_gauge     with eps_gauge = al*(T-Tref)
        alFit = [-3.9811e-0 , 2* 9.2683e-2 , 3* -2.0261e-4 , 4* 1.7127e-7]
        Tref = 20 + 273.15
        return eps_f + self.al*(T-StrainGauge.Tref) #((alFit[0] - 4.58)*(T +273.15 - Tref) + 1/2 * alFit[1]*((T + 273.15)**2 - Tref**2) + 1/3*alFit[2]*((T + 273.15)**3 - Tref**3) + 1/4*alFit[3]*((T + 273.15)**4 - Tref**4))*10**-6 #self.eps_gauge( ) #self.al*(T-StrainGauge.Tref)  # return eps_substrate 
    
    
def Ohm_to_strain(Ohm,Ohmerr,T,Terr,W):
    KF = [2.07,2.31,2.37] # K_Factor
    TR = [0.002,-0.005,-0.016]
    Ep_Par = [[-9.56, 0.89, -1.82e-2, 4.35e-5 ] ,[-44.12,2.89,-3.87e-2,5.85e-5,-5.16e-9],[-44.12, 2.89, -3.87e-2, 5.85e-5 ]]
    Values = Ohm
    Valueserr = Ohmerr
    Strain = []
    Strainerr = []
    T_Strain = T
    T_Strainerr = Terr
    i = W
    DMS = StrainGauge()
    DMS.R0 = Ohm[0]#350 # Grundwiederstand, sonst nicht bekannt
    DMS.Tref =T[0]# 20 
    DMS.K = KF[i]
    DMS.TR = TR[i]
    DMS.A = Ep_Par[i][0]
    DMS.B = Ep_Par[i][1]
    DMS.C = Ep_Par[i][2]
    DMS.D = Ep_Par[i][3]
    
    SS = []
    SR = []
    for j in range(len(Values)):
        SS = DMS.strain(Values[j], T_Strain[j])
        S1 = DMS.strain(Values[j] + Valueserr[j] , T_Strain[j])
        S2 = DMS.strain(Values[j] - Valueserr[j], T_Strain[j])            
        SR =  np.abs(np.abs(S1) - np.abs(S2))/2  
        Strain.append(SS)
        Strainerr.append(SR)

    return np.array(Strain), np.array(Strainerr)
    #return np.array(Ohm), np.array(Ohmerr)
