
import csv
import matplotlib.pyplot as plt
import os.path
import numpy as np
from scipy.signal import savgol_filter
from datetime import datetime, timedelta


def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def Get_Data(Pathi):
    length = []
    strain = []
    Zeile_Start = 0
    Zeile_Starti = 23
    Path = Pathi + ".csv"
    if os.path.exists(Path) == False:
        Path = Pathi + ".txt"   
        Zeile_Starti = 15
        if os.path.exists(Path) == False:
            print("File does not exist or wrong directory!")
                
    with open(Path, "r") as csvdatei:
        csv_reader_object = csv.reader(csvdatei,delimiter=';')#, quoting=csv.QUOTE_NONNUMERIC)
            
        zeilennummer = 0
        for row in csv_reader_object:  
            
            if (zeilennummer == 0 and len(row) == 1):
                Zeile_Start = Zeile_Starti
               # print(Zeile_Start, len(row), row, isinstance(row[0], str))
                break
                
          #  if zeilennummer == 0:
             #   print(f'{", ".join(row)}')
             #   print(Zeile_Start, len(row), row, isinstance(row[0], str))
        
            zeilennummer += 1
            if zeilennummer >= Zeile_Start:
                length.append(float(row[0]))
                strain.append(float(row[1])) 
                
                
        if Zeile_Start != 0:                
            csv_reader_object = csv.reader(csvdatei,delimiter='\t')#, quoting=csv.QUOTE_NONNUMERIC)        
            for row in csv_reader_object:  
                                
              #  if zeilennummer == 0:
               #     print(f'{", ".join(row)}')
               #     print(Zeile_Start, len(row) , row, isinstance(row[0], str))
            
                zeilennummer += 1
                if zeilennummer >= Zeile_Start:
                    length.append(float(row[0]))
                    strain.append(float(row[1]))                      
   #     print(f'Number datapoints: {zeilennummer-1}')
    return length,strain

def Median_from_Position(X,Y,p1,p2):
    m = []
    mx = []
    for i in range(len(X)):
        if X[i]>= p1 and X[i]<= p2:
            m.append(Y[i])
            mx.append(X[i])
    #return [np.mean(mx), np.mean(m), np.std(m,ddof = 1)/np.sqrt(len(m))]
    return [np.mean(mx), np.median(m), np.std(m,ddof = 1)]

def Mean_from_Position(X,Y,p1,p2):
    m = []
    mx = []
    for i in range(len(X)):
        if X[i]>= p1 and X[i]<= p2:
            m.append(Y[i])
            mx.append(X[i])
    #return [np.mean(mx), np.mean(m), np.std(m,ddof = 1)/np.sqrt(len(m))]
    return [np.mean(mx), np.mean(m), np.std(m,ddof = 1)]

def Mean_from_Position_WEM(X,Y,p1,p2):
    m = []
    mx = []
    for i in range(len(X)):
        if X[i]>= p1 and X[i]<= p2:
            m.append(Y[i])
            mx.append(X[i])
    return [np.mean(mx), np.mean(m), np.std(m,ddof = 1)/np.sqrt(len(m))]
    #return [np.mean(mx), np.mean(m), np.std(m,ddof = 1)]
    
def Weighted_mean_from_Position(X,Y,Yerr,p1,p2):
    m = []
    mx = []
    merr = []
    for i in range(len(X)):
        if X[i]>= p1 and X[i]<= p2:
            m.append(Y[i])
            mx.append(X[i])
            merr.append(Yerr[i])
    #return [np.mean(mx), np.mean(m), np.std(m,ddof = 1)/np.sqrt(len(m))]
    return [np.mean(mx), np.average(m,weights = 1/np.array(merr)**2), np.sqrt(np.std(m,ddof = 1)**2 + np.mean(merr )**2)  ]


def Get_Data_Position(X,Y,p1,p2):
    m = []
    mx = []
    for i in range(len(X)):
        if X[i]>= p1 and X[i]<= p2:
            m.append(Y[i])
            mx.append(X[i])
    return [mx, m]



def Get_Data_Temp(Pathi):
    Date = []
    Time = []
    T1 = []
    T2 = []
    T3 = []
    T4 = []
    D = []
    Zeile_Start = 2
    Zeile_Starti = 1
    Path = Pathi + ".csv"
    if os.path.exists(Path) == False:
        Path = Pathi + ".txt"   
        Zeile_Starti = 1
        if os.path.exists(Path) == False:
            print("File does not exist or wrong directory!")
                
    with open(Path, "r") as csvdatei:
        csv_reader_object = csv.reader(csvdatei,delimiter=';')#, quoting=csv.QUOTE_NONNUMERIC)
            
        zeilennummer = 0
        for row in csv_reader_object:  
            
            if (zeilennummer == 0 and len(row) == 1):
                Zeile_Start = Zeile_Starti
               # print(Zeile_Start, len(row), row, isinstance(row[0], str))
                break
                
         #   if zeilennummer == 0:
         #       Date.append(float(row[0])) 
        
            zeilennummer += 1
            if zeilennummer >= Zeile_Start:
                Time.append(float(row[0]))
                T1.append(float(row[1])+273.15) 
                T2.append(float(row[2])+273.15) 
                T3.append(float(row[3])+273.15) 
                T4.append(float(row[4])+273.15) 
                D.append(float(row[5])) 
      
        if Zeile_Start != 0:                
            csv_reader_object = csv.reader(csvdatei,delimiter=';')#, quoting=csv.QUOTE_NONNUMERIC)        
            for row in csv_reader_object:  
                                
                if zeilennummer == 0:
                     Date.append(float(row[1])) 
            
                zeilennummer += 1
                if zeilennummer >= Zeile_Start:
                    Time.append(float(row[1]))
                    T1.append(float(row[2])+273.15) 
                    T2.append(float(row[3])+273.15) 
                    T3.append(float(row[4])+273.15) 
                    T4.append(float(row[5])+273.15) 
                    D.append(float(row[6])) 

    return [Time,T1,T2,T3,T4,D]

def Get_Data_Temp_v2(Pathi):
    Date = []
    Time = []
    T1 = []
    T2 = []
    D1 = []
    D2 = []
    Zeile_Start = 2
    Zeile_Starti = 1
    Path = Pathi + ".csv"
    if os.path.exists(Path) == False:
        Path = Pathi + ".txt"   
        Zeile_Starti = 1
        if os.path.exists(Path) == False:
            print("File does not exist or wrong directory!")
                
    with open(Path, "r") as csvdatei:
        csv_reader_object = csv.reader(csvdatei,delimiter=';')#, quoting=csv.QUOTE_NONNUMERIC)
            
        zeilennummer = 0
        for row in csv_reader_object:  
            
            if (zeilennummer == 0 and len(row) == 1):
                Zeile_Start = Zeile_Starti
               # print(Zeile_Start, len(row), row, isinstance(row[0], str))
                break
            zeilennummer += 1
            if zeilennummer >= Zeile_Start:
                Time.append(float(row[0]))
                T1.append(float(row[1])) 
                D1.append(float(row[2]))
                T2.append(float(row[3])) 
                D2.append(float(row[4]))

                 
      
        if Zeile_Start != 0:                
            csv_reader_object = csv.reader(csvdatei,delimiter=';')#, quoting=csv.QUOTE_NONNUMERIC)        
            for row in csv_reader_object:  
                                
                if zeilennummer == 0:
                     Date.append(float(row[1])) 
            
                zeilennummer += 1
                if zeilennummer >= Zeile_Start:
                    Time.append(float(row[1]))
                    T1.append(float(row[1])) 
                    D1.append(float(row[2]))
                    T2.append(float(row[3])) 
                    D2.append(float(row[4]))
    
    return [np.array(Time),np.array(T1),np.array(D1),np.array(T2),np.array(D2)]


def Get_Data_Temp_Mean(PathT,t1,t2):
    T1 = []
    T2 = []
    T3 = []
    T4 = []
    DMS = []
    TList = []     
    TLists = []    
    T1s = t1[0]*60*60*24 + t1[1]*60*60 + t1[2]*60 + t1[3]
    T2s = t2[0]*60*60*24 + t2[1]*60*60 + t2[2]*60 + t2[3]
  #  print(i,"and",T/60/60, T)
    DT = Get_Data_Temp(PathT)
 #   print(DT[1])
    
 #   for i in range(len(DT)-1):
 #       DA = Mean_from_Position(DT[0],DT[1],T1s,T2s)
 #       TList.append(DA[1])
 #       TLists.append(DA[0])
  #  print(DT[0][0], DT[1][0]-273.15, DT[2][0]-273.15, DT[3][0]-273.15)   
    for j in range(len(DT[0])-1):
        if DT[0][j] >=T1s and DT[0][j] <=T2s:
            T1.append(DT[1][j])
            T2.append(DT[2][j])
            T3.append(DT[3][j])
            T4.append(DT[4][j])
            DMS.append(DT[5][j])
            
    TList = [np.mean(T1)-273.15, np.mean(T2)-273.15, np.mean(T3)-273.15, np.mean(T4)-273.15, np.mean(DMS)]
    TLists = [np.std(T1,ddof = 1),np.std(T2,ddof = 1),np.std(T3,ddof = 1),np.std(T4,ddof = 1),np.std(DMS,ddof = 1)]

                
    return [TList,TLists]

def P_Mean_Fit(X,Y,Xi):
    m = (Y[1]-Y[0])/(X[1]-X[0])
    b = Y[0] - m*X[0]
    return m*Xi + b

def P_Mean_sy( x ,xl, yerr):
    s1, s2 = yerr
    x1, x2 = xl
    dx = x2 - x1
    out  = ( s1 * ( 1 + x1 / dx ) )**2
    out += ( s2 * x1 / dx  )**2
    out += ( s1**2 + s2**2 ) / dx**2 * x**2
    out -= 2 * ( s2**2 * x1 / dx**2 + s1**2 * (  1 + x1 / dx ) / dx ) * x
    return np.sqrt( out )


def Noise_Cut(Datax,Datay,Dm,Dstd):
    Dm = Dm
    N = 3
    Ds = np.abs(Dstd)
    Thr = Ds*N
    if Dstd > 100:
        Thr = 500
    if Dstd < 10:
        Thr = 10
   # print(Thr, Ds)
    Dr = np.array([])
    Dr2 = np.array([])
    for i in range(len(Datay)):
        if Dm >0:
            if Datay[i] > Dm - Thr and  Datay[i] < Dm + Thr :
               # print(Dm - N*Ds ,Dm + N*Ds,Datay[i])
                Dr = np.append(Dr,Datay[i])
                Dr2 = np.append(Dr2,Datax[i])
        if Dm<0:
            if Datay[i] < Dm + Thr and Datay[i] > Dm - Thr :         
                #print(Dm - N*Ds ,Dm + N*Ds,Datay[i])
                Dr = np.append(Dr,Datay[i])
                Dr2 = np.append(Dr2,Datax[i])
    return [Dr2,Dr]

def IQR_Filter(Data):
    data = np.sort(Data[1])
    # First quartile (Q1)
    Q1 = np.median(data[:int(len(data)/2)])
      
    # Third quartile (Q3)
    Q3 = np.median(data[int(len(data)/2):])
      
    # Interquartile range (IQR)
    IQR = Q3 - Q1
    N = 1.5
    rex = []
    rey = []
    for i in range(len(Data[0])):
        if Data[1][i] > Q1 - N*IQR and Data[1][i] < Q3 + N*IQR:
            rex.append(Data[0][i])
            rey.append(Data[1][i])
            
        if Data[1][i] <= Q1 - (N+1)*IQR or Data[1][i] >= Q3 + (N+1)*IQR:
            rex.append(Data[0][i])
            rey.append(np.median(Data[1]))
            
    return [np.array(rex),np.array(rey)]

def Stufen_Filter(Data,Thri):
    NStufe = 100
    Thr = Thri#*10**-6
    y = Data[1]
    yerr = Data[2]
    
    for i in range(len(Data[0]) - 1 - NStufe):
        Ni = 50
        if Ni>i:
            Ni = 0
        ym = np.median(y[i-Ni:i]) #np.mean( IQR_Filter( [np.zeros(len(y[Ni:i])) , y[Ni:i]] )[1] )
        EndeStufe = NStufe
        if np.abs(ym-y[i+1]) > Thr:
            
           # print(np.abs(y[i]-y[i+1]))
            
            for j in range(1,NStufe):

                if np.abs(ym - y[i+j]) < Thr:
                    EndeStufe = j - 1
                    #print(EndeStufe)
                    break
                    
            if EndeStufe > 1:
                         
                KorreM = ym - np.mean( y[i+1:i +EndeStufe])
              #  KorreM_std = np.sqrt(np.std( y[i+1:i +EndeStufe],ddof = 1)**2 )#+ np.std(y[i-Ni:i],ddof = 1) )
                
                for j in range(1,EndeStufe):
                    y[i + j] += KorreM
                    yerr[i + j] += np.abs(KorreM)
                    
                    
           # i += EndeStufe -2
   
                
    return [Data[0],y,yerr]
 
    
def Filter(Data,IQR,SMO,Stu):
    
    
    length, strain = Data[0], Data[1]
    strain_std = np.zeros(len(Data[0]))
    if len(Data) >2:
        strain_std = Data[2]
    
    
    if Stu !=0:
        
        FilterData = Stufen_Filter([length,strain,strain_std],Stu)
        length,strain = FilterData[0],FilterData[1]      
        if len(Data) >2:
            strain_std = FilterData[2] 
        
    if SMO != 0:
        strain = smooth(strain,SMO)
        
    if IQR == True:
        if np.abs(np.mean(strain)) > 200:
            FilterData = IQR_Filter([length,strain])
            length,strain = FilterData[0],FilterData[1]
    if len(Data) >2:   
        return [length, strain, strain_std]
    else:
        return [length, strain]
                       
    
    
def Get_Data_WarmUp(Path,PathT,XP,NS,NE,Day,Th,Tm,Ts):
    Strain = []
    Strain_std = []
    T1 = []
    T2 = []
    T3 = []
    T4 = []
    DMS = []
    DMS_std = []
    for i in range(NS,NE):
        if i<10:
            D = Get_Data(Path+"/" + "LN2_WarmUp_0" + "0" + str(i))
        if i>=10 and i<100:
            D = Get_Data(Path+"/" + "LN2_WarmUp_0" + str(i))
        if i>=100:
            D = Get_Data(Path+"/" + "LN2_WarmUp_" + str(i))
            
        Dat = Mean_from_Position(D[0],D[1],XP[0],XP[1]) 
        Ti = Noise_Cut(D[0],D[1],Dat[1],Dat[2])
        
        Dat = Mean_from_Position(Ti[0],Ti[1],XP[0],XP[1])  
        Ti = Noise_Cut(D[0],D[1],Dat[1],Dat[2])
        
        Dat = Mean_from_Position(Ti[0],Ti[1],XP[0],XP[1])
        
    #    print(Path+"/" + "LN2_WarmUp_0" + str(i))
        
       # print(i, "   " , Dat[2])
        T = i*(60.435643 + 0.605) + Th*60*60 + Tm*60 + Ts + Day*60*60*24
       # print(T)
       # print(i,"and",T/60/60/24, "and" , i*(60.435643 + 0.605) , Dat[1])
        DT = Get_Data_Temp(PathT)
        #print(DT[0][100])
        if Dat[2]< 40000: #400:  # 190 für dv
            for j in range(len(DT[0])-1):
                if DT[0][j] <T and DT[0][j+1] >=T:
                    T1.append(P_Mean_Fit([DT[0][j],DT[0][j+1]],[DT[1][j],DT[1][j+1]],T) )
                    T2.append(P_Mean_Fit([DT[0][j],DT[0][j+1]],[DT[2][j],DT[2][j+1]],T) )
                    T3.append(P_Mean_Fit([DT[0][j],DT[0][j+1]],[DT[3][j],DT[3][j+1]],T) )
                    T4.append(P_Mean_Fit([DT[0][j],DT[0][j+1]],[DT[4][j],DT[4][j+1]],T) )
                    DMS.append(P_Mean_Fit([DT[0][j],DT[0][j+1]],[DT[5][j],DT[5][j+1]],T) )
                    #DMS_std.append( np.abs(DT[5][j] - DT[5][j+1])/2 ) 
                    DMS_std.append( P_Mean_sy(T,[DT[0][j],DT[0][j+1]], [1,1])  )
                    Strain.append(Dat[1])
                    
                    Strain_std.append(Dat[2])
              
                    
                    
                    break
      #  print(T1)
    return [T1,T2,T3,T4,DMS,Strain,Strain_std, DMS_std]
   
def readinDate(Pathi):
    for i in range(2):
        Path = Pathi + ".csv"
        if os.path.exists(Path) == False:
            Path = Pathi + ".txt"
            if os.path.exists(Path) == False:
                Pathi = Pathi + "_Upper"
                if i>0:
                    print("File does not exist or wrong directory!")
        
    Text = ""
    TT = ""
    with open(Path, "r") as csvdatei:
        csv_reader_object = csv.reader(csvdatei,delimiter=';')#, quoting=csv.QUOTE_NONNUMERIC)
        for row in csv_reader_object:  
            Text = row[0]
            break 
    
    if Text[13] == "/":
        TT = Text[17:21]+'-0' +Text[12:13]+ '-' + Text[14:16]   + Text[24:]
      #  print(TT)
    if Text[13] == "/" and Text[15] =="/":
        TT = Text[16:20]+'-0' +Text[12:13]+ '-0' + Text[14:15] + " "   + Text[24:]
       # print(TT)
        
    if Text[13] != "/":
        TT = Text[18:22]+'-' +Text[12:14]+ '-' + Text[15:17]   + Text[25:]
       # print(TT)
        
        
   # print(TT)
    return  TT
            
#print(readinDate("/home/kappe/Desktop/Masterarbeit/OFM/Data/2022_02_22_FOSCryo/Strain_online/Min40_cool_down_00d"))

def Get_Data_All(Pathi):
    Data = []
    Zeile_Start = 1
    Path = Pathi #+ ".csv"
    if os.path.exists(Path) == False:
        Path = Pathi #+ ".txt"   
        if os.path.exists(Path) == False:
            print("File does not exist or wrong directory!")
                
    with open(Path, "r") as csvdatei:
        csv_reader_object = csv.reader(csvdatei,delimiter=' ')#, quoting=csv.QUOTE_NONNUMERIC)
            
        zeilennummer = 0
        for row in csv_reader_object:  
            zeilennummer += 1
            if zeilennummer > Zeile_Start:
                DZ = []
                DZ.append(row[0] + " " +row[1])
                for i in range(2,len(row)):
                    if len(row[i])!=0:
                        DZ.append(row[i])
                        
                Data.append(DZ)
    return Data

def Time_corr(T):
    return T*0.99999629 + 0.00347605 
def Get_Data_WarmUp_v2(Path,PathT,PathOBR,DatName,DatName_List,XP,NS,NE):
    
    DT = Get_Data_All(PathT)
  #  print(DT[0])
    

    Test = [[] for i in range(len(DT[0]))]
    Testss = [[] for i in range(len(DT[0]))]
    D_Path = np.array([])
    for i in range(NS,NE):
        ND = ""
        if i<10:
            #ND =  DatName + "_0" + str(i)
            ND =  DatName + "_00" + str(i)
        if i>=10 and i<100:
        #ND =  DatName + "_" + str(i)
            ND =  DatName + "_0" + str(i)
        if i>=100:
            ND =  DatName + "_" + str(i)
            #ND =  DatName + "_" + str(i)
        D_Path = np.append(D_Path, ND)
        
    D_Path = np.append(np.array(DatName_List),D_Path)    
    
   # plt.figure(11)
   # plt.title(Path)
   # plt.xlim(XP[0]-0.01,XP[1]+0.01)
   # plt.ylim(-6000,-1000)
   # print(D_Path)
    for i in range(len(D_Path)):
        D = Get_Data(Path+ "/" + D_Path[i])
        
        Dat = Mean_from_Position(D[0],D[1],XP[0],XP[1]) 
        Ti = Noise_Cut(D[0],D[1],Dat[1],Dat[2])
        
        Dat = Mean_from_Position(Ti[0],Ti[1],XP[0],XP[1])  
        Ti = Noise_Cut(D[0],D[1],Dat[1],Dat[2])
        
        Dat = Mean_from_Position(Ti[0],Ti[1],XP[0],XP[1])

        
        date_time_str = readinDate(PathOBR + "/" + D_Path[i])
       # print(date_time_str)
        DateOb = datetime.strptime(date_time_str, '%Y-%m-%d %H:%M:%S')
        #print(DateOb,date_time_str)
        #print(Test[0])
        DateText = '%Y-%m-%d %H:%M:%S'
        if DatName[0:3] == "OBR":
           # DateText = '%Y-%m-%d %H:%M:%S.%f'
            print(D_Path[i])
        
        #print('\n')
        #print(DT[1])
        #print(DateOb)
        #print(datetime.strptime(DT[1][0], '%Y-%m-%d %H:%M:%S') ,datetime.strptime(DT[1+1][0], '%Y-%m-%d %H:%M:%S'))
        #print(datetime.strptime(DT[1000][0], '%Y-%m-%d %H:%M:%S') ,datetime.strptime(DT[1000+1][0], '%Y-%m-%d %H:%M:%S'))
        if Dat[2]< 40000000: #400:  # 190 für dv
            for j in range(len(DT)-1):   
                #print(datetime.strptime(DT[j][0], '%Y-%m-%d %H:%M:%S') ,DateOb)
                if datetime.strptime(DT[j][0], DateText) <DateOb and datetime.strptime(DT[j+1][0], DateText) >= DateOb:
                  #  print(DateOb)
               #     print(datetime.strptime(DT[j][0], '%Y-%m-%d %H:%M:%S') ,DateOb)
                    for k in range(1,len(DT[0])):
                        DD1 = datetime.strptime(DT[j][0], DateText)
                        DD2 = datetime.strptime(DT[j+1][0], DateText)
                        #D1 = DD1.hour*60*60 + DD1.minute*60 + DD1.second + DD1.day*60*60*24 
                        #D2 = DD2.hour*60*60 + DD2.minute*60 + DD2.second + DD2.day*60*60*24 
                        #T = DateOb.hour*60*60 + DateOb.minute*60 + DateOb.second + DD1.day*60*60*24 
                        
                        D1 = DD1.hour*60*60 + DD1.day*60*60*24 + DD1.minute*60 + DD1.second  
                        D2 = DD2.hour*60*60 + DD2.day*60*60*24  + DD2.minute*60 + DD2.second 
                        T = DateOb.hour*60*60 + DD1.day*60*60*24 + DateOb.minute*60 + DateOb.second                        
                        Test[k-1].append(P_Mean_Fit([D1,D2],[float(DT[j][k]),float(DT[j+1][k])],T)) 
                        if float(DT[j][k]) > 200:
                            Testss[k-1].append(P_Mean_sy(T,[D1,D2], [0.01,0.01]))
                        if float(DT[j][k]) < 200:
                            Testss[k-1].append(P_Mean_sy(T,[D1,D2], [0.0005,0.005]))
                    Test[-1].append(Dat[1])
                    Testss[-1].append(Dat[2])
                    if Dat[2] == 0.0:
                        Testss[-1] == 10
                        print("NEEEEEEEEEEEEEEEEEEEEEEEEEEEE")
                
                    break
                
        
    plt.xlim(XP[0]-0.01,XP[1]+0.01)
    plt.ylim(-6000,-1000)
 #   plt.show()
    return Test, Testss
        
def Get_Data_WarmUp_v3(Path,PathT,PathOBR,DatName,DatName_List,XP,NS,NE):
    
    DT = Get_Data_All(PathT)
  #  print(DT[0])
    

    Test = [[] for i in range(len(DT[0]))]
    Testss = [[] for i in range(len(DT[0]))]
    D_Path = np.array([])
    for i in range(NS,NE):
        ND = ""
        if i<10:
            #ND =  DatName + "_0" + str(i)
            ND =  DatName + "_0" + str(i)
        if i>=10 and i<100:
        #ND =  DatName + "_" + str(i)
            ND =  DatName + "_" + str(i)
        if i>=100:
            ND =  DatName + "_" + str(i)
            #ND =  DatName + "_" + str(i)
        D_Path = np.append(D_Path, ND)
        
    D_Path = np.append(np.array(DatName_List),D_Path)    
    
   # plt.figure(11)
   # plt.title(Path)
   # plt.xlim(XP[0]-0.01,XP[1]+0.01)
   # plt.ylim(-6000,-1000)
   # print(D_Path)
   
    for i in range(len(D_Path)):
        D = Get_Data(Path+ "/" + D_Path[i])               
        Da = Get_Data_Position(D[0],D[1],XP[0],XP[1])
        Ti = IQR_Filter(Da)
      #  Ti = [Da[0],smooth(Da[1],4)]
        
       # Ti = Filter(Data,True,2,20)
        
       # Dat = Mean_from_Position(D[0],D[1],XP[0],XP[1]) 
        #Ti = Noise_Cut(D[0],D[1],Dat[1],Dat[2])
        
       # Dat = Mean_from_Position(Ti[0],Ti[1],XP[0],XP[1])  
      #  Ti = Noise_Cut(D[0],D[1],Dat[1],Dat[2])
        
        Dat = Mean_from_Position(Ti[0],Ti[1],XP[0],XP[1])

        
        date_time_str = readinDate(PathOBR + "/" + D_Path[i])
       # print(date_time_str)
        DateOb = datetime.strptime(date_time_str, '%Y-%m-%d %H:%M:%S')
        #print(DateOb,date_time_str)
        #print(Test[0])
        DateText = '%Y-%m-%d %H:%M:%S'
        if DatName[0:3] == "OBR":
           # DateText = '%Y-%m-%d %H:%M:%S.%f'
            print(D_Path[i])
        
        #print('\n')
        #print(DT[1])
        #print(DateOb)
        #print(datetime.strptime(DT[1][0], '%Y-%m-%d %H:%M:%S') ,datetime.strptime(DT[1+1][0], '%Y-%m-%d %H:%M:%S'))
        #print(datetime.strptime(DT[1000][0], '%Y-%m-%d %H:%M:%S') ,datetime.strptime(DT[1000+1][0], '%Y-%m-%d %H:%M:%S'))
        if Dat[2]< 40000000: #400:  # 190 für dv
            for j in range(len(DT)-1):   
                #print(datetime.strptime(DT[j][0], '%Y-%m-%d %H:%M:%S') ,DateOb)
                if datetime.strptime(DT[j][0], DateText) <DateOb and datetime.strptime(DT[j+1][0], DateText) >= DateOb:
                  #  print(DateOb)
               #     print(datetime.strptime(DT[j][0], '%Y-%m-%d %H:%M:%S') ,DateOb)
                    for k in range(1,len(DT[0])):
                        DD1 = datetime.strptime(DT[j][0], DateText)
                        DD2 = datetime.strptime(DT[j+1][0], DateText)
                        #D1 = DD1.hour*60*60 + DD1.minute*60 + DD1.second + DD1.day*60*60*24 
                        #D2 = DD2.hour*60*60 + DD2.minute*60 + DD2.second + DD2.day*60*60*24 
                        #T = DateOb.hour*60*60 + DateOb.minute*60 + DateOb.second + DD1.day*60*60*24 
                        
                        D1 = DD1.hour*60*60 + DD1.day*60*60*24 + DD1.minute*60 + DD1.second  
                        D2 = DD2.hour*60*60 + DD2.day*60*60*24  + DD2.minute*60 + DD2.second 
                        T = DateOb.hour*60*60 + DD1.day*60*60*24 + DateOb.minute*60 + DateOb.second                        
                        Test[k-1].append(P_Mean_Fit([D1,D2],[float(DT[j][k]),float(DT[j+1][k])],T)) 
                        if float(DT[j][k]) > 200:
                            Testss[k-1].append(P_Mean_sy(T,[D1,D2], [0.01,0.01]))
                        if float(DT[j][k]) < 200:
                            Testss[k-1].append(P_Mean_sy(T,[D1,D2], [0.0005,0.005]))
                    Test[-1].append(Dat[1])
                    Testss[-1].append(Dat[2])
                    if Dat[2] == 0.0:
                        Testss[-1] == 10
                        print("NEEEEEEEEEEEEEEEEEEEEEEEEEEEE " + D_Path[i] )
                
                    break
                
        
   # plt.xlim(XP[0]-0.01,XP[1]+0.01)
    #plt.ylim(-6000,-1000)
 #   plt.show()
    return Test, Testss

# IQR - filter 
def transform_txt_to_list(file_path: str) -> list:
    """
    Transform a text file into a list.

    Parameters:
    - file_path (str): The path to the text file.

    Returns:
    - list: A list containing the transformed data from the text file.

    Note:
    - Each line in the text file is expected to represent an item in the list.
    - If a line starts with a tab, it will be appended to the previous line (indicated by the indentation).
    - Only every other line is considered (ignoring empty lines).
    """
    # Read content of the text file
    with open(file_path, 'r') as file:
        lines = file.readlines()

    result_list = []
    current_line = ""

    for line in lines:
        # Remove leading and trailing whitespaces
        line = line.strip()
        
        if line.startswith('\t'):
            # Append the indented line to the current line
            current_line += line.strip()
        else:
            # If it's not indented, add the current line to the result list
            if current_line:
                result_list.append(current_line)
            # Start a new current line
            current_line = line

    # Add the last line if it's not added
    if current_line:
        result_list.append(current_line)

    filtered_list = result_list[::2]
    return filtered_list

import pandas as pd
from typing import List, Tuple

def upper_and_lower_limit(data: Tuple[List[float], List[float]]) -> Tuple[float, float]:
    """
    Calculate the upper and lower limits based on the interquartile range (IQR).
    
    Parameters:
    - data (Tuple[List[float], List[float]]): A tuple containing lists for x and y values.
    
    Returns:
    - Tuple[float, float]: The lower and upper limits calculated using the IQR method.
    """
    y = pd.Series(data[1])  # Convert y list to a pandas Series for easy quantile calculations
    Q1 = y.quantile(0.25)
    Q3 = y.quantile(0.75)
    IQR = Q3 - Q1
    lower_limit = Q1 - 1.5 * IQR
    upper_limit = Q3 + 1.5 * IQR

    return lower_limit, upper_limit

def apply_iqr_filter_with_interpolation(data: Tuple[List[float], List[float]]) -> Tuple[List[float], List[float]]:
    """
    Apply the interquartile range (IQR) filter to replace outliers with interpolated values.
    
    Parameters:
    - data (Tuple[List[float], List[float]]): A tuple containing lists for x and y values.
    
    Returns:
    - Tuple[List[float], List[float]]: x and y lists with outliers replaced by interpolated values.
    
    Note:
    - Outliers are identified using the upper and lower limits calculated
      based on the interquartile range (IQR) method.
    - Outliers are replaced by the average of their neighboring values (i-1 and i+1).
    """
    x, y = data
    lower_limit, upper_limit = upper_and_lower_limit(data)
    print(f"Lower limit: {lower_limit}, Upper limit: {upper_limit}")
    
    # Replace outliers with interpolated values
    y_filtered = y.copy()  # Create a copy to modify
    outlier_count = 0

    for i in range(1, len(y) - 1):  # Start from 1 and end at len(y) - 1 to avoid edge cases
        if y[i] < lower_limit or y[i] > upper_limit:
            # Interpolate by averaging the neighboring values
            y_filtered[i] = (y[i - 1] + y[i + 1]) / 2
            outlier_count += 1

    print(f"{outlier_count} outliers replaced with interpolated values")
    return x, y_filtered

from typing import List

def replace_with_interpolation(lst: List[float], i: int) -> List[float]:
    """
    Replace the element at index i in the list with the mean of its neighbors (i-1 and i+1).
    
    Parameters:
    - lst (List[float]): The input list.
    - i (int): The index of the element to replace.
    
    Returns:
    - List[float]: The modified list with the element at index i replaced by interpolation.
    
    Note:
    - If i is 0 (the first element), the element is replaced by the value at i+1.
    - If i is the last index, the element is replaced by the value at i-1.
    """
    # Copy the list to avoid modifying the original list
    lst_copy = lst.copy()
    
    # Handle edge cases
    if i == 0:
        # If i is the first element, replace it with the next element
        lst_copy[i] = lst[i + 1]
    elif i == len(lst) - 1:
        # If i is the last element, replace it with the previous element
        lst_copy[i] = lst[i - 1]
    else:
        # Replace element at i with the mean of the neighbors
        lst_copy[i] = (lst[i - 1] + lst[i + 1]) / 2

    return lst_copy


def Get_Data_WarmUp_v4(Path,PathT,timestamp_path,DatName,DatName_List,XP,NS,NE,IQR,Stu,SMO=True):
    #print(f"DEBUG 1")
    DT = Get_Data_All(PathT)
    print(f"DEBUG DT: {DT}")
    #np.savetxt("Test" + '.out', (DT)) 
    timestamp_list = transform_txt_to_list(timestamp_path)

    Test = [[] for i in range(len(DT[0]))]
    Testss = [[] for i in range(len(DT[0]))]
    TS = []
    D_Path = np.array([])
    for i in range(NS,NE):
        
        ND = DatName+ "_"+ str(i)
        D_Path = np.append(D_Path, ND)
        
    D_Path = np.append(np.array(DatName_List),D_Path)    
   
    for i in range(len(D_Path)):
        D2 = Get_Data(Path+ "/" + D_Path[i])
        #D2 = apply_iqr_filter_with_interpolation(D)
        Da = Get_Data_Position(D2[0],D2[1],XP[0],XP[1])
        Da.append(np.zeros(len(Da[0])))
        Ti = Filter(Da,IQR,SMO,Stu)
        
        Dat = Mean_from_Position(Ti[0],Ti[1],XP[0],XP[1])
        Dat[2] = np.sqrt(Dat[2]**2 + Mean_from_Position_WEM(Ti[0],Ti[2],XP[0],XP[1])[1]**2  ) # if Stufenfilter korrektet data for mean, korr. is added to error
        
        date_time_str = timestamp_list[i].replace('_', ' ')
        
        DateOb = datetime.strptime(date_time_str, '%Y-%m-%d %H:%M:%S')
        
        DateText = '%Y-%m-%d %H:%M:%S'
        
        if Dat[2]< 40000000: #400:  # 190 für dv
            for j in range(len(DT)-1):
                #IQR Filter of DV before mean
                if datetime.strptime(DT[j][0], DateText) <DateOb and datetime.strptime(DT[j+1][0], DateText) >= DateOb:
                    
                    for k in range(1,len(DT[0])):
                        
                        DD1 = datetime.strptime(DT[j][0], DateText)
                        DD2 = datetime.strptime(DT[j+1][0], DateText)

                        D1 = DD1.hour*60*60 + DD1.day*60*60*24 + DD1.minute*60 + DD1.second  
                        D2 = DD2.hour*60*60 + DD2.day*60*60*24  + DD2.minute*60 + DD2.second 
                        T = DateOb.hour*60*60 + DD1.day*60*60*24 + DateOb.minute*60 + DateOb.second                        
                        Test[k-1].append(P_Mean_Fit([D1,D2],[float(DT[j][k]),float(DT[j+1][k])],T)) 
                        if float(DT[j][k]) > 200:
                            Testss[k-1].append(P_Mean_sy(T,[D1,D2], [0.01,0.01]))
                        if float(DT[j][k]) < 200:
                            Testss[k-1].append(P_Mean_sy(T,[D1,D2], [0.0005,0.005]))

                    
                    Test[-1].append(Dat[1])
                    Testss[-1].append(Dat[2])
                    TS.append(DateOb)
                    if Dat[2] == 0.0:
                        Testss[-1] == 1 #0 10 if Strain
                        print("Warning: Var = 0, only one Datapoint for mean. Error was set to 10 mu Strain or 1 Shift for " + D_Path[i] )
                
                    break
                
                        
    return Test, Testss, TS

def Get_Data_Simulation(Pathi):
    length = []
    strain = []
    Zeile_Start = 2
    Zeile_Starti = 2
    Path = Pathi + ".csv"
    if os.path.exists(Path) == False:
        Path = Pathi + ".txt"   
        Zeile_Starti = 15
        if os.path.exists(Path) == False:
            print("File does not exist or wrong directory!")
                
    with open(Path, "r") as csvdatei:
        csv_reader_object = csv.reader(csvdatei,delimiter=';')#, quoting=csv.QUOTE_NONNUMERIC)
            
        zeilennummer = 0
        for row in csv_reader_object:  
            
            if (zeilennummer == 0 and len(row) == 1):
                Zeile_Start = Zeile_Starti
               # print(Zeile_Start, len(row), row, isinstance(row[0], str))
                break
                
            if zeilennummer == 0:
                print(f'{", ".join(row)}')
                print(Zeile_Start, len(row), row, isinstance(row[8], str))
        
            zeilennummer += 1
            if zeilennummer >= Zeile_Start:
                length.append(float(row[8]))
                strain.append(float(row[7])) 
                
                
        if Zeile_Start != 0:                
            csv_reader_object = csv.reader(csvdatei,delimiter='\t')#, quoting=csv.QUOTE_NONNUMERIC)        
            for row in csv_reader_object:  
                                
                if zeilennummer == 0:
                    print(f'{", ".join(row)}')
                    print(Zeile_Start, len(row) , row, isinstance(row[8], str))
            
                zeilennummer += 1
                if zeilennummer >= Zeile_Start:
                    length.append(float(row[8]))
                    strain.append(float(row[7]))                      
        print(f'Number datapoints: {zeilennummer-1}')
    return length,strain

def Plotter(Pathi,Datei,xl,yl,FIQR = None,YLabel = "Strain*E-6",Stdvv = False,Titel_Text=None):
    k_S_C = -6.668 #-Lambda/K_S_C/c *10**6  # -6.67#muStrain/GHz # used are -6.668
    c = 299792458.0 #m/s
    Lambda = 1310.0
    St_to_dvv = Lambda/c/k_S_C*10**6 # [10**9]  ## Umrechnung von Strain in norm. Shift
    
    Path = Pathi + "/"  + Datei 
      
    length, strain = Get_Data(Path)
    
    if FIQR != None:
       # Data = [length,strain]
       # Data = IQR_Filter(Data)
       # length = Data[0]
       # strain = Data[1]
        
        l1,s1 = Get_Data_Position(length,strain,xl[0],FIQR[0])
        l2,s2 = Get_Data_Position(length,strain,FIQR[1],xl[1])
        
        length1, strain1 = Get_Data_Position(length,strain,FIQR[0],FIQR[1])
        FilterData =  IQR_Filter([length1, strain1])
        length,strain = FilterData[0],FilterData[1]
        
        length = np.concatenate((l1,length,l2))
        strain = np.concatenate((s1,strain,s2))
        
    if Stdvv == True:
        strain = np.array(strain)*St_to_dvv
        
    
    N1 = 4.73 #4.75
    N2 = 4.9
    N3 = 4.8 #5.0
    N4 = 5.03
    Start = []
    End = []
    for i in range(len(length)):
        if length[i] >= N1 and len(Start) == 0:
            Start.append(i)
        if length[i] >= N2 and len(Start) == 1:
            Start.append(i)
        if length[i] >= N3 and len(End) == 0:
            End.append(i)
        if length[i] >= N4 and len(End) == 1:
            End.append(i)
    
  #  print("Faser strain:",np.mean(strain[Start[0]:End[0]]),"Alu Strain: ",np.mean(strain[Start[1]:End[1]]))
    plt.figure(0, dpi = 180)
    plt.plot(length,strain)
    plt.xlabel("Length [m]")
    plt.ylabel(YLabel)
    TText = "Measurment: " + Datei
    plt.title(TText)
    if Titel_Text != None:
        plt.title(Titel_Text)
    #plt.legend()
  #  plt.xlim(4.5,5.2)
    plt.xlim(xl[0],xl[1])
    plt.ylim(yl[0],yl[1])
   # plt.xlim(4.0,5.2)
    plt.grid()
    SText = "Pictures_Raw_Data/" + Datei + ".png"
    plt.tight_layout()
    plt.savefig(SText)
    plt.show()
    
    
def Plotter_SGF(Pathi,Datei,Ref,a,b):
    Path1 = Pathi + "/"  + Datei 
    Path2 = Pathi + "/"  + Ref
    length,strain = Differece_calc(Path1,Path2)
    
   # Path = Pathi + "/"  + Datei     
    #length, strain = Get_Data(Path)    

    #strain_f = savgol_filter(strain, a, b)
    strain_f = smooth(strain,a)
    
    plt.figure(0)
    plt.plot(length,strain_f, label = "SGF a: " + str(a) + " b: " + str(b) )
    plt.plot(length,strain, label = "No Filter")
    plt.xlabel("length [m]")
    plt.ylabel("Strain*E-6")
    TText = "Measurment: " + Datei
    plt.title(TText)
    plt.legend()
    #plt.xlim(4.0,6.0)
    plt.xlim(4.6,5.6)
    plt.ylim(-50,300)
    plt.grid()
    SText = "Pictures_Raw_Data/" + Datei + "_Filter_" + "a_"+ str( a) + "_b_" + str(b) + ".pdf"
    plt.savefig(SText)
    plt.show()

def Differece_calc(Path1,Path2):
    length1, strain1 = Get_Data(Path1)
    length2, strain2 = Get_Data(Path2)
    length = length1
    strain = []
    NPoint = len(length1)
    if( len(length1) >= len(length2) ):
        NPoint = len(length2)
        length = length2
    Ende = length1[-1]
    if(length1[-1] >= length2[-1]):
        Ende = length2[-1]
    
    Start = length2[0]
    ST = False
    if(length1[0] >= length2[0]):
        Start = length1[0]
        ST = True
    St = [0,0] 
    for i in range(NPoint):
        if(ST == False):
            if(Start == length1[i]):
                St[0] = i
                break
        if(ST == True):
            if(Start == length2[i]):
                St[1] = i
                break        
 #   print(length1[0],length2[0],length1[-1],length2[-1])   
    for i in range(NPoint):
        if length1[i] == length2[i]:
            strain.append(strain1[i] - strain2[i])   
            
            
        if length1[0] != length2[0]:      
            if (length[i] == Start ):                
                for j in range(i,NPoint):
                    if (ST == True):
                        strain.append(strain1[j] - strain2[j + St[1]]) 
                        if i == 100:
                            print("Data",length1[j] ,length2[j + St[1]])
                    if (ST == False):
                        strain.append(strain1[j + St[0]] - strain2[j])  
                        if i == 100:
                            print("Data",length1[j + St[0]] ,length2[j])          
            if (length[i] == Ende):    
                break
    if (ST == True):
        length = length[St[1]:]   
    if (ST == False):
        length = length[St[0]:] 
    return length,strain


def Differece_calc_mean(Path1,Path2,Datei,Ref):
    length1, strain1 = Data_mean(Path1,Datei,)
    length2, strain2 = Data_mean(Path2,Ref)
    length = length1
    strain = []
    NPoint = len(length1)
    if( len(length1) >= len(length2) ):
        NPoint = len(length2)
        length = length2
    Ende = length1[-1]
    if(length1[-1] >= length2[-1]):
        Ende = length2[-1]
    
    Start = length2[0]
    ST = False
    if(length1[0] >= length2[0]):
        Start = length1[0]
        ST = True
    St = [0,0] 
    for i in range(NPoint):
        if(ST == False):
            if(Start == length1[i]):
                St[0] = i
                break
        if(ST == True):
            if(Start == length2[i]):
                St[1] = i
                break        
 #   print(length1[0],length2[0],length1[-1],length2[-1])   
    for i in range(NPoint):
        if length1[i] == length2[i]:
            strain.append(strain1[i] - strain2[i])   
            
            
        if length1[0] != length2[0]:      
            if (length[i] == Start ):                
                for j in range(i,NPoint):
                    if (ST == True):
                        strain.append(strain1[j] - strain2[j + St[1]]) 
                        if i == 100:
                            print("Data",length1[j] ,length2[j + St[1]])
                    if (ST == False):
                        strain.append(strain1[j + St[0]] - strain2[j])  
                        if i == 100:
                            print("Data",length1[j + St[0]] ,length2[j])          
            if (length[i] == Ende):    
                break
    if (ST == True):
        length = length[St[1]:]   
    if (ST == False):
        length = length[St[0]:] 
    return length,strain
    
def Plotter_several(Path_G,Pathi,Titel,Datei, xl, yl,FIQR = None,YLabel = "Strain*E-6",Stdvv = None, Legend_Text = None,Stto = None,NamePNG= None,ShiftData = None, StufenFilter = None):
    k_S_C = -6.668 #-Lambda/K_S_C/c *10**6  # -6.67#muStrain/GHz # used are -6.668
    c = 299792458.0 #m/s
    Lambda = 1310.0
    St_to_dvv = Lambda/c/k_S_C *10**6# [10**9]  ## Umrechnung von Strain in norm. Shift

    plt.figure(0, dpi =180)
    plt.xlabel("Length [m]")
    #plt.ylabel("Temperatur [K]")  #
    plt.ylabel(YLabel)
    plt.title(Titel)
    #plt.xlim(4.8,5.2)
    plt.grid()
    
    for i in range(len(Datei)):
        Path = Path_G + Pathi[i] + "/"  + Datei[i] 
        if np.size(Pathi) == 1:
            Path = Path_G + Pathi + "/"  + Datei[i] 
        length, strain = Get_Data(Path)
        
        if ShiftData != None:
            length =  np.array(length) +ShiftData[i][0]
            strain = np.array(strain) +ShiftData[i][1]
        if FIQR != None:
      #      if FIQR[i] != None:
       #         Data = [length,strain]
        #        Data = IQR_Filter(Data)
         #       length = Data[0]
          #      strain = Data[1]
                
            if FIQR[i] != None:
                l1,s1 = Get_Data_Position(length,strain,xl[0],FIQR[i][0])
                l2,s2 = Get_Data_Position(length,strain,FIQR[i][1],xl[1])
                
                length1, strain1 = Get_Data_Position(length,strain,FIQR[i][0],FIQR[i][1])
                FilterData =  IQR_Filter([length1, strain1])
                length,strain = FilterData[0],FilterData[1]
                
                length = np.concatenate((l1,length,l2))
                strain = np.concatenate((s1,strain,s2))
                
        if StufenFilter !=None:
            if len(StufenFilter[i])>0:
                VL = StufenFilter[i]
                l1,s1 = Get_Data_Position(length,strain,xl[0],VL[0])
                l2,s2 = Get_Data_Position(length,strain,VL[1],xl[1])
                
                length, strain = Get_Data_Position(length,strain,VL[0],VL[1])
                FilterData = Filter([length,strain],False,0,1)
                length,strain = FilterData[0],FilterData[1]
                
                length = np.concatenate((l1,length,l2))
                strain = np.concatenate((s1,strain,s2))


        
        if Stdvv != None:
            strain = np.array(strain)*St_to_dvv
        if Stto != None:
            strain = np.array(strain)*Stto
            
                   
            
        TextL = Datei[i]
        if Legend_Text !=None:
            TextL = Legend_Text[i]
        plt.plot(length,strain,alpha = 0.6, label = TextL)

    SText = "Pictures_Raw_Data/" + Titel 
    plt.legend()
    if Legend_Text != None:
        if len(Legend_Text)>5:
            plt.legend(fontsize = 5)
    
   # plt.xlim(4.0,5.6)
    plt.xlim(xl[0],xl[1])
    plt.ylim(yl[0],yl[1])
    plt.tight_layout()
    if NamePNG == None:
        plt.savefig(SText + ".pdf")
        plt.savefig(SText + ".png")
    if NamePNG != None:
        plt.savefig(NamePNG + ".pdf")
        plt.savefig(NamePNG + ".png")
    plt.show()


def Plotter_several_with_Name(Path_G,Pathi,Titel,Datei, DateiText, xl, yl):
    plt.figure(0, dpi =180)
    plt.xlabel("Length [m]")
    #plt.ylabel("Temperatur [K]")  #
    plt.ylabel("Strain*E-6")
    plt.title(Titel)
    #plt.xlim(4.8,5.2)
    plt.grid()
    
    for i in range(len(Datei)):
        Path = Path_G + Pathi[i] + "/"  + Datei[i] 
        if np.size(Pathi) == 1:
            Path = Path_G + Pathi + "/"  + Datei[i] 
        if np.size(Pathi) == 0:
            Path = Path_G + "/"  + Datei[i] 
            
        length, strain = Get_Data(Path)
        
        plt.plot(length,strain,alpha = 0.6, label = DateiText[i])

    SText = "Pictures_Raw_Data/" + Titel 
    plt.legend(fontsize = 8)
   # plt.xlim(4.0,5.6)
    plt.xlim(xl[0],xl[1])
    plt.ylim(yl[0],yl[1])
    plt.tight_layout()
    plt.savefig(SText + ".pdf")
    plt.savefig(SText + ".png")
    plt.show()
    
    

    #DOUBT what are these parameters and what does this function do 
def Plotter_several_with_Name_vLine_Filter(Path_G,Pathi,Titel,Datei, DateiText, xl, yl,VL,IQR,SMO,Stu):
    
    plt.figure(0, dpi =180)
    plt.xlabel("Length [m]")
    #plt.ylabel("Temperatur [K]")  #
    plt.ylabel("Strain*E-6")
    plt.title(Titel)
    #plt.xlim(4.8,5.2)
    plt.grid()
    
    for i in range(len(Datei)):
        
        if np.size(Pathi) == 1:
            Path = Path_G + Pathi + "/"  + Datei[i] 
        if np.size(Pathi) == 0:
            Path = Path_G + "/"  + Datei[i] 
        else:
            Path = Path_G + Pathi[i] + "/"  + Datei[i] 
        length, strain = Get_Data(Path)
        if len(VL) != 0:
            l1,s1 = Get_Data_Position(length,strain,xl[0],VL[0])
            l2,s2 = Get_Data_Position(length,strain,VL[1],xl[1])
            
            length, strain = Get_Data_Position(length,strain,VL[0],VL[1])
            FilterData = Filter([length,strain],IQR,SMO,Stu)
            length,strain = FilterData[0],FilterData[1]
            
            length = np.concatenate((l1,length,l2))
            strain = np.concatenate((s1,strain,s2))
        else:
            length, strain = Get_Data_Position(length,strain,xl[0],xl[1])
            FilterData = Filter([length,strain],IQR,SMO,Stu)
            length,strain = FilterData[0],FilterData[1]
        #length, strain = Get_Data_Position(length,strain,xl[0],xl[1])           
            
        
        
        
        plt.plot(length,strain,alpha = 0.6, label = DateiText[i])

    SText = "Pictures_Raw_Data/" + Titel 
    if len(VL) != 0:
        for i in range(len(VL)):
            plt.vlines(VL[i],yl[0],yl[1],linestyles="--",color="r")
    
   # plt.legend()
   # plt.xlim(4.0,5.6)
    plt.xlim(xl[0],xl[1])
    plt.ylim(yl[0],yl[1])
    plt.tight_layout()
    #plt.savefig(SText + ".pdf")
    plt.savefig(SText + ".png")
    plt.show()
    
def Plotter_several_with_Name_and_TorS(Path_G,Pathi,Titel,Datei, DateiText, xl, yl,ToS):
    Faktor = 1
    TextY = "Strain*E-6"
    if ToS == 1:
        Faktor = 0.1016003693934653
        TextY = "T - T$_{Ref. \, =\, 30}$ [°C]"  
    if ToS == 2:
        Faktor = 0.026762
        TextY =  "T - T$_{Ref. \, = \, 30}$ [°C]"  
        
    plt.figure(0, dpi =100)
    plt.xlabel("Length [m]")
    #plt.ylabel("Temperatur [K]")  #
    plt.ylabel(TextY)
    plt.title(Titel)
    #plt.xlim(4.8,5.2)
    plt.grid()
    
    for i in range(len(Datei)):
        Path = Path_G + Pathi[i] + "/"  + Datei[i] 
        if np.size(Pathi) == 1:
            Path = Path_G + Pathi + "/"  + Datei[i] 
        length, strain = Get_Data(Path)
        
        plt.plot(length,np.array(strain)*Faktor,alpha = 0.6, label = DateiText[i])

    SText = "Pictures_Raw_Data/" + Titel 
    plt.legend(title = "File Number")
   # plt.xlim(4.0,5.6)
    plt.xlim(xl[0],xl[1])
    plt.ylim(yl[0],yl[1])
    plt.tight_layout()
    plt.savefig(SText + ".pdf")
    plt.savefig(SText + ".png")
    plt.show()
 
    
def Plotter_several_with_Name_and_TorS_and_diffXmax(Path_G,Pathi,Titel,Datei, DateiText, DataCut, xl, yl,ToS):
    Faktor = 1
    TextY = "Strain*E-6"
    if ToS == 1:
        Faktor = 0.1016003693934653
        TextY = "T - T$_{Ref. \, =\, 30}$ [°C]"  
    if ToS == 2:
        Faktor = 0.026762
        TextY =  "T - T$_{Ref. \, = \, 30}$ [°C]"  
        
    plt.figure(0, dpi =100)
    plt.xlabel("Length [m]")
    #plt.ylabel("Temperatur [K]")  #
    plt.ylabel(TextY)
    plt.title(Titel)
    #plt.xlim(4.8,5.2)
    plt.grid()
    
    for i in range(len(Datei)):
        Path = Path_G + Pathi[i] + "/"  + Datei[i] 
        if np.size(Pathi) == 1:
            Path = Path_G + Pathi + "/"  + Datei[i] 
        length, strain = Get_Data(Path)
        if len(DataCut[i]):
            D = Get_Data_Position(length,strain,DataCut[i][0],DataCut[i][1])
            length = D[0]
            strain = D[1]
        
        plt.plot(length,np.array(strain)*Faktor,alpha = 0.6, label = DateiText[i])

    SText = "Pictures_Raw_Data/" + Titel 
    plt.legend(title = "File Number")
   # plt.xlim(4.0,5.6)
    plt.xlim(xl[0],xl[1])
    plt.ylim(yl[0],yl[1])
    plt.tight_layout()
    plt.savefig(SText + ".pdf")
    plt.savefig(SText + ".png")
    plt.show()
    
def Data_diff(Path_1,Path_2,Datei,Ref):
    Path1 = Path_1 + "/"  + Datei 
    Path2 = Path_2 + "/"  + Ref
    length,strain = Differece_calc(Path1,Path2)
    return length,strain

def Data_mean(Pathi,Dat):
    l = []
    s = []
    for i in range(len(Dat)):
        
        Path = Pathi + "/" + Dat[i]
        length, strain = Get_Data(Path)
        l.append(length)
        s.append(strain)
        
    lm = np.array(l[0])
    sm = np.array(s[0])
    for i in range(len(l)):
        if i > 0:             
            lm = lm + np.array(l[i])
            sm = sm + np.array(s[i])
    return lm/len(Dat),sm/len(Dat)

def Plotter_diff_mean(Path1,Path2,Dat,Ref,Name,xl,yl):

    length,strain = Differece_calc_mean(Path1,Path2,Dat,Ref)
   
    plt.figure(0)
    plt.plot(length,strain)
    plt.xlabel("length [m]")
    plt.ylabel("Strain*E-6")
    TText = Name
    plt.title(TText)
    #plt.legend()
    #plt.xlim(4.9,5.2)
    #plt.ylim(-50,6000)
    plt.xlim(xl[0],xl[1])
    
    plt.ylim(yl[0],yl[1])
    plt.grid()
    SText = "Pictures_Raw_Data/" + Name  + ".pdf"
    plt.savefig(SText)
    plt.show()
    
    
    
    
def Plotter_to_Reference(Path_1,Path_2,Datei,Ref,xl,yl):

    Path1 = Path_1 + "/"  + Datei 
    Path2 = Path_2 + "/"  + Ref
    length,strain = Differece_calc(Path1,Path2)
   
    plt.figure(0)
    plt.plot(length,strain)
    plt.xlabel("length [m]")
    plt.ylabel("Strain*E-6")
    TText = "Measurment: " + Datei + " to " + Ref
    plt.title(TText)
    #plt.legend()
    #plt.xlim(4.9,5.2)
    #plt.ylim(-50,6000)
    plt.xlim(xl[0],xl[1])

    plt.ylim(yl[0],yl[1])
    plt.grid()
    SText = "Pictures_Raw_Data/" + Datei + "_to_" + Ref + ".pdf"
    plt.savefig(SText)
    plt.show()

def Plotter_to_Ref_several(Path_G,Path_1,Path_2,Titel,Datei,Ref,a,b, xl,yl):
    plt.figure(0,dpi = 180)
    plt.xlabel("length [m]")
    plt.ylabel("Strain*E-6")
    plt.title(Titel)
    #plt.xlim(4.8,5.2)
    plt.grid()
    
    for i in range(len(Datei)):
        
        Path1 =Path_G +  Path_1[i] + "/"  + Datei[i] 
        
        if Ref[i] != 0:
            Path2 = Path_G + Path_2[i] + "/"  + Ref[i]
            length,strain = Differece_calc(Path1,Path2)
            Text = " - " + Ref[i]
        if Ref[i] == 0:
            length,strain = Get_Data(Path1)
            Text = " "
        #strain_f = savgol_filter(strain, a, b)
        #strain_f = smooth(strain, a)
        
        #plt.scatter(length,strain,alpha = 0.2,label = Datei[i])
        plt.plot(length,strain,label = Datei[i] + Text)
       # plt.plot(length,strain_f,label = Datei[i] + "SGF")
       # plt.xlim(4.9,5.55)
        plt.xlim(xl[0],xl[1])
        #plt.ylim(-20,200)
        plt.ylim(yl[0],yl[1])

    SText = "Pictures_Raw_Data/" + Titel 
    plt.legend()
    plt.savefig(SText + ".pdf")
    plt.savefig(SText + ".png")
    plt.show()

def Plotter_Roh(Pathi,Datei):
    Path = "/home/kappe/Desktop/Masterarbeit/OFM/First Measurment/" + Pathi + "/"  + Datei + ".csv"
    length = []
    strain = []
    with open(Path, "r") as csvdatei:
        csv_reader_object = csv.reader(csvdatei,delimiter='\t')#, quoting=csv.QUOTE_NONNUMERIC)
    
        zeilennummer = 0

        for row in csv_reader_object:
    
            if zeilennummer == 0:
                print(f'Spaltennamen sind: {", ".join(row)}')
            #else:
                #print(f'- Nachname: {row[0]} \t| Vorname: {row[1]} \t| Geburtstag: {row[2]}.')
            #    print(f'Spaltennamen sind: {", ".join(row)}')
            zeilennummer += 1
            if zeilennummer >= 15:
               # print(f'Spaltennamen sind: {", ".join(row)}')
                length.append(float(row[0]))
                strain.append(float(row[1]))
                
            
        print(f'Anzahl Datensätze: {zeilennummer-1}')
       # print(length[0],strain[0])
        
    plt.figure(0)
    plt.plot(length,strain)
    plt.xlabel("length [m]")
    plt.ylabel("Amplitude (dB/mm)")
    TText = "Measurment: " + Datei
    plt.title(TText)
    #plt.legend()
    plt.xlim(1.9,6)
    plt.grid()
    SText = "Pictures_Raw_Data/" + Datei + ".pdf"
    plt.savefig(SText)
    plt.show()
 
def Plotter_Roh_txt(Pathi,Datei):
    Path =  Pathi + "/"  + Datei + ".txt"
    length = []
    strain = []
    with open(Path, "r") as csvdatei:
        csv_reader_object = csv.reader(csvdatei,delimiter='\t')#, quoting=csv.QUOTE_NONNUMERIC)
    
        zeilennummer = 0

        for row in csv_reader_object:
    
            if zeilennummer == 0:
                print(f'Spaltennamen sind: {", ".join(row)}')
            #else:
                #print(f'- Nachname: {row[0]} \t| Vorname: {row[1]} \t| Geburtstag: {row[2]}.')
            #    print(f'Spaltennamen sind: {", ".join(row)}')
            zeilennummer += 1
            if zeilennummer >= 15:
               # print(f'Spaltennamen sind: {", ".join(row)}')
                length.append(float(row[0]))
                strain.append(float(row[1]))
                
            
        print(f'Anzahl Datensätze: {zeilennummer-1}')
        print(length[0],strain[0])
        #-1.06840516E-005	-3.31337585E+001
    plt.figure(0)
    plt.plot(length,strain)
    plt.xlabel("length [m]")
    plt.ylabel("Amplitude (dB/mm)")
    TText = "Measurment: " + Datei
    plt.title(TText)
    #plt.legend()
    #plt.xlim(1.9,6)
    plt.grid()
    SText = "Pictures_Raw_Data/" + Datei + ".pdf"
    plt.savefig(SText)
    plt.show()
#Plotter("-10oC_0V_Lower")
    
def Data(Pathi, Datei):
    Path = "/home/kappe/Desktop/Masterarbeit/OFM/First Measurment/"+ Pathi + "/" + Datei + ".csv"
    length = []
    strain = []
    with open(Path, "r") as csvdatei:
        csv_reader_object = csv.reader(csvdatei,delimiter='\t')#, quoting=csv.QUOTE_NONNUMERIC)
    
        zeilennummer = 0
        for row in csv_reader_object:
            zeilennummer += 1
            if zeilennummer >= 23:
               # print(f'Spaltennamen sind: {", ".join(row)}')
                length.append(float(row[0]))
                strain.append(float(row[1]))
    return [length,strain]
                
    
def Data_Strain(Pathi, Datei):
    Path = Pathi + "/"  + Datei       
    length, strain = Get_Data(Path)
    return length,strain
# Rohdaten



def Roh_Data(Pathi, Datei):
    Path = "/home/kappe/Desktop/Masterarbeit/OFM/First Measurment/"+ Pathi + "/" + Datei + ".csv"
    length = []
    strain = []
    with open(Path, "r") as csvdatei:
        csv_reader_object = csv.reader(csvdatei,delimiter='\t')#, quoting=csv.QUOTE_NONNUMERIC)
    
        zeilennummer = 0
        for row in csv_reader_object:
            zeilennummer += 1
            if zeilennummer >= 15:
               # print(f'Spaltennamen sind: {", ".join(row)}')
                length.append(float(row[0]))
                strain.append(float(row[2]))
    return [length,strain]      


def Roh_Data_txt(Pathi, Datei):
    Path = Pathi + "/" + Datei + ".txt"
    length = []
    strain = []
    with open(Path, "r") as csvdatei:
        csv_reader_object = csv.reader(csvdatei,delimiter='\t')#, quoting=csv.QUOTE_NONNUMERIC)
    
        zeilennummer = 0
        for row in csv_reader_object:
            zeilennummer += 1
            if zeilennummer >= 15:
               # print(f'Spaltennamen sind: {", ".join(row)}')
                length.append(float(row[0]))
                strain.append(float(row[2]))
    return [length,strain]   


