import Einlesen as Ein
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


############# Constant definition
c = 299792458.0 #m/s
Lambda = 1310.0  #1550.0# * 10**-9 #m

#Luna Parameter
K_S_C = 0.78 
K_S_M = 0.81
K_S_M_std = K_S_M*0.01
K_T_C =   6.45*10**(-6) #1/C
k_S_C = -6.668 #-Lambda/K_S_C/c *10**6  # -6.67#muStrain/GHz # used are -6.668
k_T_C = Lambda/K_T_C/c # -0.801 #C/GHz


def MakeData_v1(PT,timestamp_path,P2,N_SE,Text,PPosition,PTNummer,NameWU ="test1", D_Liste = [],IQR = True,StufeGB = [0,0],List_T = None, T_Ref_Differ_from_file_data = None,BlackList = None):
    D_L = []
    for i in range(N_SE[0],N_SE[1]):
        D_L.append(NameWU + "_" + str(i))
    
    if len(D_Liste) != 0:
        D_L = D_List         

    Ein.Plotter_several_with_Name_vLine_Filter(P2,[],Text + "_Plot",D_L,D_L,[PPosition[0]-0.1,PPosition[1]+0.1],[-2,5],PPosition,IQR,StufeGB[0],StufeGB[1])

    if List_T == None:
        #N_SE - 
        print(".")
        WA, WAS, TS = Ein.Get_Data_WarmUp_v4(P2,PT,timestamp_path,NameWU,D_Liste,PPosition,N_SE[0],N_SE[1],IQR,StufeGB[0],StufeGB[1])
        print(f"DEBUG - WA:{WA}")
        print(f"DEBUG - WAS:{WAS}")
        dv = np.array(WA[-1])
        dv_std = np.abs(np.array(WAS[-1]))
    
        T = np.array(WA[PTNummer]) + 273.15
        T_std = np.array(WAS[PTNummer])

    dv = np.array(dv)
    dv_std = np.array(dv_std)

    print(f"DEBUG - dv:{dv}")
    print(f"DEBUG - dv_std:{dv_std}")


    T_sys = np.zeros(len(T)) + 0.1
    
    if T_Ref_Differ_from_file_data != None:
        Ref_T = T_Ref_Differ_from_file_data[0]
        Ref_T_std = T_Ref_Differ_from_file_data[1]

    DMS_original = pd.read_csv(PT+".csv", usecols=['Strain'])['Strain'].to_numpy()
    DMS = WA[1]

    DMS_std =  np.full(len(DMS), np.std(DMS, ddof=1) / np.sqrt(len(DMS)))
    DMS_sys = np.zeros(len(DMS)) + 0.1

    np.savetxt(Text + '.out', (dv, dv_std, T, T_std, DMS, DMS_std, DMS_sys, T_sys)) 
    df_plot = pd.DataFrame({
        'Timestamp': TS,
        'dv': dv,
        'dv_std': dv_std,
        'T': T,
        'T_std': T_std,
        'T_sys': T_sys,
        'DMS': DMS
    })

    df_plot.to_csv(Text + '.csv', index=False, float_format='%.10f')

    NEnde = len(T) #200
    print(NEnde)
    for i in range(len(T)):
        if T[i]>= (-39.9 + 273.15):
            print("WU at",-39.9 + 273.15,  " :",i, "  ",T[i])
            break
    
    
    plt.figure(1234775)
    plt.plot(T_std[:NEnde], label = "stat")
    plt.plot(T_sys[:NEnde], label = "sys")
    plt.title("T_M error")
    plt.legend()
    plt.show()

    plt.figure()
    plt.plot(DMS_original, label = "stat")
    plt.title("DMS")
    plt.legend()
    plt.show()

    plt.figure()
    plt.plot(-dv, label = "stat")
    plt.title("DV")
    plt.legend()
    plt.show()

    plt.figure(12347752)
    plt.plot(dv_std[:NEnde]*10**6, label = "G stat")
    plt.title("dv error")
    plt.legend()
    plt.show()

    print(DMS)
    print(T)
    plt.figure(123543,dpi = 180)
    plt.plot((T - T[0])/np.max(np.abs(T - T[0])), label = "T")
    plt.plot(-dv[:NEnde]/np.max(np.abs(dv[:NEnde])), label = "-dv")
    plt.plot(DMS/np.max(np.abs(DMS)), label = "DMS")
    plt.grid()
    plt.legend()
    plt.xlabel("Time [m]")
    plt.ylabel("Norm. Signal")
    #plt.savefig("Plot_All_Data_OBR.png")
    plt.show()
    
    return 0

obr_csv_output_path = 'M:/Intern/OE300/E2_OE320/85_Hiwi_Ordner/gmn-js/20240703_Bachelorarbeit/Results/wednesday/joined/CSV'
timestamps = 'M:/Intern/OE300/E2_OE320/85_Hiwi_Ordner/gmn-js/20240703_Bachelorarbeit/Results/wednesday/joined/timestamps_test.txt'
Step_PT = "M:/Intern/OE300/E2_OE320/85_Hiwi_Ordner/gmn-js/20240703_Bachelorarbeit/Results/wednesday/joined/merged_without_nan"


timestamp_path = "Data_Test/new_data/timestamps_test.txt"
P2 = "Data_Test/new_data/CSV3"
PT = "Data_Test/new_data/merged_without_nan"
D_List = ["test1"]
#BOR - 2.3, 2.5
# GER - 1.6 , 1.8
P_X_A = [2.3, 2.5]
Text = "new plot- wed data - bor"

MakeData_v1(PT,timestamp_path,P2,[1,129],Text,P_X_A,0,NameWU = "test",StufeGB = [0,0])

