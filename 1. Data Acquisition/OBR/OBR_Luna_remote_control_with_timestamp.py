import time
import socket
import datetime
from pathlib import Path
from subprocess import Popen

class Instrument:
    
    def __init__(self,IP,PORT):
        self.Socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.Socket.connect((IP,PORT))  

    def write(self,CMD): 
        self.Socket.send((CMD+'\n').encode('ascii'))

    def query(self,CMD):
        self.write(CMD)
        time.sleep(1)
        
        reply = self.Socket.recv(1024).strip().decode('ascii')
        return reply


        


n_measurements=1 #measurements performed in a block 
N_measurements=5 #total amount of measurement blocks performed
sleep_time=1 # sleep time in seconds between the two blocks
exp_name = 'test' #name that appears on the OBR file and the TXT file after performing scanning
DirRoot = "C:\\REMOTE_CONTROL_TESTS"
TestName = "\\20242511DecouplingValidationBattery"
FileNameRoot = DirRoot + TestName# save directory on the LUNA laptop


#For direct connection
#LUNA = Instrument('192.168.1.4',80)

#For IPT Network connection (Shared Folder)
LUNA = Instrument('10.100.160.72', 80) 

print(LUNA.query('*IDN?'))
time.sleep(1)
start_time=time.time()




LUNA.write('CONF:XUNI 1')
time.sleep(1)

LUNA.write('CONF:INTW {}'.format(60))
time.sleep(1)
LUNA.write('CONF:INTC {}'.format(30))
time.sleep(1)


print(LUNA.query('CONF:SRBW?'))
time.sleep(1)
print(LUNA.query('CONF:XUNI?'))
time.sleep(1)
print(LUNA.query('CONF:INTW?'))
time.sleep(1)


    
Path(FileNameRoot + '/OBR_files').mkdir(parents=True, exist_ok=True)
Path(FileNameRoot + '/TXT_files').mkdir(parents=True, exist_ok=True)
Path(FileNameRoot + '/Timestamps').mkdir(parents=True, exist_ok=True)

#separate path where a folder 'time_stamp' is created 
file_path = Path(FileNameRoot + '/Timestamps/timestamps_test.txt')
file = open(file_path, 'w')
print(f"Opening {file_path} to save timestamps.")

for i in range(N_measurements):
    print(i)
    file.write('\t'+str(datetime.datetime.now().strftime("%Y-%m-%d_%H:%M:%S"))+'\n')
    
    
    for j in range(n_measurements):
        
        print(j)
        LUNA.write('SCAN')
        #LUNA.write('SYST:SPOT {}, 1'.format(2.0))
        time.sleep(1)
        
        
            
        for k in range(50):
            
            OPC = LUNA.query('*OPC?')
            print(OPC)
            time.sleep(0.2)
            if "1" in OPC: 
                
               
                #save obr file
                filename = "%s_%d.obr" % (FileNameRoot + '/OBR_files/' + exp_name, j + n_measurements * i)
                print("save file", filename)
                CMD = 'SYST:SAVE "%s",A' % (filename)
                LUNA.write(CMD)
                time.sleep(1)
                
                #save txt file
                filename = "%s_%d.txt" % (FileNameRoot + '/TXT_files/' + exp_name, j + n_measurements * i)
                print("save file", filename)
                CMD = 'SYST:SAVT "%s",A' % (filename)
                LUNA.write(CMD)
                time.sleep(1)
                
                CMD='FETC:XAXI? A, 0, 1, 0'
                Xaxis = LUNA.query(CMD)
                time.sleep(1)
                print(Xaxis)
                
                CMD = "FETC:MEAS? A, 0, 4, 0"
                print(LUNA.query(CMD))
                time.sleep(1)
                
                
                
                break
    
    
    file.write('\t'+'\t'+str(datetime.datetime.now().strftime("%Y-%m-%d_%H:%M:%S").rstrip('0'))+'\n')
    time.sleep(sleep_time)
    
    

#%% end

file.close()
print(f"Timestamps successfully saved.")
LUNA.Socket.close();


# %%
