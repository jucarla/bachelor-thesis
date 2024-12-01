import serial
import csv
import time
from datetime import datetime
import os

ser = serial.Serial('COM8', 115200, timeout=1)
time.sleep(2)  

current_date = datetime.now().strftime('%Y-%m-%d')

base_filename = f'data/temperature_data_{current_date}_50c'

file_index = 1
filename = f'{base_filename}.csv'
while os.path.exists(filename):
    filename = f'{base_filename}_{file_index}.csv'
    file_index += 1

with open(filename, 'a', newline='') as file:
    writer = csv.writer(file)
    
    if file.tell() == 0:
        writer.writerow(['timestamps', 'temp'])

    while True:
        try:
            line = ser.readline().decode('utf-8').strip()

            if line.startswith("TEMP:"):
                temperature = line.split(":")[1]

                timestamp = datetime.now().strftime("%Y-%m-%d_%H:%M:%S")

                print(f"{timestamp} - Temperature: {temperature} C")
                writer.writerow([timestamp, temperature])
                

        except Exception as e:
            print(f"Error: {e}")
            break

ser.close()
