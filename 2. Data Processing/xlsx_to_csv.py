from xlsx2csv import Xlsx2csv
import pandas as pd

def excel_to_csv_with_date_filter(excel_path, csv_path="output.csv"):
    """
    Convert an Excel file to a CSV file and filter rows by a target date.
    
    Parameters:
    - excel_path: str, path to the Excel file.
    - csv_path: str, path to save the CSV file (default is 'output.csv').
    - target_date: str, date in 'MM/DD/YYYY' format to filter rows (e.g., '11/09/2024').
    """
    try:
        temp_csv_path = "original_" + excel_path.split("/")[-1].split(".")[0] + ".csv"
        Xlsx2csv(excel_path, outputencoding="utf-8").convert(temp_csv_path)
        
        df = pd.read_csv(temp_csv_path)
        df['Time'] = pd.to_datetime(df['Time'], format='%m/%d/%Y %H:%M:%S.%f')

        df['Time_Second'] = df['Time'].dt.strftime('%Y-%m-%d %H:%M:%S')

        means_df = df.groupby('Time_Second')['Strain (Filtered)'].mean().reset_index()
        means_list = means_df['Strain (Filtered)'].tolist()
        time_list = means_df['Time_Second'].tolist()
        
        df_date = pd.DataFrame({'Time': time_list, "Strain": means_list})
        
        df_date.to_csv(csv_path, index=False, header=False, float_format='%.10f')
        
        print(f"Excel file has been successfully converted to {csv_path} with date filter applied.")
        
    except Exception as e:
        print("An error occurred:", e)
if __name__ == '__main__':
    excel_to_csv_with_date_filter("M:/Intern/OE300/E2_OE320/85_Hiwi_Ordner/gmn-js/20240703_Bachelorarbeit/Results/final validation/20242511_battery_50c.xlsx", csv_path="converted_strain_file_battery_50c.csv")


