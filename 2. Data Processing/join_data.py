from dataprocessing import OBR, PT1000, Results
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def read_and_save_data(obr_csv_output_path: str, strain_and_temp_output_path: str, timestamps: str, range: list, output_name: str):
    obr_result = OBR(obr_csv_output_path, "Save",range=range)
    obr_result.plot_animated_chart("Obr result")

    strain_and_temp_output = PT1000(strain_and_temp_output_path)
    strain_and_temp_output.plot("Temperature")

    results = Results(obr_result, strain_and_temp_output,timestamps, "Save")

    plt.figure(figsize=(12, 6))
    plt.scatter(results.df['timestep'], results.df['dv'], label='dv', alpha=0.7)

    # Customizing the plot
    plt.xlabel("Timestep")
    plt.ylabel("dv")
    plt.title("Evolution of 'dv' over timesteps")
    plt.legend()
    plt.xticks(rotation=45)
    plt.tight_layout()

    # Display the plot
    plt.show()

    plt.figure(figsize=(12, 6))
    plt.scatter(results.df['Temperature'], results.df['dv'], label='dv', alpha=0.7)

    # Customizing the plot
    plt.xlabel("Timestep")
    plt.ylabel("dv")
    plt.title("Evolution of 'dv' over temperatur")
    plt.legend()
    plt.xticks(rotation=45)
    plt.tight_layout()

    # Display the plot
    plt.show()

    # Calculate errors as the absolute difference between consecutive values
    results.df['dv_error'] = results.df['dv'].diff().abs()
    results.df['temperature_error'] = results.df['Temperature'].diff().abs()
    results.df['strain_error'] = results.df['Strain'].diff().abs()
    
    dv_std = results.df['dv'].std()
    temperature_std = results.df['Temperature'].std()
    strain_std = results.df['Strain'].std()
    dv_sys = (results.df['dv'] - results.df['dv'].mean()).abs().mean()
    temperature_sys = (results.df['Temperature'] - results.df['Temperature'].mean()).abs().mean()
    strain_sys = (results.df['Strain'] - results.df['Strain'].mean()).abs().mean()

    results.df['dv_error_std'] = dv_std
    results.df['dv_error_sys'] = dv_sys

    results.df['temperature_error_std'] = temperature_std
    results.df['temperature_error_sys'] = temperature_sys

    results.df['strain_error_std'] = strain_std
    results.df['strain_error_sys'] = strain_sys
    # Filter out rows where any of the batched errors exceed 1
    results.df = results.df[
        (results.df['dv_error'] <= 1) &
        (results.df['temperature_error'] <= 1) &
        (results.df['strain_error'] <= 1)
    ]
    # Save the filtered data to a new CSV file
    results.df.to_csv(output_name + '.csv', index=False)

    dv = results.df['dv'].to_numpy()
    dv_std = results.df['dv_error_std'].to_numpy()
    T = results.df['Temperature'].to_numpy()
    T_std = results.df['temperature_error'].to_numpy()
    DMS = results.df['Strain'].to_numpy()
    DMS_std = results.df['strain_error_std'].to_numpy()
    DMS_sys = results.df['strain_error_sys'].to_numpy()

    np.savetxt(output_name + '.out', (dv, dv_std, T, T_std, DMS, DMS_std, DMS_sys, DMS_sys)) 

    print(f"Calculated errors saved to {output_name}.csv and {output_name}.out")


if __name__ == '__main__':
    #obr_csv_output_path = 'M:/Intern/OE300/E2_OE320/85_Hiwi_Ordner/gmn-js/20240703_Bachelorarbeit/Results/saturday/20240911DecouplingTests_Heating_Fiber_1_2/CSV3'
    #timestamps = 'M:/Intern/OE300/E2_OE320/85_Hiwi_Ordner/gmn-js/20240703_Bachelorarbeit/Results/saturday/Timestamps/timestamps_test.txt'
    #strain_and_temp_output_path = 'M:/Intern/OE300/E2_OE320/85_Hiwi_Ordner/gmn-js/20240703_Bachelorarbeit/Results/saturday/merged_without_nan.csv'

    obr_csv_output_path = 'M:/Intern/OE300/E2_OE320/85_Hiwi_Ordner/gmn-js/20240703_Bachelorarbeit/Results/wednesday/joined/CSV'
    timestamps = 'M:/Intern/OE300/E2_OE320/85_Hiwi_Ordner/gmn-js/20240703_Bachelorarbeit/Results/wednesday/joined/timestamps_test.txt'
    strain_and_temp_output_path = 'M:/Intern/OE300/E2_OE320/85_Hiwi_Ordner/gmn-js/20240703_Bachelorarbeit/Results/wednesday/joined/merged_without_nan.csv'
    

    range_bor = [2.3, 2.5]
    range_ger = [1.6 , 1.8]

    output = 'step_test_time_dv_temp_strain_with_error_bor_check_error'
    read_and_save_data(obr_csv_output_path=obr_csv_output_path, 
                       strain_and_temp_output_path=strain_and_temp_output_path, 
                       timestamps=timestamps, output_name=output, range=range_bor)
    

     