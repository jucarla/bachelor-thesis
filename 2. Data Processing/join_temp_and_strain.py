import pandas as pd

def merge_files_with_tolerance(temperature_file, strain_file, output_with_nan="final validation/merged_with_nan.csv", output_without_nan="final validation/merged_without_nan.csv", tolerance_seconds=2):
    """
    Merges two files on nearest timestamp within a specified tolerance.
    
    Parameters:
    - temperature_file: str, path to the temperature CSV file.
    - strain_file: str, path to the strain CSV file.
    - output_with_nan: str, path to save the merged file with NaN values.
    - output_without_nan: str, path to save the merged file without NaN values.
    - tolerance_seconds: int, tolerance in seconds for nearest timestamp match.
    """
    temperature_data = pd.read_csv(temperature_file)
    strain_data = pd.read_csv(strain_file)

    temperature_data['temp'] = pd.to_numeric(temperature_data['temp'], errors='coerce')

    temperature_data['temp']=temperature_data['temp']+273
    temperature_data['timestamps'] = pd.to_datetime(temperature_data['timestamps'], format='%Y-%m-%d_%H:%M:%S')
    strain_data['Timestamp'] = pd.to_datetime(strain_data['Timestamp'], format='%Y-%m-%d %H:%M:%S')
    
    merged_data = pd.merge_asof(
        temperature_data.sort_values('timestamps'),
        strain_data.sort_values('Timestamp'),
        left_on='timestamps',
        right_on='Timestamp',
        direction='nearest',
        tolerance=pd.Timedelta(seconds=tolerance_seconds)
    )
    
    merged_data = merged_data.drop(columns=['Timestamp']).rename(columns={'timestamps': 'Timestamp', 'temp': 'Temperature'})
    
    merged_data['Strain'] = merged_data['Strain'].map(lambda x: '{:.10f}'.format(x) if pd.notnull(x) else x)
    
    merged_data.to_csv(output_with_nan, index=False)
    merged_data.dropna().to_csv(output_without_nan, index=False)
    
    print(f"Merged file with NaN values saved as: {output_with_nan}")
    print(f"Merged file without NaN values saved as: {output_without_nan}")

if __name__ == '__main__':
    merge_files_with_tolerance("final validation/temperature_data_2024-11-25.csv", "final validation/converted_joined_validation.csv")
