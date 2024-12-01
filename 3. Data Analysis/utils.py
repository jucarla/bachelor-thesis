import os
import pandas as pd
import plotly.express as px
import plotly.express as px
from dash import Dash, dcc, html, Input, Output

def concat_dataframes(folder_path: str) -> pd.DataFrame:
    """
    Concatenates dataframes from CSV files in the specified folder path.

    Parameters:
    - folder_path (str): Path to the folder containing CSV files.

    Returns:
    - pd.DataFrame: Concatenated dataframe containing data from all CSV files.

    Note:
    - CSV files are expected to have two columns named 'x' and 'y' separated by semicolons (';').
    - Dataframes are filtered to include rows where the 'x' column values are greater than 7.
    - The 'timestep' column is created to indicate the time step based on the filename.
    - Dataframes are concatenated, ordered by 'timestep' and 'x' columns, and null values are removed.
    """
    # Create dataframe list to allocate each of the files dataframes
    dataframes_list = []

    #Iterate through the files in the directory to crate a datframe for each one
    for file in os.listdir(folder_path):
        if file.endswith('.csv'):
            file_path = os.path.join(folder_path, file)
            df = pd.read_csv(file_path, sep = ";", names=["x","y"])
            #df = df[(df.x>7)]

            # Create a new column to reference which time step (contained in the file title) the results are
            df['timestep'] = int(os.path.splitext(file)[0].split("_")[-1:][0])
            dataframes_list.append(df)

    #Concat all dataframes in the list
    dataframe_final = pd.concat(dataframes_list, ignore_index=True)

    #Order dataframe
    dataframe_final.sort_values(by=["timestep","x"], inplace=True)
    #Remove null values
    dataframe_final.dropna(inplace=True)
    return dataframe_final

def new_concat_dataframes(folder_path: str, range: list = None) -> pd.DataFrame:
    """
    Concatenates dataframes from CSV files in the specified folder path.

    Parameters:
    - folder_path (str): Path to the folder containing CSV files.

    Returns:
    - pd.DataFrame: Concatenated dataframe containing data from all CSV files.

    Note:
    - CSV files are expected to have two columns named 'x' and 'y' separated by semicolons (';').
    - Dataframes are filtered to include rows where the 'x' column values are greater than 7.
    - The 'timestep' column is created to indicate the time step based on the filename.
    - Dataframes are concatenated, ordered by 'timestep' and 'x' columns, and null values are removed.
    """
    # Create dataframe list to allocate each of the files dataframes
    dataframes_list = []

    #Iterate through the files in the directory to crate a datframe for each one
    for file in os.listdir(folder_path):
        if file.endswith('.csv'):
            file_path = os.path.join(folder_path, file)
            df = pd.read_csv(file_path, sep = ";", names=["x","y"])

            if range:
                df = df[(df['x'] >= range[0]) & (df['x'] <= range[1])]
                df= apply_iqr_filter(df)
            #df = df[(df.x>7)]

            # Create a new column to reference which time step (contained in the file title) the results are
            df['timestep'] = int(os.path.splitext(file)[0].split("_")[-1:][0])
            dataframes_list.append(df)

    #Concat all dataframes in the list
    dataframe_final = pd.concat(dataframes_list, ignore_index=True)

    #Order dataframe
    dataframe_final.sort_values(by=["timestep","x"], inplace=True)
    #Remove null values
    dataframe_final.dropna(inplace=True)
    return dataframe_final

def filter_by_timestep(dataframe: pd.DataFrame, timestep) -> pd.DataFrame:
    """
    Filter the dataframe to select rows corresponding to the specified timestep(s).

    Parameters:
    - dataframe (pd.DataFrame): The input dataframe to be filtered.
    - timestep (int or list of int): The timestep(s) to filter the dataframe.

    Returns:
    - pd.DataFrame: Filtered dataframe containing rows corresponding to the specified timestep(s).

    Note:
    - 'dataframe' is expected to contain a 'timestep' column indicating the timestep of each row.
    - 'timestep' can be a single integer or a list of integers representing the timestep(s) to be selected.
    """
    # Filter df
    mask = dataframe['timestep'].isin(timestep)
    result = dataframe[mask]
    
    return result

def upper_and_lower_limit(df: pd.DataFrame):
    """
    Calculate the upper and lower limits based on the interquartile range (IQR).

    Parameters:
    - df (pd.DataFrame): The input dataframe.

    Returns:
    - float: Lower limit calculated using the IQR method.
    - float: Upper limit calculated using the IQR method.

    Note:
    - 'df' is expected to contain numerical data in the column 'y'.
    - The lower and upper limits are calculated based on the first quartile (Q1),
      third quartile (Q3), and interquartile range (IQR) of the 'y' column.
    - The lower limit is calculated as Q1 - 1.5*IQR, and the upper limit is calculated as Q3 + 1.5*IQR.
    """
    Q1 = df.y.quantile(0.25)
    Q3 = df.y.quantile(0.75)
    IQR = Q3 - Q1
    lower_limit = Q1 - 1.5*IQR
    upper_limit = Q3 + 1.5*IQR

    return lower_limit, upper_limit

def apply_iqr_filter(df: pd.DataFrame):
    """
    Apply the interquartile range (IQR) filter to remove outliers from the dataframe.

    Parameters:
    - df (pd.DataFrame): The input dataframe.

    Returns:
    - pd.DataFrame: The filtered dataframe with outliers removed.

    Note:
    - 'df' is expected to contain numerical data in the column 'y'.
    - Outliers are identified and removed using the upper and lower limits calculated
      based on the interquartile range (IQR) method.
    - The lower and upper limits are calculated using the 'upper_and_lower_limit' function.
    - Outliers are defined as data points outside the range [lower_limit, upper_limit].
    - The number of outliers removed is printed as part of the process.
    """
    lower_limit, upper_limit = upper_and_lower_limit(df)
    print(f"Lower limit: {lower_limit}, Upper Limit: {upper_limit}")

    df_no_outlier = df[(df.y>lower_limit)&(df.y<upper_limit)]
    print(f"{df.shape[0] - df_no_outlier.shape[0]} outliers removed")
    return df_no_outlier

def calculate_mean_per_timestep(df: pd.DataFrame):
    """
    Calculate the mean spectral shift per timestep from the dataframe.

    Parameters:
    - df (pd.DataFrame): The input dataframe containing the data.

    Returns:
    - pd.DataFrame: Dataframe containing the mean spectral shift per timestep.

    Note:
    - 'df' is expected to have columns 'x', 'y', and 'timestep'.
    - The mean spectral shift is calculated for each timestep by grouping the dataframe
      by the 'timestep' column and computing the mean of the 'y' column.
    - The resulting dataframe contains two columns: 'timestep' and 'mean_spectral_shift'.
    """
    df_mean = df.groupby('timestep').mean().rename(columns={'y': 'dv'}).reset_index().drop('x', axis=1)
    return df_mean


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

def transform_csv_to_dataframe(csv_file_path: str) ->pd.DataFrame:
    """
    Transform a CSV file into a pandas DataFrame.

    Parameters:
    - csv_file_path (str): The path to the CSV file.

    Returns:
    - pd.DataFrame: A DataFrame containing the transformed data.

    Note:
    - The CSV file is expected to have three columns: 'timestamp', 'module1_temp', and 'module2_temp'.
    - The 'timestamp' column represents the timestamp of the data.
    - The 'module1_temp' column contains temperature data for module 1 (which is irrelevant for us).
    - The 'module2_temp' column contains temperature data for module 2 (which is irrelevant for us).
    - Only lines where 'Module 3 -' is found in the 'module2_temp' column are considered.
    """
    # Define the column names for the dataframe
    columns = ['timestamps', 'temp']
    print(os.getcwd())
    # Initialize an empty list to store the data
    data = []

    # Read the CSV file line by line
    with open(csv_file_path, 'r') as file:
        for line in file:
            # Split each line by comma
            parts = line.strip().split(',')
            if len(parts) == 3:
                # Extract timestamp
                timestamp = parts[0]

                if "Module 3 -" in parts[1]:
                # Extract temperature for each module
                    module3_temp = float(parts[1].split('=')[1].strip().split()[0])
                    # Append the extracted data to the list
                    data.append([timestamp, module3_temp])

    # Create a DataFrame from the list of data
    df = pd.DataFrame(data, columns=columns)

    return df
    
def get_spectral_shift(temp:float, t0: float = 29.9)->float:
    """
    Calculate the spectral shift based on temperature difference.

    Parameters:
    - temp (float): The temperature value.
    - t0 (float, optional): Reference temperature (default: 29.9).

    Returns:
    - float: The calculated spectral shift.

    Note:
    - The spectral shift is calculated based on the temperature difference between 'temp' and 't0'.
    - The calculation formula is applied: spectral_shift = (temp - t0) * (0.1306 * -19.542).
    """
    deltat = temp - t0
    
    spectral_shift = deltat * (0.1306*-19.542)
    return spectral_shift

def get_strain(spectral_shift: float) -> float:
    """
    Calculate the strain based on the spectral shift.

    Parameters:
    - spectral_shift (float): The spectral shift value.

    Returns:
    - float: The calculated strain.

    Note:
    - The strain is calculated based on the provided spectral shift value.
    - Constants used in the calculation:
        - Speed of light (c) = 299792458 m/s
        - Center wavelength (center_wl) = 1306 nm
        - Effective Kerr constant (ke) = (0.471 + 0.4985) / 2
    - The calculation formula is applied: strain = - (spectral_shift * center_wl) / (c * ke).
    """
    c = 299792458
    center_wl = 1306
    ke = (0.471 +0.4985)/2
    strain = - (spectral_shift*center_wl) / (c*ke)
    return strain

def get_temp(self, spectral_shift:int, t0):
        #deltaat = spectral_shift / (0.1306*-19.542)
        #t0 = 29.9
        c = 299792458
        center_wl = 1306
        kt=1.9542*(10**-5)
        deltaat = - (spectral_shift*center_wl) / (c*kt)
        #print(deltaat)

        new_temp = deltaat + t0
        return new_temp