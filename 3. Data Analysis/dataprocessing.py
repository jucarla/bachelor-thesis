from utils import get_spectral_shift, get_strain, calculate_mean_per_timestep, \
    transform_txt_to_list, transform_csv_to_dataframe,\
    concat_dataframes, new_concat_dataframes, apply_iqr_filter, upper_and_lower_limit, \
    filter_by_timestep, calculate_mean_per_timestep, transform_txt_to_list, \
    transform_csv_to_dataframe


import plotly.express as px
from dash import Dash, dcc, html, Input, Output
import pandas as pd
import os
import numpy as np

class OBR:
    def __init__(self, csv_folder_path: str, results_path:str, save=None, iqr=None, range=None):
        self.saving_path = results_path
        self.df = new_concat_dataframes(csv_folder_path, range)
        if iqr:
            self.df= apply_iqr_filter(self.df)
        if save:
            self.save_plots()
        #if range:
        #    self.df = self.df[(self.df['x'] >= range[0]) & (self.df['x'] <= range[1])]
        #    self.df= apply_iqr_filter(self.df)

    def plot_animated_chart(self, title: str, delay: int = 1000):
        """
        Plot an animated line chart from the dataframe.

        Parameters:
        - df (pd.DataFrame): The input dataframe containing the data to be plotted.
        - title (str): The title of the chart.
        - delay (int, optional): Delay between frames in milliseconds (default: 1000).

        Returns:
        - None

        Note:
        - 'df' is expected to have columns 'x', 'y', and 'timestep'.
        - The chart is animated based on the 'timestep' column.
        - The chart's title is specified by the 'title' parameter.
        - The y-axis range of the chart is automatically set based on the upper and lower limits
        calculated using the 'upper_and_lower_limit' function.
        - The 'delay' parameter specifies the duration between frames in milliseconds.
        """
        lower_limit, upper_limit = upper_and_lower_limit(self.df)
        chart_range = [lower_limit, upper_limit]
        graph = px.line(self.df,x=["x"], y=self.df["y"],  title=title, animation_frame=self.df["timestep"], animation_group=self.df["timestep"], 
                        range_y = chart_range, labels={'x':'Fiber Lenght (m)',
                                        'y':'Spectral Shift'})
        graph.layout.updatemenus[0].buttons[0].args[1]["frame"]["duration"] = delay
        graph.show()

    def plot_overlaped_chart(self, title: str):
        """
        Plot an overlapped line chart using Dash.

        Parameters:
        - df (pd.DataFrame): The input dataframe containing the data to be plotted.
        - title (str): The title of the chart.

        Returns:
        - None

        Note:
        - This function uses Dash to create an interactive web application for visualizing overlapped line charts.
        - 'df' is expected to have columns 'x', 'y', and 'timestep'.
        - Users can select which timesteps to display on the chart using a dropdown menu.
        """
        app = Dash(__name__)

        app.layout = html.Div([
            dcc.Dropdown(
                id='selection',
                options=self.df["timestep"].unique(),
                value=[1,2],
                multi=True
            ),
            dcc.Loading(dcc.Graph(id="graph"), type="cube")
        ])
        @app.callback(
            Output("graph", "figure"), 
            Input("selection", "value"))
        def display_animated_graph(selection):
            delta_label = r"$-\frac{\Delta \nu}{\hat{\nu}}_{Thermal} \times 10^{-6}$"
            df_filtered = filter_by_timestep(self.df, selection)
            graph = px.line(df_filtered,x=["x"], y=df_filtered["y"],  title=title, color=df_filtered["timestep"],
                                labels={'value':'Fiber Lenght (m)',
                                        'y':'Spectral Shift'})
            
            graph.update_layout(
            template="plotly_white",
            plot_bgcolor="white",
            paper_bgcolor="white"
            )
            return graph

        app.run_server(debug=True)
    
    def save_plots(self):
        """
        Save individual plots for each timestep in the dataframe.

        Parameters:
        - df (pd.DataFrame): The input dataframe containing the data.
        - saving_path (str): The path where the plots will be saved.

        Returns:
        - None

        Note:
        - 'df' is expected to have columns 'x', 'y', and 'timestep'.
        - Individual plots are generated for each timestep, and saved as PNG files in the specified 'saving_path'.
        """
        path =self.saving_path + "/OBR_plots"
        for test in self.df["timestep"].unique():
            df_filtered = filter_by_timestep(self.df, [test])
            graph = px.line(df_filtered,x=["x"], y=df_filtered["y"],  title=f"Temperature response at instant {test}", labels={'x':'Fiber Lenght (m)',
                                        'y':'Spectral Shift'})
            if not os.path.exists(path):
                os.mkdir(path)
            graph.write_image(f"{path}/plot_{test}.png")
            print(f"plot_{test}.png successfully saved")

    def plot_mean_per_timestep(self, leftlimit: int, rigthlimit: int, title:str):
        """
        Plot the mean spectral shift per timestep within a specified range of fiber length.

        Parameters:
        - leftlimit (int): The left limit of the fiber length range.
        - rightlimit (int): The right limit of the fiber length range.
        - df (pd.DataFrame): The input dataframe containing the data.
        - title (str): The title of the chart.

        Returns:
        - None

        Note:
        - 'df' is expected to have columns 'x', 'y', and 'timestep'.
        - The dataframe is filtered to include rows where the 'x' column values are within the specified range.
        - The mean spectral shift per timestep is calculated using the 'calculate_mean_per_timestep' function.
        - The resulting mean spectral shifts are plotted against the timesteps.
        - The title of the chart includes the specified fiber length range.
        """
        df = df[(self.df.x>leftlimit)&(self.df.x<rigthlimit)]
        df_mean = calculate_mean_per_timestep(self.df)
        graph = px.line(df_mean,x=df_mean["timestep"], y=df_mean['dv'],  title=title + f" [Range {leftlimit}m to {rigthlimit}m]", labels={
                                        'dv':'Mean Spectral Shift'})
        graph.show()

    

class PT1000:
    def __init__(self, pt1000_path: str):
        df_pt1000 = pd.read_csv(pt1000_path)
        self.df = df_pt1000
        print(df_pt1000.head())

    def plot(self, title):
        graph = px.line(self.df,x=self.df["timestamps"], y=['Temperature'],  title=title)
        graph.show()

class Results(OBR, PT1000):
    def __init__(self, obr_object, pt1000_object, txt_path: str, results_path:str):

        OBR_df = obr_object.df
        PT1000_df = pt1000_object.df
        timestamps = transform_txt_to_list(txt_path)
        df_mean = calculate_mean_per_timestep(OBR_df)
        self.saving_path = results_path

        #Calculate mean spectral shift per timestep
        df_mean = OBR_df.groupby('timestep').mean().rename(columns={'y': 'dv'}).reset_index().drop('x', axis=1)
        # Create new column in the dataframe, adding the respective timestamp to each timestep
        df_mean['timestamps'] = df_mean['timestep'].map(lambda x: timestamps[x - 1] if x <= len(timestamps) else None)
        # Match timestamp formating to fit the PT1000 left join
        df_mean['timestamps'] = df_mean['timestamps'].str.replace('_', ' ')
        PT1000_df['timestamps'] = PT1000_df['timestamps'].str.replace('_', ' ')

        # For each timestamp, join PT1000 temperature data with the spectral shift data on that timestamp
        right_join_df= pd.merge(df_mean, PT1000_df, on='timestamps', how='left')
        right_join_df.interpolate(method='polynomial', order=2, inplace=True)

        self.df = right_join_df

    def plot(self, columns: list, title: str, ylabel: str, save = False):
        graph = px.scatter(self.df,x=self.df["timestep"], y=columns,  title=title, labels={
                                        'value': ylabel})
        if save:
            path = self.saving_path + "/plots"
            if not os.path.exists(path):
                os.mkdir(path)
            graph.write_image(f"{path}/{title}.png")
            print(f"{title}.png successfully saved")
        graph.show()

    def solve_delta(self, delta_lambda, Ke1, Ke2, Kt1, Kt2):
        """
        Solves for delta epsilon and delta T given delta lambda and constants Ke1, Ke2, Kt1, Kt2.

        Parameters:
        - delta_lambda (np.array): A numpy array with values [delta_lambda_1, delta_lambda_2]
        - Ke1, Ke2, Kt1, Kt2 (float): Constants for the matrix

        Returns:
        - delta_epsilon (float): Computed value of delta epsilon
        - delta_T (float): Computed value of delta T
        """
        # Define the matrix of constants
        K_matrix = np.array([[Ke1, Kt1],
                            [Ke2, Kt2]])

        # Calculate delta epsilon and delta T by inverting the matrix
        delta_values = np.linalg.inv(K_matrix).dot(delta_lambda)

        # Extract delta epsilon and delta T
        delta_epsilon, delta_T = delta_values
        return delta_epsilon, delta_T
    
    def smooth_time_series(self, df, column):
        # Tail-rolling average transform
        rolling = self.df[column].rolling(window=3)
        rolling_mean = rolling.mean()

        self.df[column+"_smooth"]=rolling_mean

    def get_temperature(self, spectral_shift:int, t0: float):
            #deltaat = spectral_shift / (0.1306*-19.542)

            c = 299792458
            center_wl = 1306
            kt=1.9542*(10**-5)
            deltaat = - (spectral_shift*center_wl) / (c*kt)
            #print(deltaat)

            new_temp = deltaat + t0
            return new_temp

    def calculate_temp(self):
        self.smooth_time_series(self.df, "dv")
        self.df['calculated_temp'] = self.df['dv'].map(lambda x: self.get_temperature(x, 29.9))
    
    def calculate_strain(self):
        self.df["calculated_spectral_shift"] = self.df['temp'].map(lambda x: get_spectral_shift(x))
        # Spectral shift from strain:
        self.df["calculated_strain"] = (self.df["dv_smooth"]-self.df["calculated_spectral_shift"]).map(lambda x: get_strain(x))
        self.smooth_time_series(self.df, "calculated_strain")

    def export_results(self, filename:str):
        path = self.saving_path + "/csv_results/"
        if not os.path.exists(path):
                os.mkdir(path)
        self.df.to_csv(path+filename)
