import numpy as np

def solve_delta(delta_lambda, Ke1, Ke2, Kt1, Kt2):
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
        delta_epsilon, delta_T = - delta_values
        return delta_epsilon, delta_T

def solve_delta_without_decoupling(delta_lambda, Kt):
    
    deltaat =  - delta_lambda / Kt
    
    return deltaat
