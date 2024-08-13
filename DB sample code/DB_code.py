# clearing the previous variables 
from IPython import get_ipython
get_ipython().magic('reset -sf')


###################################### Dynamic Budyko model function ############################################################################ 
import numpy as np
import pandas as pd

def DB(R, PET, lag):

    # W and H estimation
    len_R = len(R)
    W = np.zeros(len_R)
    H = np.zeros(len_R)
    
    for iii in range(len_R):
        if R[iii] > PET[iii]:
            W[iii] = R[iii] - PET[iii]
        else:
            W[iii] = 0
    
    for i in range(len_R):
        if PET[i] > R[i]:
            H[i] = PET[i] - R[i]
        else:
            H[i] = 0
    
    # FW and FH estimation
    a = 0.02
    FW = np.zeros(len_R)
    FW[0] = W[0]
    for i in range(1, len_R):
        FW[i] = FW[i-1] * np.exp(-a) + W[i]

    # Avoid zero in FW to prevent division by zero
    FW[FW == 0] = 0.0000001

    FH = np.zeros(len_R)
    FH[0] = H[0]
    for i in range(1, len_R):
        FH[i] = FH[i-1] * np.exp(-a) + H[i]

    # Phi and Effective rainfall
    phi = FH / FW
    ER = np.zeros(len(phi))
    for i in range(365, len(phi)-1):
        ER[i] = W[i] * (1 - ((phi[i] * np.tanh(1/phi[i])) * (1 - np.exp(-phi[i])))**0.5)

    # Calculate Discharge
    model_disch = np.zeros(len_R)
    for i in range(365, len_R-2):
        if i <= 730:
            for k in range(365, i + 1):  # k ranges from 365 to i (inclusive)
                model_disch[i + lag] += ER[k] * 0.4 / ((i - k) * 0.4 + 1) / ((i - k + 1) * 0.4 + 1)
        else:
            for k in range(i - 365, i + 1):  # k ranges from (i - 365) to i (inclusive)
                model_disch[i + lag] += ER[k] * 0.4 / ((i - k) * 0.4 + 1) / ((i - k + 1) * 0.4 + 1)

    model_disch[:366] = np.nan


    return model_disch 

############################################## NSE and B ################################################################
def calculate_nse(observed, model):
    # Exclude the first 366 values from both series
    observed_excl = observed[371:]
    model_excl = model[371:]
    
    # Calculate the mean of observed values
    mean_observed = np.mean(observed_excl)
    
    # Calculate the numerator and denominator of the NSE formula
    numerator = np.sum((observed_excl - model_excl) ** 2)
    denominator = np.sum((observed_excl - mean_observed) ** 2)
    
    # Calculate NSE
    nse = 1 - (numerator / denominator)
    
    return nse

# Define the function to calculate B
def B_metric(observed, model):
    # Calculate the mean of observed and modeled discharge
    mean_observed = np.mean(observed[371:])  # Exclude the first 366 values
    mean_model = np.mean(model[371:])        # Exclude the first 366 values
    
    # Calculate B using the given formula
    B = 1 - abs((mean_model - mean_observed) / (mean_model + mean_observed))
    
    return B
#############################################################################################################################

# loading data 
import os 
import glob
import pandas as pd
import matplotlib.pyplot as plt
os.chdir(r"E:\Prashant\DB sample\Python\DB") 
path = os.getcwd()


# Load the CSV file
data = pd.read_csv('Barman.csv')

# Extract R and PET columns
R = data['Rainfall (mm/day)'].values
PET = data['PET (mm/day)'].values
Obs=data['Observed Q (mm/day)']
# Set the lag value
lag = 1 # Adjust the value from 0 to 2 to get maximum NSE

# Run the function
model_disch = DB(R, PET, lag)

# Calculate NSE
nse = calculate_nse(Obs,model_disch)
print("NSE value is", nse)

# Calculate B
B = B_metric(Obs, model_disch)
print("B value is", B)

# Plotting the observed and modeled streamflow timeseries from 366 to 800
plt.figure(figsize=(10, 6))
plt.plot(range(366, 801), Obs[365:800], label='Observed Streamflow', color='blue')
plt.plot(range(366, 801), model_disch[365:800], label='Modeled Streamflow', color='red')
plt.xlabel('Time (days)')
plt.ylabel('Streamflow (mm/day)')
plt.title('Observed vs Modeled Streamflow')
plt.legend()
plt.grid(True)

# Save the plot as a JPEG file
plt.savefig('Observed_vs_Modeled_Streamflow.jpg', format='jpeg')

# Display the plot
plt.show()



