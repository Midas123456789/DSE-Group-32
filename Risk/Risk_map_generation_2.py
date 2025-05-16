import pandas as pd

import matplotlib.pyplot as plt
import numpy as np

# Define the path or URL to the Excel file
excel_file_path = "Risks.xlsx"  # Replace with the actual path or URL

# Read the Excel file

sheets = ["Technical risks", "cost risks", "Scheduling risks", "Programmatic risks"]
data = {sheet: pd.read_excel(excel_file_path, sheet_name=sheet) for sheet in sheets}

# Convert each sheet into a dictionary
technical_risks = data["Technical risks"].dropna().to_dict(orient="records")
cost_risks = data["cost risks"].dropna().to_dict(orient="records")
scheduling_risks = data["Scheduling risks"].dropna().to_dict(orient="records")
programmatic_risks = data["Programmatic risks"].dropna().to_dict(orient="records")

'''
# Example: Print the dictionaries (optional)
print("Technical Risks:", technical_risks)
print("Cost Risks:", cost_risks)
print("Scheduling Risks:", scheduling_risks)
print("Programmatic Risks:", programmatic_risks)
'''
# Risk map generation (example: scatter plot for risks)
def generate_risk_map(risks, title):
    boundaries = [0.2,0.5]
    # Assuming each risk has 'Likelihood' and 'Impact' columns
    
    #find all risks have the same likelihood and impact and add them to the same point
    data = {}
    for risk in risks:
    
        key = (risk['Likelihood'], risk['Impact'])
        #check if the key is numeric
        
        if key not in data:
            data[key] = f'{int(risk["ID"])}'
        else:
            data[key] += f', {int(risk["ID"])}'


    print(np.array(list(data.keys())))
    likelihood = np.array(list(data.keys()))[:, 0]
    impact = np.array(list(data.keys()))[:,1]
    labels = data.values()

    plt.figure(figsize=(10, 6))
    plt.scatter(likelihood, impact, c='red', alpha=0.6, edgecolors='w', s=100)
    for i, label in enumerate(labels):
        plt.text(likelihood[i], impact[i]-0.05, label, fontsize=9, ha='right')

    plt.title(title)
    plt.xlabel("Likelihood")
    plt.ylabel("Impact")
    plt.xlim(0, 1)  # Adjust based on your scale
    plt.ylim(0, 1)  # Adjust based on your scale
    #Add colored background for risk levels
    x=np.linspace(0, 1, 100)
    y=boundaries[1]/x
    y1=boundaries[0]/x
    plt.fill_between(x, y, 1, color='red', alpha=0.2, label='High Risk')
    plt.fill_between(x, y1, y, color='yellow', alpha=0.2, label='Medium Risk')
    plt.fill_between(x, 0, y1, color='green', alpha=0.2, label='Low Risk')

    #make the figure square
    plt.gca().set_aspect('equal', adjustable='box')
    
    #plt.legend(loc='upper right')
    #instead of a legend write the risk levels in the figure
    plt.text(0.8, 0.8, 'High Risk', fontsize=12, color='red')
    plt.text(0.5, 0.5, 'Medium Risk', fontsize=12, color='orange')
    plt.text(0.2, 0.2, 'Low Risk', fontsize=12, color='green')
    plt.grid(True)
    plt.show()

# Generate risk maps for each category
generate_risk_map(technical_risks, "Technical Risks Map")
generate_risk_map(cost_risks, "Cost Risks Map")
generate_risk_map(scheduling_risks, "Scheduling Risks Map")
generate_risk_map(programmatic_risks, "Programmatic Risks Map")