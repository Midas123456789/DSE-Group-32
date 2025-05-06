import pandas as pd

import matplotlib.pyplot as plt
import numpy as np

# Define the path or URL to the Excel file
excel_file_path = "Risks.xlsx"  # Replace with the actual path or URL

# Read the Excel file

sheets = ["Technical risks", "cost risks", "Scheduling risks", "Programmatic risks"]
data = {sheet: pd.read_excel(excel_file_path, sheet_name=sheet) for sheet in sheets}

# Convert each sheet into a dictionary
technical_risks = data["Technical risks"].to_dict(orient="records")
cost_risks = data["cost risks"].to_dict(orient="records")
scheduling_risks = data["Scheduling risks"].to_dict(orient="records")
programmatic_risks = data["Programmatic risks"].to_dict(orient="records")

'''
# Example: Print the dictionaries (optional)
print("Technical Risks:", technical_risks)
print("Cost Risks:", cost_risks)
print("Scheduling Risks:", scheduling_risks)
print("Programmatic Risks:", programmatic_risks)
'''
# Risk map generation (example: scatter plot for risks)
def generate_risk_map(risks, title):
    boundaries = [0.1,0.3]
    # Assuming each risk has 'Likelihood' and 'Impact' columns
    
    #find all risks have the same likelihood and impact and add them to the same point
    data = {}
    for risk in risks:
        if not (risk['Likelihood'] == 'TBD' or risk['Impact'] == 'TBD'):
            key = (risk['Likelihood'], risk['Impact'])

        
            if key not in data:
                data[key] = [int(risk["ID"])]
            else:
                data[key] += [int(risk["ID"])]


    likelihood = np.array(list(data.keys()))[:, 0]
    impact = np.array(list(data.keys()))[:,1]
    labels = data.values()

    plt.figure(figsize=(6, 6))
    plt.scatter(likelihood, impact, c='red', alpha=0.6, edgecolors='w', s=100)
    #add labels to the points
    #make the labels smaller and wrap them if they are too long
    #make sure they do not overlap with the points and dont go outside the figure

    for i, label in enumerate(labels):
        #plt.text(likelihood[i], impact[i], label, fontsize=9, ha='right',wrap=True, bbox=dict(facecolor='none', alpha=0., edgecolor='none', boxstyle='round,pad=0.3', lw=0.5))
        # Adjust label positions to ensure they stay within the figure bounds
        x_offset = -0.01 if likelihood[i] > 0.5 else 0.01
        x_offset=0
        y_offset = -0.05 if impact[i] > 0.5 else 0.05
        #y_offset=0
        adjusted_x = max(0.01, min(0.99, likelihood[i] + x_offset))  # Ensure labels stay within bounds
        adjusted_y = max(0.05, min(0.95, impact[i] + y_offset))  # Adjust y bounds similarly
        
        # Wrap labels within a box of width 0.1 in x
        #wrapped_label = "\n".join([label[j:j+10] for j in range(0, len(label), 10)])  # Approximate wrapping
        #wrap the label 3 entries per line
        #print(f'label: {label}')
        nl='\n'
        wrapped_label = ''.join([(f"{label[j]}{nl if (j+1)%3==0 else ', '}")  for j in range(len(label))])
        
        plt.annotate(wrapped_label, (adjusted_x, adjusted_y), fontsize=9, ha='left', wrap=True,rotation=-45)
        x_offset = -0.02 if likelihood[i] > 0.5 else 0.02
        y_offset = -0.02 if impact[i] > 0.5 else 0.02
        adjusted_x = max(0.05, min(0.95, likelihood[i] + x_offset))  # Ensure labels stay within bounds
        adjusted_y = max(0.05, min(0.95, impact[i] + y_offset))  # Adjust y bounds similarly
        
        #plt.annotate(label, (adjusted_x, adjusted_y), fontsize=9, ha='center', wrap=True)

    #plt.title(title)
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
    #save the figure
    plt.savefig(f"Figure for baseline/{title}.png", dpi=600, bbox_inches='tight')
    #plt.show()

# Generate risk maps for each category
generate_risk_map(technical_risks, "Technical riskmap")
generate_risk_map(cost_risks, "cost riskmap")
generate_risk_map(scheduling_risks, "scheduling riskmap")
generate_risk_map(programmatic_risks, "programmatic riskmap")