import pandas as pd

import matplotlib.pyplot as plt
import numpy as np

# Define the path or URL to the Excel file
excel_file_path = "Risks.xlsx"  # Replace with the actual path or URL

# Read the Excel file

sheets = ["RiskRevolution","Final Risks","Riskmap","Technical risks", "Riskmapmitigation","cost risks", "Scheduling risks", "Programmatic risks", "Mitigated Cost"]
data = {sheet: pd.read_excel(excel_file_path, sheet_name=sheet) for sheet in sheets}

# Convert each sheet into a dictionary
technical_risks = data["Riskmap"].to_dict(orient="records")
mitigated_Trisks = data["Riskmapmitigation"].to_dict(orient="records")
scheduling_risks = data["Scheduling risks"].to_dict(orient="records")
programmatic_risks = data["Programmatic risks"].to_dict(orient="records")
cost_risks = data["cost risks"].to_dict(orient="records")
mitigated_Crisks = data["Mitigated Cost"].to_dict(orient="records")

'''
# Example: Print the dictionaries (optional)

print("Cost Risks:", cost_risks)
print("Scheduling Risks:", scheduling_risks)
print("Programmatic Risks:", programmatic_risks)
'''
# Risk map generation (example: scatter plot for risks)
def generate_risk_map(risks, title):
    boundaries = [2, 10]  # Adjusted thresholds for 0–5 scale

    data = {}
    for risk in risks:
        if not (risk['Likelihood'] == 'TBD' or risk['Impact'] == 'TBD'):
            try:
                likelihood = float(risk['Likelihood'])
                impact = float(risk['Impact'])
                key = (likelihood, impact)

                if key not in data:
                    data[key] = [risk["ID"]]
                else:
                    data[key] += [risk["ID"]]
            except (ValueError, TypeError, KeyError):
                print(f"Skipping malformed risk: {risk}")
                continue

    likelihood = np.array(list(data.keys()))[:, 0]
    impact = np.array(list(data.keys()))[:, 1]
    labels = list(data.values())

    plt.figure(figsize=(7, 7))
    plt.scatter(likelihood, impact, c='red', alpha=0.6, edgecolors='w', s=100)

    for i, label in enumerate(labels):
        nl = '\n'
        if len(label) > 3:
            # 2 columns: group labels in pairs, separated by comma and newline after every 2nd label
            wrapped_label = ''
            for j in range(len(label)):
                wrapped_label += label[j]
                if (j + 1) % 2 == 0 and (j + 1) != len(label):
                    wrapped_label += '\n'  # newline every 2 labels except last
                elif (j + 1) != len(label):
                    wrapped_label += ', '  # comma between labels in a pair
        else:
            # 1 column: each label on its own line
            wrapped_label = '\n'.join(label)

        # Slight offset based on index to reduce overlap
        x_offset = 0.05 # cycles through -0.05, 0, 0.05
        y_offset = 0.05  # cycles through -0.15 to 0.15

        adjusted_x = np.clip(likelihood[i] + x_offset, 0, 5)
        adjusted_y = np.clip(impact[i] + y_offset, 0, 5)

        plt.annotate(wrapped_label, (adjusted_x, adjusted_y), fontsize=15, ha='center', rotation=-10)

    plt.xlabel("Likelihood")
    plt.ylabel("Impact")
    plt.xlim(0, 5)
    plt.ylim(0, 5)

    # Draw risk zones
    x = np.linspace(0.01, 5, 500)
    y_high = boundaries[1] / x
    y_med = boundaries[0] / x

    plt.fill_between(x, y_high, 5, color='red', alpha=0.2, label='High Risk')
    plt.fill_between(x, y_med, y_high, color='yellow', alpha=0.2, label='Medium Risk')
    plt.fill_between(x, 0, y_med, color='green', alpha=0.2, label='Low Risk')

    # Risk level labels
    plt.text(3.8, 4.05, 'High Risk', fontsize=18, color='red')
    plt.text(1.5 , 2.25, 'Medium Risk', fontsize=18, color='orange')
    plt.text(0.1, 0.2, 'Low Risk', fontsize=18, color='green')

    plt.gca().set_aspect('equal', adjustable='box')
    plt.grid(True)
    plt.savefig(f"{title}.png", dpi=600, bbox_inches='tight')
    # plt.show()

# Generate risk maps for each category
#generate_risk_map(technical_risks, "Technical riskmap")
#generate_risk_map(mitigated_Trisks, "Mitigated Technical riskmap")
generate_risk_map(scheduling_risks, "scheduling riskmap")
generate_risk_map(programmatic_risks, "programmatic riskmap")
generate_risk_map(cost_risks, "cost riskmap")
generate_risk_map(mitigated_Crisks, "Mitigated Cost")