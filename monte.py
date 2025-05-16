import numpy as np
import matplotlib.pyplot as plt

#varying_weights = False
#varying_scores = True
tittles = {(True, True): "Weights and Scores",
           (True, False): "Weights",
          (False, True): "Scores",
          (False, False): "Base Case"} 
#AHP=False
folder = 'figs_for_MTR'
def sensitivity_analysis(varying_weights, varying_scores, AHP, savefigs=False):
    """
    This function performs a sensitivity analysis on the design selection process.
    It uses Monte Carlo simulation to evaluate the impact of uncertainty in weights and scores.
    """
    # --- 0. Set Random Seed for Reproducibility ---
    np.random.seed(42)  # Set a random seed for reproducibility    

    # --- 1. Define Base Table Data ---
    criteria_names = ["Cost", "Flight Performance", "Complexity", "Availability", "Payload mass", "Power"]
    design_names = ["Tandem Aircraft", "Conventional Aircraft", "Flying Wing", "Hybrid", "Conventional Airship"]

    # Base weights (must sum to 1): normalize given weights

    #standerd weighted
    raw_weights = np.array([20, 15, 15, 20, 10, 15])

    if AHP:
        #AHP weighted
        raw_weights = np.array([0.35,0.2,0.1,0.2,0.15])

    base_weights = raw_weights / np.sum(raw_weights)

    # Base scores (Designs x Criteria)
    base_scores = np.array([
        [3, 3, 2, 4, 2, 3],  # Tandem Aircraft
        [3, 4, 2, 3, 2, 2],  # Conventional A
        [2, 3, 1, 2, 2, 2],  # Flying Wing
        [2, 2, 1, 2, 4, 1],  # Hybrid
        [3, 3, 2, 2, 4, 2]   # Conventional B
    ])
    if AHP:
        base_scores = np.array([
        [3, 3, 2, 4, 3],  # Tandem Aircraft
        [3, 4, 2, 3, 2],  # Conventional A
        [2, 3, 1, 2, 2],  # Flying Wing
        [2, 2, 1, 2, 1],  # Hybrid
        [3, 3, 2, 2, 2]   # Conventional B
        ])
        criteria_names = ["Cost", "Flight Performance", "Complexity", "Availability", "Power"]


    # --- 2. Define Uncertainty Parameters for Monte Carlo ---
    N_ITERATIONS = 10000
    MIN_SCORE = 1
    MAX_SCORE = 4

    # For weights: min, mode, max for triangular distribution
    # We'll define a relative range for weights, e.g., +/- 20% of the mode
    WEIGHT_UNCERTAINTY_FACTOR = 0.3 # e.g., 0.2 means +/- 20%
    weight_params = []
    for w in base_weights:
        min_w = max(0.01, w * (1 - WEIGHT_UNCERTAINTY_FACTOR)) # Ensure weight doesn't go to zero or negative
        max_w = min(1.0, w * (1 + WEIGHT_UNCERTAINTY_FACTOR))  # Ensure weight doesn't exceed 1
        mode_w = w
        weight_params.append((min_w, mode_w, max_w))

    # For scores: min, mode, max for triangular distribution
    # We'll assume scores can vary by +/- 1 point from their base value (within 1-10)
    score_params = [] # List of lists of tuples: designs x criteria x (min, mode, max)
    for design_scores in base_scores:
        design_score_params = []
        for s in design_scores:
            min_s = max(MIN_SCORE, s - 1)
            max_s = min(MAX_SCORE, s + 1)
            mode_s = s
            design_score_params.append((min_s, mode_s, max_s))
        score_params.append(design_score_params)

    # --- 3. Run Monte Carlo Simulation ---
    simulation_total_scores = np.zeros((N_ITERATIONS, len(design_names)))
    design_wins = np.zeros(len(design_names))

    for i in range(N_ITERATIONS):

        # Sample new weights
        sampled_weights_raw = np.array([
            np.random.triangular(low, mode, high) for low, mode, high in weight_params
        ])
        # Normalize weights to sum to 1
        sampled_weights = sampled_weights_raw / np.sum(sampled_weights_raw)
        if not varying_weights:
            sampled_weights = base_weights
        # Sample new scores
        sampled_scores = np.zeros_like(base_scores, dtype=float)
        for d_idx in range(len(design_names)):
            for c_idx in range(len(criteria_names)):
                min_s, mode_s, max_s = score_params[d_idx][c_idx]
                sampled_scores[d_idx, c_idx] = np.random.triangular(min_s, mode_s, max_s)
                # Clamp scores just in case triangular goes slightly out due to float precision
                sampled_scores[d_idx, c_idx] = np.clip(sampled_scores[d_idx, c_idx], MIN_SCORE, MAX_SCORE)

        if not varying_scores:
            sampled_scores = base_scores
        # Calculate total scores for each design in this iteration
        current_iteration_total_scores = np.dot(sampled_scores, sampled_weights) # (Designs x Criteria) . (Criteria) = (Designs)
        simulation_total_scores[i, :] = current_iteration_total_scores

        # Determine winning design for this iteration
        winner_idx = np.argmax(current_iteration_total_scores)
        design_wins[winner_idx] += 1

    # --- 4. Analyze and Display Results ---

    # Calculate base case total scores for comparison
    base_total_scores = np.dot(base_scores, base_weights)
    print("--- Base Case Results ---")
    for i, name in enumerate(design_names):
        print(f"{name}: {base_total_scores[i]:.2f}")
    print(f"Base Case Winner: {design_names[np.argmax(base_total_scores)]}\n")


    print("--- Monte Carlo Simulation Results ---")
    print(f"Number of Iterations: {N_ITERATIONS}\n")

    print("Probability of being the best design:")
    for i, name in enumerate(design_names):
        probability = (design_wins[i] / N_ITERATIONS) * 100
        print(f"{name}: {probability:.2f}%")

    print("\nDistribution of Total Scores (Mean +/- Std Dev):")
    for i, name in enumerate(design_names):
        mean_score = np.mean(simulation_total_scores[:, i])
        std_score = np.std(simulation_total_scores[:, i])
        print(f"{name}: {mean_score:.2f} +/- {std_score:.2f}")

    # Plotting the results
    # 1. Histograms of total scores for each design
    plt.figure(figsize=(12, 8))
    for i, name in enumerate(design_names):
        plt.hist(simulation_total_scores[:, i], bins=50, alpha=0.7, label=f'{name} (Mean: {np.mean(simulation_total_scores[:, i]):.2f})')
    plt.title('Distribution of Total Scores from Monte Carlo Simulation')
    plt.xlabel('Total Score')
    plt.ylabel('Frequency')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    if not savefigs:
        plt.show()
    else:
        plt.savefig(f'{folder}/{"AHP" if AHP else "Standard"} Sensitivity to {tittles[(varying_weights,varying_scores)]}1.png', dpi=600)
        plt.close()

    # 2. Bar chart of win probabilities
    plt.figure(figsize=(10, 6))
    win_probabilities = (design_wins / N_ITERATIONS) * 100
    plt.bar(design_names, win_probabilities, color=['skyblue', 'lightcoral', 'lightgreen', 'gold', 'lightsalmon'])
    plt.title('Probability of Each Design Being Optimal')
    plt.xlabel('Design')
    plt.ylabel('Probability (%)')
    for i, prob in enumerate(win_probabilities):
        plt.text(i, prob + 0.5, f'{prob:.1f}%', ha='center', va='bottom')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    if not savefigs:
        plt.show()
    else:
        plt.savefig(f'{folder}/{"AHP" if AHP else "Standard"} Sensitivity to {tittles[(varying_weights,varying_scores)]}2.png', dpi=600)
        plt.close()

if __name__ == "__main__":
    #global varying_weights, varying_scores, AHP
    combs = [(True, True, False), (True, False, False), (False, True, False), (False, False, False), (True, True, True), (True, False, True), (False, True, True), (False, False, True)]
    idk = [sensitivity_analysis(*comb,savefigs=True) for comb in combs]
    '''
    varying_weights, varying_scores = True, True
    AHP = False
    sensitivity_analysis(savefigs=True)
    AHP = True

    varying_weights, varying_scores = True, False
    sensitivity_analysis(savefigs=True)
    varying_weights, varying_scores = False, True
    sensitivity_analysis(savefigs=True)
    varying_weights, varying_scores = False, False
    sensitivity_analysis(savefigs=True)
    '''
