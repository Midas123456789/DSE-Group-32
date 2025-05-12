import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from airship import Airship

class visualise:

    def __init__(self, airship):
        self.airship = airship


        # At a certain FR, at a certain amount of lobes
        # FR = [2.5, 3, 4, 5]

    def altitude_graph(self):

        alt_step = 1000
        graph_df = [self.altitude, self.volume, self.wg, self.CL, self.CD]
        altitude = np.linspace(self.altitude + alt_step, 2 * self.altitude, alt_step)

        for altitudes in altitude:
            airship = Airship(self.FR, self.volume, self.n_lobes, self.velocity, altitudes, self.payload)
            airship.iterate_to_exact()
            row = [altitudes, airship.volume, airship.wg, airship.CL, airship.CD]

            graph_df = np.vstack((graph_df, row))

        graph_df = pd.DataFrame(graph_df, columns=['altitude', 'Volume', 'WeightG', 'CL', 'CD'])

        fig, axes = plt.subplots(2, 2, figsize=(8, 6))
        axes = axes.flatten()  # Flatten the 2x2 array of axes for easy indexing

        # Create subplots: 2 rows, 2 columns
        fig, axs = plt.subplots(2, 2, figsize=(10, 8))

        # Plot in each subplot
        axs[0, 0].plot(graph_df['Volume'], graph_df['altitude'])
        axs[0, 0].set_title('altitude vs Volume')

        axs[0, 1].plot(graph_df['WeightG'], graph_df['altitude'])
        axs[0, 1].set_title('altitude vs WeightG')

        axs[1, 0].plot(graph_df['CL'], graph_df['altitude'])
        axs[1, 0].set_title('altitude vs CL')

        axs[1, 1].plot(graph_df['CD'], graph_df['altitude'])
        axs[1, 1].set_title('altitude vs CD')

        # Adjust layout
        plt.tight_layout()
        plt.show()

        return