import matplotlib.pyplot as plt

# Create a 2x2 subplot
fig, axs = plt.subplots(2, 2)

# Create your individual plots or figures
plot1 = plt.plot([0, 1], [0, 1], label='Plot 1')
plot2 = plt.plot([0, 1], [1, 0], label='Plot 2')
plot3 = plt.plot([0, 1], [0.5, 0.5], label='Plot 3')
plot4 = plt.plot([0, 1], [0, 0], label='Plot 4')

# Assign each plot to a specific subplot position
axs[0, 0].plot(plot1)
axs[0, 1].plot(plot2)
axs[1, 0].plot(plot3)
axs[1, 1].plot(plot4)

# Optionally, you can add titles or labels to subplots
axs[0, 0].set_title('Subplot 1')
axs[0, 1].set_title('Subplot 2')
axs[1, 0].set_title('Subplot 3')
axs[1, 1].set_title('Subplot 4')

# Adjust layout to prevent overlap
plt.tight_layout()

# Show the plot
plt.show()