import numpy as np
import matplotlib.pyplot as plt

# Define focal lengths
focal_lengths = [10, 15, 20]
colors = ['blue', 'red', 'green']

# Create figure and axis
plt.figure(figsize=(10, 6))

# For each focal length
for i, f in enumerate(focal_lengths):
    # Create array of object distances (p)
    # Start slightly above focal length to avoid division by zero
    p_values = np.concatenate([
        np.linspace(f*1.01, 50, 100),  # Dense sampling near focal length
        np.linspace(50, 200, 50),      # Medium sampling
        np.linspace(200, 1000, 20)     # Sparse sampling for large distances
    ])
    
    # Calculate image distances (q) using lens equation: 1/p + 1/q = 1/f
    q_values = f * p_values / (p_values - f)
    
    # Calculate total distances (object-to-lens + lens-to-image)
    total_distances = p_values + q_values
    
    # Plot the curve
    plt.plot(total_distances, q_values, label=f'f = {f} cm', color=colors[i], linewidth=2)

# Add vertical asymptotes at 2f for each focal length
for i, f in enumerate(focal_lengths):
    # plt.axvline(x=4*f, color=colors[i], linestyle='--', alpha=0.5)
    plt.axline((f,0), slope=1, color=colors[i], linestyle='--', alpha=0.5)
    plt.axhline(y=f, color=colors[i], linestyle='--', alpha=0.5)

plt.axline((0,0), slope=1/2, color='grey', linestyle='--', alpha=0.5)

plt.rcParams['font.sans-serif']=['STHeiti'] # 用来正常显示中文标签

# Add labels and legend
plt.xlabel(u'物体到相机的距离(u+v) [cm]', fontsize=12)
plt.ylabel(u'透镜到相机的距离(v) [cm]', fontsize=12)
plt.title(u'对于不同距离的物体，透镜到相机的理想距离', fontsize=14)
plt.legend(fontsize=10)
plt.grid(True, alpha=0.3)

# Set appropriate axis limits
plt.xlim(20, 120)  # Limit x-axis for better visualization
plt.ylim(0, 60)  # Limit y-axis for better visualization
plt.xticks(np.arange(20, 125, 5))
plt.yticks(np.arange(0, 65, 5))

plt.gca().set_aspect('equal', adjustable='box')
# # Add annotations
# for i, f in enumerate(focal_lengths):
#     plt.annotate(f'2f = {2*f} cm', 
#                  xy=(2*f, 5), 
#                  xytext=(2*f+5, 10),
#                  arrowprops=dict(arrowstyle='->'),
#                  color=colors[i])

plt.tight_layout()
plt.show()
