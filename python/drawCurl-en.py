import numpy as np
import matplotlib.pyplot as plt

def calculate_curl(normal_field, res):
    """
    Calculate the discrete curl of the normal field.
    normal_field: An array with shape (res*res, 3), where each row contains (nx, ny, nz).
    """
    # Reshape the array to shape (res, res, 3)
    field = normal_field.reshape(res, res, 3)
    dx = 1.0 / res
    dy = 1.0 / res
    
    # Compute the partial derivatives in the x and y directions
    # Using central difference
    dnx_dy = np.zeros_like(field[:,:,0])
    dny_dx = np.zeros_like(field[:,:,1])
    
    # Partial derivative in the y direction
    dnx_dy[1:-1, :] = (field[2:, :, 0] - field[:-2, :, 0]) / (2 * dy)
    # Partial derivative in the x direction
    dny_dx[:, 1:-1] = (field[:, 2:, 1] - field[:, :-2, 1]) / (2 * dx)
    
    # The z-component of the curl
    # (we mainly focus on the z-component as it reflects rotation in the xy-plane)
    curl_z = dny_dx - dnx_dy
    
    return curl_z

def visualize_curl(normal_field, res):
    """
    Visualize the curl of the normal field.
    """
    curl_z = calculate_curl(normal_field, res)
    
    # Create the coordinate grid
    x = np.linspace(-0.5, 0.5, res)
    y = np.linspace(-0.5, 0.5, res)
    X, Y = np.meshgrid(x, y)
    
    # Plot the heatmap
    plt.figure(figsize=(10, 8))
    plt.pcolormesh(-X, -Y, curl_z, cmap='RdBu', shading='auto', vmin=-0.2, vmax=0.2)
    plt.colorbar(label='Curl magnitude')
    plt.title('Curl of Normal Field')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.axis('equal')
    
    # # Mark regions with large curl
    # threshold = np.std(curl_z) * 2  # Use 2 times the standard deviation as the threshold
    # high_curl = np.abs(curl_z) > threshold
    
    # # Plot points on regions with high curl
    # plt.plot(X[high_curl], Y[high_curl], 'k.', markersize=1, alpha=0.5)
    
    plt.show()

def prepare_normal_field_otmap(file, res):
    data = np.loadtxt(file)
    normal_field = np.zeros((res, res, 3))
    for i, dot in enumerate(data):
        x1, y1, z1, x2, y2, z2 = dot[1], dot[0], 0, dot[3], dot[2], 3
        normal = np.array([x2 - x1, y2 - y1, z2 - z1])
        normal = normal / np.linalg.norm(normal)
        normal_field[int(i % res), int(i / res)] = normal
    return normal_field

res = 500
normal_field = prepare_normal_field_otmap('x2.dat', res)

visualize_curl(normal_field, res)

