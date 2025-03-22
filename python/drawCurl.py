import numpy as np
import matplotlib.pyplot as plt

def calculate_curl(normal_field, res):
    """
    计算法线场的离散旋度
    normal_field: shape为(res*res, 3)的数组，每行包含(nx, ny, nz)
    """
    # 重塑数组为resxresx3的形状
    field = normal_field.reshape(res, res, 3)
    dx = 1.0/res
    dy = 1.0/res
    
    # 计算x和y方向的偏导数
    # 使用中心差分
    dnx_dy = np.zeros_like(field[:,:,0])
    dny_dx = np.zeros_like(field[:,:,1])
    
    # y方向的偏导数
    dnx_dy[1:-1, :] = (field[2:, :, 0] - field[:-2, :, 0]) / (2*dy)
    # x方向的偏导数
    dny_dx[:, 1:-1] = (field[:, 2:, 1] - field[:, :-2, 1]) / (2*dx)
    
    # 旋度的z分量 (我们主要关注z分量，因为这反映了xy平面内的旋转)
    curl_z = dny_dx - dnx_dy
    
    return curl_z

def visualize_curl(normal_field, res):
    """
    可视化法线场的旋度
    """
    curl_z = calculate_curl(normal_field, res)
    
    # 创建坐标网格
    x = np.linspace(-0.5, 0.5, res)
    y = np.linspace(-0.5, 0.5, res)
    X, Y = np.meshgrid(x, y)
    
    # 绘制热力图
    plt.figure(figsize=(10, 8))
    plt.pcolormesh(-X, -Y, curl_z, cmap='RdBu', shading='auto', vmin=-0.2, vmax=0.2)
    plt.colorbar(label='Curl magnitude')
    plt.title('Curl of Normal Field')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.axis('equal')
    
    # # 标记旋度较大的区域
    # threshold = np.std(curl_z) * 2  # 使用2倍标准差作为阈值
    # high_curl = np.abs(curl_z) > threshold
    
    # # 在旋度较大的区域上画点
    # plt.plot(X[high_curl], Y[high_curl], 'k.', markersize=1, alpha=0.5)
    
    plt.show()


res = 1200
data = np.loadtxt('pixlens-test.dat')
normal_field = np.zeros((res, res, 3))

for i, dot in enumerate(data):
    x1, y1, z1, x2, y2, z2 = dot[1], dot[0], 0, dot[3], dot[2], 3
    normal = np.array([x2-x1, y2-y1, z2-z1])
    normal = normal / np.linalg.norm(normal)
    normal_field[int(i%res), int(i/res)] = normal

visualize_curl(normal_field, res)

