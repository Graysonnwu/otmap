import numpy as np
import matplotlib.pyplot as plt

# 读取数据文件 x2.dat
data = []
with open('iq-r.dat', 'r') as file:
    for line in file:
        x1, y1, x2, y2, _ = map(float, line.split())  # 读取 3D 坐标
        source = (-y1, -x1, 0)  # source 的 z 坐标为 0
        target = (-y2, -x2, 3)  # target 的 z 坐标为 z2
        data.append((source, target))

# 设置折射率
n = 1.49

# 创建网格
grid_size = 127
x_range = (-0.5, 0.5)
y_range = (-0.5, 0.5)
x_vals = np.linspace(x_range[0], x_range[1], grid_size)
y_vals = np.linspace(y_range[0], y_range[1], grid_size)
curl_field = np.zeros((grid_size, grid_size))

# 计算每条光线的法线并填充旋度场
for (source, target) in data:
    x1, y1, z1 = source
    x2, y2, z2 = target

    # 计算出射光线方向
    transmitted = np.array([x2 - x1, y2 - y1, z2 - z1])
    transmitted = transmitted / np.linalg.norm(transmitted)  # 归一化

    # 入射光线方向
    incidentLight = np.array([0, 0, -1])

    # 计算法线方向
    normal = transmitted + incidentLight * n
    normal = normal / np.linalg.norm(normal)  # 归一化法线

    # 获取源点的 (x1, y1) 映射到网格上的位置
    ix = np.digitize(x1, x_vals) - 1
    iy = np.digitize(y1, y_vals) - 1

    # 累加法线的 z 组件（假设旋度基于法线方向的变化）
    if 0 <= ix < grid_size and 0 <= iy < grid_size:
        curl_field[iy, ix] += normal[1] - normal[0]

# 可视化旋度场
plt.figure(figsize=(8, 6))
plt.imshow(curl_field, extent=(x_range[0], x_range[1], y_range[0], y_range[1]),
           origin='lower', cmap='RdBu_r', interpolation='nearest')
plt.colorbar(label='Curl Value')
plt.title('Curl Field Visualization')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()
