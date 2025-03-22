import numpy as np

# 设置参数
SOURCE_GRID_SIZE = 800  # source点的细分度
TARGET_GRID_SIZE = 40   # target点区域的细分度
SCALE_FACTOR = 100      # target点的放大倍数（小于1为缩小，大于1为放大）
PATTERN = 'square'   # 密铺模式：'square' 或 'hex'
FOCAL_LENGTH = 0.6     # 焦距
OBJECT_DISTANCE = 60   # 物距

def generate_square_targets():
    """生成正方形密铺的target点"""
    target_x = np.linspace(1/(2*TARGET_GRID_SIZE), 1-1/(2*TARGET_GRID_SIZE), TARGET_GRID_SIZE)
    target_y = np.linspace(1/(2*TARGET_GRID_SIZE), 1-1/(2*TARGET_GRID_SIZE), TARGET_GRID_SIZE)
    target_x, target_y = np.meshgrid(target_x, target_y)
    return np.column_stack((target_x.flatten(), target_y.flatten()))

def generate_hexagonal_targets():
    """生成正六边形密铺的target点"""
    # 计算六边形网格的基本参数
    h = 1.0 / TARGET_GRID_SIZE  # 基本水平间距
    v = h * np.sqrt(3) / 2      # 垂直间距
    
    # 计算需要的行数和列数以确保完全覆盖
    rows = int(1.0 / v) + 2
    cols = int(1.0 / h) + 2
    
    points = []
    for i in range(rows):
        offset = h/2 if i % 2 else 0  # 奇数行偏移半个单位
        for j in range(cols):
            x = offset + j * h
            y = i * v
            # 只保留[0,1]区间内的点
            if -0.1 <= x <= 1.1 and -0.1 <= y <= 1.1:
                points.append([x, y])
    
    return np.array(points)

def generate_points():
    # 生成source点的坐标网格
    x = np.linspace(0, 1, SOURCE_GRID_SIZE)
    y = np.linspace(0, 1, SOURCE_GRID_SIZE)
    source_x, source_y = np.meshgrid(x, y)
    source_points = np.column_stack((source_x.flatten(), source_y.flatten()))

    # 根据选择的模式生成target点
    if PATTERN == 'square':
        target_points = generate_square_targets()
    else:  # hex
        target_points = generate_hexagonal_targets()
    
    # 为每个source点分配对应的target点
    # 使用最近邻算法找到最近的target点（基于原始位置）
    from scipy.spatial import cKDTree
    tree = cKDTree(target_points)
    _, target_indices = tree.query(source_points)
    assigned_targets = target_points[target_indices]
    
    # 计算z坐标(在缩放之前)
    # z_coords = -np.sqrt((assigned_targets[:, 0] - 0.5)**2 + (assigned_targets[:, 1] - 0.5)**2 + FOCAL_LENGTH**2)
    z_coords = np.full((source_points.shape[0],), -OBJECT_DISTANCE)
    
    # 对assigned_targets点进行缩放变换
    assigned_targets = assigned_targets - 0.5
    assigned_targets = assigned_targets * SCALE_FACTOR
    assigned_targets = assigned_targets + 0.5

    # 将结果保存到文件，包括source点、target点(x,y)和z坐标
    result = np.column_stack((source_points, assigned_targets, z_coords))
    np.savetxt('pixlens-test.dat', result, fmt='%.6f')

if __name__ == '__main__':
    generate_points()
    print(f"已生成点集文件 pixlens-test.dat")
    print(f"Source点数量: {SOURCE_GRID_SIZE * SOURCE_GRID_SIZE}")
    print(f"Target点模式: {'正方形' if PATTERN == 'square' else '正六边形'}")
    print(f"Target区域数量: {TARGET_GRID_SIZE * TARGET_GRID_SIZE if PATTERN == 'square' else int(TARGET_GRID_SIZE * TARGET_GRID_SIZE * 0.866)}")
    print(f"焦距: {FOCAL_LENGTH}， 物距: {OBJECT_DISTANCE}") 