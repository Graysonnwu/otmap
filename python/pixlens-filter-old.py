import numpy as np
import math

def generate_pixel_lens_obj(
    lens_size=15.0,        # 透镜尺寸 (cm)
    divisions=100,         # 每边细分数量
    observer_distance=15.0,# 观察点到透镜的距离 (cm)
    spread_angle=60.0,     # 扩散角度 (度)
    output_file="pixel_lens_filter.obj"  # 输出文件名
):
    """生成像素透镜的OBJ模型文件"""
    
    # 每个小方格的尺寸
    pixel_size = lens_size / divisions
    
    # 折射率（假设为玻璃，约1.5）
    refractive_index = 1.5
    
    # 准备OBJ文件的数据结构
    vertices = []
    faces = []
    
    # 计算观察点位置（位于透镜中心后方）
    observer = np.array([0, 0, -observer_distance])
    
    # 扩散角的一半（弧度）
    half_spread = math.radians(spread_angle / 2)
    
    # 为每个像素生成几何体
    for i in range(divisions):
        for j in range(divisions):
            # 计算当前像素的中心位置
            x_center = (i + 0.5) * pixel_size - lens_size / 2
            y_center = (j + 0.5) * pixel_size - lens_size / 2
            z_center = 0  # 透镜在z=0平面上
            
            # 从观察点到像素中心的向量
            to_pixel = np.array([x_center, y_center, z_center]) - observer
            to_pixel = to_pixel / np.linalg.norm(to_pixel)
            
            # 计算理想的折射方向（这里简化计算，假设我们想要的扩散方向）
            # 距离透镜中心越远，偏转角度越大
            distance_from_center = math.sqrt(x_center**2 + y_center**2)
            max_distance = math.sqrt(2) * lens_size / 2
            
            # 根据距离中心的比例确定偏转角
            deviation_factor = distance_from_center / max_distance
            deviation_angle = deviation_factor * half_spread
            
            # 计算偏转后的方向
            # 为了简化，假设偏转方向是从中心向外
            direction_from_center = np.array([x_center, y_center, 0])
            if np.linalg.norm(direction_from_center) > 0:
                direction_from_center = direction_from_center / np.linalg.norm(direction_from_center)
            else:
                direction_from_center = np.array([1, 0, 0])  # 防止除零
                
            # 构建旋转矩阵以将z轴方向偏转指定角度
            cos_dev = math.cos(deviation_angle)
            sin_dev = math.sin(deviation_angle)
            
            # 使用Rodrigues' rotation formula计算偏转后的方向
            refracted_dir = np.array([0, 0, 1]) * cos_dev + sin_dev * direction_from_center
            refracted_dir = refracted_dir / np.linalg.norm(refracted_dir)
            
            # 使用折射定律（逆向）计算表面法线
            # 简化：直接使用入射方向与折射方向的中间向量作为法线方向
            normal = to_pixel + refracted_dir
            normal = normal / np.linalg.norm(normal)
            
            # 生成四个顶点位置（平面小方格）
            # 透镜整体保持平面，但各小方格之间可以不连续
            v_base = len(vertices) + 1  # OBJ索引从1开始
            
            # 计算每个顶点的z偏移，使平面法线正确
            # 计算一个参考平面，它的法线是我们想要的法线
            plane_d = -np.dot(normal, np.array([x_center, y_center, z_center]))
            
            # 计算四个顶点
            corners = [
                [x_center - pixel_size/2, y_center - pixel_size/2],  # 左下
                [x_center + pixel_size/2, y_center - pixel_size/2],  # 右下
                [x_center + pixel_size/2, y_center + pixel_size/2],  # 右上
                [x_center - pixel_size/2, y_center + pixel_size/2],  # 左上
            ]
            
            for cx, cy in corners:
                # 根据平面方程计算z坐标
                cz = (-plane_d - normal[0]*cx - normal[1]*cy) / normal[2]
                vertices.append([cx, cy, cz])
            
            # 添加面（三角形，分成两个）
            faces.append([v_base, v_base+1, v_base+2])
            faces.append([v_base, v_base+2, v_base+3])
    
    # 写入OBJ文件
    with open(output_file, 'w') as f:
        f.write("# Pixel Lens OBJ Model\n")
        f.write(f"# Size: {lens_size}cm x {lens_size}cm\n")
        f.write(f"# Divisions: {divisions}x{divisions}\n")
        f.write(f"# Observer Distance: {observer_distance}cm\n")
        f.write(f"# Spread Angle: {spread_angle} degrees\n\n")
        
        # 写入顶点
        for v in vertices:
            f.write(f"v {v[0]:.6f} {v[1]:.6f} {v[2]:.6f}\n")
        
        # 写入面
        for face in faces:
            f.write(f"f {face[0]} {face[1]} {face[2]}\n")
    
    print(f"生成完成! 模型已保存为 {output_file}")
    print(f"总计顶点数: {len(vertices)}")
    print(f"总计面数: {len(faces)}")

# 生成模型
if __name__ == "__main__":
    generate_pixel_lens_obj()
