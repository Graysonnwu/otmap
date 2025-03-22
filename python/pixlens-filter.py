import numpy as np
import math

def generate_pixel_lens_obj(
    lens_size=15.0,        # 透镜尺寸 (cm)
    divisions=30,         # 每边细分数量
    observer_distance=30.0,# 观察点到透镜的距离 (cm)
    spread_angle=100.0,     # 扩散角度 (度)
    thickness=1.0,         # 透镜厚度 (cm)
    refractive_index=1.49, # 折射率 (亚克力)
    output_file="pixel_lens_filter.obj"  # 输出文件名
):
    """生成像素透镜的OBJ实体模型文件"""
    
    # 每个小方格的尺寸
    pixel_size = lens_size / divisions
    
    # 准备OBJ文件的数据结构
    vertices = []
    faces = []
    
    # 计算观察点位置（位于透镜中心后方）
    observer = np.array([0, 0, -observer_distance])
    
    # 扩散角的一半（弧度）
    half_spread = math.radians(spread_angle / 2)
    
    # 给每个像素分配唯一的顶点（不共享顶点，避免边界问题）
    for i in range(divisions):
        for j in range(divisions):
            # 计算当前像素的中心位置
            x_center = (i + 0.5) * pixel_size - lens_size / 2
            y_center = (j + 0.5) * pixel_size - lens_size / 2
            z_center = 0  # 透镜在z=0平面上
            
            # 从观察点到像素中心的向量（入射方向）
            to_pixel = np.array([x_center, y_center, z_center]) - observer
            incident_dir = to_pixel / np.linalg.norm(to_pixel)
            
            # 计算理想的折射方向
            distance_from_center = math.sqrt(x_center**2 + y_center**2)
            max_distance = math.sqrt(2) * lens_size / 2
            
            # 根据距离中心的比例确定偏转角
            deviation_factor = distance_from_center / max_distance
            deviation_angle = deviation_factor * half_spread
            
            # 计算偏转后的方向（从透镜中心向外偏转）
            direction_from_center = np.array([x_center, y_center, 0])
            if np.linalg.norm(direction_from_center) > 0:
                direction_from_center = direction_from_center / np.linalg.norm(direction_from_center)
            else:
                direction_from_center = np.array([1, 0, 0])  # 防止除零
                
            # 构建偏转后的方向
            cos_dev = math.cos(deviation_angle)
            sin_dev = math.sin(deviation_angle)
            
            # 计算理想的折射方向
            refracted_dir = np.array([0, 0, 1]) * cos_dev + sin_dev * direction_from_center
            refracted_dir = refracted_dir / np.linalg.norm(refracted_dir)
            
            # 使用斯涅尔定律计算所需的表面法线
            # n₁(i × n) = n₂(r × n)，其中i是入射向量，r是折射向量，n是法线
            # 我们知道i和想要的r，需要计算n
            
            # 迭代求解最佳法线方向
            # 初始猜测：入射方向和折射方向的中间向量
            normal = incident_dir + refracted_dir
            normal = normal / np.linalg.norm(normal)
            
            # 迭代多次以获得更准确的法线
            for _ in range(5):
                # 计算当前法线下的折射方向
                cos_i = -np.dot(incident_dir, normal)
                if cos_i < 0:
                    normal = -normal
                    cos_i = -cos_i
                
                # 使用斯涅尔定律计算折射角
                sin_i = math.sqrt(1 - cos_i**2)
                sin_r = (1.0 / refractive_index) * sin_i
                
                # 检查全反射条件
                if sin_r >= 1.0:
                    # 如果发生全反射，调整法线
                    normal = (normal + incident_dir) / 2
                    normal = normal / np.linalg.norm(normal)
                    continue
                
                cos_r = math.sqrt(1 - sin_r**2)
                
                # 计算折射方向
                computed_refracted = (1.0/refractive_index) * incident_dir + (1.0/refractive_index * cos_i - cos_r) * normal
                
                # 调整法线以更接近目标折射方向
                adjustment = np.cross(np.cross(computed_refracted, refracted_dir), incident_dir)
                normal = normal + 0.1 * adjustment
                normal = normal / np.linalg.norm(normal)
            
            # 生成六面体的顶点（8个顶点）
            v_base = len(vertices)  # 当前顶点索引基数
            
            # 计算前表面四个顶点
            corners = [
                [x_center - pixel_size/2, y_center - pixel_size/2],  # 左下
                [x_center + pixel_size/2, y_center - pixel_size/2],  # 右下
                [x_center + pixel_size/2, y_center + pixel_size/2],  # 右上
                [x_center - pixel_size/2, y_center + pixel_size/2],  # 左上
            ]
            
            front_vertices = []
            for cx, cy in corners:
                # 根据平面方程计算z坐标
                plane_d = -np.dot(normal, np.array([x_center, y_center, z_center]))
                cz = (-plane_d - normal[0]*cx - normal[1]*cy) / normal[2]
                front_vertices.append([cx, cy, cz])
                vertices.append([cx, cy, cz])
            
            # 添加后表面的顶点（平行于xy平面，z = -thickness）
            # 确保背面是完全平的，z坐标固定为-thickness
            for cx, cy, _ in front_vertices:  # 忽略前表面的z值，固定为-thickness
                vertices.append([cx, cy, -thickness])
            
            # 添加六个面（前、后、四个侧面）
            # 前面（两个三角形）
            faces.append([v_base + 0, v_base + 1, v_base + 2])
            faces.append([v_base + 0, v_base + 2, v_base + 3])
            
            # 后面（两个三角形）
            v_back_base = v_base + 4  # 后面顶点的基索引
            faces.append([v_back_base + 0, v_back_base + 2, v_back_base + 1])
            faces.append([v_back_base + 0, v_back_base + 3, v_back_base + 2])
            
            # 四个侧面（每个侧面两个三角形）
            # 左侧面
            faces.append([v_base + 0, v_base + 3, v_back_base + 3])
            faces.append([v_base + 0, v_back_base + 3, v_back_base + 0])
            
            # 右侧面
            faces.append([v_base + 1, v_back_base + 1, v_back_base + 2])
            faces.append([v_base + 1, v_base + 2, v_back_base + 2])
            
            # 上侧面
            faces.append([v_base + 3, v_base + 2, v_back_base + 2])
            faces.append([v_base + 3, v_back_base + 2, v_back_base + 3])
            
            # 下侧面
            faces.append([v_base + 0, v_back_base + 0, v_back_base + 1])
            faces.append([v_base + 0, v_base + 1, v_back_base + 1])
    
    # 写入OBJ文件
    with open(output_file, 'w') as f:
        f.write("# 像素透镜OBJ实体模型\n")
        f.write(f"# 尺寸: {lens_size}cm x {lens_size}cm\n")
        f.write(f"# 细分: {divisions}x{divisions}\n")
        f.write(f"# 观察点距离: {observer_distance}cm\n")
        f.write(f"# 扩散角度: {spread_angle}度\n")
        f.write(f"# 厚度: {thickness}cm\n")
        f.write(f"# 折射率: {refractive_index} (亚克力)\n\n")
        
        # 写入顶点
        for v in vertices:
            f.write(f"v {v[0]:.6f} {v[1]:.6f} {v[2]:.6f}\n")
        
        # 写入面
        for face in faces:
            f.write(f"f {face[0]+1} {face[1]+1} {face[2]+1}\n")  # OBJ索引从1开始
    
    print(f"生成完成! 实体模型已保存为 {output_file}")
    print(f"总计顶点数: {len(vertices)}")
    print(f"总计面数: {len(faces)}")

# 生成模型
if __name__ == "__main__":
    generate_pixel_lens_obj()
