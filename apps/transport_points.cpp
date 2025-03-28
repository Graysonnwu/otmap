// This file is part of otmap, an optimal transport solver.
//
// Copyright (C) 2017-2018 Gael Guennebaud <gael.guennebaud@inria.fr>
// Copyright (C) 2017 Georges Nader
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <iostream>
#include "otsolver_2dgrid.h"
#include "common/otsolver_options.h"
#include "utils/eigen_addons.h"
#include "common/image_utils.h"
#include "common/generic_tasks.h"
#include "utils/BenchTimer.h"
#include <surface_mesh/Surface_mesh.h>

#include "normal_integration/normal_integration.h"
#include "normal_integration/mesh.h"

using namespace Eigen;
using namespace surface_mesh;
using namespace otmap;

void output_usage()
{
  std::cout << "___________________________________\n";
  std::cout << "Real usage:\n";
  std::cout << "./transport_points -in_src <src_image> -in_trg <trg_image> -res <grid_size> -focal <focal_length> -o <output_obj>\n\n";
  std::cout << "Example:\n";
  std::cout << "./transport_points -in_src ../data/source.png -in_trg ../data/mita.png -res 500 -o ../output.obj\n\n";
}

struct CLIopts : CLI_OTSolverOptions
{
  std::string filename_src;
  std::string filename_trg;
  std::string filename_out;

  VectorXi ores;
  double pt_scale;
  std::string pattern;
  bool inv_mode;
  bool export_maps;

  bool reflective_caustics;
  bool point_light;
  std::vector<double> lightPosition;
  std::vector<double> trgPlanePosition;
  std::vector<double> trgPlaneRotation;
  std::vector<double> trgPlaneScale;

  uint resolution;
  uint resolution_height;
  double focal_length;
  double thickness;

  std::string out_prefix;

  void set_default()
  {
    filename_src = "";
    filename_trg = "";
    filename_out = "./output.obj";

    ores.resize(1); ores.setZero();
    ores(0) = 1;

    out_prefix = "";

    pt_scale = 1;
    export_maps = 0;
    pattern = "";

    reflective_caustics = false;
    point_light = false;
    lightPosition = {0, 0, 1};
    trgPlanePosition = {0, 0, -3};
    trgPlaneRotation = {0, 0, 0};
    trgPlaneScale = {1, 1, 1};

    resolution = 100;
    resolution_height = 100;
    focal_length = 3.0;
    thickness = 0.2;

    CLI_OTSolverOptions::set_default();
  }

  static bool parseTriple(std::vector<std::string>& value, std::vector<double>& target, const std::string& errorMsg) {
    if (!value.empty() && value.size() == 3) {
        try {
          if (value[0][0] == '_') value[0][0] = '-'; // 用下划线代替负号，防止被解析器误读
          if (value[1][0] == '_') value[1][0] = '-';
          if (value[2][0] == '_') value[2][0] = '-';
          target = {std::stod(value[0]), std::stod(value[1]), std::stod(value[2])}; return true;
        }
        catch (const std::exception& e) {
          std::cerr << errorMsg << ": " << e.what() << std::endl;
        }
    } else {
        std::cerr << errorMsg << std::endl;
    }
    return false;
  }

  bool load(const InputParser &args)
  {
    set_default();

    CLI_OTSolverOptions::load(args);

    std::vector<std::string> value;

    if(args.getCmdOption("-in_src", value))
      filename_src = value[0];
    else
      return false;

    if(args.getCmdOption("-in_trg", value))
      filename_trg = value[0];
    else
      return false;

    if(args.getCmdOption("-o", value))
      filename_out = value[0];

    /*if(args.getCmdOption("-points", value))
      pattern = value[0];
    else
      return false;*/

    if(args.getCmdOption("-res", value))
      resolution = std::atof(value[0].c_str());
    
    if(args.getCmdOption("-res_h", value))
      resolution_height = std::atof(value[0].c_str());
    else
      resolution_height = resolution;

    if(args.getCmdOption("-ptscale", value))
      pt_scale = std::atof(value[0].c_str());
    
    if(args.cmdOptionExists("-export_maps"))
      export_maps = true;
    
    if(args.cmdOptionExists("-reflect")){
      reflective_caustics = true;
      lightPosition = {-1, 0, -1}; // 默认直角反射
      trgPlanePosition = {1, 0, -1};
      trgPlaneRotation = {0, -90, 0};
      trgPlaneScale = {1, 1, 1};
    }

    if(args.cmdOptionExists("-point_light")){
      point_light = true;
      //...
    }

    if(args.getCmdOption("-light_pos", value)){
      if(!parseTriple(value, lightPosition, "light position must contain 3 valid values")){return false;}
    }

    if (args.getCmdOption("-trg_pos", value)) {
      if(!parseTriple(value, trgPlanePosition, "target plane position must contain 3 valid values")){return false;}
    }

    if (args.getCmdOption("-trg_rot", value)) {
      if(!parseTriple(value, trgPlaneRotation, "target plane rotation must contain 3 valid values")){return false;}
    }

    if (args.getCmdOption("-trg_scale", value)) {
      if(!parseTriple(value, trgPlaneScale, "target plane scale must contain 3 valid values")){return false;}
    }

    if(args.getCmdOption("-focal", value)){
      focal_length = std::atof(value[0].c_str());
      trgPlanePosition[2] = -focal_length;
    }

    if(args.getCmdOption("-thick", value)){
      thickness = std::atof(value[0].c_str());
    }

    return true;
  }
};

std::vector<double> normalize_vec(std::vector<double> p1) {
    std::vector<double> vec(3);
    double squared_len = 0;
    for (int i=0; i<p1.size(); i++) {
        squared_len += p1[i] * p1[i];
    }

    double len = std::sqrt(squared_len);

    for (int i=0; i<p1.size(); i++) {
        vec[i] = p1[i] / len;
    }

    return vec;
}

// 计算旋转矩阵
std::vector<std::vector<double>> computeRotationMatrix(double x, double y, double z) {
    double cx = cos(x), cy = cos(y), cz = cos(z);
    double sx = sin(x), sy = sin(y), sz = sin(z);

    // ZYX 顺序的欧拉角旋转矩阵
    return {
        {cy * cz, cz * sx * sy - cx * sz, cx * cz * sy + sx * sz},
        {cy * sz, cx * cz + sx * sy * sz, -cz * sx + cx * sy * sz},
        {-sy, cy * sx, cx * cy}
    };
}

// 3D 向量矩阵乘法
std::vector<double> applyRotation(const std::vector<double>& vec, const std::vector<std::vector<double>>& rotMat) {
    return {
        vec[0] * rotMat[0][0] + vec[1] * rotMat[0][1] + vec[2] * rotMat[0][2],
        vec[0] * rotMat[1][0] + vec[1] * rotMat[1][1] + vec[2] * rotMat[1][2],
        vec[0] * rotMat[2][0] + vec[1] * rotMat[2][1] + vec[2] * rotMat[2][2]
    };
}

// 对 target_pts 应用变换
std::vector<std::vector<double>> applyTargetPlaneTransform(
    std::vector<std::vector<double>> target_pts,
    const std::vector<double>& targetPlanePosition,  // {x, y, z}
    const std::vector<double>& targetPlaneRotation, // {x, y, z} in degrees
    const std::vector<double>& targetPlaneScale     // {x, y, z}
) {
    // 计算旋转矩阵
    auto rotationMatrix = computeRotationMatrix(
        targetPlaneRotation[0] * M_PI / 180.0,
        targetPlaneRotation[1] * M_PI / 180.0,
        targetPlaneRotation[2] * M_PI / 180.0);

    for (auto& pt : target_pts) {
        // 缩放
        pt[0] *= targetPlaneScale[0];
        pt[1] *= targetPlaneScale[1];
        pt[2] *= targetPlaneScale[2]; //z轴缩放无用

        // 旋转
        pt = applyRotation(pt, rotationMatrix);

        // 平移
        pt[0] += targetPlanePosition[0];
        pt[1] += targetPlanePosition[1];
        pt[2] += targetPlanePosition[2]; // -focal_len
    }

    return target_pts;
}

void generate_mask(
    const MatrixXd &image,
    const std::vector<std::vector<double>> &source_points,
    double ratio,
    std::vector<int> &mask_indices)
{
    int h = image.rows();
    int w = image.cols();
    mask_indices.resize(source_points.size(), -1); // default -1 (not within shape)

    for (uint i = 0; i < source_points.size(); i++)
    {
        // source_points index -> image coordinate
        int x = i % h;
        int y = i / h;

        // boundary check, ensure the mapped point is in the image range
        if (x < 0 || x >= w || y < 0 || y >= h){
            printf("x: %d, y: %d\n", x, y);
            continue;
        }

        // if within shape, mask = 1
        if (image(y, x) > 0.0001) // threshold
        {
            mask_indices[i] = 1;
        }
    }
}

// 计算前表面折射点
std::vector<double> computeFrontRefractionPoint(
    const std::vector<double>& vertexPosition,
    const std::vector<double>& pointLightPosition,
    double refractive_index,
    double front_surface_z
) {
    const double n1 = 1.0; // 空气折射率
    const double n2 = refractive_index; // 材料折射率
    const double CONVERGENCE_THRESHOLD = 1e-6;
    const int MAX_ITERATIONS = 1000;
    const double LEARNING_RATE = 0.05;

    // 计算从光源指向顶点的光线方向，S -> T
    std::vector<double> rayDir = {
        vertexPosition[0] - pointLightPosition[0],
        vertexPosition[1] - pointLightPosition[1],
        vertexPosition[2] - pointLightPosition[2]
    };
    rayDir = normalize_vec(rayDir);

    // 计算射线与平面的交点
    double t = (front_surface_z - pointLightPosition[2]) / rayDir[2];
    if (std::abs(rayDir[2]) < 1e-6) {
        std::cerr << "Warning: No intersect point for vertex at (" << vertexPosition[0] << ", " << vertexPosition[1] << ", " << vertexPosition[2] << ")\n";
        return pointLightPosition;
    }

    std::vector<double> intersectionPoint = {
        pointLightPosition[0] + t * rayDir[0],
        pointLightPosition[1] + t * rayDir[1],
        pointLightPosition[2] + t * rayDir[2]
    };

    double x_p = intersectionPoint[0];
    double y_p = intersectionPoint[1];

    // 使用梯度下降法最小化光程，S(光源) -> P(前表面) -> T(后表面)
    bool converged = false;
    for(int iter = 0; iter < MAX_ITERATIONS; ++iter) {
        std::vector<double> P = {x_p, y_p, front_surface_z};
        std::vector<double> S_to_P = {
            P[0] - pointLightPosition[0],
            P[1] - pointLightPosition[1],
            P[2] - pointLightPosition[2]
        };
        std::vector<double> P_to_T = {
            vertexPosition[0] - P[0],
            vertexPosition[1] - P[1],
            vertexPosition[2] - P[2]
        };

        double dist_S_P = std::sqrt(S_to_P[0] * S_to_P[0] + S_to_P[1] * S_to_P[1] + S_to_P[2] * S_to_P[2]);
        double dist_P_T = std::sqrt(P_to_T[0] * P_to_T[0] + P_to_T[1] * P_to_T[1] + P_to_T[2] * P_to_T[2]);

        // 防止除零
        if (dist_S_P < 1e-6 || dist_P_T < 1e-6) {
            std::cerr << "Warning: Lens too thin for vertex at (" << vertexPosition[0] << ", " << vertexPosition[1] << ", " << vertexPosition[2] << ")\n";
            break;
        }

        // 计算光程
        // double L = dist_S_P / n2 + dist_P_T / n1;

        // 计算梯度
        double dL_dx = ((P[0] - pointLightPosition[0]) / (n2 * dist_S_P)) - 
                      ((vertexPosition[0] - P[0]) / (n1 * dist_P_T));
        double dL_dy = ((P[1] - pointLightPosition[1]) / (n2 * dist_S_P)) - 
                      ((vertexPosition[1] - P[1]) / (n1 * dist_P_T));

        // 更新x_p和y_p
        x_p -= LEARNING_RATE * dL_dx;
        y_p -= LEARNING_RATE * dL_dy;

        // 检查收敛性
        if (std::abs(dL_dx) < CONVERGENCE_THRESHOLD && 
            std::abs(dL_dy) < CONVERGENCE_THRESHOLD) {
            converged = true;
            // std::cout << "Converged after " << iter << " iterations. Final refraction point at (" << x_p << ", " << y_p << ", " << front_surface_z << ") for vertex at (" << vertexPosition[0] << ", " << vertexPosition[1] << ", " << vertexPosition[2] << ")\n";
            break;
        }
    }

    if (!converged) {
        std::cerr << "Warning: Gradient descent did not converge. Final refraction point at (" << x_p << ", " << y_p << ", " << front_surface_z << ") for vertex at (" << vertexPosition[0] << ", " << vertexPosition[1] << ", " << vertexPosition[2] << ")\n";
    }

    return {x_p, y_p, front_surface_z};
}

//compute the desired normals
std::vector<std::vector<double>> fresnelMapping(
    std::vector<std::vector<double>> &vertices,
    std::vector<std::vector<double>> &target_pts,
    double refractive_index,
    const std::vector<double>& lightPosition,
    bool reflective_caustics = false,
    bool point_light = false,
    double front_surface_z = 0.2
    ) 
{

    std::vector<std::vector<double>> desiredNormals;

    for(int i = 0; i < vertices.size(); i++) {
        // 光线方向。光源在下方，vertex -> light
        std::vector<double> incidentLight;
        if (point_light) { // 点光源，要考虑前表面折射，不能简单两者相减
            std::vector<double> frontRefrationPoint = computeFrontRefractionPoint(vertices[i], lightPosition, refractive_index, front_surface_z);
            incidentLight = {
                frontRefrationPoint[0] - vertices[i][0],
                frontRefrationPoint[1] - vertices[i][1],
                frontRefrationPoint[2] - vertices[i][2]
            };
        } else { // 平行光
            incidentLight = lightPosition;
        }
        incidentLight = normalize_vec(incidentLight);

        // 投射点方向。投射面在上方，vertex -> target
        std::vector<double> transmitted = {
            target_pts[i][0] - vertices[i][0],
            target_pts[i][1] - vertices[i][1],
            target_pts[i][2] - vertices[i][2]
        };
        transmitted = normalize_vec(transmitted);

        // 计算法线
        std::vector<double> normal;
        if (reflective_caustics) { // 反射焦散
            normal = {
                transmitted[0] + incidentLight[0],
                transmitted[1] + incidentLight[1],
                transmitted[2] + incidentLight[2]
            };
        } else { // 折射焦散
            normal = {
                (transmitted[0] + incidentLight[0] * refractive_index) * -1.0f,
                (transmitted[1] + incidentLight[1] * refractive_index) * -1.0f,
                (transmitted[2] + incidentLight[2] * refractive_index) * -1.0f,
            };
        }
        normal = normalize_vec(normal);

        desiredNormals.push_back(normal);
    }

    return desiredNormals;
}

int main(int argc, char** argv)
{
  setlocale(LC_ALL,"C");

  InputParser input(argc, argv);

  if(input.cmdOptionExists("-help") || input.cmdOptionExists("-h")){
    output_usage();
    return 0;
  }

  CLIopts opts;
  if(!opts.load(input)){
    std::cerr << "invalid input" << std::endl;
    output_usage();
    return EXIT_FAILURE;
  }

  double ratio = static_cast<double>(opts.resolution_height) / opts.resolution; // ratio必须大于1，即h>w，res_h>res
  printf("ratio = %5.2f\r\n", ratio);
  // Mesh mesh(1.0, 1.0/2, opts.resolution, opts.resolution/2);
  // Mesh mesh(1.0, opts.height, opts.resolution, static_cast<unsigned int>(opts.resolution*opts.height));
  // Mesh mesh(1.0, opts.height, opts.resolution, (unsigned int)(opts.resolution*opts.height));
  // Mesh mesh(1.0, ratio, opts.resolution, opts.resolution_height);
  Mesh mesh(ratio, 1.0, opts.resolution_height, opts.resolution); // 生成模型行列是反的，所以在此故意取反 ratio 1
  
  std::vector<Eigen::Vector2d> tile;
  normal_integration normal_int;
  /*if(!load_point_cloud_dat(opts.pattern, tile))
  {
    std::cerr << "Error loading tile \"" << opts.pattern << "\"\n";
    return EXIT_FAILURE;
  }*/

  for (int i=0; i<mesh.source_points.size(); i++)
  {
    Eigen::Vector2d point = {mesh.source_points[i][0]/ratio, mesh.source_points[i][1]/ratio}; // 1 1/ratio
    tile.push_back(point);
  }

  MatrixXd src_image;
  load_resize_image(opts.filename_src.c_str(), src_image, opts.resolution_height, opts.resolution_height); //将图像缩放为和trg_pts一样大小

  std::vector<int> mask_indices;
  generate_mask(src_image, mesh.source_points, ratio, mask_indices);

  normal_int.initialize_data(mesh, mask_indices);

  GridBasedTransportSolver otsolver;
  otsolver.set_verbose_level(opts.verbose_level-1);

  if(opts.verbose_level>=1)
    std::cout << "Generate transport map...\n";

  MatrixXd density_src;
  MatrixXd density_trg;
  if(!load_input_density(opts.filename_src, density_src))
  {
    std::cout << "Failed to load input \"" << opts.filename_src << "\" -> abort.";
    exit(EXIT_FAILURE);
  }

  if(!load_input_density(opts.filename_trg, density_trg))
  {
    std::cout << "Failed to load input \"" << opts.filename_trg << "\" -> abort.";
    exit(EXIT_FAILURE);
  }

  if(density_src.maxCoeff()>1.)
    density_src = density_src / density_src.maxCoeff(); //normalize

  if(density_trg.maxCoeff()>1.)
    density_trg = density_trg / density_trg.maxCoeff(); //normalize

  //density_src = 1. - density_src.array();
  //density_trg = 1. - density_trg.array();

  //save_image((opts.out_prefix + "_target.png").c_str(), 1.-density_src.array());

  BenchTimer t_solver_init, t_solver_compute, t_generate_uniform, t_inverse;

  t_solver_init.start();
  otsolver.init(density_src.rows());
  t_solver_init.stop();

  t_solver_compute.start();
  TransportMap tmap_src = otsolver.solve(vec(density_src), opts.solver_opt);
  t_solver_compute.stop();

  std::cout << "transport map for source image found" << std::endl;

  t_solver_init.start();
  otsolver.init(density_trg.rows());
  t_solver_init.stop();

  t_solver_compute.start();
  TransportMap tmap_trg = otsolver.solve(vec(density_trg), opts.solver_opt);
  t_solver_compute.stop();

  std::cout << "transport map for target image found" << std::endl;

  std::cout << "STATS solver -- init: " << t_solver_init.value(REAL_TIMER) << "s  solve: " << t_solver_compute.value(REAL_TIMER) << "s\n";

  for(unsigned int i=0; i<opts.ores.size(); ++i){
    
    std::vector<Eigen::Vector2d> points;

    t_inverse.start();
    apply_forward_map(tmap_src, tile, opts.verbose_level-1); // 该行即为插值步骤
    apply_inverse_map(tmap_trg, tile, opts.verbose_level-1);
    t_inverse.stop();

    std::vector<std::vector<double>> trg_pts;
    for (int i=0; i<mesh.source_points.size(); i++)
    {
      std::vector<double> point = {tile[i].x()*ratio, tile[i].y()*ratio, 0}; // ratio 1
      trg_pts.push_back(point);
    }
    
    // read points from file (migrate points from Caustic_Design repo)
    if(false) {
      std::vector<std::vector<double>> src_pts_temp;
      std::vector<std::vector<double>> trg_pts_temp;
      std::ifstream file("pixlens-test.dat");
      // std::ifstream file("head-grey-otmap.dat");
      if (file.is_open()) {
        std::string line;
        while (std::getline(file, line)) {
          std::istringstream iss(line);
          double src_x, src_y, trg_x, trg_y;
          // 读取每一行的数据
          if (iss >> src_x >> src_y >> trg_x >> trg_y) {
            std::vector<double> src_pt = {src_x, src_y, 0};
            std::vector<double> trg_pt = {trg_x, trg_y, 0};
            src_pts_temp.push_back(src_pt);
            trg_pts_temp.push_back(trg_pt);
          }
        }
        file.close();
      } else {std::cerr << "Cannot open file" << std::endl;}

      if (src_pts_temp.size() != mesh.source_points.size()) {
        std::cerr << "Point number mismatch" << std::endl;
      }

      for (int i=0; i<mesh.source_points.size(); i++){
        mesh.source_points[i] = src_pts_temp[i];
        trg_pts[i] = trg_pts_temp[i];
      }
    }

    // if not within shape, move points to source location
    for (uint i = 0; i < trg_pts.size(); i++)
    {
        if (mask_indices[i] == -1)
        {
            trg_pts[i] = mesh.source_points[i];
        }
    }
    
    // size(sp) = size(tp) = max(res,res_h)^2
    // spx[0,ratio], spy[0,1], tpx[0,ratio], tpy[0,1]
    // 将mesh.source_points和trg_pts的xy坐标中心点都从(ratio/2,0.5)移到原点
    for (int i=0; i<mesh.source_points.size(); i++)
    {
      mesh.source_points[i][0] -= ratio/2;
      mesh.source_points[i][1] -= 0.5;
      trg_pts[i][0] -= ratio/2;
      trg_pts[i][1] -= 0.5;
    }
    // spx[-ratio/2,ratio/2], spy[-0.5,0.5], tpx[-ratio/2,ratio/2], tpy[-0.5,0.5]

    double tp_max_x = -1000000.0f, tp_max_y = -1000000.0f, tp_max_z = -1000000.0f;
    double tp_min_x = 1000000.0f, tp_min_y = 1000000.0f, tp_min_z = 1000000.0f;
    for (int i=0; i<trg_pts.size(); i++)
    {
      if (mask_indices[i] == 1) {
        if (tp_max_x < trg_pts[i][0]) {tp_max_x = trg_pts[i][0];}
        if (tp_max_y < trg_pts[i][1]) {tp_max_y = trg_pts[i][1];}
        if (tp_max_z < trg_pts[i][2]) {tp_max_z = trg_pts[i][2];}
        if (tp_min_x > trg_pts[i][0]) {tp_min_x = trg_pts[i][0];}
        if (tp_min_y > trg_pts[i][1]) {tp_min_y = trg_pts[i][1];}
        if (tp_min_z > trg_pts[i][2]) {tp_min_z = trg_pts[i][2];}
      }
    }
    printf("tp_max_x  = %5.2f, tp_max_y  = %5.2f, tp_max_z  = %5.2f\r\n", tp_max_x, tp_max_y, tp_max_z);
    printf("tp_min_x  = %5.2f, tp_min_y  = %5.2f, tp_min_z  = %5.2f\r\n", tp_min_x, tp_min_y, tp_min_z);
    
    std::vector<double> lightPosition = opts.lightPosition;
    std::vector<double> targetPlanePosition = opts.trgPlanePosition;
    std::vector<double> targetPlaneRotation = opts.trgPlaneRotation;
    std::vector<double> targetPlaneScale = opts.trgPlaneScale;
    
    printf("lightPosition       = %.2f, %.2f, %.2f\r\n", lightPosition[0], lightPosition[1], lightPosition[2]);
    printf("targetPlanePosition = %.2f, %.2f, %.2f\r\n", targetPlanePosition[0], targetPlanePosition[1], targetPlanePosition[2]);
    printf("targetPlaneRotation = %.2f, %.2f, %.2f\r\n", targetPlaneRotation[0], targetPlaneRotation[1], targetPlaneRotation[2]);
    printf("targetPlaneScale    = %.2f, %.2f, %.2f\r\n", targetPlaneScale[0], targetPlaneScale[1], targetPlaneScale[2]);

    trg_pts = applyTargetPlaneTransform(trg_pts, targetPlanePosition, targetPlaneRotation, targetPlaneScale); // 默认仅平移到 (0,0,focal_len)

    double tp_max_x2 = -1000000.0f, tp_max_y2 = -1000000.0f, tp_max_z2 = -1000000.0f;
    double tp_min_x2 = 1000000.0f, tp_min_y2 = 1000000.0f, tp_min_z2 = 1000000.0f;
    for (int i=0; i<trg_pts.size(); i++)
    {
      if (mask_indices[i] == 1) {
        if (tp_max_x2 < trg_pts[i][0]) {tp_max_x2 = trg_pts[i][0];}
        if (tp_max_y2 < trg_pts[i][1]) {tp_max_y2 = trg_pts[i][1];}
        if (tp_max_z2 < trg_pts[i][2]) {tp_max_z2 = trg_pts[i][2];}
        if (tp_min_x2 > trg_pts[i][0]) {tp_min_x2 = trg_pts[i][0];}
        if (tp_min_y2 > trg_pts[i][1]) {tp_min_y2 = trg_pts[i][1];}
        if (tp_min_z2 > trg_pts[i][2]) {tp_min_z2 = trg_pts[i][2];}
      }
    }
    printf("tp_max_x2 = %5.2f, tp_max_y2 = %5.2f, tp_max_z2 = %5.2f\r\n", tp_max_x2, tp_max_y2, tp_max_z2);
    printf("tp_min_x2 = %5.2f, tp_min_y2 = %5.2f, tp_min_z2 = %5.2f\r\n", tp_min_x2, tp_min_y2, tp_min_z2);

    // 将mesh.source_points和trg_pts都输出到.dat文件中
    if(false){
      std::string filename_dat = opts.filename_out.substr(0, opts.filename_out.find_last_of(".")) + ".dat";
      std::ofstream file(filename_dat);
      for (int i=0; i<mesh.source_points.size(); i++)
      {
        file << mesh.source_points[i][0] << " " << mesh.source_points[i][1] << " " << trg_pts[i][0] << " " << trg_pts[i][1] << " " << mask_indices[i] << "\n";
      }
      file.close();
    }

    // move black points backward by 0.05
    for (int i=0; i<mesh.source_points.size(); i++)
    {
      if (mask_indices[i] == -1){
        mesh.source_points[i][2] += 0.1;
      }
    }

    std::vector<std::vector<double>> desired_normals;
    double front_surface_z = opts.thickness;

    for (int i=0; i<(opts.point_light ? 10 : 5); i++)
    {
        double src_max_z = 0;
        for (int j=0; j<mesh.source_points.size(); j++) {if (mesh.source_points[j][2] > src_max_z) {src_max_z = mesh.source_points[j][2];}}
        if(i<3){
          front_surface_z = src_max_z + opts.thickness;
          std::cout << "front_surface_z = src_max_z + thickness = " << src_max_z << " + " << opts.thickness << " = " << front_surface_z << std::endl;
        }
        else{ //迭代超过3次后，front_surface_z保持不变，加快收敛速度
          std::cout << "thickness = front_surface_z - src_max_z = " << front_surface_z << " - " << src_max_z << " = " << front_surface_z - src_max_z << std::endl;
        }
        std::vector<std::vector<double>> normals = fresnelMapping(mesh.source_points, trg_pts, 1.49, lightPosition, opts.reflective_caustics, opts.point_light, front_surface_z);

        //desired_normals.clear();

        // make a copy of the original positions of the vertices
        /*for (int i = 0; i < mesh.source_points.size(); i++) {
            std::vector<double> trg_normal = {normals[0][i], normals[1][i], normals[2][i]};
            desired_normals.push_back(trg_normal);
        }*/

        normal_int.perform_normal_integration(mesh, normals, mask_indices);

        std::cout << "Outer Loop: " << i+1 << "/" << (opts.point_light ? 10 : 5) << " done." << std::endl;
        std::cout << "--------------------------------------------------ʕ•ᴥ•ʔっ♡--------------------------------------------------" << std::endl;
    }

    /*std::string filename = opts.out_prefix + "_" + std::to_string(opts.ores[i]);
    save_point_cloud_dat(filename + ".dat", tile);
    save_point_cloud_eps(filename + ".eps", tile, opts.pt_scale);

    std::cout << " # " << opts.ores[i] << "/" << tile.size()
                << "  ;  gen: " << t_generate_uniform.value(REAL_TIMER)
                << "s  ;  bvh+inverse: " << t_inverse.value(REAL_TIMER) << "s\n";*/

    // mesh.save_solid_obj_source(0.2, "../output.obj"); // 0.2 为厚度
    mesh.save_solid_obj_source(front_surface_z, opts.filename_out);
    std::cout << "Obj file saved to " << opts.filename_out << std::endl;
  }
}
