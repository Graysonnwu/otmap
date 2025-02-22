#include "normal_integration.h"
#include <thread>

normal_integration::normal_integration(/* args */)
{
}

normal_integration::~normal_integration()
{
}


void normal_integration::addResidualBlocks(Mesh &mesh, Problem *problem, uint vertexIndex, vector<int> &neighbors, vector<int> &neighborMap, vector<int> &gridNeighbors, double *vertices, std::vector<double> &trg_normal)
{
    float weightMult = 1.0;
    if(mesh.is_border(vertexIndex)) // we have an edge. Set weight for edir extremely high
        weightMult = 10000;

    // EDir depends on the original position
    CostFunction* cost_function_edir =
            new AutoDiffCostFunction<CostFunctorEdir2, 3, 3>(new CostFunctorEdir2(&x_sources[vertexIndex], weightMult));

    problem->AddResidualBlock(
                cost_function_edir,
                NULL, // no loss function
                &vertices[vertexIndex*3]
                );


    if(!mesh.is_border(vertexIndex)){ //not an edge we optimize the normals
        // For eint we have several functors, each for a different amount of neighbors
        switch(neighbors.size()){
        case 1:
        {

            CostFunction* cost_function_eint1 =
                new AutoDiffCostFunction<CostFunctorEint1Neighbors, 3, 3, 3>(new CostFunctorEint1Neighbors(trg_normal, neighborMap));

            problem->AddResidualBlock(cost_function_eint1, NULL,
                                       &vertices[vertexIndex*3], // vertex
                                       &vertices[neighbors[0]*3]); // and the neighbors..

            break;
            }
        
        case 2:
        {

            CostFunction* cost_function_eint2 =
                new AutoDiffCostFunction<CostFunctorEint2Neighbors, 3, 3, 3, 3>(new CostFunctorEint2Neighbors(trg_normal, neighborMap));

            problem->AddResidualBlock(cost_function_eint2, NULL,
                                       &vertices[vertexIndex*3], // vertex
                                       &vertices[neighbors[0]*3], // and the neighbors..
                                       &vertices[neighbors[1]*3]);


            break;
            }

        case 3:
        {

            CostFunction* cost_function_eint3 =
                new AutoDiffCostFunction<CostFunctorEint3Neighbors, 3, 3, 3, 3, 3>(new CostFunctorEint3Neighbors(trg_normal, neighborMap));

            problem->AddResidualBlock(cost_function_eint3, NULL,
                                       &vertices[vertexIndex*3], // vertex
                                       &vertices[neighbors[0]*3], // and the neighbors..
                                       &vertices[neighbors[1]*3],
                                       &vertices[neighbors[2]*3]);


            break;
            }

        case 4:
        {

            CostFunction* cost_function_eint4 =
                new AutoDiffCostFunction<CostFunctorEint4Neighbors, 3, 3, 3, 3, 3, 3>(new CostFunctorEint4Neighbors(trg_normal, neighborMap));

            problem->AddResidualBlock(cost_function_eint4, NULL,
                                       &vertices[vertexIndex*3], // vertex
                                       &vertices[neighbors[0]*3], // and the neighbors..
                                       &vertices[neighbors[1]*3],
                                       &vertices[neighbors[2]*3],
                                       &vertices[neighbors[3]*3]);


            break;
            }

        case 5:
        {
            CostFunction* cost_function_eint5 =
                new AutoDiffCostFunction<CostFunctorEint5Neighbors, 3, 3, 3, 3, 3, 3, 3>(new CostFunctorEint5Neighbors(trg_normal, neighborMap));

            problem->AddResidualBlock(cost_function_eint5, NULL,
                           &vertices[vertexIndex*3], // vertex
                           &vertices[neighbors[0]*3], // and the neighbors..
                           &vertices[neighbors[1]*3],
                           &vertices[neighbors[2]*3],
                           &vertices[neighbors[3]*3],
                           &vertices[neighbors[4]*3]);



            break;
        }

        case 6:
        {
            CostFunction* cost_function_eint6 =
                new AutoDiffCostFunction<CostFunctorEint6Neighbors, 3, 3, 3, 3, 3, 3, 3, 3>(new CostFunctorEint6Neighbors(trg_normal, neighborMap));

            problem->AddResidualBlock( cost_function_eint6, NULL,
                       &vertices[vertexIndex*3], // vertex
                       &vertices[neighbors[0]*3], // and the neighbors..
                       &vertices[neighbors[1]*3],
                       &vertices[neighbors[2]*3],
                       &vertices[neighbors[3]*3],
                       &vertices[neighbors[4]*3],
                       &vertices[neighbors[5]*3]);



            break;
        }

        case 7:
        {
            CostFunction* cost_function_eint7 =
                new AutoDiffCostFunction<CostFunctorEint7Neighbors, 3, 3, 3, 3, 3, 3, 3, 3, 3>(new CostFunctorEint7Neighbors(trg_normal, neighborMap));

            problem->AddResidualBlock( cost_function_eint7, NULL,
                       &vertices[vertexIndex*3], // vertex
                       &vertices[neighbors[0]*3], // and the neighbors..
                       &vertices[neighbors[1]*3],
                       &vertices[neighbors[2]*3],
                       &vertices[neighbors[3]*3],
                       &vertices[neighbors[4]*3],
                       &vertices[neighbors[5]*3],
                       &vertices[neighbors[6]*3]);


            break;
        }

        case 8:
        {
            CostFunction* cost_function_eint8 =
                new AutoDiffCostFunction<CostFunctorEint8Neighbors, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3>(new CostFunctorEint8Neighbors(trg_normal, neighborMap));

            problem->AddResidualBlock( cost_function_eint8, NULL,
                           &vertices[vertexIndex*3], // vertex
                           &vertices[neighbors[0]*3], // and the neighbors..
                           &vertices[neighbors[1]*3],
                           &vertices[neighbors[2]*3],
                           &vertices[neighbors[3]*3],
                           &vertices[neighbors[4]*3],
                           &vertices[neighbors[5]*3],
                           &vertices[neighbors[6]*3],
                           &vertices[neighbors[7]*3]);

            break;

        }

        } // switch end
    }

    // regularization term with 4-neighborhood
    CostFunction* ereg;
    switch(gridNeighbors.size())
    {
        case 1:
        {
            ereg = new AutoDiffCostFunction<CostFunctorEreg1Neighbor, 3, 3, 3>(new CostFunctorEreg1Neighbor(gridNeighbors));
            problem->AddResidualBlock(
                        ereg,
                        NULL,
                        &vertices[vertexIndex*3],
                        &vertices[gridNeighbors[0]*3]
                        );
            break;
        }
        case 2:
        {
            ereg = new AutoDiffCostFunction<CostFunctorEreg2Neighbors, 3, 3, 3, 3>(new CostFunctorEreg2Neighbors(gridNeighbors));
            problem->AddResidualBlock(
                        ereg,
                        NULL,
                        &vertices[vertexIndex*3],
                        &vertices[gridNeighbors[0]*3],
                        &vertices[gridNeighbors[1]*3]
                        );
            break;
        }
        case 3:
        {
            ereg = new AutoDiffCostFunction<CostFunctorEreg3Neighbors, 3, 3, 3, 3, 3>(new CostFunctorEreg3Neighbors(gridNeighbors));
            problem->AddResidualBlock(
                        ereg,
                        NULL,
                        &vertices[vertexIndex*3],
                        &vertices[gridNeighbors[0]*3],
                        &vertices[gridNeighbors[1]*3],
                        &vertices[gridNeighbors[2]*3]
                        );
            break;
        }
        case 4:
        {
            ereg = new AutoDiffCostFunction<CostFunctorEreg4Neighbors, 3, 3, 3, 3, 3, 3>(new CostFunctorEreg4Neighbors(gridNeighbors));
            problem->AddResidualBlock(
                        ereg,
                        NULL,
                        &vertices[vertexIndex*3],
                        &vertices[gridNeighbors[0]*3],
                        &vertices[gridNeighbors[1]*3],
                        &vertices[gridNeighbors[2]*3],
                        &vertices[gridNeighbors[3]*3]
                        );
            break;
        }
    }
}

/**
 * @brief TargetOptimization::gatherVertexInformation Gathers information about neighbors of current vertex.
 * @param vertex The current vertex
 * @param vertexIndex The index of the current vertex
 * @param neighborList [out] (Eint) A list of indices to the neighbors of the current vertex. Each neighbor is only included once
 * @param neighborMap [out] (Eint) A list of indices to indices to neighbors. The value is an index to the neighborList. Is used to handle several references to one neighbor (e.g. for neighboring faces)
 * @param gridNeighbors [out] (Ereg) A list indices to of horizontal and vertical neighbors
 * @param mask_indices A list of indices of whether the vertex is within the shape (1) or not (-1)
 */
void normal_integration::gatherVertexInformation(Mesh &mesh, uint vertexIndex, std::vector<int> &neighborList, std::vector<int> &neighborMap, std::vector<int> & gridNeighbors, const std::vector<int>& mask_indices)
{
    mesh.find_vertex_connectivity(vertexIndex, neighborList, neighborMap, mask_indices);

    int left_vtx = 0;
    int right_vtx = 0;
    int top_vtx = 0;
    int bot_vtx = 0;
    
    mesh.get_vertex_neighbor_ids(vertexIndex, left_vtx, right_vtx, top_vtx, bot_vtx);

    if (left_vtx != -1 && mask_indices[left_vtx] == 1) {
        gridNeighbors.push_back(left_vtx);
    }

    if (right_vtx != -1 && mask_indices[right_vtx] == 1) {
        gridNeighbors.push_back(right_vtx);
    }

    if (top_vtx != -1 && mask_indices[top_vtx] == 1) {
        gridNeighbors.push_back(top_vtx);
    }

    if (bot_vtx != -1 && mask_indices[bot_vtx] == 1) {
        gridNeighbors.push_back(bot_vtx);
    }
}

void normal_integration::initialize_data(Mesh &mesh, const std::vector<int>& mask_indices) {
    neighborsPerVertex.resize(mesh.source_points.size());
    neighborMapPerVertex.resize(mesh.source_points.size());
    eightNeighborsPerVertex.resize(mesh.source_points.size());

        // gather information for each vertex to optimize
    for(uint i = 0; i < mesh.source_points.size(); i++)
    {
        std::vector<int> neighbors;
        std::vector<int> neighborMap;
        std::vector<int> eightNeighbors;

        gatherVertexInformation(mesh, i, neighbors, neighborMap, eightNeighbors, mask_indices);

        neighborsPerVertex[i] = neighbors;
        neighborMapPerVertex[i] = neighborMap;
        eightNeighborsPerVertex[i] = eightNeighbors;
    }

    vertices = new double[3*mesh.source_points.size()];
}

void normal_integration::perform_normal_integration(Mesh &mesh, std::vector<std::vector<double>> &desired_normals, const std::vector<int>& mask_indices) {
    // make a copy of the original positions of the vertices
    for (int i = 0; i < mesh.source_points.size(); i++) {
        x_sources.push_back(mesh.source_points[i]);
    }

    // put all positions in one big list that we access later
    for(uint i=0; i<mesh.source_points.size(); i++)
    {
        vertices[3*i + 0] = mesh.source_points[i][0];
        vertices[3*i + 1] = mesh.source_points[i][1];
        vertices[3*i + 2] = mesh.source_points[i][2];
    }

    Problem prob;

    // iterate over all vertices and add the corresponding residual blocks
    for(uint i=0; i<mesh.source_points.size(); i++)
    {
      if (mask_indices[i] == 1) // only optimize vertices within shape
        addResidualBlocks(mesh, &prob, i, neighborsPerVertex[i], neighborMapPerVertex[i], eightNeighborsPerVertex[i], vertices, desired_normals[i]);
    }

    int threads_num = std::thread::hardware_concurrency();
    if (threads_num == 0) threads_num = 4;
    std::cout << "Using " << threads_num << " threads." << std::endl;

    Solver::Options options;
    options.minimizer_progress_to_stdout = true;
    options.linear_solver_type = ceres::ITERATIVE_SCHUR;
    options.max_num_iterations = 200;
    options.dense_linear_algebra_library_type = ceres::LAPACK;
    options.num_threads = threads_num;

    std::string error;
    if(!options.IsValid(&error))
    {
        std::cout << "Options not valid: " << error << std::endl;
    }

    Solver::Summary summary;
    Solve(options, &prob, &summary);

    // std::cout << summary.FullReport() << std::endl;
    std::cout << summary.BriefReport() << std::endl;

    std::vector<double> * pos;
    for(uint i=0; i<mesh.source_points.size(); i++)
    {
        mesh.source_points[i][0] = vertices[3*i + 0];
        mesh.source_points[i][1] = vertices[3*i + 1];
        mesh.source_points[i][2] = vertices[3*i + 2];
    }
}