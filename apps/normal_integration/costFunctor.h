#ifndef COSTFUNCTOR_H
#define COSTFUNCTOR_H

#include "ceres/ceres.h"
#include "glog/logging.h"

#define EBAR_DETH 39
#define EINT_WEIGHT 0.1
#define EBAR_WEIGHT 0.0
#define EDIR_WEIGHT 5.0
#define EREG_WEIGHT 20.0


/******************************************************/
/*************** Ceres-Helper-Methods *****************/
/******************************************************/

using namespace std;

template<typename T> void cross(T* v1, T* v2, T* result);
template<typename T> void calcFaceNormal(const T* const v1, const T* const v2, const T* const v3, T* result);
template<typename T> T angle(T* v1, T* v2);
template<typename T> void normalize(T* v);
template<typename T> T evaluateInt(const T* const vertex, const T** neighbors, uint nNeighbors, const vector<int> & neighborMap, const std::vector<double> &desiredNormal, T* result);
template<typename T> void calcVertexNormal(const T* vertex, T* result, T** faceNormals, const T** neighbors, const vector<int> & neighborMap);
template<typename T> void evaluateReg(const T** const allVertices, const float* L, uint nVertices, T* res);

/******************************************************/
/*************** Ceres-Helper-Methods *****************/
/******************************************************/

// For EInt
template<typename T> void cross(T* v1, T* v2, T* result)
{
    result[0] = T(v1[1]*v2[2] - v1[2]*v2[1]);
    result[1] = T(v1[2]*v2[0] - v1[0]*v2[2]);
    result[2] = T(v1[0]*v2[1] - v1[1]*v2[0]);
}

template<typename T> void calcFaceNormal(const T* const v1, const T* const v2, const T* const v3, T* result)
{
    T vertex1[3];
    T vertex2[3];

    vertex1[0] = v2[0] - v1[0];
    vertex1[1] = v2[1] - v1[1];
    vertex1[2] = v2[2] - v1[2];

    vertex2[0] = v3[0] - v1[0];
    vertex2[1] = v3[1] - v1[1];
    vertex2[2] = v3[2] - v1[2];

    cross(vertex1, vertex2, result);
}

template<typename T> T angle(T* v1, T* v2)
{
    return ceres::acos(v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]);
}

template<typename T> void normalize(T* v)
{
    T sum = T(0);
    for (uint i=0; i<3; i++)
    {
        sum += v[i] * v[i];
    }

    sum = ceres::sqrt(sum);

    //if (sum > T(10e-12)) {
    for (uint i=0; i<3; i++)
    {
        v[i] = v[i]/sum;
    }
    //}
}

// Eint，即顶点法线与OTM计算出的法线的残差，允许部分偏离，该项由Edir限制
template<typename T> T evaluateInt(const T* const vertex, const T** neighbors, uint nNeighbors, const vector<int> & neighborMap, const std::vector<double> &desiredNormal, T* result)
{

    T** faceNormals = new T*[neighborMap.size()/2];

    // -- face normals
    for(uint i=0; i<neighborMap.size(); i+=2)
    {
        int faceIndex = i/2;
        faceNormals[faceIndex] = new T[3];
        T* faceNormal = faceNormals[faceIndex];
        const T* n1 = neighbors[neighborMap[i]];
        const T* n2 = neighbors[neighborMap[i+1]];

        calcFaceNormal(vertex, n1, n2, faceNormal); // 传入三个点，用叉积算出面法线

        if (faceNormal[2] < 0) {
            //faceNormal[0] = -faceNormal[0];
            //faceNormal[1] = -faceNormal[1];
            //faceNormal[2] = -faceNormal[2];
        }
    }

    // -- vertex normal
    T* vertexNormal = new T[3];
    calcVertexNormal(vertex, vertexNormal, faceNormals, neighbors, neighborMap); // 对相邻面的法线加权平均，算出顶点法线，权重为两边夹角大小


    // evaluation（顶点法线 - OTM预期法线）
    T x = vertexNormal[0] - T(desiredNormal[0]);
    T y = vertexNormal[1] - T(desiredNormal[1]);
    T z = vertexNormal[2] - T(desiredNormal[2]);

    T res = T(EINT_WEIGHT) * ceres::sqrt(x*x + y*y + z*z); // 这个才是二范数，与论文一致，但返回了没用到
    result[0] = x*T(EINT_WEIGHT); // 这个返回的是残差，而不是二范数，与论文有差别
    result[1] = y*T(EINT_WEIGHT);
    result[2] = z*T(EINT_WEIGHT);

    // -- clean
    delete[] vertexNormal;

    for(uint i=0; i<neighborMap.size(); i+=2)
    {
        delete[] faceNormals[i/2];
    }

    delete[] faceNormals;

    // -- done

    return res;
}

// 计算顶点法线，方法：对相邻面的法线加权平均，权重为(相邻面上与顶点相连)两条边的夹角大小
template<typename T> void calcVertexNormal(const T* vertex, T* result, T** faceNormals, const T** neighbors, const vector<int> & neighborMap)
{
    result[0] = result[1] = result[2] = T(0);

    T edge1Res[3];
    T edge2Res[3];
    T edge1[3];
    T edge2[3];
    T edge1Sub[3];
    T edge2Sub[3];

    for(uint i=0; i<neighborMap.size(); i+=2) // 遍历相邻面
    {

        for (uint j=0; j<3; j++) // 构造两邻边向量
        {
            edge1[j] = vertex[j];
            edge2[j] = vertex[j];
            edge1Sub[j] = neighbors[neighborMap[i]][j];
            edge2Sub[j] = neighbors[neighborMap[i+1]][j];
            
            edge1Res[j] = edge1[j] - edge1Sub[j];
            edge2Res[j] = edge2[j] - edge2Sub[j];
        }

        normalize(edge1Res);
        normalize(edge2Res);

        T incidentAngle = angle(edge1Res, edge2Res); // 两邻边夹角

        T* faceNormal = faceNormals[i/2];

        // 将夹角作为权重
        for (uint j=0; j<3; j++)
        {
            T val = faceNormal[j];
            result[j] += val * incidentAngle;
        }

    }

    normalize(result);
}


// Ereg（影响最大的一项），说白了就是计算邻居点坐标平均值和当前点坐标的差值，最小化该残差以确保网格表面平滑（孙宇欧论文中的E_lap，但他允许部分尖边）
template<typename T> void evaluateReg(const T** const allVertices, const float* L, uint nVertices, T* res)
{

    for (uint i=0; i<3; i++)
        res[i] = T(0);

    if(nVertices == uint(5)) // 标准情况，包括自身在内一共5个点，xyz坐标分别计算 -4*中+上+下+左+右
    {
        for(int i=0; i<nVertices; i++){
            for (int j=0; j<3; j++){
                res[j] += T(L[i]) * allVertices[i][j];  // L为简化的拉普拉斯矩阵，形式为[-邻居数, 1, 1, 1, 1...]，有多少邻居就有多少1
            }
        }
    }
    else // 非标准情况（边界），只叠加z分量，x分量和y分量强制设为0，乘对应拉普拉斯算子（也就是只看z轴上的偏离情况）
    {
        for(uint i=0; i<nVertices; i++)
        {
            res[0] = T(0);
            res[1] = T(0);
            res[2] += T(L[i]) * allVertices[i][2];
        }
    }

    for (uint i=0; i<3; i++)
        res[i] = res[i] * T(EREG_WEIGHT);
        // res[i] = atan(res[i]*1e2 * T(EREG_WEIGHT))/1e2; // 限制大值影响
}


/**********************************************/
/************** Cost-Functors *****************/
/**********************************************/

/********* EInt *********/

class CostFunctorEint8Neighbors{
public:
    CostFunctorEint8Neighbors(std::vector<double> &desiredNormal, vector<int> & neighborMap): desiredNormal(desiredNormal), neighborMap(neighborMap)
    {}

    template <typename T> bool operator()(const T* const vertex,
                                          const T* const neighbor1,
                                          const T* const neighbor2,
                                          const T* const neighbor3,
                                          const T* const neighbor4,
                                          const T* const neighbor5,
                                          const T* const neighbor6,
                                          const T* const neighbor7,
                                          const T* const neighbor8,
                                          T* residual) const
    {

        // -- preparation
        const T* neighbors[8];
        neighbors[0] = neighbor1;
        neighbors[1] = neighbor2;
        neighbors[2] = neighbor3;
        neighbors[3] = neighbor4;
        neighbors[4] = neighbor5;
        neighbors[5] = neighbor6;
        neighbors[6] = neighbor7;
        neighbors[7] = neighbor8;


        evaluateInt(vertex, neighbors, 8, neighborMap, desiredNormal, residual);

        return true;
    }

private:
    std::vector<double> desiredNormal;
    vector<int> neighborMap;

};


class CostFunctorEint7Neighbors{
public:
    CostFunctorEint7Neighbors(std::vector<double> &desiredNormal, vector<int> & neighborMap): desiredNormal(desiredNormal), neighborMap(neighborMap)
    {}

    template <typename T> bool operator()(const T* const vertex,
                                          const T* const neighbor1,
                                          const T* const neighbor2,
                                          const T* const neighbor3,
                                          const T* const neighbor4,
                                          const T* const neighbor5,
                                          const T* const neighbor6,
                                          const T* const neighbor7,
                                          T* residual) const
    {

        // -- preparation
        const T* neighbors[7];
        neighbors[0] = neighbor1;
        neighbors[1] = neighbor2;
        neighbors[2] = neighbor3;
        neighbors[3] = neighbor4;
        neighbors[4] = neighbor5;
        neighbors[5] = neighbor6;
        neighbors[6] = neighbor7;


        evaluateInt(vertex, neighbors, 7, neighborMap, desiredNormal, residual);

        return true;
    }

private:
    std::vector<double> desiredNormal;
    vector<int> neighborMap;

};


class CostFunctorEint6Neighbors{
public:
    CostFunctorEint6Neighbors(std::vector<double> &desiredNormal, vector<int> & neighborMap): desiredNormal(desiredNormal), neighborMap(neighborMap)
    {}

    template <typename T> bool operator()(const T* const vertex,
                                          const T* const neighbor1,
                                          const T* const neighbor2,
                                          const T* const neighbor3,
                                          const T* const neighbor4,
                                          const T* const neighbor5,
                                          const T* const neighbor6,
                                          T* residual) const
    {

        // -- preparation
        const T* neighbors[6];
        neighbors[0] = neighbor1;
        neighbors[1] = neighbor2;
        neighbors[2] = neighbor3;
        neighbors[3] = neighbor4;
        neighbors[4] = neighbor5;
        neighbors[5] = neighbor6;


        evaluateInt(vertex, neighbors, 6, neighborMap, desiredNormal, residual);

        return true;
    }

private:
    std::vector<double> desiredNormal;
    vector<int> neighborMap;

};

class CostFunctorEint5Neighbors{
public:
    CostFunctorEint5Neighbors(std::vector<double> &desiredNormal, vector<int> & neighborMap): desiredNormal(desiredNormal), neighborMap(neighborMap)
    {}

    template <typename T> bool operator()(const T* const vertex,
                                          const T* const neighbor1,
                                          const T* const neighbor2,
                                          const T* const neighbor3,
                                          const T* const neighbor4,
                                          const T* const neighbor5,
                                          T* residual) const
    {

        // -- preparation
        const T* neighbors[5];
        neighbors[0] = neighbor1;
        neighbors[1] = neighbor2;
        neighbors[2] = neighbor3;
        neighbors[3] = neighbor4;
        neighbors[4] = neighbor5;


        evaluateInt(vertex, neighbors, 5, neighborMap, desiredNormal, residual);

        return true;
    }

private:
    std::vector<double> desiredNormal;
    vector<int> neighborMap;

};


class CostFunctorEint4Neighbors{
public:
    CostFunctorEint4Neighbors(std::vector<double> &desiredNormal, vector<int> & neighborMap): desiredNormal(desiredNormal), neighborMap(neighborMap)
    {}

    template <typename T> bool operator()(const T* const vertex,
                                          const T* const neighbor1,
                                          const T* const neighbor2,
                                          const T* const neighbor3,
                                          const T* const neighbor4,
                                          T* residual) const
    {

        // -- preparation
        const T* neighbors[4];
        neighbors[0] = neighbor1;
        neighbors[1] = neighbor2;
        neighbors[2] = neighbor3;
        neighbors[3] = neighbor4;


        evaluateInt(vertex, neighbors, 4, neighborMap, desiredNormal, residual);

        return true;
    }

private:
    std::vector<double> desiredNormal;
    vector<int> neighborMap;

};


class CostFunctorEint3Neighbors{
public:
    CostFunctorEint3Neighbors(std::vector<double> &desiredNormal, vector<int> & neighborMap): desiredNormal(desiredNormal), neighborMap(neighborMap)
    {}

    template <typename T> bool operator()(const T* const vertex,
                                          const T* const neighbor1,
                                          const T* const neighbor2,
                                          const T* const neighbor3,
                                          T* residual) const
    {

        // -- preparation
        const T* neighbors[3];
        neighbors[0] = neighbor1;
        neighbors[1] = neighbor2;
        neighbors[2] = neighbor3;


        evaluateInt(vertex, neighbors, 3, neighborMap, desiredNormal, residual);

        return true;
    }

private:
    std::vector<double> desiredNormal;
    vector<int> neighborMap;

};


class CostFunctorEint2Neighbors{
public:
    CostFunctorEint2Neighbors(std::vector<double> &desiredNormal, vector<int> & neighborMap): desiredNormal(desiredNormal), neighborMap(neighborMap)
    {}

    template <typename T> bool operator()(const T* const vertex,
                                          const T* const neighbor1,
                                          const T* const neighbor2,
                                          T* residual) const
    {

        // -- preparation
        const T* neighbors[2];
        neighbors[0] = neighbor1;
        neighbors[1] = neighbor2;

        evaluateInt(vertex, neighbors, 2, neighborMap, desiredNormal, residual);

        return true;
    }

private:
    std::vector<double> desiredNormal;
    vector<int> neighborMap;

};

class CostFunctorEint1Neighbors{
public:
    CostFunctorEint1Neighbors(std::vector<double> &desiredNormal, vector<int> & neighborMap): desiredNormal(desiredNormal), neighborMap(neighborMap)
    {}

    template <typename T> bool operator()(const T* const vertex,
                                          const T* const neighbor1,
                                          T* residual) const
    {
        // -- preparation
        const T* neighbors[1];
        neighbors[0] = neighbor1;

        evaluateInt(vertex, neighbors, 1, neighborMap, desiredNormal, residual);

        return true;
    }

private:
    std::vector<double> desiredNormal;
    vector<int> neighborMap;
};



/****/

inline void printLaplacien(float* L, int size){
    std::cout<<"|"<<std::endl;
    for (uint i=0; i<size; i++){
        std::cout<<L[i]<<std::endl;
    }
    std::cout<<"|"<<std::endl;
};

/********* EReg *********/ // 邻居数大于4的都被移除了，因为不可能（网格点邻居无斜交邻居），而且从未被引用
class CostFunctorEreg4Neighbors{
public:
    CostFunctorEreg4Neighbors(vector<int> neighbors){

        uint size = neighbors.size() + 1;
        L = new float[size];
        L[0] = - float(neighbors.size());
        for (uint i=1; i<size; i++)
        {
            L[i] = 1;
        }// L = [-4 1 1 1 1]
    }

    ~CostFunctorEreg4Neighbors()
    {
        delete[] L;
    }

    template <typename T> bool operator()(const T* const vertex,
                                          const T* const neighbor1,
                                          const T* const neighbor2,
                                          const T* const neighbor3,
                                          const T* const neighbor4,
                                          T* residual) const
    {

        const T* allVertices[5];
        allVertices[0] = vertex;
        allVertices[1] = neighbor1;
        allVertices[2] = neighbor2;
        allVertices[3] = neighbor3;
        allVertices[4] = neighbor4;


        evaluateReg(allVertices, L, 5, residual);

        return true;
    }

private:
    float* L;
};

class CostFunctorEreg3Neighbors{
public:
    CostFunctorEreg3Neighbors(vector<int> neighbors){

        uint size = neighbors.size() + 1;
        L = new float[size];
        L[0] = - float(neighbors.size());
        for (uint i=1; i<size; i++)
        {
            L[i] = 1;
        } // L = [-3, 1, 1, 1]
    }

    ~CostFunctorEreg3Neighbors()
    {
        delete[] L;
    }

    template <typename T> bool operator()(const T* const vertex,
                                          const T* const neighbor1,
                                          const T* const neighbor2,
                                          const T* const neighbor3,
                                          T* residual) const
    {

        const T* allVertices[4];
        allVertices[0] = vertex;
        allVertices[1] = neighbor1;
        allVertices[2] = neighbor2;
        allVertices[3] = neighbor3;


        evaluateReg(allVertices, L, 4, residual);

        return true;
    }

private:
    float* L;
};

class CostFunctorEreg2Neighbors{
public:
    CostFunctorEreg2Neighbors(vector<int> neighbors){

        uint size = neighbors.size() + 1;
        L = new float[size];
        L[0] = - float(neighbors.size());
        for (uint i=1; i<size; i++)
        {
            L[i] = 1;
        } // L = [-2, 1, 1]
    }

    ~CostFunctorEreg2Neighbors()
    {
        delete[] L;
    }

    template <typename T> bool operator()(const T* const vertex,
                                          const T* const neighbor1,
                                          const T* const neighbor2,
                                          T* residual) const
    {

        const T* allVertices[3];
        allVertices[0] = vertex;
        allVertices[1] = neighbor1;
        allVertices[2] = neighbor2;


        evaluateReg(allVertices, L, 3, residual);

        return true;
    }

private:
    float* L;
};

class CostFunctorEreg1Neighbor{
public:
    CostFunctorEreg1Neighbor(vector<int> neighbors){

        uint size = neighbors.size() + 1;
        L = new float[size];
        L[0] = - float(neighbors.size());
        for (uint i=1; i<size; i++)
        {
            L[i] = 1;
        } // L = [-1, 1]
    }

    ~CostFunctorEreg1Neighbor()
    {
        delete[] L;
    }

    template <typename T> bool operator()(const T* const vertex,
                                          const T* const neighbor1,
                                          T* residual) const
    {

        const T* allVertices[2];
        allVertices[0] = vertex;
        allVertices[1] = neighbor1;


        evaluateReg(allVertices, L, 2, residual);

        return true;
    }

private:
    float* L;
};

/********* EBar *********/

/*class CostFunctorEbar2 {
public:
    CostFunctorEbar2 (std::vector<double>* receiverPosition): receiverPosition(receiverPosition)
    {}

    template <typename T> bool operator()(const T* const vertex, T* e) const{

        T dth = T(EBAR_DETH);

        T distance = vertex[0] - T(receiverPosition->x);
        e[0] = T(EBAR_WEIGHT) * ceres::fmax(T(0), - ceres::log( (T(1)-distance) + dth) );

        return true;
    }

private:
    std::vector<double>* receiverPosition;

};*/



/********* EDir *********/
// 即顶点偏移的距离（xy平面上，不计算z轴）
class CostFunctorEdir2{
public:
    CostFunctorEdir2(std::vector<double> * s, float weightMultiplicator): source(s), weightMultiplicator(weightMultiplicator) {
    }

    template <typename T> bool operator()(const T* const vertex, T* e) const{
        // z - proj(..) will only return the difference in x- and y-coordinates. we may ignore the rest
        T x = vertex[0] - T((*source)[0]);
        T y = vertex[1] - T((*source)[1]);
        e[0] = T(EDIR_WEIGHT*weightMultiplicator) * x;
        e[1] = T(EDIR_WEIGHT*weightMultiplicator) * y;
        e[2] = T(0);
        return true;
    }

private:
    std::vector<double> * source;
    float weightMultiplicator;
};


#endif // COSTFUNCTOR_H
