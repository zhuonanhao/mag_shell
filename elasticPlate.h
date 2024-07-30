#ifndef ELASTICPLATE_H
#define ELASTICPLATE_H

#include "eigenIncludes.h"
#include <fstream>

struct edgeElement
{
	int nv_1;
	int nv_2;

	Vector3d x_1;
	Vector3d x_2;

	double refLength;
	double edgeLength;

	VectorXi arrayNum;
};

struct triangularElement
{
	int nv_1;
	int nv_2;
	int nv_3;

	Vector3d x_1;
	Vector3d x_2;
	Vector3d x_3;

	Vector3d v_c;

	Vector3d a1, a2, a3;
	Vector3d norm;

	Matrix3d abar;
	Matrix3d abarinv;

	double area;
};

struct bendingElement
{
	int triangular_1;
	int triangular_2;

	int nv_1;
	int nv_2;
	int nv_3;
	int nv_4;

	VectorXi arrayNum;

	// x_1, x_2, and x_3 for triangular_1;
	// x_4, x_2, and x_3 for triangular_2;
	Vector3d x_1;
	Vector3d x_2;
	Vector3d x_3;
	Vector3d x_4;

	Vector3d e_1;
	Vector3d e_2;
	Vector3d e_3;
	Vector3d e_4;

	Vector3d n_1;
	Vector3d n_2;
	Vector3d norm_1;
	Vector3d norm_2;

	Vector3d nBar;
};

class elasticPlate
{
	public:
	elasticPlate(double m_YoungM, double m_density, double m_thickness, double m_Possion, double m_dt);
	~elasticPlate();

	double YoungM;
	double thickness;
	double Possion;
	double dt;
	double density;

	double ratio;

	Vector3d getVertex(int i);
	Vector3d getVertexOld(int i);
	Vector3d getVelocity(int i);

	void setOneVertexBoundaryCondition(double position, int i, int k);

	VectorXd x;
	VectorXd x0;
	VectorXd u;

	std::vector<Vector3d> v_nodes;
	std::vector<edgeElement> v_edgeElement;
	std::vector<triangularElement> v_triangularElement;
	std::vector<bendingElement> v_bendingElement;
	std::vector<Vector2i> constraint;
	std::vector<Vector3d> br_vec;
	std::vector<double> timeSeries;
    std::vector<Vector3d> ba_vec;

	int nv;
	int edgeNum;
	int triangularNum;
	int bendingNum;

	int ndof;
	int uncons;
	int ncons;

	double deltaArea;
	double deltaLength;

	void setupGeometry();

	void setVertexBoundaryCondition(Vector3d position, int k);

	void readInputNodes();
	void readInputEdge();
	void readInputTriangular();
	void computeBendingPair();

	void updateEdgePair();
	void updateBendingPair();

	// boundary conditions
	int* isConstrained;
	int getIfConstrained(int k);
	int* unconstrainedMap;
	int* fullToUnconsMap;
	void setup();
	void setupMap();

	void updateTimeStep();
	void updateGuess();
	void updateNewtonMethod(VectorXd m_motion);
	void prepareForIteration();

	VectorXd massArray;
	void setupMass();
	VectorXi boundaryIndex;

	double EA;
	double EI;

	double length;
    double width;
    double deltaEdge;
    double height;

    int xNodeNum;
    int yNodeNum;

    std::vector<Vector3d> nodes;
    std::vector<int> startnode;
    std::vector<Vector2i> edge;
    std::vector<Vector3i> triangular;

    std::vector<int> boundaryArray_1;
    std::vector<int> boundaryArray_2;

    int temp;

	private:
};

#endif
