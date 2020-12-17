//#include "tools.h"

// sofcla.h
// headers to define all the class use in sofea.cpp
int EndOfFile = 0;
//#define PI 3.14159265359
using namespace std;
void NextWord(FILE *file, char *ptr) // what is the next word/Number
{
	int Letter, inc = 0;
	Letter = fgetc(file);
	while ((Letter != (int)'\n' && Letter != '\t' && Letter != ',' && Letter != ' ' || inc == 0) && Letter != EOF)
	{
		if ((Letter == (int)'\n' || Letter == '\t' || Letter == ',' || Letter == ' ') && inc == 0)
			inc = 0;
		else if (Letter != (int)'\n' && Letter != '\t' && Letter != ',' && Letter != ' ' && Letter != EOF)
		{
			*ptr = (char)Letter;
			ptr++;
			inc++;
		}
		Letter = fgetc(file);
		if (Letter == EOF)
			::EndOfFile = 1;
	}
	*ptr = '\0';
}

void SortNumber(long *NbLink, int index)
{
	int Temp;
	if ((index - 1) != 0) // break the recursivity
	{
		SortNumber(&NbLink[1], index - 1);
		for (int inc = 0; inc < index - 1; inc++)
			if (NbLink[inc] > NbLink[inc + 1])
			{
				Temp = NbLink[inc];
				NbLink[inc] = NbLink[inc + 1];
				NbLink[inc + 1] = Temp;
			}
			else
				inc = index;
	}
}

//-----------NODE--------------------------------------------------------
class Node
{
public:			   // Caution No private because couldn't use some other function!!
				   // So to change
	long number;   // What's the node's name ?
	double nx, ny; // Nodes Coordonates
	long DOF;	   // Number of DOF for this element
	long CumSum;   // "Adress" in the vectors.
	void Init(FILE *file);
};

void Node::Init(FILE *file) // function to initialize data during the file reading
{
	char Command[100];
	sscanf(Command, "%ld", &number);
	NextWord(file, &Command[0]);
	sscanf(Command, "%le", &nx);
	NextWord(file, &Command[0]);
	sscanf(Command, "%le", &ny);
	DOF = 0;
	CumSum = 0;
}

//----------NODE----------------------------------------------------------
class Material
{
public:
	long number;		// What's the Material's name ?
	double ex, poisson; // Material properties...
						// ex = Young modulus on x
						// poisson = Poisson Coef

	// put the constructor to create the class
};

//-----------LINK1---------------------------------------------------------
class Link1 // Ansys code for 2D-spar
{
private:
	long number;		   // What's the element's name ?
	long NbLink[2];		   // define how many connection  And see alla bout the array!!
	long Mat;			   // Material
	long R;				   // Real
	double Angle;		   // Angle between the x-axis and the bar
	double Length;		   // lenght of the bar
	double StiffMat[4][4]; // Stiffness Matrix dim = 2
public:
	void CalculateLength(Node *DNodes, Material *DMat);
	Link1(FILE *file, long Num, Node *DNodes);
	void Display(int m);
	void AssembleMatrix(double **Global, int MaxDOF, Node *DNodes);
	void DefineNodeDOF(Node *DNodes);
};

void Link1::CalculateLength(Node *DNodes, Material *DMat)
{

	double CosA, SinA;
	Length = sqrt(pow(((DNodes[(NbLink[0] - 1)].nx) - (DNodes[(NbLink[1] - 1)].nx)), 2) + pow((DNodes[(NbLink[0] - 1)].ny) - (DNodes[(NbLink[1] - 1)].ny), 2));
	CosA = ((DNodes[(NbLink[0] - 1)].nx) - (DNodes[(NbLink[1] - 1)].nx)) / Length;
	SinA = ((DNodes[(NbLink[0] - 1)].ny) - (DNodes[(NbLink[1] - 1)].ny)) / Length;
	Angle = acos(CosA) * 180 / PI;
	StiffMat[0][0] = DMat[(Mat - 1)].ex / Length * pow(CosA, 2);
	StiffMat[0][1] = DMat[(Mat - 1)].ex / Length * SinA * CosA;
	StiffMat[0][2] = -StiffMat[0][0];
	StiffMat[0][3] = -StiffMat[0][1];
	StiffMat[1][0] = StiffMat[0][1];
	StiffMat[1][1] = DMat[(Mat - 1)].ex / Length * pow(SinA, 2);
	StiffMat[1][2] = -StiffMat[0][1];
	StiffMat[1][3] = -StiffMat[1][1];
	StiffMat[2][0] = -StiffMat[0][0];
	StiffMat[2][1] = -StiffMat[0][1];
	StiffMat[2][2] = StiffMat[0][0];
	StiffMat[2][3] = StiffMat[0][1];
	StiffMat[3][0] = -StiffMat[0][1];
	StiffMat[3][1] = -StiffMat[1][1];
	StiffMat[3][2] = StiffMat[0][1];
	StiffMat[3][3] = StiffMat[1][1];
}

Link1::Link1(FILE *file, long Num, Node *DNodes) // function to initialize data during the file reading
{
	char Command[100];

	number = Num;
	NextWord(file, &Command[0]);
	sscanf(Command, "%ld", &Mat);
	NextWord(file, &Command[0]);
	sscanf(Command, "%ld", &R);
	NextWord(file, &Command[0]);
	sscanf(Command, "%ld", &NbLink[0]);
	NextWord(file, &Command[0]);
	sscanf(Command, "%ld", &NbLink[1]);

	SortNumber(&NbLink[0], 2);
	DefineNodeDOF(DNodes);
}

void Link1::Display(int m)
{
	cout << "\n Element : "
		 << m + 1
		 << "\n";
	for (int inc = 0; inc < 4; inc++)
	{

		for (int in = 0; in < 4; in++)
			cout << StiffMat[inc][in] << "      ";
		cout << "\n";
	}
}

void Link1::AssembleMatrix(double **Global, int MaxDOF, Node *DNodes)
{
	// Function which assemble the Matrix !

	for (int inc = 0; inc < 2; inc++)			 //step for column (matrix 2x2)
		for (int inc2 = 0; inc2 < 2; inc2++)	 // step for row (matrix 2x2)
			for (int inc3 = 0; inc3 < 2; inc3++) // step for line
				for (int inc4 = 0; inc4 < 2; inc4++)
					Global[(DNodes[(NbLink[inc] - 1)].CumSum + inc3)][(DNodes[NbLink[inc2] - 1].CumSum + inc4)] += StiffMat[(inc * 2) + inc3][(inc2 * 2) + inc4];
}

void Link1::DefineNodeDOF(Node *DNodes)
{
	for (int inc = 0; inc < 2; inc++)
		if (!(DNodes[(NbLink[inc] - 1)].DOF > 2))
			DNodes[(NbLink[inc] - 1)].DOF = 2;
}

//--------Beam3-----------------------------------------------------------------------
class Beam3
{
private:
	long number;		   // What's the element's name ?
	long NbLink[2];		   // define how many connection  And see alla bout the array!!
	long Mat;			   // Material
	long R;				   // Real
	double EL;			   // EI/L
	double Area;		   // Cross area
	double Izz;			   // inertia moment about z-axis
	double Angle;		   // Angle between the x-axis and the bar
	double Length;		   // lenght of the bar
	double StiffMat[6][6]; // Stiffness Matrix dim = 3 (3*2Nodes)
public:
	void CalculateLength(Node *DNodes, Material *DMat);
	Beam3(FILE *file, long Num, Node *DNodes);
	void Display(int m);
	void AssembleMatrix(double **Global, int MaxDOF, Node *DNodes);
	void DefineNodeDOF(Node *DNodes);
};

void Beam3::CalculateLength(Node *DNodes, Material *DMat)
{
	double CosA, SinA;
	double A, B, K, M, F, G, P, H, Q, Z; // some variables
	Length = sqrt(pow(((DNodes[(NbLink[0] - 1)].nx) - (DNodes[(NbLink[1] - 1)].nx)), 2) + pow((DNodes[(NbLink[0] - 1)].ny) - (DNodes[(NbLink[1] - 1)].ny), 2));
	CosA = ((DNodes[(NbLink[0] - 1)].nx) - (DNodes[(NbLink[1] - 1)].nx)) / Length;
	SinA = ((DNodes[(NbLink[0] - 1)].ny) - (DNodes[(NbLink[1] - 1)].ny)) / Length;
	Angle = acos(CosA) * 180 / PI;
	EL = DMat[(Mat - 1)].ex / Length; // E/L ratio Young modulus over the length
	//Previous Calculations : Check them 'cos I'm not sure of the result...
	A = 4 * EL * Izz; //A=Area Pense que ce n'est pas bon du tout;
	Z = Area * EL;	  // Verify this is the area use instead of A
	B = 2 * EL * Izz;
	K = 12 * EL * Izz / Length / Length;
	M = 6 * EL * Izz / Length;
	F = Z * pow(CosA, 2) + K * pow(SinA, 2);
	G = (Z - K) * CosA * SinA;
	P = Z * pow(SinA, 2) + K * pow(CosA, 2);
	H = -M * SinA;
	Q = M * CosA;
	StiffMat[0][0] = F;
	StiffMat[0][1] = G;
	StiffMat[0][2] = H;
	StiffMat[0][3] = -F;
	StiffMat[0][4] = -G;
	StiffMat[0][5] = H;
	StiffMat[1][1] = P;
	StiffMat[1][2] = Q;
	StiffMat[1][3] = -G;
	StiffMat[1][4] = -P;
	StiffMat[1][5] = Q;
	StiffMat[2][2] = A;
	StiffMat[2][3] = -H;
	StiffMat[2][4] = -Q;
	StiffMat[2][5] = B;
	StiffMat[3][3] = F;
	StiffMat[3][4] = G;
	StiffMat[3][5] = -H;
	StiffMat[4][4] = P;
	StiffMat[4][5] = -Q;
	StiffMat[5][5] = A;
	// Matrix Symmetry
	for (int inc = 0; inc < 6; inc++)
		for (int inc2 = inc + 1; inc2 < 6; inc2++)
			StiffMat[inc2][inc] = StiffMat[inc][inc2]; // coordinate with the definition
}

Beam3::Beam3(FILE *file, long Num, Node *DNodes) // function to initialize data during the file reading
{
	char Command[100];
	number = Num;
	NextWord(file, &Command[0]);
	sscanf(Command, "%ld", &Mat); // read the material number
	NextWord(file, &Command[0]);
	sscanf(Command, "%le", &Area); // read the area
	NextWord(file, &Command[0]);
	sscanf(Command, "%le", &Izz); // read the moment of inertia
	NextWord(file, &Command[0]);
	sscanf(Command, "%ld", &NbLink[0]); // Stiffness matrix
	NextWord(file, &Command[0]);
	sscanf(Command, "%ld", &NbLink[1]);
	SortNumber(&NbLink[0], 2);

	DefineNodeDOF(DNodes);
}

void Beam3::Display(int m)
{
	cout << "\n Element : "
		 << m + 1
		 << "\n";
	for (int inc = 0; inc < 6; inc++)
	{
		for (int in = 0; in < 6; in++)
			cout << StiffMat[inc][in] << "      ";
		cout << "\n";
	}
}

void Beam3::AssembleMatrix(double **Global, int MaxDOF, Node *DNodes)
{
	// Function which assemble the Matrix !
	for (int inc = 0; inc < 2; inc++)				 //step for column (4 matrices 3x3)
		for (int inc2 = 0; inc2 < 2; inc2++)		 // step for row (4 matrices 3x3)
			for (int inc3 = 0; inc3 < 3; inc3++)	 // step for line (row)
				for (int inc4 = 0; inc4 < 3; inc4++) // step for line (column)
					Global[(DNodes[(NbLink[inc] - 1)].CumSum + inc3)][(DNodes[NbLink[inc2] - 1].CumSum + inc4)] += StiffMat[(inc * 3) + inc3][(inc2 * 3) + inc4];
}

void Beam3::DefineNodeDOF(Node *DNodes)
{
	for (int inc = 0; inc < 2; inc++)
		if (!(DNodes[(NbLink[inc] - 1)].DOF > 3))
			DNodes[(NbLink[inc] - 1)].DOF = 3;
}

//--------Plane2----------------------------------------------------------------------
//Comment on the Matrices used later on
//
//
//
//
//
//
//
//
//
//-Plane 2---------------------------------------------------------------------------

class Plane2
{
private:
	long number;
	long NbLink[3];	   // 3 connections
	long Mat;		   // Material
	long R;			   // Real
	double EL;		   // EI/L
	double Area;	   // Cross area
	double Izz;		   // inertia moment about z-axis
	double Angle;	   // Angle between the x-axis and the bar
	double Length;	   // lenght of the bar
	double Thickness;  //Thickness of the element
	double **StiffMat; //
	double ElemArea;   //Element's Area..
	void DefineNodeDOF(Node *DNodes);

public:
	void CalculateLength(Node *DNodes, Material *DMat);
	Plane2(FILE *file, long Num, Node *DNodes);
	void Display(int m);
	void AssembleMatrix(double **Global, int MaxDOF, Node *DNodes);
};

void Plane2::CalculateLength(Node *DNodes, Material *DMat)
{

	double *col, *d, twodet;
	int *indx, inc2, inc, indxx1, indxx2, indxy1, indxy2;
	double DMatrixConst;
	double **DMatrix = new double *[3];

	for (inc = 0; inc < 3; inc++) // 3*3
		DMatrix[inc] = new double[3];

	// Calculation of D
	//DMatrixConst
	DMatrixConst = DMat[(Mat - 1)].ex / (1 - pow(DMat[(Mat - 1)].poisson, 2));
	for (inc = 0; inc < 3; inc++)
		for (inc2 = 0; inc2 < 3; inc2++)
			DMatrix[inc][inc2] = 0;
	DMatrix[0][0] = DMatrixConst;
	DMatrix[0][1] = DMat[(Mat - 1)].poisson * DMatrixConst;
	DMatrix[1][0] = DMat[(Mat - 1)].poisson * DMatrixConst;
	DMatrix[1][1] = DMatrixConst;
	DMatrix[2][2] = 0.5 * (1 - DMat[(Mat - 1)].poisson) * DMatrixConst;

	//Calculation of the determinant

	twodet = DNodes[(NbLink[1] - 1)].nx * DNodes[(NbLink[2] - 1)].ny - DNodes[(NbLink[2] - 1)].nx * DNodes[(NbLink[1] - 1)].ny - DNodes[(NbLink[0] - 1)].nx * DNodes[(NbLink[2] - 1)].ny + DNodes[(NbLink[0] - 1)].nx * DNodes[(NbLink[1] - 1)].ny + DNodes[(NbLink[2] - 1)].nx * DNodes[(NbLink[0] - 1)].ny - DNodes[(NbLink[1] - 1)].nx * DNodes[(NbLink[0] - 1)].ny;

	// calculation of the local Matrix using the algorithms of the log book

	for (inc = 0; inc < 3; inc++)
		for (inc2 = 0; inc2 < 3; inc2++)
		{
			if (inc == 0)
			{
				indxx1 = 3 - 1;
				indxx2 = 2 - 1;
			}
			if (inc == 1)
			{
				indxx1 = 1 - 1;
				indxx2 = 3 - 1;
			}
			if (inc == 2)
			{
				indxx1 = 2 - 1;
				indxx2 = 1 - 1;
			}

			if (inc2 == 0)
			{
				indxy1 = 3 - 1;
				indxy2 = 2 - 1;
			}
			if (inc2 == 1)
			{
				indxy1 = 1 - 1;
				indxy2 = 3 - 1;
			}
			if (inc2 == 2)
			{
				indxy1 = 2 - 1;
				indxy2 = 1 - 1;
			}

			for (int inc3 = 0; inc3 < 2; inc3++)
				for (int inc4 = 0; inc4 < 2; inc4++)
				{
					if (inc3 == 0 && inc4 == 0)
						StiffMat[inc * 2 + inc3][inc2 * 2 + inc4] = Thickness / (2 * twodet) * (DMatrix[0][0] * ((DNodes[(NbLink[indxx2] - 1)].ny) - (DNodes[(NbLink[indxx1] - 1)].ny)) * ((DNodes[(NbLink[indxy2] - 1)].ny) - (DNodes[(NbLink[indxy1] - 1)].ny)) + DMatrix[2][2] * ((DNodes[(NbLink[indxy1] - 1)].nx) - (DNodes[(NbLink[indxy2] - 1)].nx)) * ((DNodes[(NbLink[indxx1] - 1)].nx) - (DNodes[(NbLink[indxx2] - 1)].nx)));

					if (inc3 == 1 && inc4 == 0)
						StiffMat[inc * 2 + inc3][inc2 * 2 + inc4] = Thickness / (2 * twodet) * (DMatrix[0][1] * ((DNodes[(NbLink[indxy2] - 1)].ny) - (DNodes[(NbLink[indxy1] - 1)].ny)) * ((DNodes[(NbLink[indxx1] - 1)].nx) - (DNodes[(NbLink[indxx2] - 1)].nx)) + DMatrix[2][2] * ((DNodes[(NbLink[indxy1] - 1)].nx) - (DNodes[(NbLink[indxy2] - 1)].nx)) * ((DNodes[(NbLink[indxx2] - 1)].ny) - (DNodes[(NbLink[indxx1] - 1)].ny)));

					if (inc3 == 0 && inc4 == 1)
						StiffMat[inc * 2 + inc3][inc2 * 2 + inc4] = Thickness / (2 * twodet) * (DMatrix[0][1] * ((DNodes[(NbLink[indxx1] - 1)].nx) - (DNodes[(NbLink[indxx2] - 1)].nx)) * ((DNodes[(NbLink[indxy2] - 1)].ny) - (DNodes[(NbLink[indxy1] - 1)].ny)) + DMatrix[2][2] * ((DNodes[(NbLink[indxx1] - 1)].nx) - (DNodes[(NbLink[indxx2] - 1)].nx)) * ((DNodes[(NbLink[indxy2] - 1)].ny) - (DNodes[(NbLink[indxy1] - 1)].ny)));

					if (inc3 == 1 && inc4 == 1)
						StiffMat[inc * 2 + inc3][inc2 * 2 + inc4] = Thickness / (2 * twodet) * (DMatrix[1][1] * ((DNodes[(NbLink[indxx1] - 1)].nx) - (DNodes[(NbLink[indxx2] - 1)].nx)) * ((DNodes[(NbLink[indxy1] - 1)].nx) - (DNodes[(NbLink[indxy2] - 1)].nx)) + DMatrix[2][2] * ((DNodes[(NbLink[indxy2] - 1)].ny) - (DNodes[(NbLink[indxy1] - 1)].ny)) * ((DNodes[(NbLink[indxx2] - 1)].ny) - (DNodes[(NbLink[indxx1] - 1)].ny)));
				}
		}
}

Plane2::Plane2(FILE *file, long Num, Node *DNodes) // function to initialize data during the file reading
{
	char Command[100];
	number = Num;
	NextWord(file, &Command[0]);
	sscanf(Command, "%ld", &Mat); // read the material number
	NextWord(file, &Command[0]);
	sscanf(Command, "%le", &Area); // read the area
	NextWord(file, &Command[0]);
	sscanf(Command, "%le", &Izz); // read the moment of inertia
	NextWord(file, &Command[0]);
	sscanf(Command, "%le", &Thickness); // read the thickness
	for (int inc = 0; inc < 3; inc++)
	{
		NextWord(file, &Command[0]);
		sscanf(Command, "%ld", &NbLink[inc]); // Stiffness matrix
	}

	StiffMat = new double *[6];
	for (int inc = 0; inc < 6; inc++) // 6*6
		StiffMat[inc] = new double[6];

	// Sort the number of Node
	SortNumber(&NbLink[0], 3);
	DefineNodeDOF(DNodes); // To read the boundary conditions
						   // as forces or displacements
}

void Plane2::Display(int m)
{
	cout << "\n Element : "
		 << m + 1
		 << "\n";
	for (int inc = 0; inc < 6; inc++)
	{
		for (int in = 0; in < 6; in++)
			cout << StiffMat[inc][in] << "      ";
		cout << "\n";
	}
}

void Plane2::AssembleMatrix(double **Global, int MaxDOF, Node *DNodes)
{
	// Function which assemble the Matrix !
	for (int inc = 0; inc < 3; inc++)				 //step for column (3 Nodes 2x2)
		for (int inc2 = 0; inc2 < 3; inc2++)		 // step for row (3 Nodes 2x2)
			for (int inc3 = 0; inc3 < 2; inc3++)	 // step for line (row)
				for (int inc4 = 0; inc4 < 2; inc4++) // step for line (column)
				{
					Global[(DNodes[(NbLink[inc] - 1)].CumSum + inc3)][(DNodes[NbLink[inc2] - 1].CumSum + inc4)] += StiffMat[(inc * 2) + inc3][(inc2 * 2) + inc4];
					inc4 = inc4;
				}
}
void Plane2::DefineNodeDOF(Node *DNodes)
{
	for (int inc = 0; inc < 3; inc++)
		if (!(DNodes[(NbLink[inc] - 1)].DOF > 2))
			DNodes[(NbLink[inc] - 1)].DOF = 2;
}

//--------Plane 42--------------------------------------------------------------------
//Comment on the Matrices used later on
//
//
//
//
//
//
//
//
//
//-Plane 42--------------------------------------------------------------------------

class Plane42
{
private:
	long number;
	long NbLink[4];	   // 4 connections
	long Mat;		   // Material
	double Thickness;  //Thickness of the element
	double **StiffMat; //
	void DefineNodeDOF(Node *DNodes);
	void Sort(Node *DNodes);

public:
	void CalculateMatrix(Node *DNodes, Material *DMat);
	Plane42(FILE *file, long Num, Node *DNodes);
	void AssembleMatrix(double **Global, int MaxDOF, Node *DNodes);
};

void Plane42::CalculateMatrix(Node *DNodes, Material *DMat)
{
	double Coordonate[2], Jacobian[2][2], JacDet, Temp;
	double EMatrixConst;
	double RValue[4] = {-.57735, .57735, .57735, -.57735}; // for 0,1,2 and 3
	double SValue[4] = {-.57735, -.57735, .57735, .57735}; // for 0,1,2 and 3
	double NShFu[4];									   // shape function
	double NDerR[4];									   // derivative relative to R
	double NDerS[4];									   // derivative relative to S
	int inc, inc2, inc3, inc4, inc5;
	double DMatrixConst;
	double **EMatrix = new double *[3];
	double **BMatrix = new double *[3];	 //matrix B
	double **BtMatrix = new double *[8]; // Transpose Matrix
	double **TempMatrix = new double *[3];
	double temp02[8][8];

	Sort(DNodes); // Position node in the right order
	for (inc = 0; inc < 8; inc++)
		for (inc2 = 0; inc2 < 8; inc2++)
			temp02[inc][inc2] = 0;

	for (inc = 0; inc < 3; inc++) // 3*8
	{
		BMatrix[inc] = new double[8];
		TempMatrix[inc] = new double[8];
	}
	for (inc = 0; inc < 3; inc++) // 3*3
		EMatrix[inc] = new double[3];
	for (inc = 0; inc < 8; inc++) // 8*3
		BtMatrix[inc] = new double[3];

	//quadrature order is 2
	Coordonate[0] = 1;
	Coordonate[1] = -1;
	//Calculatoin de EMatrix
	EMatrixConst = DMat[(Mat - 1)].ex / (1 - pow(DMat[(Mat - 1)].poisson, 2));
	for (inc = 0; inc < 3; inc++)
		for (inc2 = 0; inc2 < 3; inc2++)
			EMatrix[inc][inc2] = 0;
	EMatrix[0][0] = EMatrixConst;
	EMatrix[0][1] = (DMat[(Mat - 1)].poisson * EMatrixConst);
	EMatrix[1][0] = (DMat[(Mat - 1)].poisson * EMatrixConst);
	EMatrix[1][1] = EMatrixConst;
	EMatrix[2][2] = (0.5 * (1 - DMat[(Mat - 1)].poisson) * EMatrixConst);

	StiffMat = new double *[8];
	for (inc = 0; inc < 8; inc++) // 8*8
		StiffMat[inc] = new double[8];

	for (inc = 0; inc < 4; inc++) // where 2 is the order of the integration
	{
		// Definition of Ni, and their derivatives
		// faux

		NShFu[0] = 0.25 * (1 - RValue[inc]) * (1 - SValue[inc]);
		NShFu[1] = 0.25 * (1 + RValue[inc]) * (1 - SValue[inc]);
		NShFu[2] = 0.25 * (1 + RValue[inc]) * (1 + SValue[inc]);
		NShFu[3] = 0.25 * (1 - RValue[inc]) * (1 + SValue[inc]);

		NDerR[0] = -0.25 * (1 - SValue[inc]);
		NDerR[1] = 0.25 * (1 - SValue[inc]);
		NDerR[2] = 0.25 * (1 + SValue[inc]);
		NDerR[3] = -0.25 * (1 + SValue[inc]);

		NDerS[0] = -0.25 * (1 - RValue[inc]);
		NDerS[1] = -0.25 * (1 + RValue[inc]);
		NDerS[2] = 0.25 * (1 + RValue[inc]);
		NDerS[3] = 0.25 * (1 - RValue[inc]);

		// initialisation of the Jacobian Matrix
		for (inc4 = 0; inc4 < 2; inc4++)
			for (inc5 = 0; inc5 < 2; inc5++)
				Jacobian[inc4][inc5] = 0;

		for (inc3 = 0; inc3 < 4; inc3++)
		{
			Jacobian[0][0] += NDerR[inc3] * DNodes[(NbLink[inc3] - 1)].nx;
			Jacobian[0][1] += NDerS[inc3] * DNodes[(NbLink[inc3] - 1)].nx;
			Jacobian[1][0] += NDerR[inc3] * DNodes[(NbLink[inc3] - 1)].ny;
			Jacobian[1][1] += NDerS[inc3] * DNodes[(NbLink[inc3] - 1)].ny;
		}

		JacDet = Jacobian[0][0] * Jacobian[1][1] - Jacobian[1][0] * Jacobian[0][1];

		// Replacement Of the jacobian matrix by its inverse
		Temp = Jacobian[0][0] / JacDet;
		Jacobian[0][0] = Jacobian[1][1] / JacDet;
		Jacobian[1][1] = Temp;
		Temp = Jacobian[1][0] / (-JacDet);
		Jacobian[1][0] = Jacobian[0][1] / (-JacDet);
		Jacobian[0][1] = Temp;
		// Calculation of BMatrix[][]
		for (inc3 = 0; inc3 < 3; inc3++)
			for (inc4 = 0; inc4 < 8; inc4++)
				BMatrix[inc3][inc4] = 0;

		for (inc3 = 0; inc3 < 4; inc3++)
		{
			BMatrix[0][inc3 * 2] = Jacobian[0][0] * NDerR[inc3] + Jacobian[1][0] * NDerS[inc3];
			BMatrix[1][inc3 * 2 + 1] = Jacobian[0][1] * NDerR[inc3] + Jacobian[1][1] * NDerS[inc3];
			BMatrix[2][inc3 * 2] = BMatrix[1][inc3 * 2 + 1];
			BMatrix[2][inc3 * 2 + 1] = BMatrix[0][inc3 * 2];
		}
		// Calculation of the Transpose of BMatrix[][]
		matrix_transpose(BtMatrix, BMatrix, 8, 3);
		for (inc3 = 0; inc3 < 3; inc3++)
			for (inc4 = 0; inc4 < 8; inc4++)
				BMatrix[inc3][inc4] = BMatrix[inc3][inc4] * JacDet * 1; // 1 for the weight

		matrix_multiply(EMatrix, BMatrix, TempMatrix, 3, 3, 8);
		matrix_multiply(BtMatrix, TempMatrix, StiffMat, 8, 3, 8);
		for (inc3 = 0; inc3 < 8; inc3++)
			for (inc4 = 0; inc4 < 8; inc4++)
			{
				temp02[inc3][inc4] += StiffMat[inc3][inc4];
				inc5 = 0;
			}
	}
	for (inc3 = 0; inc3 < 8; inc3++)
		for (inc4 = 0; inc4 < 8; inc4++)
			StiffMat[inc3][inc4] = temp02[inc3][inc4];
}

Plane42::Plane42(FILE *file, long Num, Node *DNodes) // function to initialize data during the file reading
{

	char Command[100];

	number = Num;
	NextWord(file, &Command[0]);
	sscanf(Command, "%ld", &Mat); // read the material number
								  //	NextWord(file,&Command[0]);
								  //	sscanf(Command,"%le", &Thickness);		// read the thickness
	for (int inc = 0; inc < 4; inc++)
	{
		NextWord(file, &Command[0]);
		sscanf(Command, "%ld", &NbLink[inc]); // Stiffness matrix
	}

	// Sort the number of Node
	//SortNumber (&NbLink[0],6);
	DefineNodeDOF(DNodes); // To read the boundary conditions
						   // as forces or displacements
}

void Plane42::AssembleMatrix(double **Global, int MaxDOF, Node *DNodes)
{
	// Function which assemble the Matrix !
	for (int inc = 0; inc < 4; inc++)				 //step for column (3 Nodes 2x2)
		for (int inc2 = 0; inc2 < 4; inc2++)		 // step for row (3 Nodes 2x2)
			for (int inc3 = 0; inc3 < 2; inc3++)	 // step for line (row) (2 DOF)
				for (int inc4 = 0; inc4 < 2; inc4++) // step for line (column) (2 DOF)
				{
					Global[(DNodes[(NbLink[inc] - 1)].CumSum + inc3)][(DNodes[NbLink[inc2] - 1].CumSum + inc4)] += StiffMat[(inc * 2) + inc3][(inc2 * 2) + inc4];
					inc4 = inc4;
				}
}
void Plane42::DefineNodeDOF(Node *DNodes)
{
	for (int inc = 0; inc < 4; inc++)
		if (!(DNodes[(NbLink[inc] - 1)].DOF > 2))
			DNodes[(NbLink[inc] - 1)].DOF = 2;
}

void Plane42::Sort(Node *DNodes)
{
	long Temp;
	double vectorAB[2], vectorBC[2], vectorBD[2], vectorDA[2], vectorCD[2];
	/*	A equiv 0
	B equiv 1
	C equiv 2
	D equiv 3 */
	// vector definition
	vectorAB[0] = DNodes[(NbLink[0] - 1)].nx - DNodes[(NbLink[1] - 1)].nx;
	vectorAB[1] = DNodes[(NbLink[0] - 1)].ny - DNodes[(NbLink[1] - 1)].ny;
	vectorBC[0] = DNodes[(NbLink[1] - 1)].nx - DNodes[(NbLink[2] - 1)].nx;
	vectorBC[1] = DNodes[(NbLink[1] - 1)].ny - DNodes[(NbLink[2] - 1)].ny;
	//Is C at the right place
	// Exchange C and B
	if ((vectorAB[0] * vectorBC[1] - vectorBC[0] * vectorAB[1]) < 0) // if =0 if ABC is straight
	{
		Temp = NbLink[1];
		NbLink[1] = NbLink[2];
		NbLink[2] = Temp;
	}

	//Is D at the right place
	//CD^DA

	vectorAB[0] = DNodes[(NbLink[0] - 1)].nx - DNodes[(NbLink[1] - 1)].nx;
	vectorAB[1] = DNodes[(NbLink[0] - 1)].ny - DNodes[(NbLink[1] - 1)].ny;
	vectorBD[0] = DNodes[(NbLink[1] - 1)].nx - DNodes[(NbLink[3] - 1)].nx;
	vectorBD[1] = DNodes[(NbLink[1] - 1)].ny - DNodes[(NbLink[3] - 1)].ny;
	vectorCD[0] = DNodes[(NbLink[2] - 1)].nx - DNodes[(NbLink[3] - 1)].nx;
	vectorCD[1] = DNodes[(NbLink[2] - 1)].ny - DNodes[(NbLink[3] - 1)].ny;
	vectorDA[0] = DNodes[(NbLink[3] - 1)].nx - DNodes[(NbLink[0] - 1)].nx;
	vectorDA[1] = DNodes[(NbLink[3] - 1)].ny - DNodes[(NbLink[0] - 1)].ny;

	if ((vectorCD[0] * vectorDA[1] - vectorDA[0] * vectorCD[1]) < 0) // if =0 if ABC is straight
	{
		//AB^BD
		if ((vectorAB[0] * vectorBD[1] - vectorBD[0] * vectorAB[1]) < 0) // if =0 if ABC is straight
		{
			Temp = NbLink[1];	   // temp takes B value
			NbLink[1] = NbLink[3]; // B takes D value
			NbLink[3] = NbLink[2]; // D takes C value
			NbLink[2] = Temp;	   // C takes B value
		}
		else
		{
			//exchange C and D
			Temp = NbLink[2];	   // temp takes C value
			NbLink[2] = NbLink[3]; // C takes D value
			NbLink[3] = Temp;	   // D takes C value
		}
	}
}