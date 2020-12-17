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
    long NbLink[3];    // 3 connections
    long Mat;          // Material
    long R;            // Real
    double EL;         // EI/L
    double Area;       // Cross area
    double Izz;        // inertia moment about z-axis
    double Angle;      // Angle between the x-axis and the bar
    double Length;     // lenght of the bar
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
    // NextWord(file, &Command[0]);
    sscanf(Command, "%ld", &Mat); // read the material number
    // NextWord(file, &Command[0]);
    sscanf(Command, "%le", &Area); // read the area
    // NextWord(file, &Command[0]);
    sscanf(Command, "%le", &Izz); // read the moment of inertia
    // NextWord(file, &Command[0]);
    sscanf(Command, "%le", &Thickness); // read the thickness
    for (int inc = 0; inc < 3; inc++)
    {
        // NextWord(file, &Command[0]);
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
    std::cout << "\n Element : "
              << m + 1
              << "\n";
    for (int inc = 0; inc < 6; inc++)
    {
        for (int in = 0; in < 6; in++)
            std::cout << StiffMat[inc][in] << "      ";
        std::cout << "\n";
    }
}

void Plane2::AssembleMatrix(double **Global, int MaxDOF, Node *DNodes)
{
    // Function which assemble the Matrix !
    for (int inc = 0; inc < 3; inc++)                //step for column (3 Nodes 2x2)
        for (int inc2 = 0; inc2 < 3; inc2++)         // step for row (3 Nodes 2x2)
            for (int inc3 = 0; inc3 < 2; inc3++)     // step for line (row)
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