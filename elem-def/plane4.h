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
    long NbLink[4];    // 4 connections
    long Mat;          // Material
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
    double NShFu[4];                                       // shape function
    double NDerR[4];                                       // derivative relative to R
    double NDerS[4];                                       // derivative relative to S
    int inc, inc2, inc3, inc4, inc5;
    double DMatrixConst;
    double **EMatrix = new double *[3];
    double **BMatrix = new double *[3];  //matrix B
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
    for (int inc = 0; inc < 4; inc++)                //step for column (3 Nodes 2x2)
        for (int inc2 = 0; inc2 < 4; inc2++)         // step for row (3 Nodes 2x2)
            for (int inc3 = 0; inc3 < 2; inc3++)     // step for line (row) (2 DOF)
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
            Temp = NbLink[1];      // temp takes B value
            NbLink[1] = NbLink[3]; // B takes D value
            NbLink[3] = NbLink[2]; // D takes C value
            NbLink[2] = Temp;      // C takes B value
        }
        else
        {
            //exchange C and D
            Temp = NbLink[2];      // temp takes C value
            NbLink[2] = NbLink[3]; // C takes D value
            NbLink[3] = Temp;      // D takes C value
        }
    }
}