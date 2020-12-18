/** 
 * Plane 4
 * 4 Nodes
 * 2 dimension
 * 
*/

class Plane4
{
private:
    long id;
    int nb_nodes = 4;  // the element is defined by 'ends' node
    long node_ids[4];  // define which node the element he is connected to
    long Mat;          // Material
    double Thickness;  //Thickness of the element
    double **StiffMat; //
    void DefineNodeDOF(Node *DNodes);
    void Sort(Node *DNodes);

public:
    Plane4(long element_id, const std::string &values, Node *DNodes);
    void CreateStiffMatrix(Node *DNodes, Material *DMat);
    void AssembleMatrix(double **main_matrix, int DOF_max, Node *DNodes);
};

Plane4::Plane4(long element_id, const std::string &values, Node *DNodes) // function to initialize data during the file reading
{
    id = element_id;

    int delimiterPos = values.find(",");
    std::string mat = values.substr(0, delimiterPos);
    Mat = std::stol(mat, nullptr, 10);
    delimiterPos = values.find(",", delimiterPos + 1);

    std::string z = values.substr(delimiterPos + 1);
    Thickness = std::atol(z.c_str());
    delimiterPos = values.find(",", delimiterPos + 1);

    for (int index = 0; index < nb_nodes; index++)
    {
        std::string node_id = values.substr(delimiterPos + 1);
        node_ids[index] = std::stol(node_id, nullptr, 10);
        delimiterPos = values.find(",", delimiterPos + 1);
    }

    StiffMat = new double *[6];
    for (int inc = 0; inc < 6; inc++) // 6*6
        StiffMat[inc] = new double[6];

    DefineNodeDOF(DNodes); // To read the boundary conditions
                           // as forces or displacements
}

void Plane4::CreateStiffMatrix(Node *DNodes, Material *DMat)
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
            Jacobian[0][0] += NDerR[inc3] * DNodes[(node_ids[inc3] - 1)].nx;
            Jacobian[0][1] += NDerS[inc3] * DNodes[(node_ids[inc3] - 1)].nx;
            Jacobian[1][0] += NDerR[inc3] * DNodes[(node_ids[inc3] - 1)].ny;
            Jacobian[1][1] += NDerS[inc3] * DNodes[(node_ids[inc3] - 1)].ny;
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

void Plane4::AssembleMatrix(double **main_matrix, int DOF_max, Node *DNodes)
{

    for (int inode1 = 0; inode1 < 4; inode1++)          //step for column (matrix 2x2)
        for (int inode2 = 0; inode2 < 4; inode2++)      // step for row (matrix 2x2)
            for (int iDOF1 = 0; iDOF1 < 2; iDOF1++)     // for all 2 DOF
                for (int iDOF2 = 0; iDOF2 < 2; iDOF2++) // for all 2 DOF
                    main_matrix[(DNodes[(node_ids[inode1] - 1)].absolute_DOF_addr + iDOF1)][(DNodes[node_ids[inode2] - 1].absolute_DOF_addr + iDOF2)] += StiffMat[(inode1 * 2) + iDOF1][(inode2 * 2) + iDOF2];
}

void Plane4::DefineNodeDOF(Node *DNodes)
{
    for (int inc = 0; inc < 4; inc++)
        if (!(DNodes[(node_ids[inc] - 1)].DOF > 2))
            DNodes[(node_ids[inc] - 1)].DOF = 2;
}

void Plane4::Sort(Node *DNodes)
{
    long Temp;
    double vectorAB[2], vectorBC[2], vectorBD[2], vectorDA[2], vectorCD[2];
    /*	A equiv 0
	B equiv 1
	C equiv 2
	D equiv 3 */
    // vector definition
    vectorAB[0] = DNodes[(node_ids[0] - 1)].nx - DNodes[(node_ids[1] - 1)].nx;
    vectorAB[1] = DNodes[(node_ids[0] - 1)].ny - DNodes[(node_ids[1] - 1)].ny;
    vectorBC[0] = DNodes[(node_ids[1] - 1)].nx - DNodes[(node_ids[2] - 1)].nx;
    vectorBC[1] = DNodes[(node_ids[1] - 1)].ny - DNodes[(node_ids[2] - 1)].ny;
    //Is C at the right place
    // Exchange C and B
    if ((vectorAB[0] * vectorBC[1] - vectorBC[0] * vectorAB[1]) < 0) // if =0 if ABC is straight
    {
        Temp = node_ids[1];
        node_ids[1] = node_ids[2];
        node_ids[2] = Temp;
    }

    //Is D at the right place
    //CD^DA

    vectorAB[0] = DNodes[(node_ids[0] - 1)].nx - DNodes[(node_ids[1] - 1)].nx;
    vectorAB[1] = DNodes[(node_ids[0] - 1)].ny - DNodes[(node_ids[1] - 1)].ny;
    vectorBD[0] = DNodes[(node_ids[1] - 1)].nx - DNodes[(node_ids[3] - 1)].nx;
    vectorBD[1] = DNodes[(node_ids[1] - 1)].ny - DNodes[(node_ids[3] - 1)].ny;
    vectorCD[0] = DNodes[(node_ids[2] - 1)].nx - DNodes[(node_ids[3] - 1)].nx;
    vectorCD[1] = DNodes[(node_ids[2] - 1)].ny - DNodes[(node_ids[3] - 1)].ny;
    vectorDA[0] = DNodes[(node_ids[3] - 1)].nx - DNodes[(node_ids[0] - 1)].nx;
    vectorDA[1] = DNodes[(node_ids[3] - 1)].ny - DNodes[(node_ids[0] - 1)].ny;

    if ((vectorCD[0] * vectorDA[1] - vectorDA[0] * vectorCD[1]) < 0) // if =0 if ABC is straight
    {
        //AB^BD
        if ((vectorAB[0] * vectorBD[1] - vectorBD[0] * vectorAB[1]) < 0) // if =0 if ABC is straight
        {
            Temp = node_ids[1];        // temp takes B value
            node_ids[1] = node_ids[3]; // B takes D value
            node_ids[3] = node_ids[2]; // D takes C value
            node_ids[2] = Temp;        // C takes B value
        }
        else
        {
            //exchange C and D
            Temp = node_ids[2];        // temp takes C value
            node_ids[2] = node_ids[3]; // C takes D value
            node_ids[3] = Temp;        // D takes C value
        }
    }
}