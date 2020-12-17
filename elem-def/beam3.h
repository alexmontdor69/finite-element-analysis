//--------Beam3-----------------------------------------------------------------------
class Beam3
{
private:
    long number;           // What's the element's name ?
    long NbLink[2];        // define how many connection  And see alla bout the array!!
    long Mat;              // Material
    long R;                // Real
    double EL;             // EI/L
    double Area;           // Cross area
    double Izz;            // inertia moment about z-axis
    double Angle;          // Angle between the x-axis and the bar
    double Length;         // lenght of the bar
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
    Z = Area * EL;    // Verify this is the area use instead of A
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
    /*NextWord(file, &Command[0]);
    sscanf(Command, "%ld", &Mat); // read the material number
    NextWord(file, &Command[0]);
    sscanf(Command, "%le", &Area); // read the area
    NextWord(file, &Command[0]);
    sscanf(Command, "%le", &Izz); // read the moment of inertia
    NextWord(file, &Command[0]);
    sscanf(Command, "%ld", &NbLink[0]); // Stiffness matrix
    NextWord(file, &Command[0]);
    sscanf(Command, "%ld", &NbLink[1]);
    SortNumber(&NbLink[0], 2);*/

    DefineNodeDOF(DNodes);
}

void Beam3::Display(int m)
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

void Beam3::AssembleMatrix(double **Global, int MaxDOF, Node *DNodes)
{
    // Function which assemble the Matrix !
    for (int inc = 0; inc < 2; inc++)                //step for column (4 matrices 3x3)
        for (int inc2 = 0; inc2 < 2; inc2++)         // step for row (4 matrices 3x3)
            for (int inc3 = 0; inc3 < 3; inc3++)     // step for line (row)
                for (int inc4 = 0; inc4 < 3; inc4++) // step for line (column)
                    Global[(DNodes[(NbLink[inc] - 1)].CumSum + inc3)][(DNodes[NbLink[inc2] - 1].CumSum + inc4)] += StiffMat[(inc * 3) + inc3][(inc2 * 3) + inc4];
}

void Beam3::DefineNodeDOF(Node *DNodes)
{
    for (int inc = 0; inc < 2; inc++)
        if (!(DNodes[(NbLink[inc] - 1)].DOF > 3))
            DNodes[(NbLink[inc] - 1)].DOF = 3;
}