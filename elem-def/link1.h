//-----------LINK1---------------------------------------------------------
class Link1 // Ansys code for 2D-spar
{
private:
    long number;           // What's the element's name ?
    long NbLink[2];        // define how many connection  And see alla bout the array!!
    long Mat;              // Material
    long R;                // Real
    double Angle;          // Angle between the x-axis and the bar
    double Length;         // lenght of the bar
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

    for (int inc = 0; inc < 2; inc++)            //step for column (matrix 2x2)
        for (int inc2 = 0; inc2 < 2; inc2++)     // step for row (matrix 2x2)
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