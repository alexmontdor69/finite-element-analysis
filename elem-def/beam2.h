/** 
 * Beam2
 * 
*/

class Beam2
{
private:
    long id;
    int nb_nodes = 2;      // What's the element's name ?
    long node_ids[2];      // define which node the element he is connected to
    long Mat;              // Material
    long R;                // Real
    double EL;             // EI/L
    double Area;           // Cross area
    double Izz;            // inertia moment about z-axis
    double Angle;          // Angle between the x-axis and the bar
    double Length;         // lenght of the bar
    double StiffMat[6][6]; // Stiffness Matrix dim = 3 (3*2Nodes)
public:
    Beam2(long element_id, const std::string &values, Node *DNodes);
    void CreateStiffMatrix(Node *DNodes, Material *DMat);
    void Display(int m);
    void AssembleMatrix(double **main_matrix, int DOF_max, Node *DNodes);
    void DefineNodeDOF(Node *DNodes);
};
Beam2::Beam2(long element_id, const std::string &values, Node *DNodes) // function to initialize data during the file reading
{
    id = element_id;
    std::string data = "";

    int delimiterPos = values.find(",");
    std::string mat = values.substr(0, delimiterPos);
    data = values.substr(delimiterPos + 1);
    Mat = std::stol(mat, nullptr, 10);
    delimiterPos = data.find(",");

    std::string cross_area = data.substr(0, delimiterPos);
    data = data.substr(delimiterPos + 1);
    Area = std::atol(cross_area.c_str());
    delimiterPos = data.find(",");

    std::string moment_inertia_z = data.substr(0, delimiterPos);
    data = data.substr(delimiterPos + 1);
    Izz = std::atol(moment_inertia_z.c_str());
    delimiterPos = data.find(",");

    for (int index = 0; index < nb_nodes; index++)
    {
        std::string node_id = data.substr(0, delimiterPos);
        data = data.substr(delimiterPos + 1);
        node_ids[index] = std::stol(node_id, nullptr, 10);
        delimiterPos = data.find(",");
    }

    DefineNodeDOF(DNodes);
}
void Beam2::CreateStiffMatrix(Node *DNodes, Material *DMat)
{
    double CosA, SinA;
    double A, B, K, M, F, G, P, H, Q, Z; // Temp Variable to make calculation
    Length = sqrt(pow(((DNodes[(node_ids[0] - 1)].nx) - (DNodes[(node_ids[1] - 1)].nx)), 2) + pow((DNodes[(node_ids[0] - 1)].ny) - (DNodes[(node_ids[1] - 1)].ny), 2));
    CosA = ((DNodes[(node_ids[0] - 1)].nx) - (DNodes[(node_ids[1] - 1)].nx)) / Length;
    SinA = ((DNodes[(node_ids[0] - 1)].ny) - (DNodes[(node_ids[1] - 1)].ny)) / Length;
    Angle = acos(CosA) * 180 / PI;

    EL = DMat[(Mat - 1)].ex / Length; // E/L ratio Young modulus over the length

    // Previous Calculations : Check them 'cos I'm not sure of the result...
    A = 4 * EL * Izz;
    Z = Area * EL;
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

void Beam2::Display(int m)
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

void Beam2::AssembleMatrix(double **main_matrix, int DOF_max, Node *DNodes)
{

    for (int inode1 = 0; inode1 < 2; inode1++)          //step for 2 Nodes
        for (int inode2 = 0; inode2 < 2; inode2++)      //step for 2 Nodes
            for (int iDOF1 = 0; iDOF1 < 3; iDOF1++)     // for all 3 DOF
                for (int iDOF2 = 0; iDOF2 < 3; iDOF2++) // for all 3 DOF
                    main_matrix[(DNodes[(node_ids[inode1] - 1)].absolute_DOF_addr + iDOF1)][(DNodes[node_ids[inode2] - 1].absolute_DOF_addr + iDOF2)] += StiffMat[(inode1 * 3) + iDOF1][(inode2 * 3) + iDOF2];
}

void Beam2::DefineNodeDOF(Node *DNodes)
{
    for (int inc = 0; inc < 2; inc++)
        if (!(DNodes[(node_ids[inc] - 1)].DOF > 3))
            DNodes[(node_ids[inc] - 1)].DOF = 3;
}