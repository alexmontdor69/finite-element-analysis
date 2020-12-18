/** 
 * link 1
 * 
*/

class Link1 // Ansys code for 2D-spar
{
private:
    long id; // What's the element's name ?

    int nb_nodes = 2; // the element is defined by 'ends' node
    long node_ids[2]; // define which node the element he is connected to

    long Mat; // Material
    long R;   // Real (for cross section area)

    double Area;           //equivalent of R ? Check on the calculation
    double Angle;          // Angle between the x-axis and the bar
    double Length;         // lenght of the bar
    double StiffMat[4][4]; // Stiffness Matrix dim = 2
public:
    Link1(long element_id, const std::string &values, Node *DNodes);
    void CreateStiffMatrix(Node *DNodes, Material *DMat);
    void Display(int m);
    void AssembleMatrix(double **main_matrix, int DOF_max, Node *DNodes);
    void DefineNodeDOF(Node *DNodes);
};

/** 
 * Initialization of the Link
 * 
 * element id  : id of the element in the model
 * values : series of values obtained in the model file
 * Node : Array of existing node
*/
Link1::Link1(long element_id, const std::string &values, Node *DNodes)
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

    for (int index = 0; index < nb_nodes; index++)
    {
        std::string node_id = data.substr(0, delimiterPos);
        data = data.substr(delimiterPos + 1);
        node_ids[index] = std::stol(node_id, nullptr, 10);
        delimiterPos = data.find(",");
    }

    SortNumber(&node_ids[0], nb_nodes);
    DefineNodeDOF(DNodes);
}
void Link1::CreateStiffMatrix(Node *DNodes, Material *DMat)
{

    double CosA, SinA;

    Length = sqrt(pow(((DNodes[(node_ids[0] - 1)].nx) - (DNodes[(node_ids[1] - 1)].nx)), 2) + pow((DNodes[(node_ids[0] - 1)].ny) - (DNodes[(node_ids[1] - 1)].ny), 2));
    CosA = ((DNodes[(node_ids[0] - 1)].nx) - (DNodes[(node_ids[1] - 1)].nx)) / Length;
    SinA = ((DNodes[(node_ids[0] - 1)].ny) - (DNodes[(node_ids[1] - 1)].ny)) / Length;
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
void Link1::Display(int m)
{
    std::cout << "\n Element : "
              << m + 1
              << "\n";
    for (int inc = 0; inc < 4; inc++)
    {

        for (int in = 0; in < 4; in++)
            std::cout << StiffMat[inc][in] << "      ";
        std::cout << "\n";
    }
}

void Link1::AssembleMatrix(double **main_matrix, int DOF_max, Node *DNodes)
{
    // Function which assemble the Matrix !

    for (int inode1 = 0; inode1 < 2; inode1++)          //step for 2 Nodes
        for (int inode2 = 0; inode2 < 2; inode2++)      // step for 2 Nodes
            for (int iDOF1 = 0; iDOF1 < 2; iDOF1++)     // for all 2 DOF
                for (int iDOF2 = 0; iDOF2 < 2; iDOF2++) // for all 2 DOF
                    main_matrix[(DNodes[(node_ids[inode1] - 1)].absolute_DOF_addr + iDOF1)][(DNodes[node_ids[inode2] - 1].absolute_DOF_addr + iDOF2)] += StiffMat[(inode1 * 2) + iDOF1][(inode2 * 2) + iDOF2];
}

void Link1::DefineNodeDOF(Node *DNodes)
{
    for (int index = 0; index < nb_nodes; index++)
        if (!(DNodes[(node_ids[index] - 1)].DOF > 2))
            DNodes[(node_ids[index] - 1)].DOF = 2;
}