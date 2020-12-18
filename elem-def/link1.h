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

    double Angle;          // Angle between the x-axis and the bar
    double Length;         // lenght of the bar
    double StiffMat[4][4]; // Stiffness Matrix dim = 2
public:
    Link1(long element_id, const std::string &values, Node *DNodes);
    void CalculateLength(Node *DNodes, Material *DMat);
    void Display(int m);
    void AssembleMatrix(double **Global, int MaxDOF, Node *DNodes);
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

    int delimiterPos = values.find(",");
    std::string mat = values.substr(0, delimiterPos);
    Mat = std::stol(mat, nullptr, 10);
    delimiterPos = values.find(",", delimiterPos + 1);

    std::string cross_area = values.substr(delimiterPos + 1);
    R = std::stol(cross_area, nullptr, 10);
    delimiterPos = values.find(",", delimiterPos + 1);

    for (int index = 0; index < nb_nodes; index++)
    {
        std::string node_id = values.substr(delimiterPos + 1);
        node_ids[index] = std::stol(node_id, nullptr, 10);
        delimiterPos = values.find(",", delimiterPos + 1);
    }

    SortNumber(&node_ids[0], nb_nodes);
    DefineNodeDOF(DNodes);
}
void Link1::CalculateLength(Node *DNodes, Material *DMat)
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

void Link1::AssembleMatrix(double **Global, int MaxDOF, Node *DNodes)
{
    // Function which assemble the Matrix !

    for (int inc = 0; inc < 2; inc++)            //step for column (matrix 2x2)
        for (int inc2 = 0; inc2 < 2; inc2++)     // step for row (matrix 2x2)
            for (int inc3 = 0; inc3 < 2; inc3++) // step for line
                for (int inc4 = 0; inc4 < 2; inc4++)
                    Global[(DNodes[(node_ids[inc] - 1)].CumSum + inc3)][(DNodes[node_ids[inc2] - 1].CumSum + inc4)] += StiffMat[(inc * 2) + inc3][(inc2 * 2) + inc4];
}

void Link1::DefineNodeDOF(Node *DNodes)
{
    for (int index = 0; index < nb_nodes; index++)
        if (!(DNodes[(node_ids[index] - 1)].DOF > 2))
            DNodes[(node_ids[index] - 1)].DOF = 2;
}