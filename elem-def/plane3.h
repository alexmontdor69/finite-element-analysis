/** 
 * Plane 3
 * 3 Nodes
 * 2 Dimensions
 * 
*/

class Plane3
{
private:
    long id;
    int nb_nodes = 3;  // the element is defined by 'ends' node
    long node_ids[3];  // define which node the element he is connected to
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
    Plane3(long element_id, const std::string &values, Node *DNodes);
    void CalculateLength(Node *DNodes, Material *DMat);
    void Display(int m);
    void AssembleMatrix(double **Global, int MaxDOF, Node *DNodes);
};
Plane3::Plane3(long element_id, const std::string &values, Node *DNodes) // function to initialize data during the file reading
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

    /*  No Thickness in the data file
   std::string z = data.substr(0, delimiterPos + 1);
    data = data.substr(delimiterPos + 1);
    Thickness = std::atol(z.c_str());
    delimiterPos = data.find(","); */

    for (int index = 0; index < nb_nodes; index++)
    {
        std::string node_id = data.substr(0, delimiterPos);
        data = data.substr(delimiterPos + 1);
        node_ids[index] = std::stol(node_id, nullptr, 10);
        delimiterPos = data.find(",");
    }

    StiffMat = new double *[6];
    for (int inc = 0; inc < 6; inc++) // 6*6
        StiffMat[inc] = new double[6];

    SortNumber(&node_ids[0], nb_nodes);
    DefineNodeDOF(DNodes);
}

void Plane3::CalculateLength(Node *DNodes, Material *DMat)
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

    twodet = DNodes[(node_ids[1] - 1)].nx * DNodes[(node_ids[2] - 1)].ny - DNodes[(node_ids[2] - 1)].nx * DNodes[(node_ids[1] - 1)].ny - DNodes[(node_ids[0] - 1)].nx * DNodes[(node_ids[2] - 1)].ny + DNodes[(node_ids[0] - 1)].nx * DNodes[(node_ids[1] - 1)].ny + DNodes[(node_ids[2] - 1)].nx * DNodes[(node_ids[0] - 1)].ny - DNodes[(node_ids[1] - 1)].nx * DNodes[(node_ids[0] - 1)].ny;

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
                        StiffMat[inc * 2 + inc3][inc2 * 2 + inc4] = Thickness / (2 * twodet) * (DMatrix[0][0] * ((DNodes[(node_ids[indxx2] - 1)].ny) - (DNodes[(node_ids[indxx1] - 1)].ny)) * ((DNodes[(node_ids[indxy2] - 1)].ny) - (DNodes[(node_ids[indxy1] - 1)].ny)) + DMatrix[2][2] * ((DNodes[(node_ids[indxy1] - 1)].nx) - (DNodes[(node_ids[indxy2] - 1)].nx)) * ((DNodes[(node_ids[indxx1] - 1)].nx) - (DNodes[(node_ids[indxx2] - 1)].nx)));

                    if (inc3 == 1 && inc4 == 0)
                        StiffMat[inc * 2 + inc3][inc2 * 2 + inc4] = Thickness / (2 * twodet) * (DMatrix[0][1] * ((DNodes[(node_ids[indxy2] - 1)].ny) - (DNodes[(node_ids[indxy1] - 1)].ny)) * ((DNodes[(node_ids[indxx1] - 1)].nx) - (DNodes[(node_ids[indxx2] - 1)].nx)) + DMatrix[2][2] * ((DNodes[(node_ids[indxy1] - 1)].nx) - (DNodes[(node_ids[indxy2] - 1)].nx)) * ((DNodes[(node_ids[indxx2] - 1)].ny) - (DNodes[(node_ids[indxx1] - 1)].ny)));

                    if (inc3 == 0 && inc4 == 1)
                        StiffMat[inc * 2 + inc3][inc2 * 2 + inc4] = Thickness / (2 * twodet) * (DMatrix[0][1] * ((DNodes[(node_ids[indxx1] - 1)].nx) - (DNodes[(node_ids[indxx2] - 1)].nx)) * ((DNodes[(node_ids[indxy2] - 1)].ny) - (DNodes[(node_ids[indxy1] - 1)].ny)) + DMatrix[2][2] * ((DNodes[(node_ids[indxx1] - 1)].nx) - (DNodes[(node_ids[indxx2] - 1)].nx)) * ((DNodes[(node_ids[indxy2] - 1)].ny) - (DNodes[(node_ids[indxy1] - 1)].ny)));

                    if (inc3 == 1 && inc4 == 1)
                        StiffMat[inc * 2 + inc3][inc2 * 2 + inc4] = Thickness / (2 * twodet) * (DMatrix[1][1] * ((DNodes[(node_ids[indxx1] - 1)].nx) - (DNodes[(node_ids[indxx2] - 1)].nx)) * ((DNodes[(node_ids[indxy1] - 1)].nx) - (DNodes[(node_ids[indxy2] - 1)].nx)) + DMatrix[2][2] * ((DNodes[(node_ids[indxy2] - 1)].ny) - (DNodes[(node_ids[indxy1] - 1)].ny)) * ((DNodes[(node_ids[indxx2] - 1)].ny) - (DNodes[(node_ids[indxx1] - 1)].ny)));
                }
        }
}

void Plane3::Display(int m)
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

void Plane3::AssembleMatrix(double **Global, int MaxDOF, Node *DNodes)
{
    // Function which assemble the Matrix !
    for (int inc = 0; inc < 3; inc++)                //step for column (3 Nodes 2x2)
        for (int inc2 = 0; inc2 < 3; inc2++)         // step for row (3 Nodes 2x2)
            for (int inc3 = 0; inc3 < 2; inc3++)     // step for line (row)
                for (int inc4 = 0; inc4 < 2; inc4++) // step for line (column)
                {
                    Global[(DNodes[(node_ids[inc] - 1)].absolute_DOF_addr + inc3)][(DNodes[node_ids[inc2] - 1].absolute_DOF_addr + inc4)] += StiffMat[(inc * 2) + inc3][(inc2 * 2) + inc4];
                    inc4 = inc4;
                }
}
void Plane3::DefineNodeDOF(Node *DNodes)
{
    for (int inc = 0; inc < 3; inc++)
        if (!(DNodes[(node_ids[inc] - 1)].DOF > 2))
            DNodes[(node_ids[inc] - 1)].DOF = 2;
}