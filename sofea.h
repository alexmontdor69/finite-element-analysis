// sofcla.h
// headers to define all the class use in sofea.cpp

int EndOfFile = 0;

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

/** 
 * 
 * 
*/

class Node
{
public:			   // Caution No private because couldn't use some other function!!
				   // So to change
	long number;   // What's the node's name ?
	double nx, ny; // Nodes Coordonates
	long DOF;	   // Number of DOF for this node
	long CumSum;   // "Adress" in the vectors.
	void Init(const char *node_id, const std::string &pos_x, const std::string &pos_y);
};

void Node::Init(const char *node_id, const std::string &pos_x, const std::string &pos_y) // function to initialize data during the file reading
{
	number = std::atol(node_id);
	nx = std::stol(pos_x, nullptr, 10);
	ny = std::stol(pos_y, nullptr, 10);
	DOF = 0;
	CumSum = 0;
}

//----------Material----------------------------------------------------------
class Material
{
public:
	long number;		// What's the Material's name ?
	double ex, poisson; // Material properties...
						// ex = Young modulus on x
						// poisson = Poisson Coef

	// put the constructor to create the class
};
