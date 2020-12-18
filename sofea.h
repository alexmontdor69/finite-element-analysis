// sofcla.h
// headers to define all the class use in sofea.cpp
#include <string>
#include <sstream>

namespace patch
{
template <typename T>
std::string to_string(const T &n)
{
	std::ostringstream stm;
	stm << n;
	return stm.str();
}
} // namespace patch

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
public:						// Caution No private because couldn't use some other function!!
							// So to change
	long id;				// What's the node's name ?
	double nx, ny;			// Nodes Coordonates
	long DOF = 2;			// Number of DOF for this node
	long absolute_DOF_addr; // Absolute address fo the DOF

	void Init(const char *node_id, const std::string &pos_x, const std::string &pos_y);
};

void Node::Init(const char *node_id, const std::string &pos_x, const std::string &pos_y) // function to initialize data during the file reading
{
	id = std::stol(node_id, nullptr, 10);
	nx = std::atol(pos_x.c_str());
	ny = std::atol(pos_y.c_str());
	DOF = 0;
	absolute_DOF_addr = 0;
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
