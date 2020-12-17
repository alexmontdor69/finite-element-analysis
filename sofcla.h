//#include "tools.h"
//#define PI 3.14159265359

// sofcla.h
// headers to define all the class use in sofea.cpp
int EndOfFile = 0;

using namespace std;
void NextWord(FILE *file, char *ptr) // what is the next word/Number
{
	int Letter, inc = 0;
	Letter = fgetc(file);
	while ((Letter != (int)'\n' && Letter != '\t' && Letter != ',' && Letter != ' ' || inc == 0) && Letter != EOF)
	{
		if ((Letter == (int)'\n' || Letter == '\t' || Letter == ',' || Letter == ' ') && inc == 0)
			inc = 0;
		else if (Letter != (int)'\n' && Letter != '\t' && Letter != ',' && Letter != ' ' && Letter != EOF)
		{
			*ptr = (char)Letter;
			ptr++;
			inc++;
		}
		Letter = fgetc(file);
		if (Letter == EOF)
			::EndOfFile = 1;
	}
	*ptr = '\0';
}

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

//-----------NODE--------------------------------------------------------
class Node
{
public:			   // Caution No private because couldn't use some other function!!
				   // So to change
	long number;   // What's the node's name ?
	double nx, ny; // Nno fileodes Coordonates
	long DOF;	   // Number of DOF for this element
	long CumSum;   // "Adress" in the vectors.
	void Init(FILE *file);
};

void Node::Init(FILE *file) // function to initialize data during the file reading
{
	char Command[100];
	sscanf(Command, "%ld", &number);
	NextWord(file, &Command[0]);
	sscanf(Command, "%le", &nx);
	NextWord(file, &Command[0]);
	sscanf(Command, "%le", &ny);
	DOF = 0;
	CumSum = 0;
}

//----------NODE----------------------------------------------------------
class Material
{
public:
	long number;		// What's the Material's name ?
	double ex, poisson; // Material properties...
						// ex = Young modulus on x
						// poisson = Poisson Coef

	// put the constructor to create the class
};

