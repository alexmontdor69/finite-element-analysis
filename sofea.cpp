#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <iostream>
#include <fstream>
#include <algorithm>

#include "tools.h"

#include "sofea.h"
#include "elem-def/link1.h"
#include "elem-def/beam2.h"
#include "elem-def/plane3.h"
#include "elem-def/plane4.h"

//functions declaration
void display_greetings(void);

int main(int argc, char **argv)
{
	display_greetings();

	char *filename = "data.txt";
	// Take the first argument to read the file
	if (argc > 1)
		filename = argv[1];
	std::cout << "\n Model file : " << filename << " ...";

	// Opening model file for Parsing
	std::ifstream model_file(filename);

	// Init Side information
	Node *DNodes;
	Material *DMat;
	DMat = new Material[2];

	// Init elements
	Link1 *DLink1[300];
	Beam2 *DBeam2[3000];
	Plane3 *DPlane3[3000];
	Plane4 *DPlane4[3000];

	char Command[100], *ptr;

	long NbElems = 0, NbNodes = 0, Number, NbObjectCreated = 0;

	int inc, *indx, MaxNode = 2;
	long count[256] = {0}, num; // to know how objects created per kind of element
	double *Forces, max;
	ptr = &Command[0];

	// Function to read the file
	// Function to build the matrix
	// Function to solve the model
	// Setting the model filename

	if (model_file.is_open())
	{
		std::cout << "\nParsing File\n"
				  << filename;
		//parsing_file(model_file);
		std::string line;
		while (getline(model_file, line))
		{
			line.erase(std::remove_if(line.begin(), line.end(), [](char c) {
						   return std::isspace(static_cast<unsigned char>(c));
					   }),
					   line.end());
			if (line[0] == '#' || line.empty())
				continue;
			auto delimiterPos = line.find(",");

			//Need To parse better To get all information
			auto instruction = line.substr(0, delimiterPos);
			auto value = line.substr(delimiterPos + 1);

			// Give the number of Nodes
			if (instruction == "NODE")
			{
				NbNodes = std::stol(value, nullptr, 10);
				DNodes = new Node[NbNodes + 1];
				std::cout << "\n - Total nodes number : " << NbNodes;
			}

			// Give the number of Elements
			if (instruction == "ELEM")
			{
				NbElems = std::stol(value, nullptr, 10);
				std::cout << "\n - Total elements number : " << NbElems;
			}

			// Enter in the description of Coordonates of Nodes -> X, Y

			if (instruction == "NBLOCK")
			{
				int inc = 0;
				std::string node_number = "";
				while (getline(model_file, line) && node_number != "ENDNODEBLOCK")
				{
					line.erase(std::remove_if(line.begin(), line.end(), [](char c) {
								   return std::isspace(static_cast<unsigned char>(c));
							   }),
							   line.end());
					if (line[0] == '#' || line.empty())
						continue;
					delimiterPos = line.find(",");
					//Need To parse better To get all information
					node_number = line.substr(0, delimiterPos);
					auto coordonates = line.substr(delimiterPos + 1);
					if (node_number != "ENDNODEBLOCK")
					{
						delimiterPos = coordonates.find(",");
						std::string x = coordonates.substr(0, delimiterPos);
						std::string y = coordonates.substr(delimiterPos + 1);
						DNodes[inc].Init(node_number.c_str(), x, y);
						std::cout << "\nNode at x:" << x << " y:" << y;
					}
					inc++;
				}
			}
			// Enter in the description of Elements -> type of element & elementId, material, cross area, second moment of area / inertia, End Node A, End Node B

			if (instruction == "EBLOCK")
			{
				int inc = 0, element_type = 0;

				std::string instruction = "";
				while (getline(model_file, line) && instruction != "ENDEBLOCK")
				{
					line.erase(std::remove_if(line.begin(), line.end(), [](char c) {
								   return std::isspace(static_cast<unsigned char>(c));
							   }),
							   line.end());
					if (line[0] == '#' || line.empty())
						continue;
					delimiterPos = line.find(",");
					//Need To parse better To get all information
					instruction = line.substr(0, delimiterPos);
					std::string data = line.substr(delimiterPos + 1);
					if (instruction == "ET")
						element_type = std::stoi(data);
					// Management of "no element_type value" case
					if (element_type == 0)
						std::cout
							<< "\nERROR : No Element Type written in the model file\n";

					if (instruction != "ENDEBLOCK" && instruction != "ET")
					{
						long element_id = std::stol(instruction, nullptr, 10);

						switch (element_type)
						{
						case 1: // Link

							DLink1[count[element_type]] = new Link1(element_id, data, DNodes);
							count[element_type]++;
							NbObjectCreated++;
							break;
						case 3: // Plane 3
							DPlane3[count[element_type]] = new Plane3(element_id, data, DNodes);
							count[element_type]++;
							NbObjectCreated++;
							break;
						case 2: // beam 2
							DBeam2[count[element_type]] = new Beam2(element_id, data, DNodes);
							count[element_type]++;
							NbObjectCreated++;
							break;
						case 4: // place 4
							DPlane4[count[element_type]] = new Plane4(element_id, data, DNodes);
							count[element_type]++;
							NbObjectCreated++;
							break;
						}
					}
					inc++;
				}
			}
		}
		/* 
			case ("EBLOCK"): 
				break;

			case "ET":
				break;
			case "ENDEBLOCK":
				break;



// Description of Material K,v
			case "MAT": 
				break;
			case "K":
				break;
			case "v":
				break;
			case "EMAT":
				break;
// Description of Boundary forces
			case "BForces": 
				break;
			case "ENDBForces":
				break;
// Description of Boundary Displacements (constraints)
			case "BDis": 
				break;
			case "ENDBDis":
				break;
			 */
	}
	else
	{
		std::cout << "\nCan't read " << filename;
	}
	/*
	

	//Coordonates of ELEM
	if (!strcmp("EBLOCK", Command) && !::EndOfFile)
	{
		NextWord(file, &Command[0]);
		inc = 0;
		countlev1 = 0; // First use no element defined
		while (strcmp("ENDEBLOCK", Command) && !::EndOfFile)
		{
			if (!strcmp("ET", Command))
			{
				NextWord(file, &Command[0]);
				sscanf(Command, "%i", &countlev1);
				NextWord(file, &Command[0]);
			}

			//Utilisation de CASE perhaps !...
			// The object relating to element must be created and then deleted
			// A simple bar : Link1 on Ansys
			if (countlev1 == 1 && strcmp("ENDEBLOCK", Command) && !::EndOfFile)
			{

			}

			if (countlev1 == 2 && strcmp("ENDEBLOCK", Command) && !::EndOfFile)
			{
				DPlane3[count[countlev1]] = new Plane3(file, Number, DNodes);
				count[countlev1]++;
				NbObjectCreated++;
			}

			// Beam2 !
			if (countlev1 == 3 && strcmp("ENDEBLOCK", Command) && !::EndOfFile)
			{
				DBeam2[count[countlev1]] = new Beam2(file, Number, DNodes);
				count[countlev1]++;
				NbObjectCreated++;
			}
			// Plane4 !
			if (countlev1 == 42 && strcmp("ENDEBLOCK", Command) && !::EndOfFile)
			{
				DPlane4[count[countlev1]] = new Plane4(file, Number, DNodes);
				count[countlev1]++;
				NbObjectCreated++;
			}
			NextWord(file, &Command[0]);
		}
	}

	// Material data reading process
	if (!strcmp("MAT", Command) && !::EndOfFile)
	{
		inc = 0;
		NextWord(file, &Command[0]);
		while (strcmp("EMAT", Command) && !::EndOfFile)
		{

			if (!strcmp("K", Command))
			{
				NextWord(file, &Command[0]);
				sscanf(Command, "%le", &DMat[0].ex);
			}
			if (!strcmp("v", Command))
			{
				NextWord(file, &Command[0]);
				sscanf(Command, "%le", &DMat[0].poisson);
			}
			NextWord(file, &Command[0]);
		}
	}

	// Forces data reading process
	if (!strcmp("BForces", Command) && !::EndOfFile)
	{
		MaxNode = 0;
		//Calculation of the position and the maximum DOF
		for (inc = 0; inc < NbNodes; inc++)
		{
			DNodes[inc].CumSum = MaxNode;
			MaxNode += DNodes[inc].DOF;
		}
		Forces = new double[MaxNode];
		for (inc = 0; inc < MaxNode; inc++)
			Forces[inc] = 0;
		inc = 0;
		NextWord(file, &Command[0]);
		while (strcmp("ENDBForces", Command) && !::EndOfFile)
		{
			sscanf(Command, "%ld", &Number);
			for (inc = 0; inc < DNodes[Number - 1].DOF; inc++)
			{
				NextWord(file, &Command[0]);
				sscanf(Command, "%le", &Forces[DNodes[Number - 1].CumSum + inc]);
			}
			NextWord(file, &Command[0]);
		}
	}

	// Displacement data reading process
	if (!strcmp("BDis", Command) && !::EndOfFile)
	{
		inc = 0;
		indx = new int[MaxNode];
		for (inc = 0; inc < MaxNode; inc++)
			indx[inc] = 1;
		NextWord(file, &Command[0]);
		while (strcmp("ENDDis", Command) && !::EndOfFile)
		{
			sscanf(Command, "%ld", &Number);
			for (inc = 0; inc < DNodes[Number - 1].DOF; inc++)
			{
				NextWord(file, &Command[0]);
				sscanf(Command, "%i", &indx[DNodes[Number - 1].CumSum + inc]);
			}
			NextWord(file, &Command[0]);
		}
	}
}
printf(" completed!\n");
fclose(file);
}
* /
		std::cout
	<< "Global Matrix building...";

// Global matrix creation ! square matrix for the moment
double **Global = new double *[MaxNode];
for (inc = 0; inc < MaxNode; inc++)
	Global[inc] = new double[MaxNode];

// Initialization of Global Matrix
for (inc = 0; inc < MaxNode; inc++)
	for (int inc2 = 0; inc2 < MaxNode; inc2++)
		Global[inc][inc2] = 0;

// Assemble Matrix
for (inc = 0; inc < count[1]; inc++)
{
	DLink1[inc]->CalculateLength(DNodes, DMat);
	DLink1[inc]->AssembleMatrix(Global, MaxNode, DNodes);
}
for (inc = 0; inc < count[2]; inc++)
{
	DPlane3[inc]->CalculateLength(DNodes, DMat);
	DPlane3[inc]->AssembleMatrix(Global, MaxNode, DNodes);
}

for (inc = 0; inc < count[3]; inc++)
{
	DBeam2[inc]->CalculateLength(DNodes, DMat);
	DBeam2[inc]->AssembleMatrix(Global, MaxNode, DNodes);
}
for (inc = 0; inc < count[42]; inc++)
{
	DPlane4[inc]->CalculateMatrix(DNodes, DMat);
	DPlane4[inc]->AssembleMatrix(Global, MaxNode, DNodes);
}

// It is a "diagonal matrix" so to implement Ks only on the diagonal so far
// implementation of Ks function of the Maximum degree of freedom to locate the Node
/* for (inc = 0; inc < MaxNode; inc++)
		if (indx[inc] == 0)
			Global[inc][inc] = 1e200; */

	/*
std::cout<<"\n"<<"Global matrix"<<"\n";
	// Globlal Matrix display
for (inc=0; inc<MaxNode;inc++)
{
	for (int inc2=0; inc2 <MaxNode ; inc2++)
		std::cout<<Global[inc][inc2]
			<<"\t";
	std::cout<<"\n";
}
*/
	printf("completed!\n");
	printf("Solving by LU decomposition...");
	// Solving part______________________________________________________________________________
	//display the Result the backsubstition put all the result in the Forces vector..
	//(unfortunatly deleting them)
	// First step use the tool header and the Lower Upper Decomposition

	//ludcmp(Global, MaxNode, indx, Forces);

	// Then let's do a back substitution to get the result
	//lubksb(Global, MaxNode, indx, Forces);

	printf("completed!\nDisplacements sorting process...");
	max = 0;
	for (inc = 0; inc < NbNodes; inc++)
	{
		//	std::cout<<"\n Node "
		//		<< inc+1;
		//	for (int inc2=0;inc2 < DNodes[inc].DOF;inc2++)
		//	{
		//		std::cout<<"\n"
		//			<<Forces[DNodes[inc].CumSum+ inc2];
		//	}
		if ((sqrt(pow(Forces[DNodes[inc].CumSum], 2) + pow(Forces[DNodes[inc].CumSum + 1], 2))) > max)
		{
			max = sqrt(pow(Forces[DNodes[inc].CumSum], 2) + pow(Forces[DNodes[inc].CumSum + 1], 2));
			num = inc + 1;
		}
	}
	printf("completed!\n");
	std::cout << "\nThe Greater displacement is : " << max
			  << " and the Node number of the greater element is : " << num;

	std::cout << "\nux = " << Forces[DNodes[num - 1].CumSum]
			  << " and uy = " << Forces[DNodes[num - 1].CumSum + 1];

	std::cout << "\nThis problem includes " << NbNodes << " nodes and " << NbElems << " elements.\n";

	return 0;
}

void display_greetings(void)
{
	// Greetings
	std::cout << "\nWELCOME TO SOFEA";
	std::cout << "\n----------------\n\n";
}