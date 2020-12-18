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

	long NbElems = 0, NbNodes = 0, tmp_id, NbObjectCreated = 0;

	int inc, *indx, DOF_max = 2;
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
				int element_type = 0;

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
				}
			}

			if (instruction == "MAT")
			{
				std::string instruction = "";
				while (getline(model_file, line) && instruction != "EMAT")
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
					if (instruction == "K")
						DMat[0].ex = std::atol(data.c_str());
					if (instruction == "v")
						DMat[0]
							.poisson = std::atol(data.c_str());
				}
			}

			// Description of Boundary forces
			if (instruction == "BForces")
			{
				DOF_max = 0;
				//Calculation of the position and the maximum DOF
				for (int inc = 0; inc < NbNodes; inc++)
				{
					DNodes[inc].absolute_DOF_addr = DOF_max;
					DOF_max += DNodes[inc].DOF;
				}

				// Initialization of forces
				Forces = new double[DOF_max];
				for (int DOF_id = 0; DOF_id < DOF_max; DOF_id++)
					Forces[DOF_id] = 0;

				std::string instruction = "";
				while (getline(model_file, line) && instruction != "ENDBForces")
				{
					line.erase(std::remove_if(line.begin(), line.end(), [](char c) {
								   return std::isspace(static_cast<unsigned char>(c));
							   }),
							   line.end());
					if (line[0] == '#' || line.empty())
						continue;
					int delimiterPos = line.find(",");
					//Need To parse better To get all information
					instruction = line.substr(0, delimiterPos);
					std::string data = line.substr(delimiterPos + 1);

					if (instruction != "ENDBForces")
					{
						long node_id = std::stol(instruction, nullptr, 10);
						long node_index = node_id - 1;
						for (int DOF_Node_index = 0; DOF_Node_index < DNodes[node_index].DOF; DOF_Node_index++)
						{
							std::string value = data.substr(0, delimiterPos);
							data = data.substr(delimiterPos + 1);
							Forces[DNodes[node_index].absolute_DOF_addr + DOF_Node_index] = std::atol(value.c_str());
							delimiterPos = data.find(",");
						}
					}
				}
			}
			// Description of Boundary Displacements (constraints)

			if (instruction == "BDis")
			{
				indx = new int[DOF_max];
				for (int DOF_index = 0; DOF_index < DOF_max; DOF_index++)
					indx[DOF_index] = 1;
				std::string instruction = "";
				while (getline(model_file, line) && instruction != "ENDBDis")
				{

					line.erase(std::remove_if(line.begin(), line.end(), [](char c) {
								   return std::isspace(static_cast<unsigned char>(c));
							   }),
							   line.end());
					if (line[0] == '#' || line.empty())
						continue;
					int delimiterPos = line.find(",");
					//Need To parse better To get all information
					instruction = line.substr(0, delimiterPos);
					std::string data = line.substr(delimiterPos + 1);
					if (instruction != "ENDBDis")
					{
						long node_id = std::stol(instruction, nullptr, 10);
						long node_index = node_id - 1;

						for (int DOF_Node_index = 0; DOF_Node_index < DNodes[node_index - 1].DOF; DOF_Node_index++)
						{
							std::string value = data.substr(0, delimiterPos);
							data = data.substr(delimiterPos + 1);
							indx[DNodes[node_index].absolute_DOF_addr + DOF_Node_index] = std::stoi(value, nullptr, 10); // indx cest les contraintes ?
							delimiterPos = data.find(",");
						}
					}
				}
			}

			/* Qu'est ce que Max Node?
que font les sous-fonctions?
Peux-t-on simplifier?
Debug  et simplify the structure 
- rename + put more function

Clearly check  the result even display result
Display the process

and sort stuff in specific header 
(one for parsing maybe)
*/
		}
	}
	else
	{
		std::cout << "\nCan't read " << filename;
	}

	//Display nodes

	for (int index = 0; index < NbNodes; index++)
	{
		std::cout << '\n'
				  << DNodes[index].id << " (" << DNodes[index].nx << ", " << DNodes[index].nx << ") " << DNodes[index].DOF;
		std::cout << '\n';
	}
	/*
	

}
printf(" completed!\n");
fclose(file);
}
* /
		std::cout
	<< "Global Matrix building...";

// Global matrix creation ! square matrix for the moment
double **Global = new double *[DOF_max];
for (inc = 0; inc < DOF_max; inc++)
	Global[inc] = new double[DOF_max];

// Initialization of Global Matrix
for (inc = 0; inc < DOF_max; inc++)
	for (int inc2 = 0; inc2 < DOF_max; inc2++)
		Global[inc][inc2] = 0;

// Assemble Matrix
for (inc = 0; inc < count[1]; inc++)
{
	DLink1[inc]->CalculateLength(DNodes, DMat);
	DLink1[inc]->AssembleMatrix(Global, DOF_max, DNodes);
}
for (inc = 0; inc < count[2]; inc++)
{
	DPlane3[inc]->CalculateLength(DNodes, DMat);
	DPlane3[inc]->AssembleMatrix(Global, DOF_max, DNodes);
}

for (inc = 0; inc < count[3]; inc++)
{
	DBeam2[inc]->CalculateLength(DNodes, DMat);
	DBeam2[inc]->AssembleMatrix(Global, DOF_max, DNodes);
}
for (inc = 0; inc < count[42]; inc++)
{
	DPlane4[inc]->CalculateMatrix(DNodes, DMat);
	DPlane4[inc]->AssembleMatrix(Global, DOF_max, DNodes);
}

// It is a "diagonal matrix" so to implement Ks only on the diagonal so far
// implementation of Ks function of the Maximum degree of freedom to locate the Node
/* for (inc = 0; inc < DOF_max; inc++)
		if (indx[inc] == 0)
			Global[inc][inc] = 1e200; */

	/*
std::cout<<"\n"<<"Global matrix"<<"\n";
	// Globlal Matrix display
for (inc=0; inc<DOF_max;inc++)
{
	for (int inc2=0; inc2 <DOF_max ; inc2++)
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

	//ludcmp(Global, DOF_max, indx, Forces);

	// Then let's do a back substitution to get the result
	//lubksb(Global, DOF_max, indx, Forces);

	printf("completed!\nDisplacements sorting process...");
	max = 0;
	for (inc = 0; inc < NbNodes; inc++)
	{
		//	std::cout<<"\n Node "
		//		<< inc+1;
		//	for (int inc2=0;inc2 < DNodes[inc].DOF;inc2++)
		//	{
		//		std::cout<<"\n"
		//			<<Forces[DNodes[inc].absolute_DOF_addr+ inc2];
		//	}
		if ((sqrt(pow(Forces[DNodes[inc].absolute_DOF_addr], 2) + pow(Forces[DNodes[inc].absolute_DOF_addr + 1], 2))) > max)
		{
			max = sqrt(pow(Forces[DNodes[inc].absolute_DOF_addr], 2) + pow(Forces[DNodes[inc].absolute_DOF_addr + 1], 2));
			num = inc + 1;
		}
	}
	printf("completed!\n");
	std::cout << "\nThe Greater displacement is : " << max
			  << " and the Node number of the greater element is : " << num;

	std::cout << "\nux = " << Forces[DNodes[num - 1].absolute_DOF_addr]
			  << " and uy = " << Forces[DNodes[num - 1].absolute_DOF_addr + 1];

	std::cout << "\nThis problem includes " << NbNodes << " nodes and " << NbElems << " elements.\n";

	return 0;
}

void display_greetings(void)
{
	// Greetings
	std::cout << "\nWELCOME TO SOFEA";
	std::cout << "\n----------------\n\n";
}