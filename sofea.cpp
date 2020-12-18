#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <iostream>
#include <fstream>
#include <algorithm>

#include "tools.h"

#include "sofea.h"
#include "parsing.h"
#include "elem-def/link1.h"
#include "elem-def/beam2.h"
#include "elem-def/plane3.h"
#include "elem-def/plane4.h"

//functions declaration
void display_greetings(void);

// Parsing functions

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
	/*
	* Optimization possible if we pre-read the file and create the exact amount of element
	*/
	Link1 *DLink1[300];
	Beam2 *DBeam2[3000];
	Plane3 *DPlane3[3000];
	Plane4 *DPlane4[3000];

	long NbElems = 0, NbNodes = 0;

	int *BDis, DOF_max = 0;
	long count[256] = {0}, num; // to know how objects created per kind of element
	double *Forces, *Displacements;

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
			line = clear_line(line);
			if (line_not_readable(line))
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
				int node_index = 0;
				std::string node_number = "";
				while (getline(model_file, line) && node_number != "ENDNODEBLOCK")
				{
					line = clear_line(line);
					if (line_not_readable(line))
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
						DNodes[node_index].Init(node_number.c_str(), x, y);
						std::cout << "\nNode at x:" << x << " y:" << y;
					}
					node_index++;
				}
			}
			// Enter in the description of Elements -> type of element & elementId, material, cross area, second moment of area / inertia, End Node A, End Node B

			if (instruction == "EBLOCK")
			{
				int element_type = 0;

				std::string instruction = "";
				while (getline(model_file, line) && instruction != "ENDEBLOCK")
				{
					line = clear_line(line);
					if (line_not_readable(line))
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
							break;
						case 3: // Plane 3
							DPlane3[count[element_type]] = new Plane3(element_id, data, DNodes);
							count[element_type]++;
							break;
						case 2: // beam 2
							DBeam2[count[element_type]] = new Beam2(element_id, data, DNodes);
							count[element_type]++;
							break;
						case 4: // place 4
							DPlane4[count[element_type]] = new Plane4(element_id, data, DNodes);
							count[element_type]++;
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
					line = clear_line(line);
					if (line_not_readable(line))
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
				for (int node_index = 0; node_index < NbNodes; node_index++)
				{
					DNodes[node_index].absolute_DOF_addr = DOF_max;
					DOF_max += DNodes[node_index].DOF;
				}

				// Initialization of forces
				Forces = new double[DOF_max];
				for (int DOF_id = 0; DOF_id < DOF_max; DOF_id++)
					Forces[DOF_id] = 0;

				std::string instruction = "";
				while (getline(model_file, line) && instruction != "ENDBForces")
				{
					line = clear_line(line);
					if (line_not_readable(line))
						continue;

					int delimiterPos = line.find(",");
					//Need To parse better To get all information
					instruction = line.substr(0, delimiterPos);
					std::string data = line.substr(delimiterPos + 1);

					if (instruction != "ENDBForces")
					{
						long node_id = std::stol(instruction, nullptr, 10);
						long node_index = node_id - 1;

						create_boundary_forces(data, Forces, DNodes, node_index);
					}
				}
			}
			// Description of Boundary Displacements (constraints)

			if (instruction == "BDis")
			{
				BDis = new int[DOF_max];
				for (int DOF_index = 0; DOF_index < DOF_max; DOF_index++)
					BDis[DOF_index] = 1;
				std::string instruction = "";
				while (getline(model_file, line) && instruction != "ENDBDis")
				{
					line = clear_line(line);
					if (line_not_readable(line))
						continue;

					int delimiterPos = line.find(",");
					//Need To parse better To get all information
					instruction = line.substr(0, delimiterPos);
					std::string data = line.substr(delimiterPos + 1);
					if (instruction != "ENDBDis")
					{
						long node_id = std::stol(instruction, nullptr, 10);
						long node_index = node_id - 1;

						create_boundary_displacement(data, BDis, DNodes, node_index);
					}
				}
			}

			/*
que font les sous-fonctions?
Peux-t-on simplifier?
Debug  et simplify the structure index
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

	/* 	std::cout << "\n\n---BDis content----\n\n";
	for (int index = 0; index < DOF_max; index++)
	{
		std::cout << "\n DOF : " << patch::to_string(index + 1) << " : " << patch::to_string(*(Forces + index)).c_str();
		std::cout << "\n";
	}
	for (int index = 0; index < DOF_max; index++)
	{
		std::cout << "\n DOF : " << patch::to_string(index + 1) << " : " << patch::to_string(*(BDis + index)).c_str();
		std::cout << "\n";
	} */

	std::cout << "\nParsing completed!\n";

	std::cout << "\n----------------------\n\n";

	std::cout << "\nAssemble Matrices...\n";
	std::cout << "\nBuilding main Matrix...\n";

	// main_matrix : square matrix with DOF_max side length
	double **main_matrix = new double *[DOF_max];
	for (int DOF_index = 0; DOF_index < DOF_max; DOF_index++)
		main_matrix[DOF_index] = new double[DOF_max];

	// Initialization of main_matrix Matrix
	for (int index_line = 0; index_line < DOF_max; index_line++)
		for (int index_row = 0; index_row < DOF_max; index_row++)
			main_matrix[index_line][index_row] = 0;

	//Filling main matrix
	// Assemble Matrix
	for (int element_index = 0; element_index < count[1]; element_index++)
	{
		DLink1[element_index]->CreateStiffMatrix(DNodes, DMat);
		DLink1[element_index]->AssembleMatrix(main_matrix, DOF_max, DNodes);
	}
	for (int element_index = 0; element_index < count[2]; element_index++)
	{
		DBeam2[element_index]->CreateStiffMatrix(DNodes, DMat);
		DBeam2[element_index]->AssembleMatrix(main_matrix, DOF_max, DNodes);
	}
	for (int element_index = 0; element_index < count[3]; element_index++)
	{
		DPlane3[element_index]->CreateStiffMatrix(DNodes, DMat);
		DPlane3[element_index]->AssembleMatrix(main_matrix, DOF_max, DNodes);
	}

	for (int element_index = 0; element_index < count[4]; element_index++)
	{
		DPlane4[element_index]->CreateStiffMatrix(DNodes, DMat);
		DPlane4[element_index]->AssembleMatrix(main_matrix, DOF_max, DNodes);
	}

	/**
 * Positionning the Boundary displacement constraints
 * 
 * It is a "diagonal matrix" 
 * Node stiffness in the diagonal..
 * If constraints exists, stiffness is considered infinite (value 1e200)
 * 
 * 
 * */
	std::cout << "\nAdding Boundary Displacement constraints main Matrix...\n";
	for (int DOF_index = 0; DOF_index < DOF_max; DOF_index++)
		if (BDis[DOF_index] == 0)
			main_matrix[DOF_index][DOF_index] = 1e200;

	std::cout << "\nMain Matric Creation completed!\n";
	std::cout << "\n----------------------\n\n";
	std::cout << "\nSolving by LU decomposition...";

	printf("Solving by LU decomposition...");
	// Solving part______________________________________________________________________________
	//display the Result the backsubstition put all the result in the Forces vector..
	//(unfortunatly deleting them)
	// First step use the tool header and the Lower Upper Decomposition
	ludcmp(main_matrix, DOF_max, BDis, Forces);

	// Then let's do a back substitution to get the result
	lubksb(main_matrix, DOF_max, BDis, Forces);

	std::cout << "completed!\n";
	std::cout << "\n----------------------\n\n";
	std::cout << "RESULTS\n";
	std::cout << "Sorting Displacements ...";

	//Forces becomes the array of displacement - Could be fixed

	double max_displacement = 0;
	long max_displaced_node_id;
	for (long node_index = 0; node_index < NbNodes; node_index++)
	{
		if ((sqrt(pow(Forces[DNodes[node_index].absolute_DOF_addr], 2) + pow(Forces[DNodes[node_index].absolute_DOF_addr + 1], 2))) > max_displacement)
		{
			max_displacement = sqrt(pow(Forces[DNodes[node_index].absolute_DOF_addr], 2) + pow(Forces[DNodes[node_index].absolute_DOF_addr + 1], 2));
			max_displaced_node_id = node_index + 1;
		}
	}

	std::cout << "completed!\n";
	std::cout << "\n----------------------\n\n";
	std::cout << "RESULTS\n";
	std::cout << "\nThe Greater displacement is : " << max_displacement
			  << " and the Node number of the greater element is : " << max_displaced_node_id;

	std::cout << "\nux = " << Forces[DNodes[max_displaced_node_id - 1].absolute_DOF_addr]
			  << " and uy = " << Forces[DNodes[max_displaced_node_id - 1].absolute_DOF_addr + 1];

	return 0;
}

void display_greetings(void)
{
	// Greetings
	std::cout << "\nWELCOME TO SOFEA";
	std::cout << "\n----------------\n\n";
}
