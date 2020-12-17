#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>

#include "tools.h"

#include "sofcla.h"
#include "elem-def/link1.h"
#include "elem-def/beam3.h"
#include "elem-def/plane2.h"
#include "elem-def/plane4.h"

using namespace std;

//functions declaration
void display_greetings(void);

int main(int argc, char **argv)
{
	display_greetings();

	char *filename = "data.txt";
	// Take the first argument to read the file
	if (argc > 1)
		filename = argv[1];
	std::cout << "\n Parsing model file " << filename << " ...";

	// Objects Creation
	std::ifstream model_file(filename);
	Node *DNodes;

	// Object the material
	Material *DMat;
	DMat = new Material[2];

	// Object for the elements
	Link1 *DLink1[300];
	Plane2 *DPlane2[3000];
	Beam3 *DBeam3[3000];
	Plane42 *DPlane42[3000];

	char Command[100], *ptr;

	long NbElems = 0, NbNodes = 0, Number, NbObjectCreated = 0;
	int inc, countlev1, *indx, MaxNode = 2;
	long count[256] = {0}, num; // to know how objects created per kind of element
	double *Forces, max;
	ptr = &Command[0];

	// Function to read the file
	// Function to build the matrix
	// Function to solve the model
	// Setting the model filename

	if (model_file.is_open())
	{
		std::cout << "\ncan read " << filename;
	}
	else
	{
		std::cout << "\nCan't read " << filename;
	}

	/*
	if ((model_file = fopen(filename, "r")) != NULL)
	{
		while (!::EndOfFile)
		{
			NextWord(file, &Command[0]);
			//To know how nodes is there
			if (!strcmp("NODE", Command) && !::EndOfFile)
			{
				NextWord(file, &Command[0]);
				sscanf(Command, "%ld", &NbNodes);
				DNodes = new Node[NbNodes + 1];
			}

			//To know how nodes there is
			if (!strcmp("ELEM", Command) && !::EndOfFile)
			{
				NextWord(file, &Command[0]);
				sscanf(Command, "%d", &NbElems);
			}

			// Coordonates of NODE
			if (!strcmp("NBLOCK", Command) && !::EndOfFile)
			{
				NextWord(file, &Command[0]);
				inc = 0;
				while (strcmp("ENDNODEBLOCK", Command) && !::EndOfFile)
				{
					DNodes[inc].Init(file);
					NextWord(file, &Command[0]);
					inc++;
				}
			}

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
						DLink1[count[countlev1]] = new Link1(file, Number, DNodes); // Replace by a constructor
					cout << "\n Reading file " << filename << " ...";

	model_file.open(model_file, iosif (countlev1 == 2 && strcmp("ENDEBLOCK", Command) && !::EndOfFile)
					{
						DPlane2[count[countlev1]] = new Plane2(file, Number, DNodes);
						count[countlev1]++;
						NbObjectCreated++;
					}

					// beam3 !
					if (countlev1 == 3 && strcmp("ENDEBLOCK", Command) && !::EndOfFile)
					{
						DBeam3[count[countlev1]] = new Beam3(file, Number, DNodes);
						count[countlev1]++;
						NbObjectCreated++;
					}
					if (countlev1 == 42 && strcmp("ENDEBLOCK", Command) && !::EndOfFile)
					{

						DPlane42[count[countlev1]] = new Plane42(file, Number, DNodes);
						count[countlev1]++;
						NbObjectCreated++;
					}
					NextWord(file, &Command[0]);
				}
			}
cout << "\n Reading file " << filename << " ...";

	model_file.open(model_file, ios
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
*/
	std::cout << "Global Matrix building...";

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
		DPlane2[inc]->CalculateLength(DNodes, DMat);
		DPlane2[inc]->AssembleMatrix(Global, MaxNode, DNodes);
	}

	for (inc = 0; inc < count[3]; inc++)
	{
		DBeam3[inc]->CalculateLength(DNodes, DMat);
		DBeam3[inc]->AssembleMatrix(Global, MaxNode, DNodes);
	}
	for (inc = 0; inc < count[42]; inc++)
	{
		DPlane42[inc]->CalculateMatrix(DNodes, DMat);
		DPlane42[inc]->AssembleMatrix(Global, MaxNode, DNodes);
	}

	// It is a "diagonal matrix" so to implement Ks only on the diagonal so far
	// implementation of Ks function of the Maximum degree of freedom to locate the Node
	/* for (inc = 0; inc < MaxNode; inc++)
		if (indx[inc] == 0)
			Global[inc][inc] = 1e200; */

	/*
cout<<"\n"<<"Global matrix"<<"\n";
	// Globlal Matrix display
for (inc=0; inc<MaxNode;inc++)
{
	for (int inc2=0; inc2 <MaxNode ; inc2++)
		cout<<Global[inc][inc2]
			<<"\t";
	cout<<"\n";
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
		//	cout<<"\n Node "
		//		<< inc+1;
		//	for (int inc2=0;inc2 < DNodes[inc].DOF;inc2++)
		//	{
		//		cout<<"\n"
		//			<<Forces[DNodes[inc].CumSum+ inc2];
		//	}
		if ((sqrt(pow(Forces[DNodes[inc].CumSum], 2) + pow(Forces[DNodes[inc].CumSum + 1], 2))) > max)
		{
			max = sqrt(pow(Forces[DNodes[inc].CumSum], 2) + pow(Forces[DNodes[inc].CumSum + 1], 2));
			num = inc + 1;
		}
	}
	printf("completed!\n");
	cout << "\nThe Greater displacement is : " << max
		 << " and the Node number of the greater element is : " << num;

	cout << "\nux = " << Forces[DNodes[num - 1].CumSum]
		 << " and uy = " << Forces[DNodes[num - 1].CumSum + 1];

	cout << "\nThis problem includes " << NbNodes << " nodes and " << NbElems << " elements.\n";

	return 0;
}

void display_greetings(void)
{
	// Greetings
	std::cout << "\nWELCOME TO SOFEA";
	std::cout << "\n----------------\n\n";
}