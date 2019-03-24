#pragma once
#include <vector>
#include <random>
#include <iostream>
#include "Node.h"

using namespace std;

struct Network
{
	Network(vector<int> networkTopology);
	Network(vector<int> networkTopology, vector<int> networkTopology2);
	vector<vector<Node>> nodes;
	vector<vector<Node>> nodes2;

	vector<double> results;
	vector<double> expectedValues;
	vector<double> sumUpError;

	void printInformation();
	void feedForward(vector<double> input);
	void feedForwardMultiplex(vector<double> input);
	void backPropagation();
	void backPropagationMultiplex();
	void printError();
	void adjustWeights();
	void adjustWeightsMultiplex();

	void addError();
	void resetError();

	double eta = 0.15;
	double alpha = 0.5;


	//Random generator 
	double randomReal(const double lowerBoundary, const double upperBoundary);
	int randomInt(const int lowerBoundary, const int upperBoundary);
	static random_device seed_generator;
	static unsigned seed;
	static mt19937 mersenne_generator;
	//----------------------------------------------------------------------


};


Network::Network(vector<int> networkTopology)
{
	for (int i = 0; i < networkTopology.size(); i++)
	{
		vector<Node> tmp = vector<Node>(networkTopology[i] + 1);
		tmp[tmp.size() - 1].value = 1.0;
		if (i == networkTopology.size() - 1)
			tmp = vector<Node>(networkTopology[i]);
		nodes.push_back(tmp);

		results = vector<double>(networkTopology[networkTopology.size() - 1]);
		expectedValues = vector<double>(networkTopology[networkTopology.size() - 1]);
		sumUpError = vector<double>(networkTopology[networkTopology.size() - 1]);
	}

	for (int i = 0; i < networkTopology.size() - 1; i++)
	{
		for (int j = 0; j < networkTopology[i] + 1; j++)
		{
			vector<double> tmp = vector<double>(networkTopology[i + 1]);
			if (i == networkTopology.size() - 2)
				tmp = vector<double>(networkTopology[i + 1]);
			for (int k = 0; k < tmp.size(); k++)
			{
				tmp[k] = randomReal(0.0, 1.0);
			}
			nodes[i][j].weights = tmp;
		}
	}
}


Network::Network(vector<int> networkTopology, vector<int> networkTopology2)
{

	//Network 1
	for (int i = 0; i < networkTopology.size(); i++)
	{
		vector<Node> tmp = vector<Node>(networkTopology[i] + 1);
		tmp[tmp.size() - 1].value = 1.0;
		if (i == networkTopology.size() - 1)
			tmp = vector<Node>(networkTopology[i]);
		nodes.push_back(tmp);

		results = vector<double>(networkTopology[networkTopology.size() - 1]);
		expectedValues = vector<double>(networkTopology[networkTopology.size() - 1]);
		sumUpError = vector<double>(networkTopology[networkTopology.size() - 1]);
	}

	for (int i = 0; i < networkTopology.size() - 1; i++)
	{
		for (int j = 0; j < networkTopology[i] + 1; j++)
		{
			vector<double> tmp = vector<double>(networkTopology[i + 1]);
			if (i == networkTopology.size() - 2)
				tmp = vector<double>(networkTopology[i + 1]);
			for (int k = 0; k < tmp.size(); k++)
			{
				tmp[k] = randomReal(0.0, 1.0);
			}
			nodes[i][j].weights = tmp;
		}
	}

	//Network 2
	for (int i = 0; i < networkTopology2.size(); i++)
	{
		vector<Node> tmp = vector<Node>(networkTopology2[i] + 1);
		tmp[tmp.size() - 1].value = 1.0;
		nodes2.push_back(tmp);
	}

	for (int i = 0; i < networkTopology2.size() - 1; i++)
	{
		for (int j = 0; j < networkTopology2[i]; j++)
		{
			vector<double> tmp = vector<double>(networkTopology[i + 1]);
			for (int k = 0; k < tmp.size(); k++)
			{
				tmp[k] = randomReal(0.0, 1.0);
			}
			nodes2[i][j].weights = tmp;
		}
	}

	//Fill interconnect weights
	for (int i = 0; i < networkTopology2.size(); i++)
	{
		for (int j = 0; j < networkTopology2[i] + 1; j++)
		{
			vector<vector<double>> tmp2 = vector<vector<double>>(nodes[i].size());
			nodes2[i][j].weightsInterconnect = tmp2;
			for (int k = 0; k < nodes[i].size(); k++)
			{
				vector<double> tmp = vector<double>(nodes[i][k].weights.size());
				for (int l = 0; l < nodes[i][k].weights.size(); l++)
				{
					tmp[l] = randomReal(0.0, 1.0);
				}
				nodes2[i][j].weightsInterconnect[k] = tmp;
			}
		}
	}
}


void Network::printInformation()
{
	//cout << "Printing Node numbers.\n";
	//for (int i = 0; i < nodes.size(); i++)
	//{
	//	cout << "Row: " << i << " , Nodes: " << nodes[i].size() << endl;
	//}

	cout << "Printing Outputs Network 1.\n";
	for (int i = 0; i < nodes.size(); i++)
	{
		for (int j = 0; j < nodes[i].size(); j++)
		{
			cout << "Row: " << i << " , Node: " << j << " , output: " << nodes[i][j].value << endl;
		}
	}

	cout << "Printing Outputs Network 2.\n";
	for (int i = 0; i < nodes2.size(); i++)
	{
		for (int j = 0; j < nodes2[i].size(); j++)
		{
			cout << "Row: " << i << " , Node: " << j << " , output: " << nodes2[i][j].value << endl;
		}
	}

	//cout << "\n\nPrinting weights\n\n.";
	//for (int i = 0; i < nodes.size(); i++)
	//{
	//	for (int j = 0; j < nodes[i].size(); j++)
	//	{
	//		for (int k = 0; k < nodes[i][j].weights.size(); k++)
	//		{
	//			cout << nodes[i][j].weights[k] << " ";
	//		}
	//		cout << "\t\t";
	//	}
	//	cout << "\n\n";
	//}

	cout << "\n\n";

	printError();
}


void Network::feedForward(vector<double> input)
{
	if (input.size() - nodes[nodes.size() - 1].size() > nodes[0].size() - 1 )
		cout << "Input size too large!\n";

	if (input.size() - nodes[nodes.size() - 1].size() < nodes[0].size() - 1 )
		cout << "Input size too small!\n";

	for (int i = 0; i < results.size(); i++)
	{
		expectedValues[i] = input[input.size() - nodes[nodes.size() - 1].size() + i];
	}

	//Give input to first row
	for (int i = 0; i < nodes[0].size() - 1; i++)
	{
		nodes[0][i].value = input[i];
	}

	//Reset network
	for (int i = 1; i < nodes.size(); i++)
	{
		for (int j = 0; j < nodes[i].size(); j++)
		{
			nodes[i][j].input = 0.0;
		}
	}

	//First row
	for (int i = 0; i < nodes[0].size(); i++)
	{
		for (int j = 0; j < nodes[0][i].weights.size(); j++)
		{
			nodes[1][j].input += nodes[0][i].value * nodes[0][i].weights[j];
		}
	}

	//Feed input through deep layers
	for (int i = 1; i < nodes.size() - 1; i++)
	{
		for (int j = 0; j < nodes[i].size() - 1; j++)
		{
			nodes[i][j].value = nodes[i][j].tanhFunc(nodes[i][j].input);
			for (int k = 0; k < nodes[i][j].weights.size(); k++)
			{
				nodes[i + 1][k].input += nodes[i][j].value * nodes[i][j].weights[k];
			}
		}
	}

	//Output row
	for (int i = 0; i < nodes[nodes.size() - 1].size(); i++)
	{
		nodes[nodes.size() - 1][i].value = nodes[nodes.size() - 1][i].tanhFunc(nodes[nodes.size() - 1][i].input);
	}

	//Save results
	for (int i = 0; i < nodes[nodes.size() - 1].size(); i++)
	{
		results[i] = nodes[nodes.size() - 1][i].value;
	}
}




void Network::feedForwardMultiplex(vector<double> input)
{
	if (input.size() - nodes[nodes.size() - 1].size() > nodes[0].size() - 1)
		cout << "Input size too large!\n";

	if (input.size() - nodes[nodes.size() - 1].size() < nodes[0].size() - 1)
		cout << "Input size too small!\n";

	for (int i = 0; i < results.size(); i++)
	{
		expectedValues[i] = input[input.size() - nodes[nodes.size() - 1].size() + i];
	}

	//Give input to first row of first and second network
	for (int i = 0; i < nodes[0].size() - 1; i++)
	{
		nodes[0][i].value = input[i];
		nodes2[0][i].value = input[i];
	}

	//Reset network
	for (int i = 1; i < nodes.size(); i++)
	{
		for (int j = 0; j < nodes[i].size(); j++)
		{
			nodes[i][j].input = 0.0;
		}
	}
	for (int i = 1; i < nodes2.size(); i++)
	{
		for (int j = 0; j < nodes2[i].size(); j++)
		{
			nodes2[i][j].input = 0.0;
		}
	}

	//First row
	for (int i = 0; i < nodes2[0].size(); i++)
	{
		for (int j = 0; j < nodes2[0][i].weights.size(); j++)
		{
			nodes2[1][j].input += nodes2[0][i].value * nodes2[0][i].weights[j];
		}
	}

	for (int i = 0; i < nodes[0].size(); i++)
	{
		for (int j = 0; j < nodes[0][i].weights.size(); j++)
		{
			//Get sum from other network
			double sum = 0.0;
			for (int k = 0; k < nodes2[0].size(); k++)
			{
				sum += nodes2[0][k].value * nodes2[0][k].weightsInterconnect[i][j];
			}
			nodes[1][j].input += nodes[0][i].value * (nodes[0][i].weights[j] + sum);
		}
	}


	//Feed input through deep layers
	for (int i = 1; i < nodes2.size() - 1; i++)
	{
		for (int j = 0; j < nodes2[i].size() - 1; j++)
		{
			nodes2[i][j].value = nodes2[i][j].tanhFunc(nodes2[i][j].input);
			for (int k = 0; k < nodes2[i][j].weights.size(); k++)
			{
				nodes2[i + 1][k].input += nodes2[i][j].value * nodes2[i][j].weights[k];
			}
		}
	}

	//Doing last layer of network 2
	for (int i = 1; i < nodes2[nodes2.size() - 1].size() - 1; i++)
	{
		nodes2[nodes2.size() - 1][i].value = nodes2[nodes2.size() - 1][i].tanhFunc(nodes2[nodes2.size() - 1][i].input);
	}


	for (int i = 1; i < nodes.size() - 1; i++)
	{
		for (int j = 0; j < nodes[i].size() - 1; j++)
		{
			nodes[i][j].value = nodes[i][j].tanhFunc(nodes[i][j].input);
			for (int k = 0; k < nodes[i][j].weights.size(); k++)
			{
				//Get sum from other network
				double sum = 0.0;
				for (int l = 0; l < nodes2[i].size(); l++)
				{
					sum += nodes2[i][l].value * nodes2[i][l].weightsInterconnect[j][k];
				}

				nodes[i + 1][k].input += nodes[i][j].value * (nodes[i][j].weights[k] + sum);
			}
		}
	}

	//Output row
	for (int i = 0; i < nodes[nodes.size() - 1].size(); i++)
	{
		nodes[nodes.size() - 1][i].value = nodes[nodes.size() - 1][i].tanhFunc(nodes[nodes.size() - 1][i].input);
	}

	//Save results
	for (int i = 0; i < nodes[nodes.size() - 1].size(); i++)
	{
		results[i] = nodes[nodes.size() - 1][i].value;
	}
}


void Network::addError()
{
	for (int i = 0; i < nodes[nodes.size() - 1].size(); i++)
	{
		sumUpError[i] += 0.5 * (results[i] - expectedValues[i]) * (results[i] - expectedValues[i]);
	}
}

void Network::resetError()
{
	for (int i = 0; i < nodes[nodes.size() - 1].size(); i++)
	{
		sumUpError[i] = 0.0;
	}
}


void Network::printError()
{
	double sum = 0.0;
	for (int i = 0; i < results.size(); i++)
	{
		sum += sumUpError[i];
	}
	cout << "Error: " << sum << endl;
}


void Network::backPropagation()
{
	//Derivative of last row
	for (int i = 0; i < nodes[nodes.size() - 1].size(); i++)
	{
		nodes[nodes.size() - 1][i].derivative = nodes[nodes.size() - 1][i].tanhderf(nodes[nodes.size() - 1][i].value) * sumUpError[i];
	}

	//Derivative deep layers
	for (int i = nodes.size() - 2; i > 0; i--)
	{
		for (int j = 0; j < nodes[i].size(); j++)
		{
			//Get Sum of Nodes that are connected to
			double sum = 0.0;
			for (int k = 0; k < nodes[i][j].weights.size(); k++)
			{
				sum += nodes[i][j].weights[k] * nodes[i + 1][k].derivative;
			}
			sum *= nodes[i][j].tanhderf(nodes[i][j].value);
			nodes[i][j].derivative = sum;
		}
	}

	adjustWeights();
}


void Network::adjustWeights()
{
	for (int i = 0; i < nodes.size() - 1; i++)
	{
		for (int j = 0; j < nodes[i].size(); j++)
		{
			for (int k = 0; k < nodes[i][j].weights.size(); k++)
			{
				nodes[i][j].weights[k] -= eta * nodes[i][j].value * nodes[i + 1][k].derivative + nodes[i][j].oldDerivative * alpha;
				nodes[i][j].oldDerivative = nodes[i][j].derivative;
			}
		}
	}
}


void Network::backPropagationMultiplex()
{
	//Doing first network
	//Derivative of last row
	for (int i = 0; i < nodes[nodes.size() - 1].size(); i++)
	{
		nodes[nodes.size() - 1][i].derivative = nodes[nodes.size() - 1][i].tanhderf(nodes[nodes.size() - 1][i].value) * (results[i] - expectedValues[i]);
	}

	//Derivative deep layers
	for (int i = nodes.size() - 2; i > 0; i--)
	{
		for (int j = 0; j < nodes[i].size(); j++)
		{
			//Get Sum of Nodes that are connected to
			double sum = 0.0;
			for (int k = 0; k < nodes[i][j].weights.size(); k++)
			{
				sum += nodes[i][j].weights[k] * nodes[i + 1][k].derivative;
			}
			sum *= nodes[i][j].tanhderf(nodes[i][j].value);
			nodes[i][j].derivative = sum;
		}
	}

	//Doing second network
	for (int i = nodes2.size() - 1; i > 0; i--)
	{
		for (int j = 0; j < nodes2[i].size(); j++)
		{
			//Get Sum of Nodes that are connected to
			double sum = 0.0;
			for (int k = 0; k < nodes2[i][j].weightsInterconnect.size(); k++)
			{
				for (int l = 0; l < nodes2[i][j].weightsInterconnect[k].size(); l++)
				{
					sum += nodes2[i][j].weightsInterconnect[k][l] * nodes[i + 1][l].derivative;
				}
			}

			if (i <= nodes2.size() - 2)
			{
				for (int k = 0; k < nodes2[i][j].weights.size(); k++)
				{
					sum += nodes2[i][j].weights[k] * nodes2[i + 1][k].derivative;
				}
			}

			sum *= nodes2[i][j].tanhderf(nodes2[i][j].value);
			nodes2[i][j].derivative = sum;
		}
	}

	adjustWeightsMultiplex();
}


void Network::adjustWeightsMultiplex()
{
	for (int i = 0; i < nodes.size() - 1; i++)
	{
		for (int j = 0; j < nodes[i].size(); j++)
		{
			for (int k = 0; k < nodes[i][j].weights.size(); k++)
			{
				nodes[i][j].weights[k] -= eta * nodes[i][j].value * nodes[i + 1][k].derivative + nodes[i][j].oldDerivative * alpha;
				nodes[i][j].oldDerivative = nodes[i][j].derivative;
			}
		}
	}

	for (int i = 0; i < nodes2.size(); i++)
	{
		for (int j = 0; j < nodes2[i].size(); j++)
		{
			for (int k = 0; k < nodes2[i][j].weightsInterconnect.size(); k++)
			{
				for (int l = 0; l < nodes2[i][j].weightsInterconnect[k].size(); l++)
				{
					nodes2[i][j].weightsInterconnect[k][l] -= eta * nodes2[i][j].value * nodes[i][k].value * nodes[i + 1][l].derivative + nodes2[i][j].oldDerivative * alpha;
					nodes2[i][j].oldDerivative = nodes2[i][j].derivative;
				}
			}
		}
	}

	for (int i = 0; i < nodes2.size() - 1; i++)
	{
		for (int j = 0; j < nodes2[i].size(); j++)
		{
			for (int k = 0; k < nodes2[i][j].weights.size(); k++)
			{
				nodes2[i][j].weights[k] -= eta * nodes2[i][j].value * nodes2[i + 1][k].derivative + nodes2[i][j].oldDerivative * alpha;
				nodes2[i][j].oldDerivative = nodes2[i][j].derivative;
			}
		}
	}

}









//----------------------------------------------------------------------
//Random engine call
double Network::randomReal(const double lowerBoundary, const double upperBoundary)
{
	uniform_real_distribution<double> distribution_real(lowerBoundary, upperBoundary);
	return distribution_real(mersenne_generator);
}


int Network::randomInt(const int lowerBoundary, const int upperBoundary)
{
	uniform_int_distribution<int> distribution_int(lowerBoundary, upperBoundary);
	return distribution_int(mersenne_generator);
}


//Random engine initialisation
random_device Network::seed_generator;
unsigned Network::seed = seed_generator();
mt19937 Network::mersenne_generator(Network::seed);