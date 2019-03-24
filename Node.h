#pragma once
#include <math.h>

using namespace std;


struct Node
{
	double input = 0.0;
	double value = 0.0;
	double sigmoid(double input);
	void relu(double input);
	double sigmoidDerv(double input);
	double tanhFunc(double input);
	double tanhderf(double input);
	vector<double> weights;
	vector<vector<double>> weightsInterconnect;
	double derivative = 0.0;
	double oldDerivative = 0.0;
};


double Node::sigmoid(double input)
{
	return (1.0 / (1.0 + exp(-input)));
}


double Node::sigmoidDerv(double input)
{
	return sigmoid(input) * (1.0 - sigmoid(input));
}


double Node::tanhFunc(double input)
{
	return tanh(input);
}

double Node::tanhderf(double input)
{
	return (1.0 - input * input);
}


void Node::relu(double input)
{
	value = input;
	if (input < 0)
		value = 0;

}
