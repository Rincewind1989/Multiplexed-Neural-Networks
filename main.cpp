#include <iostream>
#include <fstream>
#include <vector>
#include "Network.h"
#include <stdlib.h>
#include <string>
#include <sstream>

vector<vector<double>> createXORTest();
vector<vector<double>> createSumTest();
vector<vector<double>> loadLetterRecoigniton();

int main()
{
	Network network({ 16,10, 10 ,1 }, {16, 10, 10});

	//Loaded data
	vector<vector<double>> testArray = loadLetterRecoigniton();
	cout << "Loaded data...start training.\n";
	for (int j = 0; j < 100000; j++)
	{
		random_shuffle(testArray.begin(), testArray.end());
		cout << "Doing Iteration " << j << endl;
		network.resetError();
		for (int i = 0; i < 10; i++)
		{
			network.feedForwardMultiplex(testArray[i]);
			network.addError();
		}
		network.backPropagationMultiplex();
		system("CLS");
		network.printInformation();
	}
	int a;
	cout << "Finished\n";
	cin >> a;

	return 0;
}


vector<vector<double>> createXORTest()
{
	//----------------------------------------------------------------------
	//Random engine initialisation
	random_device seed_generator;
	unsigned seed = seed_generator();
	mt19937 mersenne_generator(seed);

	//Random engine call
	uniform_int_distribution<int> distribution_int(0, 1);

	int number = 100000;
	vector<vector<double>> tests = vector<vector<double>>(number);
	for (int i = 0; i < number; i++)
	{
		vector<double> tmp = { double(distribution_int(mersenne_generator)) , double(distribution_int(mersenne_generator)) , 0};

		tmp[2] = 0;
		if (tmp[0] != tmp[1])
			tmp[2] = 1;

		tests[i] = tmp;
	}

	return tests;
}



vector<vector<double>> createSumTest()
{
	//----------------------------------------------------------------------
	//Random engine initialisation
	random_device seed_generator;
	unsigned seed = seed_generator();
	mt19937 mersenne_generator(seed);

	//Random engine call
	uniform_int_distribution<int> distribution_int(0, 1);

	int number = 100000;
	vector<vector<double>> tests = vector<vector<double>>(number);
	for (int i = 0; i < number; i++)
	{
		vector<double> tmp = vector<double>(11);
		for (int j = 0; j < 10; j++)
			tmp[j] = double(distribution_int(mersenne_generator));

		tmp[10] = 0;
		double tmp2 = 0.0;
		for (int j = 0; j < 10; j++)
			tmp2 += tmp[j];
		if (tmp2 > 2)
			tmp[10] = 1;

		tests[i] = tmp;
	}

	return tests;
}


vector<vector<double>> loadLetterRecoigniton()
{
	vector<vector<double>> data;
	string line;
	ifstream myfile("letter-recognition.data");
	if (myfile.is_open())
	{
		while (getline(myfile, line))
		{
			vector<double> tmp;
			std::stringstream   linestream(line);
			std::string         value;

			while (getline(linestream, value, ','))
			{
				if (int(value[0]) >= 65 && int(value[0]) <= 90)
					tmp.push_back( (int(value[0]) - 65) / 26);
				else
					tmp.push_back(std::stoi(value));
			}
			data.push_back(tmp);
		}
		myfile.close();
	}

	return data;
}