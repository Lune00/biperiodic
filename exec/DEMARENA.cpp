#include<iostream>
#include<string>
#include<cmath>
#include<fstream>
#include<vector>

#include"Particle.hpp"
#include"Cell.hpp"
#include"Config.hpp"
#include"Algo.hpp"
#include"Sample.hpp"


using namespace std;

int main (int argc,char **argv)
{

	cout<<"DEM arena."<<endl;

	//Config est interface pour init tout
	//a partir fichier configuration
	Config config;
	Algo algo;
	Cell cell;
	Sample spl;

	//Initialisation simulation
	ifstream is(argv[1]);
	int init = config.init(is,algo,cell,spl);
	if(init!=0) return 0;

	//Lancement simulation
	cout<<"Lancement."<<endl;
	algo.run();

	return 0;

}
