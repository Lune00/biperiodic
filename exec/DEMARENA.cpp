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
#include"Interactions.hpp"
#include"Analyse.hpp"

using namespace std;

int main (int argc,char **argv)
{

	cout<<"DEM arena test."<<endl;


	//double theta=30.;
	//theta = theta * M_PI/180. ;
	//Vecteur fnt(0.5,0.5);
	//Vecteur n(cos(theta),sin(theta));
	//Vecteur t(-sin(theta),cos(theta));
	//Vecteur fxy = n * fnt.getx() + t * fnt.gety();
	//cerr<<"fxy (x) = "<<fxy.getx()<<" "<<fxy.gety()<<endl;


	//Config est interface pour init tout
	//a partir fichier configuration
	Config config;
	Analyse analyse;
	Algo algo;
	Cell cell;
	Sample spl;
	Interactions Int;

	//Initialisation simulation
	ifstream is(argv[1]);
	int init_status = config.init(is,algo,cell,spl,Int,analyse);
	if(init_status != 0) return 0;

	//Lancement simulation
	algo.run();

	return 0;

}
