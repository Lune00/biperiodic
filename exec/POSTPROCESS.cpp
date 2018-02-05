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
	cout<<"DEM post-processing."<<endl;

	//Eventually move the directory to a variable of the
	//input configuration file of post-processing

	//Default:
	string folder_post_process = "postprocess";

	//Config est interface pour init tout
	//a partir fichier configuration
	Config config;
	Analyse analyse;
	Cell cell;
	Sample spl;

	//On utilise la liste de verlet pour recreer la liste de contacts actifs
	Interactions Int;

	config.initfolders(cell,spl,Int,analyse);

	unsigned int nstart=0;
	unsigned int nend=0;
	unsigned int nperiod=1;
	string filetime = string();

	//Initialisation post-process data range
	ifstream is(argv[1]);

	if(!is){
		cerr<<"Post-process : cannot open file."<<endl;
		return 1;
	}

	string token;
	is >> token;

	while(is){
		if(token=="nstart") is >> nstart;
		if(token=="nend") is >> nend;
		if(token=="period") is >> nperiod;
		if(token=="time") is >> filetime;
		if(token=="dir") is >> folder_post_process;
		if(token=="Analyse{") analyse.init(is);
		is >> token;
	}

	is.close();

	cerr<<"Résumé:"<<endl;
	cerr<<"nstart : "<<nstart<<endl;
	cerr<<"nend : "<<nend<<endl;
	cerr<<"period : "<<nperiod<<endl;
	cerr<<"time : "<<filetime<<endl;
	cerr<<"directory : "<<folder_post_process<<endl;

	//Build dir:
	string makefolder_pp = "mkdir -p " + folder_post_process;
	system(makefolder_pp.c_str());

	//Reset folder for analyse -> distinguish between runanalyse
	analyse.initfolder(folder_post_process);

	//Load time: correspondance between tics and time (s)
	vector<double> t;
	ifstream time(filetime.c_str());

	if(!time.is_open()){
		cerr<<"Impossible to open "<<filetime<<"."<<endl;
		return 1;
	}
	else{
		double temp1, temp2 ;
		while(time){
			time >> temp1 >> temp2 ;
			if(time.eof()) break;
			t.push_back(temp2);
		}
	}
	cerr<<"Number of total steps : "<<t.size()<<endl;

	//Read simu-setup file for : h0 and rho
	double density = 0.;
	
	is.open("simusetup.txt");
	if(!is){
		cerr<<"Impossible to read simusetup.txt"<<endl;
		return 1;
	}
	else{
		string token;
		is >> token;
		while(is){
			if(token=="density") is >> density;
			//if(token=="h0") cell.readh0(is);
			is >> token;
		}
	}

	if(density < 1e-20) {
		cerr<<"Density was not found."<<endl;
		return 1;
	}

	for(unsigned int i = nstart ; i != nend ; i+=nperiod){

		//Load sample:
		spl.setfiletoload(i);
		spl.setrho(density);
		spl.attributeMass();
		spl.setminmax();

		//Load network:

		//Load cell:
		cell.load(i);

		//Analyse:
	}

	return 0;
}

