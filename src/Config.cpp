#include"Config.hpp"
#include"Cell.hpp"
#include"Algo.hpp"
#include"Sample.hpp"
#include"Tenseur.hpp"
#include"Interactions.hpp"
#include"Analyse.hpp"


using namespace std;

Config::Config(){

	//Defaults folder names:
	string makemain = "mkdir -p data";
	system(makemain.c_str());

	folder_spl_ = "data/sample" ;
	folder_analyse_ = "analyse";
	folder_cell_ = "data/cell" ;
	folder_Interactions_ = "data/network";
}

Config::~Config(){

}


void Config::initfolders(Cell& cell, Sample& spl, Interactions& Int, Analyse& ana){

	//Writing paths initialisation:
	spl.initfolder(folder_spl_);
	cell.initfolder(folder_cell_);
	Int.initfolder(folder_Interactions_);
	ana.initfolder(folder_analyse_);
}

void Config::initfolders(Cell& cell, Sample& spl, Interactions& Int){
	spl.initfolder(folder_spl_);
	cell.initfolder(folder_cell_);
	Int.initfolder(folder_Interactions_);
}

int Config::init(ifstream& is, Algo& algo, Cell& cell, Sample& spl, Interactions& Int, Analyse& ana){

	//Building folders if they do not exist
	string makefolder_spl = "mkdir -p " + folder_spl_;
	string makefolder_cell = "mkdir -p " + folder_cell_;
	string makefolder_Interactions = "mkdir -p " + folder_Interactions_;
	string makefolder_analyse = "mkdir -p " + folder_analyse_;

	system(makefolder_spl.c_str());
	system(makefolder_Interactions.c_str());
	system(makefolder_cell.c_str());
	system(makefolder_analyse.c_str());
	
	//Assign a folder to each one
	initfolders(cell,spl,Int,ana);

	if(!is){
		cerr<< "Config::init : cannot open file."<<endl;
		return 1;
	}
	string token;
	is >> token;

	while(is){
		if(token=="Algo{")
		{
			algo.init(is);
		}
		if(token=="Sample{")
		{
			spl.init(is);
		}
		if(token=="Cell{")
		{
			cell.init(is);
		}
		if(token=="Interactions{")
		{
			Int.init(is);
		}
		if(token=="Analyse{")
		{
			ana.init(is);
		}
		is >> token;
	}

	//Talk Cell with sample for how to init:
	cell.talkinit(spl);
	//Plugs:
	spl.plugtoCell(cell);
	Int.plug(spl,cell);
	algo.plug(cell,spl,Int,ana);
	ana.plug(spl,cell,Int);

	//Compute reduced coordinates from cell geometry
	//if build, no if loaded
	if(!spl.loaded()) spl.initReducedCoordinates();

	//Interactions talk to sample to know if need to load contact data
	Int.build();

	//Global check:

	//Check interactions parameters:
	bool checkInteractions = Int.initcheck();
	//Check algo parameters using DEM interactions parameters
	bool checkAlgo = algo.initcheck();
	//Check cell parameters:
	bool checkCell = cell.initcheck();
	//Check sample parameters:
	bool checkSample = spl.initcheck();

	if(!checkSample){
		cerr<<"Sample::initcheck() problem."<<endl;
		return 1;
	}
	if(!checkAlgo){
		cerr<<"Algo::initcheck() problem."<<endl;
		return 1;
	}
	if(!checkCell){
		cerr<<"Cell::initcheck() problem."<<endl;
		return 1;
	}
	if(!checkInteractions){
		cerr<<"Interactions::initcheck() problem."<<endl;
		return 1;
	}

	return 0;

}
