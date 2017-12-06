#include"Config.hpp"

using namespace std;

Config::Config(){
	folder_spl_ = "spl" ;
	folder_analyse_ = "analyse";
	folder_cell_ = "cell" ;
	sample_ = NULL;
	algo_ = NULL;
	cell_ = NULL;
}

Config::~Config(){
}

//Pourra etre utile plus tard
void Config::plug(Algo& algo, Cell& cell, Sample& spl){
	sample_ = &spl;
	algo_ = &algo;
	cell_ = &cell;
}

int Config::init(ifstream& is, Algo& algo, Cell& cell, Sample& spl){

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
		is >> token;
	}


	//Global check:
	bool checkSample = spl.initcheck();
	//Initialise Cell avec sample:
	if(checkSample && cell.needSample()) cell.initFromSample(spl);

	bool checkCell = cell.initcheck();
	bool checkAlgo = algo.initcheck();

	if(!checkSample){
		cerr<<" Sample::initcheck() problem."<<endl;
		return 1;
	}
	if(!checkAlgo){
		cerr<<" Algo::initcheck() problem."<<endl;
		return 1;
	}
	if(!checkCell){
		cerr<<" Cell::initcheck() problem."<<endl;
		return 1;
	}

	plug(algo,cell,spl);

	//La où sera ecrit les fichiers des particules
	sample_->initfolder(folder_spl_);

	return 0;

}

