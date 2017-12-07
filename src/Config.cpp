#include"Config.hpp"
#include"Cell.hpp"
#include"Algo.hpp"
#include"Sample.hpp"
#include"Tenseur.hpp"
#include"Interactions.hpp"

using namespace std;

Config::Config(){
	folder_spl_ = "spl" ;
	folder_analyse_ = "analyse";
	folder_cell_ = "cell" ;
//	sample_ = NULL;
//	algo_ = NULL;
//	cell_ = NULL;
}

Config::~Config(){
}

//Pourra etre utile plus tard
//void Config::plug(Algo& algo, Cell& cell, Sample& spl){
//	sample_ = &spl;
//	algo_ = &algo;
//	cell_ = &cell;
//}

int Config::init(ifstream& is, Algo& algo, Cell& cell, Sample& spl, Interactions& Int){

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
	//Compute reduced coordinates from cell geometry
	if(checkCell) spl.initReducedCoordinates(cell);

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

	//Plugs:
	spl.plugtoCell(cell);
	Int.plug(spl);
	algo.plug(cell,spl,Int);

	//Writing paths:
	spl.initfolder(folder_spl_);

	////Tests:
	//ofstream s("reduced.txt");
	//ofstream r("absolute.txt");
	//sample_->write(s);
	//sample_->writeAbsolute(r);
	

	return 0;

}

