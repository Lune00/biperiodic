#include"Config.hpp"

using namespace std;

//Essai mecanique: 2 fichiers input parametres pour Config, echantillon initial (xmin,xmax...)
//Prendra input eventuellement dans un fichier
Config::Config(){
	//A lire dans un fichier a terme
	//Pour le moment tout en vitesse (cinematique)
	for(int i=0;i<4;i++){
		BCU[i] = 'v';
	}
	//BCU[3]='f';
	//Par souci de convention, les DL non controles sont init a 0
	//Par definition, le taux de cisaillement pur est egale a la moitie de LdUser.yx(ou xy)
	LdUser_.set(0.,1.,0.,0.);
	StressUser_.set(0.,0.,0.,0.);
}

Config::~Config(){
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
		if(token=="Cell{")
		{
			cell.init(is);
		}
		if(token=="Sample{")
		{
			spl.init(is);
		}
		is >> token;
	}

	//Global check:
	bool checkSample = spl.initcheck();
	bool checkAlgo = algo.initcheck();
	bool checkCell = cell.initcheck();

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

	return 0;

}

