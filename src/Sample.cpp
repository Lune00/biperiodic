#include"Sample.hpp"


typedef std::vector<Particle>::iterator spit;

using namespace std;

Sample::Sample(){

	sampleIsLoaded_ = false ;
	sampleIsFilled_ = false ;
	rhodefined_ = false ;
}

Sample::~Sample(){
}


void Sample::init(ifstream& is){

	string token;
	is >> token;
	while(is){
		if(token=="includeFile") is >> fichsample_;
		if(token=="rho") {
			rhodefined_= true;
			is >> rho_;
		}
		is >> token ;
	}

	loadSample();
	attributeMass();
}

void Sample::loadSample(){

	ifstream is(fichsample_.c_str());
	if(!is){
		cerr<<"Sample::loadSample() : can not open file."<<endl;
		return;
	}
	else{
		if(isEmptySampleFile(is)) {
			cerr<<fichsample_<<" est vide!"<<endl;
			return ;
		}
		else
		{
			string token;
			while(is){
				Particle P(is);
				spl_.push_back(P);
				is >> token;
			}

			cout<<"Nombre de particules: "<<spl_.size()<<endl;
			sampleIsLoaded_ = true ;
		}

	}
}

void Sample::write(ofstream& os){
	for(std::vector<Particle>::const_iterator it = spl_.begin(); it!= spl_.end(); it++){
		it->write(os);
	}
}

//Fill particles with mass and inertia
void Sample::attributeMass(){
	if(rhodefined_) {
		for(std::vector<Particle>::iterator it = spl_.begin(); it!= spl_.end(); it++){
			double s = it->getR() * it->getR() * M_PI;
			double m = s * rho_ ;
			it->setMasse(m);
		}


		sampleIsFilled_ = true ;

	}
}

//Calcul rmin,rmax,xmin,xmax,ymin,ymax de l'echantillon
//Sert a initialiser geometrie cellule au debut simulation
void Sample::setminmax(){

	if(!sampleIsLoaded_){
		cerr<<"Sample::setminmax() Charger l'echantillon avant."<<endl;
	}
	double xmin = spl_[0].getx();
	double xmax = spl_[0].getx();
	double ymin = spl_[0].gety();
	double ymax = spl_[0].gety();
	double rmin = spl_[0].getRadius();
	double rmax = spl_[0].getRadius();

	for(spit it= spl_.begin(); it != spl_.end(); it++){
		rmax = max(rmax, it->getRadius());
		rmin = min(rmin, it->getRadius());
	}

}

bool Sample::initcheck() {
	if( sampleIsLoaded_ && sampleIsFilled_ ) return true;
	else return false;
}

bool Sample::isEmptySampleFile(ifstream& is){
	return is.peek() == std::ifstream::traits_type::eof();
}
