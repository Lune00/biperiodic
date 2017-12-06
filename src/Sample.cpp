#include"Sample.hpp"
#include"Cell.hpp"

typedef std::vector<Particle>::iterator spit;

using namespace std;

Sample::Sample(){

	sampleIsLoaded_ = false ;
	sampleIsFilled_ = false ;
	rhodefined_ = false ;
	M_ = 0.;
	cell_ = NULL;
}

Sample::~Sample(){
}


void Sample::plugtoCell(Cell& cell){
	cell_ = &cell;
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
	setminmax();
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
			//Pb pour la derniere lecture...
			string token;
			while(is){
				if(is.eof()) break;
				Particle P(is);
				spl_.push_back(P);
			}
			//Solution provisoire:
			spl_.pop_back();

			cout<<"Nombre de particules: "<<spl_.size()<<endl;
			sampleIsLoaded_ = true ;
		}

	}
}

//Coordonnees reduites
void Sample::write(ofstream& os){
	for(std::vector<Particle>::const_iterator it = spl_.begin(); it!= spl_.end(); it++){
		it->write(os);
	}
}

//Coordonnees absolues
void Sample::writeAbsolute(ofstream& os){
	Tensor2x2 h = cell_->geth();
	for(std::vector<Particle>::const_iterator it = spl_.begin(); it!= spl_.end(); it++){
		it->write(os,h);
	}
}


//Fill particles with mass and inertia
void Sample::attributeMass(){
	if(rhodefined_) {
		for(std::vector<Particle>::iterator it = spl_.begin(); it!= spl_.end(); it++){
			double s = it->getR() * it->getR() * M_PI;
			double m = s * rho_ ;
			it->setMasse(m);
			M_ += m;
		}
		sampleIsFilled_ = true ;
	}
}

//Calcul rmin,rmax,xmin,xmax,ymin,ymax de l'echantillon
//Sert a initialiser geometrie cellule au debut simulation
void Sample::setminmax(){

	if(!sampleIsLoaded_ || spl_.size() == 0){
		cerr<<"Sample::setminmax() Charger l'echantillon avant."<<endl;
		return;
	}

	xmin_ = spl_[0].getx() - spl_[0].getRadius();
	xmax_ = spl_[0].getx() + spl_[0].getRadius();
	ymin_ = spl_[0].gety() - spl_[0].getRadius();
	ymax_ = spl_[0].gety() + spl_[0].getRadius();

	rmin_ = spl_[0].getRadius();
	rmax_ = spl_[0].getRadius();

	for(spit it= spl_.begin(); it != spl_.end(); it++){
		rmax_ = max(rmax_, it->getRadius());
		rmin_ = min(rmin_, it->getRadius());
		xmin_ = min(xmin_, it->getx() - it->getRadius());
		xmax_ = max(xmax_, it->getx() + it->getRadius());
		ymin_ = min(ymin_, it->gety() - it->getRadius());
		ymax_ = max(ymax_, it->gety() + it->getRadius());
	}
	//cout<<"xmin = "<<xmin_<<endl;
	//cout<<"xmax = "<<xmax_<<endl;
	//cout<<"ymin = "<<ymin_<<endl;
	//cout<<"ymax = "<<ymax_<<endl;
}

bool Sample::initcheck() {
	if( sampleIsLoaded_ && sampleIsFilled_ ) return true;
	else return false;
}

bool Sample::isEmptySampleFile(ifstream& is){
	return is.peek() == std::ifstream::traits_type::eof();
}

//Shift and reduce
void Sample::initReducedCoordinates(Cell& cell){
	double Lx = cell.getLx();
	double Ly = cell.getLy();
	cout<<"Lx = "<<Lx<<endl;
	cout<<"Ly = "<<Ly<<endl;
	for(spit it= spl_.begin(); it != spl_.end(); it++){
		it->setr( (it->getx()-xmin_)/Lx, (it->gety()-ymin_)/Ly);
	}
}

