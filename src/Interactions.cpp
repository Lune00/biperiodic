#include"Interactions.hpp"
#include"Sample.hpp"
#include"Cell.hpp"
#include"Particle.hpp"

using namespace std;

Interactions::Interactions(){
	//Temp
	//Distances a mettre en Rmax
	dv_ = 0. ;
	dsv_ = 0. ;
	nv_ = 0 ;
	nsv_ = 0;
	scale_=string();
	initScale_ = false;
	initdv_ = false;
	initdsv_ = false;
	checkInteractions_ = false;
}

Interactions::~Interactions(){

}


void Interactions::init(ifstream& is){


	string token;
	is >> token;
	while(is){

		if(token=="setUnit") is >> scale_;

		if(token=="dv"){
			is >> dv_;
			initdv_ = true;
		}
		if(token=="dsv"){
			is >> dsv_;
			initdsv_ = true;
		}
		if(token=="niterv")  is >> nv_;
		if(token=="nitersv") is >> nsv_;
		if(token=="}") break;
		is >> token;
	}


}


void Interactions::initScale(){

	double scale;
	double epsilon=0.00001;

	if(scale_ == "Rmax" ) {
		scale = spl_->getrmax();
	}
	else if(scale_ == "Rmin") {
		scale = spl_->getrmin();
	}
	else{
		cerr<<"Interactions::initScale() : no internal scale found! Specify Rmin or Rmax"<<endl;
		return ;
	}

	//Check if dv et dsv have been initialised "correctly"
	if(dv_ < epsilon || dsv_ < epsilon)
	{
		initScale_ = false;
		cerr<<"Interactions::initScale() : distances are not well set! Too small or zero."<<endl;
		return;
	}

	dv_ *= scale;
	dsv_ *= scale;
	initScale_ = true ;

	return ;
}

void Interactions::plug(Sample& spl){
	spl_ = & spl;
	initScale();
}


bool Interactions::initcheck(){
	bool initNiter = false;
	if( nv_ != 0 && nsv_ != 0 ) initNiter = true; 
	checkInteractions_ = initNiter && initScale_ && initdv_ && initdsv_;


	cout<<"d verlet (scale)  = "<<dv_<<endl;
	cout<<"d sverlet (scale)  = "<<dsv_<<endl;
	return checkInteractions_;
}
	

void Interactions::updateverlet(const int tic){
	if( tic % nsv_ == 0 ) updatesvlist();
	if( tic % nv_ == 0 ) updatevlist();
}

void Interactions::updatesvlist(){

	svlist_.clear();

	//If less than two particles, no interaction possible
	if(spl_->getsize() < 2 ) return ;
	

	vector<Particle>* ps = spl_->getSample();
	//Get h:
	Tensor2x2 h = spl_->getCell()->geth();

	for(std::vector<Particle>::const_iterator iti = ps->begin(); iti!=ps->end();iti++){
		for(std::vector<Particle>::const_iterator itj = iti+1; itj!=ps->end();itj++){
			if( near( *(iti), *(itj), h , dsv_) ) {
				cout<<"Near!"<<endl;
			}
		}
	}
}

//True if distance between "surface" of particle i and j are lower than d
bool Interactions::near(const Particle& i, const Particle& j,const Tensor2x2& h,const double d) const{

	double sijx = j.getx() - i.getx();
	double sijy = j.gety() - i.gety();
	//Shortest branch through periodicity:
	sijx -= floor(sijx + 0.5);
	sijy -= floor(sijy + 0.5);

	Vecteur sij(sijx,sijy);
	//Absolute vector branch:
	sij = h * sij ;
	//Test distance compared to dsv
	if( sij.getNorme() - d < j.getRadius() + i.getRadius() ) return true;
	else return false;

}

void Interactions::updatevlist(){


}


