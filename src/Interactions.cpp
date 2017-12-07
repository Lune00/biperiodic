#include"Interactions.hpp"
#include"Sample.hpp"
#include"Cell.hpp"
#include"Particle.hpp"

using namespace std;

Interactions::Interactions(){
	//Temp
	//Distances a mettre en Rmax
	dv_ = 2. ;
	dsv_ = 3. ;
	nv_ = 1 ;
	nsv_ = 1;
	scale_=string();
	checkInteractions_ = false;
	initDistances_ = false;
}

Interactions::~Interactions(){

}


void Interactions::init(ifstream& is){




}


void Interactions::initScale(){

	double scale;

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

	dv_ *= scale;
	dsv_ *= scale;
	initDistances_ = true ;

}

void Interactions::plug(Sample& spl){
	spl_ = & spl;
	initScale();
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


