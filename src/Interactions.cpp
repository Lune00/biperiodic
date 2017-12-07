#include"Interactions.hpp"
#include"Sample.hpp"
#include"Cell.hpp"

using namespace std;

Interactions::Interactions(){
	//Temp
	dv_ = 2. ;
	dsv_ = 3. ;
	nv_ = 1 ;
	nsv_ = 1;
}

Interactions::~Interactions(){

}

void Interactions::plug(Sample& spl){
	spl_ = & spl;
}


void Interactions::updateverlet(const int tic){
	if( tic % nsv_ == 0 ) updatesvlist();
	if( tic % nv_ == 0 ) updatevlist();
}

void Interactions::updatesvlist(){

	//If less than two particles, no interaction possible
	if(spl_->getsize() < 2 ) return ;

	vector<Particle>* ps = spl_->getSample();
	//Get h:
	Tensor2x2 h = spl_->getCell()->geth();

	for(std::vector<Particle>::const_iterator iti = ps->begin(); iti!=ps->end();iti++){
		for(std::vector<Particle>::const_iterator itj = iti+1; itj!=ps->end();itj++){

			double sijx = itj->getx() - iti->getx();
			double sijy = itj->gety() - itj->gety();
			//Shortest branch through periodicity:
			sijx -= floor(sijx + 0.5);
			sijy -= floor(sijy + 0.5);

			Vecteur sij(sijx,sijy);
			sij.print();
			//Absolute vector branch:
			sij = h * sij ;
			//Test distance compared to dsv


		}
	}



}

void Interactions::updatevlist(){


}


