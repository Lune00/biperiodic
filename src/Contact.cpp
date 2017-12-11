#include"Contact.hpp"
#include"Particle.hpp"
#include"Tenseur.hpp"

using namespace std;

//Tolerance for interpenetration
//if fabs(dn_) > tolerance_, considered to be non zero
//then check if negative or not
const double Contact::tolerance_ = 1e-14 ;

Contact::Contact(Particle* i, Particle* j){
	i_ = i ;
	j_ = j ;
	isActif_ = false;
}

//Appears twice in the code
//Made this computation a value returned function of Interactions
//But contact need a pointer to Interactions then
//Or absolute postion a member of Particle... but i do not like the idea
void Contact::Frame(Tensor2x2& h){
	double sijx = j_->getx() - i_->getx();
	double sijy = j_->gety() - i_->gety();
	//Shortest branch through periodicity:
	sijx -= floor(sijx + 0.5);
	sijy -= floor(sijy + 0.5);

	Vecteur sij(sijx,sijy);
	//Absolute vector branch:
	sij = h * sij ;
	//Norm:
	double l = sij.getNorme();
	double invl = 1./l ;
	//Build local frame:
	n_ = sij * invl;
	t_.set( - n_.gety() , n_.getx() );

	Vecteur riAbs = h * i_->getR();
	r_ = riAbs + n_ * i_->getRadius();
	//Interpenetration:
	dn_ = l - (i_->getRadius() + j_->getRadius());

	
	//Erreurs d'arrondis (ordre 10e-16) peuvent apparaitre
	//Preferer un test de tolerance a inferieur a < 0
	//Si superieur a tolerance, non zero
	//if(dn_ < 0. ){
	//A test
	if(fabs(dn_) > tolerance_ && dn_ < 0.){
		cout<<dn_<<endl;
		isActif_ = true;
	}
	else {
		isActif_ = false ;
		dt_ = 0. ;
	}
}
//Temporar: debug use
void Contact::write(ofstream& os) const{
	os<<r_.getx()<<" "<<r_.gety()<<" "<<n_.getx()<<" "<<n_.gety()<<endl;
}


void Contact::updateRelativeVelocities(){



	return ;
}

void Contact::computeForce(){


	return ;
}

void Contact::updateAccelerations(){


	return;
}

