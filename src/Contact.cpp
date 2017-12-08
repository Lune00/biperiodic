#include"Contact.hpp"
#include"Particle.hpp"
#include"Tenseur.hpp"

using namespace std;


Contact::Contact(Particle* i, Particle* j){
	i_ = i ;
	j_ = j ;
	isActif_ = false;
}

void Contact::Frame(Tensor2x2& h){


	double l = (i_->getR() - j_->getR()).getNorme();


}
