#include"Particle.hpp"

using namespace std;
//sx et sy doivent etre compris entre 0 et 1 (initialiement dans la cellule)
void Particle::init(){
	//Les vecteurs appellent constructeur, tout egal a 0
	m_=1.;
	R_=0.1;
	vrot_=0.;
	rot_=0.;
	t_=0.;
}

void Particle::write(ofstream& of){
	of<< r_.getx()<<" "<<r_.gety()<<" "<<v_.getx()<<" "<<v_.gety()<<endl;
}

void Particle::affiche(){
	cout<<"Coordonnees reduites: "<<r_.getx()<<" "<<r_.gety()<<endl;
}	
