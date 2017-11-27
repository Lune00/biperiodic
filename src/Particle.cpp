#include"Particle.hpp"

using namespace std;
//sx et sy doivent etre compris entre 0 et 1 (initialiement dans la cellule)
void Particle::init(double L){
	//Les vecteurs appellent constructeur, tout egal a 0
	m_=1.;
	R_=0.1;
	vrot_=0.;
	rot_=0.;
	t_=0.;
	s_.set(r_.getx() / L,r_.gety() / L);
}

void Particle::write(ofstream& of){
	of<< r_.getx()<<" "<<r_.gety()<<" "<<v_.getx()<<" "<<v_.gety()<<endl;
}

void Particle::affiche(){
	cout<<"Coordonnees absolues: "<<r_.getx()<<" "<<r_.gety()<<endl;
	cout<<"Coordonnees reduites: "<<s_.getx()<<" "<<s_.gety()<<endl;
}	
