#include"Particle.hpp"

using namespace std;
//sx et sy doivent etre compris entre 0 et 1 (initialiement dans la cellule)
//Read from file
Particle::Particle(ifstream& is){
	double x = 0. ;
	double y = 0. ;
	double vx = 0. ;
	double vy = 0. ;
	is >> id_ >> R_ ;
	r_.load(is);
	v_.load(is);
	a_.load(is);
	is >> rot_ >> vrot_ >> arot_;
}

//Never used in practice, only for debugging
void Particle::init(){
	//Les vecteurs appellent constructeur, tout egal a 0
	m_=1.;
	R_=0.1;
	I_ = 0.5 * m_ * R_ * R_;
	arot_ = 0.;
	vrot_=0.;
	rot_=0.;
}

//To be re written for symetry to load

//Reduced coordinates: read for simulations
void Particle::write(ofstream& os) const{
	os<< id_<<" "<<R_<<" ";
	r_.write(os);
	v_.write(os);
	a_.write(os);
	os << rot_<<" "<<vrot_<<" "<<arot_<<endl;
}

//Only suitable for analysis(not used to load sample, or start
//simuation)
void Particle::write(ofstream& of, const Tensor2x2& h) const{
	Vecteur rabs = h * r_ ;
	//TODO : bring hd here!
	//Vecteur vabs = hd * r_ + h * v_;
	of<<id_<<" "<<R_<<" "<<rabs.getx()<<" "<<rabs.gety()<<" "<<v_.getx()<<" "<<v_.gety()<<" "<<rot_<<" "<<vrot_<<endl;
}

void Particle::affiche() const{
	cout<<"Coordonnees reduites: "<<r_.getx()<<" "<<r_.gety()<<endl;
	cout<<"Vitesses reduites: "<<v_.getx()<<" "<<v_.gety()<<endl;
}	

void Particle::updateR(const double dt){
	r_ = r_ + v_ * dt + a_ * dt * dt * 0.5 ;
}

//Verlet algo: half each time
void Particle::updateV(const double dt_2){
	v_ = v_ + a_ * dt_2 ;
	return;
}

void Particle::updateRot(const double dt){
	rot_ = rot_ + vrot_ * dt + arot_ * dt * dt * 0.5;
	return;
}

void Particle::updateVrot(const double dt_2){
	vrot_ = vrot_ + arot_ * dt_2 ;
}
