#include"Particle.hpp"

using namespace std;
//Read from file
Particle::Particle(ifstream& is){
	load(is);
}
//Would be better to have a load function independant of the constructor
void Particle::load(ifstream& is){
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
	R_=1.;
	I_ = 0.5 * m_ * R_ * R_;
	arot_ = 0.;
	vrot_=0.;
	rot_=0.;
}
void Particle::setInertia(const double m){
	m_ = m ;
	I_ = m_ * R_ * R_ * 0.5 ; 
}

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
void Particle::write(ofstream& of, const Tensor2x2& h, const Tensor2x2& hd) const{
	Vecteur rabs = h * r_ ;
	Vecteur vabs = hd * r_ + h * v_;
	of<<id_<<" "<<R_<<" "<<rabs.getx()<<" "<<rabs.gety()<<" "<<vabs.getx()<<" "<<vabs.gety()<<" "<<rot_<<" "<<vrot_<<endl;
}

void Particle::print() const{
	cerr<<"------------------------------------"<<endl;
	cerr<<"Particule "<<id_<<endl;
	cerr<<"positions reduites: "<<r_.getx()<<" "<<r_.gety()<<endl;
	cerr<<"vitesses reduites: "<<v_.getx()<<" "<<v_.gety()<<endl;
	cerr<<"accelerations reduites: "<<a_.getx()<<" "<<a_.gety()<<endl;
	cerr<<"------------------------------------"<<endl;
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

void Particle::resetA() {
	a_.set(0.,0.); 
}
