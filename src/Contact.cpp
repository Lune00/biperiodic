#include"Contact.hpp"
#include"Particle.hpp"
#include"Tenseur.hpp"
#include"Cell.hpp"

using namespace std;

//Tolerance for interpenetration
//if fabs(dn_) > tolerance_, considered to be non zero
//then check if negative or not
const double Contact::tolerance_ = 1e-20 ;

Contact::Contact(Particle* i, Particle* j,Cell& cell){
	i_ = i ;
	j_ = j ;
	isActif_ = false;
	cell_ = &cell;
}

//Appears twice in the code
//Made this computation a value returned function of Interactions
//But contact need a pointer to Interactions then
//Or absolute postion a member of Particle... but i do not like the idea
void Contact::Frame(){

	double sijx = j_->getx() - i_->getx();
	double sijy = j_->gety() - i_->gety();
	//Shortest branch through periodicity:
	sijx -= floor(sijx + 0.5);
	sijy -= floor(sijy + 0.5);

	Vecteur sij(sijx,sijy);
	//Absolute vector branch:
	sij = cell_->geth() * sij ;
	//Norm:
	double l = sij.getNorme();
	double invl = 1./l ;
	//Build local frame:
	n_ = sij * invl;
	t_.set( - n_.gety() , n_.getx() );

	Vecteur riAbs = cell_->geth() * i_->getR();
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
	os<<r_.getx()<<" "<<r_.gety()<<" "<<n_.getx()<<" "<<n_.gety()<<" "<<f_.getx()<<" "<<f_.gety()<<endl;
}


void Contact::updateRelativeVelocity(){
	//Real velocities
	Vecteur vj = cell_->gethd() * j_->getR() +cell_->geth() * j_->getV();
	Vecteur vi = cell_->gethd() * i_->getR() +cell_->geth() * i_->getV();
	 v_ = vj - vi;
	 //Components in the contact frame:
	 v_.setx( v_ * n_ );
	 v_.sety( v_ * t_ );
	 //cout<<"Vitesse relative "<<j_->getId()<<" "<<i_->getId()<<endl;
	 //v_.print();
	 rvrot_ = j_->getVrot() - i_->getVrot();

	 //Rotational contribution added to the relative tangential componant
	 double vtr = -j_->getRadius() * j_->getVrot() -i_->getRadius() * i_->getVrot();
	 v_.addy(vtr); 
	return ;
}

//A checker:
//Temp form for DEM parameters
void Contact::computeForce(const double kn, const double kt, const double gn, const double gt, const double mus, const double dt){
	
	double fn = - kn * dn_ - gn * v_.getx();
	if(fn < 0.) fn = 0.;

	double ft = - kt * dt_ - gt * v_.gety();
	//cout<<"vrel.n = "<<v_.getx()<<" vrel.t = "<<v_.gety()<<endl;
	//cout<<"gn * v_t = "<< gt * v_.gety()<<endl;
	const double ftmax = fabs( fn * mus);
	if(fabs(ft) > ftmax){
		ft = sign(ft) * ftmax;
		dt_= ft/kt;
	}
	else dt_ += v_.gety() * dt ;
	cout<<"ftmax = "<<mus * fn<<endl;
	cout<<"dt_ = "<<dt_<<endl;
	cout<<"ft = "<<ft<<" "<<"ft/kt="<<ft/kt<<endl;
	cout<<"fn = "<<fn<<endl;

	f_.set(fn,ft);

	return ;
}

//Notice that f_.gety() refers here to ft_ (tangential force in contact frame)
//Acclerations computed in the absolute sense (lab frame)
void Contact::updateAccelerations(){
	//Expression force vector in the lab frames:
	Vecteur fxy = getfxy();
	//Turns acceleration vector in reduced coordinates:
	fxy = cell_->geth() * fxy;
	//Update linear acceleration
	j_->updateA(fxy);
	i_->updateA(-fxy);
	//Update rotational acceleration
	j_->updateArot(-f_.gety() * j_->getRadius());
	i_->updateArot(-f_.gety() * i_->getRadius());
	return;
}


void Contact::print() const{


	cerr<<"Contact entre la particule "<<i_->getId()<<" et "<<j_->getId()<<endl;
	cerr<<"interpenetration = "<<dn_<<endl;
	cerr<<"fn = "<<f_.getx()<<endl;
	cerr<<"ft = "<<f_.gety()<<endl;
	cerr<<"vx particule "<<i_->getId()<<" = "<<i_->getV().getx()<<endl;
	cerr<<"vx partjcule "<<j_->getId()<<" = "<<j_->getV().getx()<<endl;


}
