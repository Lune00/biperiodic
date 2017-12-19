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

//BEFORE WAS IN FRAME
//Return the shortest branch between the two particles in contact
//The branch is the absolute distance betwwen center of particles
//in the LAB frame
Vecteur Contact::getbranch() const {
	double sijx = j_->getx() - i_->getx();
	double sijy = j_->gety() - i_->gety();


	sijx -= floor(sijx + 0.5);
	sijy -= floor(sijy + 0.5);

	Vecteur sij(sijx,sijy);
	//Absolute vector branch:
	sij = cell_->geth() * sij ;
	return sij;
}


//WIPPP
//Only j can be image, to check
void Contact::detectImage(){

	//cerr<<" "<<endl;
	//cerr<<"In interaction "<<i_->getId()<<" "<<j_->getId()<<endl;
	double sijx = j_->getx() - i_->getx();
	double sijy = j_->gety() - i_->gety();

	//WIP: indice for the add term for image particle
	ix_ = floor(sijx + 0.5);
	iy_ = floor(sijy + 0.5);

	//cerr<<"sijx = "<<sijx<<" sijy = "<<sijy<<endl;
	//Shortest branch through periodicity:
	double short_sijx = sijx -floor(sijx + 0.5);
	double short_sijy = sijy -floor(sijy + 0.5);

	//cerr<<"short_sijx = "<<short_sijx<<" sijy = "<<short_sijy<<endl;

	//cerr<<abs(ix_)<<" "<<abs(iy_)<<endl;
	Vecteur shortest(short_sijx,short_sijy);
	Vecteur realbranch(sijx,sijy);

	shortest = cell_->geth() * shortest ;
	realbranch = cell_->geth() * realbranch ;

	//cerr<<"real branch :";
	//realbranch.print();
	//cerr<<"shortest branch : ";
	//shortest.print();
	
	if(realbranch.getNorme() > shortest.getNorme()) {

		//cerr<<"Particle "<<j_->getId()<<" is an image"<<endl;
		j_is_image = true ;

	}
}

void Contact::Frame(){
	//Norm:
	//double sijx = j_->getx() - i_->getx();
	//double sijy = j_->gety() - i_->gety();
	////Shortest branch through periodicity:
	//sijx -= floor(sijx + 0.5);
	//sijy -= floor(sijy + 0.5);
	//Vecteur sij(sijx,sijy);
	////Absolute vector branch:
	//sij = cell_->geth() * sij ;
	detectImage();

	Vecteur sij = getbranch();
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
		//cout<<dn_<<endl;
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


//Take into account case of contact with images
void Contact::updateRelativeVelocity(){

	//Real velocities
	Vecteur vj = cell_->gethd() * j_->getR() +cell_->geth() * j_->getV();
	Vecteur vi = cell_->gethd() * i_->getR() +cell_->geth() * i_->getV();

	//WIP
	//ToTEST
	//cerr<<"Relative velo term : "<<ix_<<" "<<iy_<<endl;

	//Add affine term transformation for image/real particle inter
	//If real, ix and iy have been set to zero in previous called function
	vj.addx( - cell_->gethd().getxx() * ix_ - cell_->gethd().getxy() * iy_);
	vj.addy( - cell_->gethd().getyx() * ix_ - cell_->gethd().getyy() * iy_);

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
	//	cout<<"ftmax = "<<mus * fn<<endl;
	//	cout<<"dt_ = "<<dt_<<endl;
	//	cout<<"ft = "<<ft<<" "<<"ft/kt="<<ft/kt<<endl;
	//	cout<<"fn = "<<fn<<endl;

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
