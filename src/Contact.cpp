#include"Contact.hpp"
#include"Particle.hpp"
#include"Tenseur.hpp"
#include"Cell.hpp"

using namespace std;

//Tolerance for interpenetration
//if fabs(dn_) > tolerance_, considered to be non zero
//then check if negative or not
//NOT used
const double Contact::tolerance_ = 1e-20 ;

Contact::Contact(Particle* i, Particle* j,Cell* cell){
	i_ = i ;
	j_ = j ;
	isActif_ = false;
	cell_ = cell;
	branch_.setx(0.);
	branch_.sety(0.);
	dt_ = 0. ;
}

Contact::Contact(Particle* i, Particle* j,Cell* cell,pair<int,int> indexes){
	i_ = i ;
	j_ = j ;
	isActif_ = false;
	cell_ = cell;
	branch_.setx(0.);
	branch_.sety(0.);
	dt_ = 0. ;
	indexes_ = indexes;
}

Contact::Contact(ifstream& is){
}

//On suppose que les indexes ne changent pas
//entre 2 renouvellements de la liste de verlet
//on ne change pas d'image
void Contact::computeBranch(){
	Vecteur u(indexes_.first,indexes_.second); 
	branch_ = cell_->geth() * ( j_->getR() + u  - i_->getR()); 
}

//On recalcule les indexes ici
//Call for postprocessing
void Contact::computeShortestBranch(){

	double sijx = j_->getx() - i_->getx() ;
	double sijy = j_->gety() - i_->gety() ;

	indexes_.first = -(int)floor(sijx + 0.5) ;
	indexes_.second = -(int)floor(sijy + 0.5) ;

	sijx += indexes_.first;
	sijy += indexes_.second;

	Vecteur sij(sijx,sijy);
	branch_ = cell_->geth() * sij;
}

void Contact::Frame(){

	computeBranch();

	double l = branch_.getNorme();
	double invl = 1./l ;

	//Build local frame:
	n_ = branch_ * invl;
	t_.set( - n_.gety() , n_.getx() );

	Vecteur riAbs = cell_->geth() * i_->getR();
	r_ = riAbs + n_ * i_->getRadius();
	//Interpenetration:
	dn_ = l - (i_->getRadius() + j_->getRadius());

	if(dn_ < 0.){
		isActif_ = true;
	}
	else {
		isActif_ = false ;
		dt_ = 0. ;
	}
}

void Contact::write(ofstream& os) const{
	os<<i_->getId()<<" "<<j_->getId()<<" "<<v_.getx()<<" "<<v_.gety()<<" "<<n_.getx()<<" "<<n_.gety()<<" "<<f_.getx()<<" "<<f_.gety()<<" "<<dt_<<endl;
}

//Take into account case of contact with images
void Contact::updateRelativeVelocity(){

	//Add affine term transformation for image/real particle inter

	Tensor2x2 h = cell_->geth();
	Tensor2x2 hd = cell_->gethd();

	//Real velocities
	Vecteur u(indexes_.first,indexes_.second); 
	Vecteur sj = j_->getR() + u;

	Vecteur vj = hd * sj + h * j_->getV();
	Vecteur vi = hd * i_->getR() + h * i_->getV();

	Vecteur vxy = vj - vi;

	//Components in the contact frame:
	//v_n = vx nx + vy ny

	v_.setx( vxy * n_ );
	v_.sety( vxy * t_ );

	//Rotational contribution added to the relative tangential componant
	double vtr = -j_->getRadius() * j_->getVrot() -i_->getRadius() * i_->getVrot();
	v_.addy(vtr); 

	return ;
}

void Contact::computeForce(const double kn, const double kt, const double gn, const double gt, const double mus, const double dt){

	double fn = - kn * dn_ - gn * v_.getx();
	if(fn < 0.) fn = 0.;

	double ft = - kt * dt_ - gt * v_.gety();

	const double ftmax = fabs( fn * mus);

	if(fabs(ft) > ftmax){
		//Si glissant on fixe dt_, s'il est negatif on le laisse negatif sinon positif
		ft = sign(ft) * ftmax;
		dt_ =  ft/kt;
	}
	else dt_ += v_.gety() * dt ;

	f_.set(fn,ft);

	return ;
}


void Contact::updateAccelerations(){

	//Expression force vector in the lab frames:
	Vecteur fxy = getfxy();

	const double mi = i_->getMasse();
	const double mj = j_->getMasse();

	Vecteur ai = -fxy / mi ;
	Vecteur aj = fxy / mj ;

	//Update linear acceleration
	j_->add_acc(aj);
	i_->add_acc(ai);

	const double ft = - f_.gety();
	//Update rotational acceleration
	j_->updateArot(ft);
	i_->updateArot(ft);

	return;

}

void Contact::print() const{
	cerr<<"Contact entre la particule "<<i_->getId()<<" et "<<j_->getId()<<endl;
	cerr<<"interpenetration = "<<dn_<<endl;
	cerr<<"fn = "<<f_.getx()<<endl;
	cerr<<"ft = "<<f_.gety()<<endl;
}
