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

//Compute the shortest branch between the two particles in contact
//The branch is the absolute distance betwwen center of particles
//Takes into account for periodicity no matter the cell shape
//More general algo to find first image (calculation performed in absolute space instead of reduced space)
//branch_ is a member of Contact because its value is needed again for computing stress or any  further analysis...
//Can be optimized (surely but see that later...)
//Compute also indexes_ (i and j) for knowing which particle image if j is an image in the interaction
//These indexes are needed to take into account the affine term interaction between particles in contact at the edges of the cell
void Contact::computeShortestBranch() {

	cerr<<"Called from Frame()"<<endl;
	//Need the branch vector (in absolute units)
	double sijx = j_->getx() - i_->getx();
	double sijy = j_->gety() - i_->gety();
	Vecteur sij(sijx,sijy);
	//Branch vector
	Vecteur d = cell_->geth() * sij;
	//Cell basis vectors: maybe should be member functions of cell
	Vecteur a0(cell_->geth().getxx(), cell_->geth().getyx());
	Vecteur a1(cell_->geth().getxy(), cell_->geth().getyy());

	//Test for indices that minimize the distance
	//Interaction can only be with original particles (0,0)
	//or first cell (-1,-1), (1,1) etc...
	vector<pair<int,int> > pairs;
	vector<double> l_dcarre;
	for (int i = -1 ; i != 2 ; i++){
		for(int j = -1; j != 2; j++){
			double dcarre = d * d + d * a0 * 2 * i + a0 * a0 * i * i +  d * a1 * 2 * j + a0 * a1 * 2 * i * j + a1 * a1 * j * j ;
			l_dcarre.push_back(dcarre);
			pairs.push_back(std::make_pair(i,j));
		}
	}
	//Find minimum:
	std::vector<double>::iterator it = std::min_element(l_dcarre.begin(),l_dcarre.end());
	double dmin = * it ;
	int k = distance( l_dcarre.begin(), it);
	//Get matching pairs of indexes
	indexes_ = pairs[k];
	//cerr<<"dmin="<<sqrt(dmin)<<" for "<<indexes.first<<" "<<indexes.second<<endl;
	if(i_->getId()==2 && j_->getId()==3) cerr<<"dmin="<<sqrt(dmin)<<" for "<<indexes_.first<<" "<<indexes_.second<<endl;

	//Return shortest vector branch:
	branch_ = d + a0 * indexes_.first + a1 * indexes_.second ;

	return;
}

void Contact::Frame(){

	computeShortestBranch();

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
	//os<<r_.getx()<<" "<<r_.gety()<<" "<<n_.getx()<<" "<<n_.gety()<<" "<<f_.getx()<<" "<<f_.gety()<<endl;
	os<<r_.getx()<<" "<<r_.gety()<<" "<<getbranch().getx()<<" "<<getbranch().gety()<<endl;
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
	//If j is real, ix and iy have been set to zero in previous called function
	vj.addx( - cell_->gethd().getxx() * indexes_.first - cell_->gethd().getxy() * indexes_.second);
	vj.addy( - cell_->gethd().getyx() * indexes_.first - cell_->gethd().getyy() * indexes_.second);

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
