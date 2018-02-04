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

//Compute the shortest branch between the two particles in contact
//The branch is the absolute distance betwwen center of particles
//Takes into account for periodicity no matter the cell shape
//More general algo to find first image (calculation performed in absolute space instead of reduced space)
//branch_ is a member of Contact because its value is needed again for computing stress or any  further analysis...
//Can be optimized (surely but see that later...)
//Compute also indexes_ (i and j) for knowing which particle image if j is an image in the interaction
//These indexes are needed to take into account the affine term interaction between particles in contact at the edges of the cell
void Contact::computeShortestBranch() {
	//Branch vector
	//Cell basis vectors: maybe should be member functions of cell
	Vecteur a0(cell_->geth().getxx(), cell_->geth().getyx());
	Vecteur a1(cell_->geth().getxy(), cell_->geth().getyy());

	//Branch vector in absolute frame:
	Vecteur rj = cell_->geth() * j_->getR();
	Vecteur ri = cell_->geth() * i_->getR();
	Vecteur d = rj - ri;
	//Test for indices that minimize the distance
	//Interaction can only be with original particles (0,0)
	//or first cell (-1,-1), (1,1) etc...
	vector<pair<int,int> > pairs;
	vector<double> l_dcarre;
	for (int i = -1 ; i != 2 ; i++){
		for(int j = -1; j != 2; j++){
			Vecteur u = d + a0 * i + a1 * j;
                        double dcarre = u * u ;
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
	//Return shortest vector branch:
	branch_ = d + a0 * indexes_.first + a1 * indexes_.second ;
	return;
}

void Contact::Frame(){

	computeShortestBranch();

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
	}
}

//TODO
void Contact::write(ofstream& os) const{
	os<<i_->getId()<<" "<<j_->getId()<<" "<<r_.getx()<<" "<<r_.gety()<<" "<<n_.getx()<<" "<<n_.gety()<<" "<<f_.getx()<<" "<<f_.gety()<<" "<<dt_<<endl;
}

//Take into account case of contact with images
void Contact::updateRelativeVelocity(){

  //Add affine term transformation for image/real particle inter
  //If j is real, ix and iy have been set to zero in previous called function

  //Real velocities
  Vecteur u(indexes_.first,indexes_.second); 
  Vecteur sj = j_->getR() + u;

  Vecteur vj = cell_->gethd() * sj + cell_->geth() * j_->getV();
  Vecteur vi = cell_->gethd() * i_->getR() + cell_->geth() * i_->getV();

  v_ = vj - vi;

  //Components in the contact frame:
  //v_n = vx nx + vy ny
  Vecteur vtmp = v_ ;

  v_.setx( vtmp * n_ );
  v_.sety( vtmp * t_ );

  //Useless rvrot_?
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
  j_->update_a(aj);
  i_->update_a(ai);

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
