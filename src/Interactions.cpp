#include"Interactions.hpp"
#include"Sample.hpp"
#include"Cell.hpp"
#include"Particle.hpp"

using namespace std;

Interactions::Interactions(){
	dv_ = 0. ;
	nv_ = 0 ;
	scale_=string();
	folder_ = string();
	initScale_ = false;
	fInteractions_ = "inter.txt";
	initdv_ = false;
	checkInteractions_ = false;
	initkn_ = false;
	initkt_ = false;
	initgn_ = false;
	initgt_ = false;
	initmus_ = false;
	setgnmax_ = false;

	ofstream debug("debugInteractions.txt");
	debug.close();
	debug.open("tmp.txt");
	debug.close();
}

Interactions::~Interactions(){

}


void Interactions::init(ifstream& is){

	string token;
	is >> token;
	while(is){

		if(token=="setUnit") is >> scale_;

		if(token=="dv"){
			is >> dv_;
			initdv_ = true;
		}
		if(token=="niterv"){
		  is >> nv_;
		}

		if(token=="kn") {
			is >> kn_;
			initkn_ = true;
		}
		if(token=="kt"){
			is >> kt_;
			initkt_ = true;
		}
		if(token=="gn"){
			is >> gn_;
			initgn_ = true;
		}
		if(token=="gnmax"){
			setgnmax_ = true;
			initgn_ = true;
		}
		if(token=="gt"){
			is >> gt_;
			initgt_ = true;
		}
		if(token=="mu"){
			is >> mus_;
			initmus_ = true;
		}

		if(token=="}") break;
		is >> token;
	}
}


void Interactions::initScale(){

	double scale;
	double epsilon=0.001;

	if(scale_ == "Rmax" ) {
		scale = spl_->getrmax();
	}
	else if(scale_ == "Rmin") {
		scale = spl_->getrmin();
	}
	else{
		cerr<<"Interactions::initScale() : no internal scale found! Specify Rmin or Rmax"<<endl;
		return ;
	}

	//Check if dv has been initialised "decently"
	if(dv_ < epsilon )
	{
		initScale_ = false;
		cerr<<"Interactions::initScale() : distances are not well set! Too small or zero."<<endl;
		return;
	}

	dv_ *= scale;
	initScale_ = true ;

	return ;
}

bool Interactions::initcheck() {

  bool initNiter = false;

  if( nv_ != 0 ) initNiter = true; 

  checkInteractions_ = initNiter && initScale_ && initdv_ && checkDEMparameters();

  cout<<"distance Verlet  = "<<dv_<<endl;

  if(!checkDEMparameters()) cout<<"Interactions::initcheck() : a DEM parameter is not initialised in config file. stop"<<endl;

  return (checkInteractions_);
}

void Interactions::plug(Sample& spl,Cell& cell){
	spl_ = & spl;
	cell_ = &cell;
}

//called in config.cpp after plug. Si array exists already
//Ask sample if loaded sample or new one (using build)
//If loaded, then load dt for tangential forces
void Interactions::talkinit() {
  if(spl_->loaded()){
    unsigned int filetoload = spl_->filetoload();
    load(filetoload);
  }
  return ;
}


//Build pairs_ list et setscale to verlet cutoff. Load existing contact network if a sample is loaded(todo)
void Interactions::build(){

  //On cree la liste de toutes les interactions possibles:
  //pairs_:
  unsigned int N = spl_->getsize();

  for(std::vector<Particle>::iterator iti = spl_->getSample()->begin();iti!=spl_->getSample()->end();iti++)
  {
    for(std::vector<Particle>::iterator itj = iti + 1; itj!=spl_->getSample()->end();itj++)
    {
      Contact c(*iti,*itj,*cell_);
      pairs_.push_back(c);
    }
  }

  cerr<<"Nombre d'interactions possibles: "<<pairs_.size()<<endl;
  initScale();

  //TODO:Est ce qu'on load? on lit les contacts dans un fichier
  talkinit();

}

//Load the dt for the contacts
void Interactions::load(const int k){

  string filename = formatfile( folder_, fInteractions_, k );
  ifstream is(filename.c_str());
  if(!is){
    cerr<<"Interactions::load "<<filename<<" fail."<<endl;
    return ;
  }
  else if(is.peek() == std::ifstream::traits_type::eof()){
    cerr<<"Pas de contacts a charger."<<endl;
    return ;
  }
  else{
    //Build contacts list to fill array_dt:
    string token;
    while(is){
      if(is.eof()) break;
      //TODO
    }

  }
  is.close();
}

void Interactions::read(ifstream& is){

}

bool Interactions::checkDEMparameters() const{
  return (initkn_ && initkt_ && initgn_ && initgt_ && initmus_);
}



void Interactions::updateverlet(const int tic){
  if( tic % nv_ == 0 ) updatevlist();
}

//Return smallest vector branch between particle i and particle j
Vecteur Interactions::getShortestBranch(const Particle* i, const Particle* j) const{
  //Need the branch vector (in absolute units)
  double sijx = j->getx() - i->getx();
  double sijy = j->gety() - i->gety();
  Vecteur sij(sijx,sijy);
  //Branch vector
  Vecteur d = cell_->geth() * sij;
  //Cell basis vectors
  Vecteur a0(cell_->geth().getxx(), cell_->geth().getyx());
  Vecteur a1(cell_->geth().getxy(), cell_->geth().getyy());

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
  pair<int,int> indexes = pairs[k];
  //Return shortest vector branch:
  return (d + a0 * indexes.first + a1 * indexes.second) ; 
}


//True if distance between "surface" of particle i and j are lower than d
bool Interactions::near(const Particle* i, const Particle* j,const double d) const{

  Vecteur shortest_branch = getShortestBranch(i,j);

  //Test distance compared to dsv

  if( shortest_branch.getNorme() - d < j->getRadius() + i->getRadius() ) return true;

  else return false;
}

void Interactions::updatevlist(){

  vlist_.clear();
  Tensor2x2 h = cell_->geth();

  for(vector<Contact>::iterator it = pairs_.begin() ; it != pairs_.end(); it ++ )
  {
    if( near( it->geti(), it->getj(), dv_ ) ){
      Contact * pc  = &(*it) ;
      vlist_.push_back(pc);
    }
  }

  cerr<<"Taille de la liste de Verlet: "<<vlist_.size()<<endl;
}


//Build contact list (activated interactions)
void Interactions::detectContacts(){

  clist_.clear();
  for( vector<Contact*>::iterator it = vlist_.begin() ; it != vlist_.end(); it++){

    (*it)->Frame();

    if( (*it)->isActif() ) {
      Contact * pc = (*it) ;
      clist_.push_back(pc);
    }
  }

  cerr<<"Nombre de contacts : "<< clist_.size()<<endl;
}


void Interactions::computeForces(const double dt){


  for(vector<Contact*>::iterator it = clist_.begin() ; it != clist_.end() ; it++){

    (*it)->updateRelativeVelocity();
    (*it)->computeForce(kn_,kt_,gn_,gt_,mus_,dt);
    (*it)->updateAccelerations();

  //    //TODO: METTRE CALCUL DU TENSEUR CONTRAINTES STATIQUE ICI

  }

  // Transform acceleration into reduced coordinates:
  Tensor2x2 hdd = cell_->gethdd();
  Tensor2x2 hd = cell_->gethd();
  Tensor2x2 h = cell_->geth();
  Tensor2x2 hinv = h.getInverse();


  for(std::vector<Particle>::iterator it = spl_->getSample()->begin();it!=spl_->getSample()->end();it++)
  {

    Vecteur a_red = hinv * ( it->getA() - hd * (it)->getV() * 2. - hdd * (it)->getR());
    //Vecteur a_red = hinv * (it->getA());
    (it)->setAcceleration(a_red);

    //TODO CALCUL DU TENSEUR CONTRAINTES CINEMATIQUES
  }

}


void Interactions::askNumberOfContacts() const{
  cerr<<"Nombre de contacts: "<<clist_.size()<<endl;
}


//Write contact network
void Interactions::writeContacts(int k) const {

  //if(clist_.size()==0) cerr<<"step "<<k<<" Il n'y a pas de contact"<<endl;
  string filename = formatfile(folder_, fInteractions_, k);
  ofstream file(filename.c_str());

  //  for(vector<int>::const_iterator it = clist_.begin(); it != clist_.end(); it++){
  //TODO
  // }
}

double Interactions::getElasticEnergy() const {

  double E = 0.;
  //for(vector<int>::const_iterator it = clist_.begin(); it != clist_.end(); it++){
  //  //	double dn = vlist_[*it].getdn();
  //  //	double dt = vlist_[*it].getdt();
  //  //	E += 0.5 * (dn) * (dn) * kn_ ;
  //}
  return E;
}

//TODO : Cette loop pourrait etre margee avec celle du calcul des forces dans computeForces pour gagner en efficacite
void Interactions::computeInternalStress(){
  //Static
  double sxx_s = 0. ;
  double sxy_s = 0. ;
  double syx_s = 0. ;
  double syy_s = 0. ;
  //Kinetic:
  double sxx_c = 0. ;
  double sxy_c = 0. ;
  double syy_c = 0. ;

  //Static stress : Loop over contacts:
  //for(vector<int>::const_iterator it = clist_.begin(); it != clist_.end(); it++){ 
  //  //	Vecteur branch = vlist_[*it].getbranch();
  //  //	Vecteur force = vlist_[*it].getfxy();
  //  //	sxx_s += branch.getx() * force.getx();
  //  //	sxy_s += branch.gety() * force.getx();
  //  //	syx_s += branch.getx() * force.gety();
  //  //	syy_s += branch.gety() * force.gety();
  //}

  ////Kinetic stress : Loop over particles:
  ////	for(std::vector<Particle>::const_iterator it = spl_->inspectSample().begin();it!=spl_->inspectSample().end();it++)
  ////	{
  ////		//Partie fluctuante : v = h * sdot
  ////		Vecteur v = cell_->geth() * it->getV();
  ////		double m = it->getMasse();
  ////		sxx_c += m * v.getx() * v.getx();
  ////		sxy_c += m * v.getx() * v.gety();
  ////		syy_c += m * v.gety() * v.gety();
  ////	}

  ////Update:
  //stress_s.set(sxx_s,sxy_s,syx_s,syy_s);
  //stress_c.set(sxx_c,sxy_c,sxy_c,syy_c);
  ////Overload division by double for Tensor2x2
  ////Divide or not by volume here???
  //stress_s = stress_s * (1. / cell_->getVolume());
  //stress_c = stress_c * (1. / cell_->getVolume());

  ////Total stress:
  //stress_ = stress_s + stress_c;
}


vector<double> Interactions::getAverageMaxPenetration()const{

  double dn_average = 0. ;
  double dn_max = 0. ;

  for(vector<Contact*>::const_iterator it = clist_.begin() ; it != clist_.end() ; it++){
    double dn = (*it)->getdn();
    dn_average += fabs(dn);
    dn_max = max(fabs(dn),fabs(dn_max));
  }
  if(clist_.size()!=0) dn_average /= (double)clist_.size();

  vector<double> dns;
  dns.push_back(dn_average);
  dns.push_back(dn_max);
  return dns;
}


void Interactions::debug(const int k) const{

}



