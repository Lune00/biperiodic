#include"Interactions.hpp"
#include"Sample.hpp"
#include"Cell.hpp"
#include"Particle.hpp"

using namespace std;

Interactions::Interactions(){

	dv_ = 0. ;
	dsv_ = 0.;

	nsv_ = 0;
	nv_ = 0 ;

	scale_=string();
	folder_ = string();

	initScale_ = false;
	fInteractions_ = "inter.txt";

	initdv_ = false;
	initdsv_ = false;

	checkInteractions_ = false;

	initkn_ = false;
	initkt_ = false;
	initgn_ = false;
	initgt_ = false;
	initmus_ = false;
	setgnmax_ = false;
	setgtmax_ = false;

	ofstream debug("tmp.txt");
	debug.close();
}

Interactions::~Interactions(){

  delete [] dts_ ;
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
		if(token=="dsv"){
			is >> dsv_;
			initdsv_ = true;
		}
		if(token=="niterv"){
		  is >> nv_;
		}
		if(token=="nitersv"){
		  is >> nsv_;
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
		if(token=="gtmax"){
		  setgtmax_ = true ;
		  initgt_ = true ;
		}
		if(token=="mu"){
			is >> mus_;
			initmus_ = true;
		}

		if(token=="}") break;
		is >> token;
	}
}

//For post-processing only:
void Interactions::setparameters(double kn,double kt,double gn,double gt, double mu){
	kn_ = kn ;
	kt_ = kt ;
	gn_ = gn ;
	gt_ = gt ;
	mus_ = mu ;
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
	dsv_ *= scale;
	initScale_ = true ;

	return ;
}

bool Interactions::initcheck() {

  bool initNiter = false;

  if( nv_ != 0 && nsv_ != 0 ) initNiter = true; 

  checkInteractions_ = initNiter && initScale_ && initdv_ && initdv_ && (dsv_>dv_) && checkDEMparameters();

  cout<<"distance SVerlet  = "<<dsv_<<endl;
  cout<<"distance Verlet  = "<<dv_<<endl;

  if(!checkDEMparameters()) cout<<"Interactions::initcheck() : a DEM parameter is not initialised in config file. stop"<<endl;

  return (checkInteractions_);
}

void Interactions::plug(Sample& spl,Cell& cell){
	spl_ = & spl;
	cell_ = &cell;
}

void Interactions::build(){

  init_array_dt();

  //Est ce qu'on load? on lit les contacts dans un fichier
  if(spl_->loaded()){
    unsigned int filetoload = spl_->filetoload();
    load(filetoload);
  }

  initScale();
}

//On suppose que les Id des particules vont de 0 à N-1
void Interactions::init_array_dt(){

	N_ = spl_->getsize();

	dts_ = new double [ N_ * N_ ] ;

	for(unsigned int i = 0 ; i < N_ ; i++){
		for(unsigned int j = 0 ; j < N_; j++){
			dts_[i * N_ + j ] = 0. ;
		}
	}
}

//Load the network for continuing simulation
//Only needs to restore dt (tangential spring)
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
  else  {
    read_dt(is);
  }
  is.close();
}

//Read contact network for continuing simulation: only need dt
void Interactions::read_dt(ifstream& is){

    string token;
    while(is){
      if(is.eof()) break;

      int idi=0;
      int idj=0;
      double vn;
      double vt;
      double nx;
      double ny;
      double fn;
      double ft;
      double dt;

      is >> idi >> idj >> vn >> vt >> nx >> ny >> fn >> ft>> dt;

      if( idi != idj){
      //particle_pair i__j = { spl_->getP(idi), spl_->getP(idj)};
      dts_[ idi * N_ + idj ] = dt;
      //svlist_.push_back(i__j);
      }

    }
}


//Write contact network
void Interactions::writeContacts(int k) const {

	string filename = formatfile(folder_, fInteractions_, k);
	ofstream file(filename.c_str());
	file.precision(12);

	for(vector<int>::const_iterator it = clist_.begin(); it != clist_.end();it++){
		int k = *it;
		vlist_[k].write(file);
	}
	file.close();
}

bool Interactions::checkDEMparameters() const{
	return (initkn_ && initkt_ && initgn_ && initgt_ && initmus_);
}

void Interactions::updateverlet(const int tic){
	if( tic % nsv_ == 0 ) updatesvlist();
	if( tic % nv_ == 0 ) updatevlist();
}

//Return smallest vector branch between particle i and particle j
Vecteur Interactions::getShortestBranch(const Particle& i, const Particle& j) const{

	//Need the branch vector (in absolute units)
	double sijx = j.getx() - i.getx();
	double sijy = j.gety() - i.gety();

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
	//double dmin = * it ;
	int k = distance( l_dcarre.begin(), it);
	//Get matching pairs of indexes
	pair<int,int> indexes = pairs[k];
	//Return shortest vector branch:
	return (d + a0 * indexes.first + a1 * indexes.second) ; 
}


//True if distance between "surface" of particle i and j are lower than d
bool Interactions::near(const Particle& i, const Particle& j,const double d) const{

	Vecteur shortest_branch = getShortestBranch(i,j);

	//Test distance compared to d

	if( shortest_branch.getNorme() - d < j.getRadius() + i.getRadius() ) return true;

	else return false;
}

void Interactions::updatesvlist(){

	svlist_.clear();

	for(std::vector<Particle>::iterator iti = spl_->getSample()->begin();iti!=spl_->getSample()->end();iti++)
	{
		for(std::vector<Particle>::iterator itj = iti + 1; itj!=spl_->getSample()->end();itj++)
		{
			if( near( *iti , *itj , dsv_ ) ){
				particle_pair i__j = { &(*iti), &(*itj)};
				svlist_.push_back(i__j);
			}
		}
	}

}

void Interactions::updatevlist(){

	vlist_.clear();

	for(vector<particle_pair>::iterator it = svlist_.begin() ; it != svlist_.end(); it ++ )
	{
		if( near( *(it->i), *(it->j), dv_ ) ){
			Contact c(it->i, it->j,cell_);
			//Attribue dt
			c.setdt(get_dt(c));
			vlist_.push_back(c);
		}
	}
}


//Build contact list (activated interactions)
void Interactions::detectContacts(){

	clist_.clear();

	for( vector<Contact>::iterator it = vlist_.begin() ; it != vlist_.end(); it++){

		it->Frame();
		//cerr<<"Frame : "<<it->geti()->getId()<<" "<<it->getj()->getId()<<endl;

		int k = distance (vlist_.begin(), it );

		if( it->isActif() ) {
			clist_.push_back(k);
		}
		else{
			//dt_ is set to zero in the contact
			set_dt(*it);
		}

	}
	return ;

}

double Interactions::get_dt(Contact& c) const{
	int i = c.geti()->getId();
	int j = c.getj()->getId();
	return dts_[ i * N_ + j];
}

void Interactions::set_dt(Contact& c){
	int i = c.geti()->getId();
	int j = c.getj()->getId();
	dts_[ i * N_ + j ] = c.getdt() ;
}

//Not used
void Interactions::reset_dt(Contact& c){
	int i = c.geti()->getId();
	int j = c.getj()->getId();
	dts_[ i * N_ + j ] = 0. ;
}

//Compute force at each contact and compute procedurally internal stress on the fly
void Interactions::computeForces(const double dt){

	//Static stress:
	double sxx_s = 0. ;
	double sxy_s = 0. ;
	double syx_s = 0. ;
	double syy_s = 0. ;

	//Kinetic stress:
	double sxx_c = 0. ;
	double sxy_c = 0. ;
	double syy_c = 0. ;

	for(vector<int>::iterator it = clist_.begin() ; it != clist_.end() ; it++){

		vlist_[*it].updateRelativeVelocity();
		vlist_[*it].computeForce(kn_,kt_,gn_,gt_,mus_,dt);
		//Update array table:
		set_dt(vlist_[*it]);

		vlist_[*it].updateAccelerations();

		Vecteur branch = vlist_[*it].getbranch();
		Vecteur force  = vlist_[*it].getfxy();

		sxx_s += branch.getx() * force.getx();
		sxy_s += branch.gety() * force.getx();
		syx_s += branch.getx() * force.gety();
		syy_s += branch.gety() * force.gety();
	}


	// Transform acceleration into reduced coordinates:
	Tensor2x2 hdd = cell_->gethdd();
	Tensor2x2 hd = cell_->gethd();
	Tensor2x2 h = cell_->geth();
	Tensor2x2 hinv = h.getInverse();

	for(std::vector<Particle>::iterator it = spl_->getSample()->begin();it!=spl_->getSample()->end();it++)
	{
		//Acceleration from contact force:
		Vecteur a_red = hinv * ( it->getA() - hd * (it)->getV() * 2. - hdd * (it)->getR());

		it->setAcceleration(a_red);

		//Acceleration from external drive:
		if(cell_->imposeForce()) addForce(*it);

		//CALCUL DU TENSEUR CONTRAINTES CINEMATIQUES
		//MAIS pour ca on a besoin de la vitesse au debut, au milieu ou a la fin du pas de temps?
		//On peut prendre celle au milieu, ca ne devrait pas changer grand chose...
		//Kinetic stress : Loop over particles:

		//for(std::vector<Particle>::const_iterator it = spl_->inspectSample().begin();it!=spl_->inspectSample().end();it++)
		//{
		//Partie fluctuante : v = h * sdot
		Vecteur v = cell_->geth() * it->getV();
		double m = it->getMasse();
		sxx_c += m * v.getx() * v.getx();
		sxy_c += m * v.getx() * v.gety();
		syy_c += m * v.gety() * v.gety();
		//	}

	}

	//Update:
	stress_s.set(sxx_s,sxy_s,syx_s,syy_s);
	stress_c.set(sxx_c,sxy_c,sxy_c,syy_c);

	//cerr<<"sxx_s = "<<sxx_s<<" sxx_c = "<<sxx_c<<endl;

	stress_s = stress_s * (1. / cell_->getVolume());
	stress_c = stress_c * (1. / cell_->getVolume());

	//Total stress:
	stress_ = stress_s + stress_c;

}

//Add a force to each particule in the horizontal direction
//fx = A * sin (2pi y/Ly * mode)
void Interactions::addForce(Particle& p){

	double A = cell_->getAmplitudeForce();
	int mode = cell_->getModeForce();
	//double Ly = cell_->geth().getyy();

	double yLy = p.getR().gety() ;
	double fx = A * sin ( 2. * M_PI * yLy * (double)mode);
	Vecteur f(fx,0.);
	//Warning: j'impose une acceleration direct en coord reduites
	//Peut etre imposer plutot en coord absolu et transformer
	//pour avoir une meilleur idee de l'echelle de force a mettre
	p.updateA(f);
	//cerr<<"Particule "<<p.getId()<<": "<<f.getx()<<" y/Ly = "<<yLy<<endl;
}

void Interactions::askNumberOfContacts() const{
	cerr<<"Nombre de contacts: "<<clist_.size()<<endl;
}



double Interactions::getElasticEnergy() const {

	double E = 0.;

	//Tangential elastic component kt*dt*dt?
	for(vector<int>::const_iterator it = clist_.begin() ; it != clist_.end() ; it++){

		double dn = vlist_[*it].getdn();
		E += 0.5 * (dn) * (dn) * kn_ ;
	}

	return E;
}

//Only call for post-processing
void Interactions::computeInternalStress(){

	//Static stress:
	double sxx_s = 0. ;
	double sxy_s = 0. ;
	double syx_s = 0. ;
	double syy_s = 0. ;

	//Kinetic stress:
	double sxx_c = 0. ;
	double sxy_c = 0. ;
	double syy_c = 0. ;

	for(vector<int>::iterator it = clist_.begin() ; it != clist_.end() ; it++){
		Vecteur branch = vlist_[*it].getbranch();
		Vecteur force  = vlist_[*it].getfxy();

		sxx_s += branch.getx() * force.getx();
		sxy_s += branch.gety() * force.getx();
		syx_s += branch.getx() * force.gety();
		syy_s += branch.gety() * force.gety();
	}

	for(std::vector<Particle>::const_iterator it = spl_->inspectSample().begin();it!=spl_->inspectSample().end();it++)
	{
		//Partie fluctuante : v = h * sdot
		Vecteur v = cell_->geth() * it->getV();
		double m = it->getMasse();
		sxx_c += m * v.getx() * v.getx();
		sxy_c += m * v.getx() * v.gety();
		syy_c += m * v.gety() * v.gety();
	}


	//Update:
	stress_s.set(sxx_s,sxy_s,syx_s,syy_s);
	stress_c.set(sxx_c,sxy_c,sxy_c,syy_c);


	stress_s = stress_s * (1. / cell_->getVolume());
	stress_c = stress_c * (1. / cell_->getVolume());

	//Total stress:
	stress_ = stress_s + stress_c;
}


vector<double> Interactions::getAverageMaxPenetration()const{

	double dn_average = 0. ;
	double dn_max = 0. ;

	for(vector<int>::const_iterator it = clist_.begin() ; it != clist_.end() ; it++){
		double dn = vlist_[*it].getdn();
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
	ofstream tmp("tmp.txt",ios::app);

	int idi = 20 ;
	int idj = 24 ;
	int a = idi * N_ + idj ;

	for(vector<int>::const_iterator it = clist_.begin(); it != clist_.end(); it++){
		if( vlist_[*it].geti()->getId() == idi && vlist_[*it].getj()->getId() == idj ) {
			double dt = vlist_[*it].getdt();
			tmp<<k<<" "<<dts_[a]<<" "<<dt<<" "<<vlist_[*it].getdn()<<endl;
			break;
		}
	}

	tmp.close();

}

void Interactions::print() const {
	cout<<"Nombre d'interactions (Super Verlet): "<<svlist_.size()<<endl;
	cout<<"Nombre d'interactions (Verlet): "<<vlist_.size()<<endl;
	cout<<"Nombre de contacts actifs: "<<clist_.size()<<endl;
}


//Load the network for postprocessing only
void Interactions::loadnetwork(const int k){


	vlist_.clear();
	clist_.clear();

	string filename = formatfile( folder_, fInteractions_, k );
	ifstream is(filename.c_str());
	if(!is){
		cerr<<"Interactions::load "<<filename<<" fail."<<endl;
		return ;
	}
	else if(is.peek() == std::ifstream::traits_type::eof()){
		cerr<<"There are no contacts to load."<<endl;
		return ;
	}
	else  {
		read_contact(is);
	}
	is.close();
}


//Read contact network: build the contact list ONLY for post-processing
void Interactions::read_contact(ifstream& is){

	string token;
	while(is){
		if(is.eof()) break;

		int idi=0;
		int idj=0;
		double vn;
		double vt;
		double nx;
		double ny;
		double fn;
		double ft;
		double dt;

		is >> idi >> idj >> vn >> vt >> nx >> ny >> fn >> ft>> dt;

		Vecteur v(vn,vt);
		Vecteur f(fn,ft);
		Vecteur n(nx,ny);

		//Build verlet list
		//Si trop long d'appeler Frame() on peut ecrire tout a la place
		if( idi != idj){
			//Set properties of the contact:
			Contact c( spl_->getP(idi), spl_->getP(idj),cell_ );
			//Set branch, n, t, dn,dt
			c.Frame();
			c.setrv(v);
			c.setf(f);
			c.setdt(dt);
			vlist_.push_back(c);
		}

	}
	//Build clist:
	for( vector<Contact>::iterator it = vlist_.begin() ; it != vlist_.end(); it++){
		int k = distance (vlist_.begin(), it );
		clist_.push_back(k);
		//cerr<<vlist_[k].getfn()<<" "<<vlist_[k].getdn()<<endl;
	}

	cerr<<"Nombre de contacts : "<<clist_.size()<<endl;
}

