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
			dts_[ idi * N_ + idj ] = dt;
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

void Interactions::updatesvlist(){

	svlist_.clear();
	Tensor2x2 h = cell_->geth();

	for(std::vector<Particle>::iterator iti = spl_->getSample()->begin();iti!=spl_->getSample()->end();iti++)
	{
		double six = iti->getx();
		double siy = iti->gety();
		double ri = iti->getRadius();

		for(std::vector<Particle>::iterator itj = iti + 1; itj!=spl_->getSample()->end();itj++)
		{
			double sijx =itj->getx() - six;
			double sijy =itj->gety() - siy;

			sijx -= floor(sijx + 0.5);
			sijy -= floor(sijy + 0.5);

			Vecteur sij(sijx,sijy);
			Vecteur d = h * sij;

			if( d.getNorme() - dsv_ < ri + itj->getRadius() ){
				particle_pair i__j = { &(*iti), &(*itj)};
				svlist_.push_back(i__j);
			}
		}
	}

}

void Interactions::updatevlist(){

	vlist_.clear();

	Tensor2x2 h = cell_->geth();

	for(vector<particle_pair>::iterator it = svlist_.begin() ; it != svlist_.end(); it ++ )
	{
		double sijx = it->j->getx() - it->i->getx() ;
		double sijy = it->j->gety() - it->i->gety() ;

		pair<int,int> indexes;
		indexes.first = -(int)floor(sijx + 0.5) ;
		indexes.second = -(int)floor(sijy + 0.5) ;
		sijx += indexes.first;
		sijy += indexes.second;

		Vecteur sij(sijx,sijy);
		Vecteur d = h * sij;

		if( d.getNorme() - dv_ < it->i->getRadius() + it->j->getRadius() ){
			Contact c(it->i, it->j,cell_,indexes);
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

//Compute force at each contact and compute procedurally internal stress on the fly
void Interactions::computeForces(const double dt,const double t){

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
		//Impose a real force
		if(cell_->imposeForce()){
			addForce(*it,t);
		}
		//Acceleration from contact force:
		//Vecteur a_red = hinv * ( it->getA() - hd * (it)->getV() * 2. - hdd * (it)->getR());
		Vecteur a_red = hinv * ( it->getA()) ; // - hd * (it)->getV() * 2. - hdd * (it)->getR());
		//double a = (cell_->getStressExt()+stress_s).getyy()*cell_->geth().getxx()/cell_->getMasse();
		//double b = (hdd * (it)->getR()).gety();
		//double c = b/a;
		//ofstream os("terms.txt",ios::app);
		//os<<it->getA().getNorme()<<" "<<a<<" "<<b<<" "<<c<<endl;
		//os.close();

		it->setAcceleration(a_red);

		//CALCUL DU TENSEUR CONTRAINTES CINEMATIQUES
		//MAIS pour ca on a besoin de la vitesse au debut, au milieu ou a la fin du pas de temps?
		//On peut prendre celle au milieu, ca ne devrait pas changer grand chose...
		//Kinetic stress : Loop over particles:

		//TMP: A AJOUTER

		//Partie fluctuante : v = h * sdot
		//Vecteur v = cell_->geth() * it->getV();
		//double m = it->getMasse();
		//sxx_c += m * v.getx() * v.getx();
		//sxy_c += m * v.getx() * v.gety();
		//syy_c += m * v.gety() * v.gety();

	}

	//Update:
	stress_s.set(sxx_s,sxy_s,syx_s,syy_s);
	stress_c.set(sxx_c,sxy_c,sxy_c,syy_c);

	stress_s = stress_s * (1. / cell_->getVolume());
	//stress_c = stress_c * (1. / cell_->getVolume());

	//Total stress:
	stress_ = stress_s + stress_c;

}

//Add a force to each particule in the horizontal direction
//fx = A * sin (2pi y/Ly * mode)
void Interactions::addForce(Particle& p,const double t){

	double A = cell_->getAmplitudeForce();
	int mode = cell_->getModeForce();
	double yLy = p.getR().gety() ;
	Vecteur v = cell_->geth() * p.getV();
	double vx = v.getx();
	//draging force prop to a fictuous sinusoidal velocity field
	double Vxy = A * sin ( 2. * M_PI * yLy * (double)mode);
	double fx = 0.5 * (Vxy-vx);
	Vecteur f(fx,0.);
	p.add_force(f);
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
			c.computeShortestBranch();
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
	}

	cerr<<"Nombre de contacts : "<<clist_.size()<<endl;
}

