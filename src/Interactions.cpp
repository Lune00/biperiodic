#include"Interactions.hpp"
#include"Sample.hpp"
#include"Cell.hpp"
#include"Particle.hpp"

using namespace std;

Interactions::Interactions(){
	dv_ = 0. ;
	dsv_ = 0. ;
	nv_ = 0 ;
	nsv_ = 0;
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

	ofstream debug("debugInteractions.txt");
	debug.close();
}

Interactions::~Interactions(){

	delete [] array_dt ;
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
		if(token=="niterv")  is >> nv_;
		if(token=="nitersv") is >> nsv_;

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
	//Default:
	//kn_ = 1.e7 ;
	//kt_ = 1.e7 ;
	//gn_ = 0. ;
	//gt_ = 20. ;
	//mus_ = 0.5 ;
}


void Interactions::initScale(){

	double scale;
	double epsilon=0.00001;

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

	//Check if dv et dsv have been initialised "correctly"
	if(dv_ < epsilon || dsv_ < epsilon)
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

void Interactions::plug(Sample& spl,Cell& cell){
	spl_ = & spl;
	//Allow accessing h at any time for force computations
	cell_ = &cell;
	initScale();
	//Store tangantial displacement;
	init_array_dt();
}

//Should be N * (N - 1) but more complicated no?
//WARNING: on suppose que l'id des particules va de 0 a N-1
void Interactions::init_array_dt(){
	N_ = spl_->getsize();
	array_dt = new double [ N_ * N_ ] ;
}


bool Interactions::checkDEMparameters() const{
	return (initkn_ && initkt_ && initgn_ && initgt_ && initmus_);
}

bool Interactions::initcheck() {
	bool initNiter = false;
	if( nv_ != 0 && nsv_ != 0 ) initNiter = true; 
	checkInteractions_ = initNiter && initScale_ && initdv_ && initdsv_ && checkDEMparameters();
	cout<<"d verlet (scale)  = "<<dv_<<endl;
	cout<<"d sverlet (scale)  = "<<dsv_<<endl;
	if(!checkDEMparameters()) cout<<"Interactions::initcheck() : a DEM parameter is not initialised in config file. stop"<<endl;
	return (checkInteractions_);
}
	

void Interactions::updateverlet(const int tic){
	if( tic % nsv_ == 0 ) updatesvlist();
	if( tic % nv_ == 0 ) updatevlist();
}


//On pourrait stocker la distance interbranche dans la liste de superverlet
//au lieu de la recalculer, et dans verlet tester la distance au cut off plus court, optimisation legere...
void Interactions::updatesvlist(){

//	cout<<"UPDATE VLIST"<<endl;
	svlist_.clear();

	//If less than two particles, no interaction possible
	if(spl_->getsize() < 2 ) return ;

	//TMP, a terme la vlist sera passée a sample plutot
	vector<Particle> * ps = spl_->getSample();
	//Get h:
	Tensor2x2 h = cell_->geth();

	//Const iterator? But how, after we need these pointers to modify
	for(std::vector<Particle>::iterator iti = ps->begin(); iti!=ps->end();iti++){
		for(std::vector<Particle>::iterator itj = iti+1; itj!=ps->end();itj++){
			if( near( *(iti), *(itj), h , dsv_) ) {
				particle_pair O__O = { &(*iti), &(*itj) };
				svlist_.push_back(O__O);
			}
		}
	}

	//cout<<"Super verlet list size: "<<svlist_.size()<<endl;
}

//Return the pair of indices for which the branch vector between particle i and particle j is minimum (first image algo for any cell shape)
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
			//double dcarre = d * d + d * a0 * 2 * i + a0 * a0 * i * i +  d * a1 * 2 * j + a0 * a1 * 2 * i * j + a1 * a1 * j * j ;
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
bool Interactions::near(const Particle& i, const Particle& j,const Tensor2x2& h,const double d) const{
	Vecteur shortest_branch = getShortestBranch(i,j);
	//Test distance compared to dsv
	if( shortest_branch.getNorme() - d < j.getRadius() + i.getRadius() ) return true;
	else return false;
}

void Interactions::updatevlist(){

	vlist_.clear();
	Tensor2x2 h = cell_->geth();

	for(vector<particle_pair>::iterator it = svlist_.begin(); it != svlist_.end(); it++){
		//if near dverlet
		if( near(*(it->i),*(it->j),h,dv_) )
		{
		//Create a Contact
		Contact c(it->i,it->j,*cell_);	
		//Push contact to vlist_;
		vlist_.push_back(c);
		}
	}

	//cout<<"Verlet list size: "<<vlist_.size()<<endl;
}

//Original, no memory procedure
//void Interactions::detectContacts(){
//	clist_.clear();
//	for(vector<Contact>::iterator it = vlist_.begin(); it != vlist_.end(); it++){
//		it->Frame();
//		int k = distance( vlist_.begin(), it);
//		if(it->isActif()){
//			clist_.push_back(k);
//		}
//	}
//}

//Build contact list (activated interactions)
//On doit renouveler la liste en gardant en memoire les contacts deja actifs (pas suprrimer et refaire comme avant)...
//Ici on essaye en gardant juste le dt stocker dans un grand tableau
void Interactions::detectContacts(){

	clist_.clear();

	for(vector<Contact>::iterator it = vlist_.begin(); it != vlist_.end(); it++){

		it->Frame();
		int k = distance( vlist_.begin(), it);

		if(it->isActif()){
			//On set le dt
	//		cerr<<"Contact entre particule "<<it->geti()->getId()<<" et "<<it->getj()->getId()<<" dt = "<<get_dt(*it)<<endl;
			it->set_dt(get_dt(*it));
			clist_.push_back(k);
		}
		else{
			//On me a zero le dt
			//Deja fait si le contact n'est pas actif
			//set_dt(0.,*it);
		}
	}
}

double Interactions::get_dt(Contact& c) const{
	int i = c.geti()->getId();
	int j = c.getj()->getId();
	return array_dt[ i * N_ + j];
}

void Interactions::set_dt(Contact& c){
	int i = c.geti()->getId();
	int j = c.getj()->getId();
	array_dt[ i * N_ + j ] = c.getdt() ;
}

void Interactions::computeForces(const double dt){

	for(vector<int>::iterator it = clist_.begin(); it != clist_.end();it++){
		vlist_[*it].updateRelativeVelocity();

		vlist_[*it].computeForce(kn_,kt_,gn_,gt_,mus_,dt);
		vlist_[*it].updateAccelerations();
		//Update dt table
		set_dt(vlist_[*it]);
		//ofstream fo("debugInteractions.txt",ios::app);
		//fo<< vlist_[*it].getfn()<<" "<<vlist_[*it].getft()<<" "<<vlist_[*it].getdn()<<" "<<vlist_[*it].getdt()<<" "<<vlist_[*it].getrv().gety()<<endl;
		//fo.close();
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

	for(vector<int>::const_iterator it = clist_.begin(); it != clist_.end(); it++){
		vlist_[*it].write(file);
	}
}

double Interactions::getElasticEnergy() const {

	double E = 0.;
	for(vector<int>::const_iterator it = clist_.begin(); it != clist_.end(); it++){
		double dn = vlist_[*it].getdn();
		double mj = vlist_[*it].getj()->getMasse();
		double mi = vlist_[*it].geti()->getMasse();
		double dt = vlist_[*it].getdt();
		E += 0.5 * (dn) * (dn) * kn_ ;
			//+ 0.5 * (dt) * (dt) * kt_  ;
	}
	return E;
}

//Compute static stress & kinetic stress
//Add them to compute totalstress stress_
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
	for(vector<int>::const_iterator it = clist_.begin(); it != clist_.end(); it++){ 
		Vecteur branch = vlist_[*it].getbranch();
		Vecteur force = vlist_[*it].getfxy();
		sxx_s += branch.getx() * force.getx();
		sxy_s += branch.gety() * force.getx();
		syx_s += branch.getx() * force.gety();
		syy_s += branch.gety() * force.gety();
	}

	//Kinetic stress : Loop over particles:
	for(std::vector<Particle>::const_iterator it = spl_->inspectSample().begin();it!=spl_->inspectSample().end();it++)
	{
		Vecteur v = it->getV();
		double m = it->getMasse();
		sxx_c += m * v.getx() * v.getx();
		sxy_c += m * v.getx() * v.gety();
		syy_c += m * v.gety() * v.gety();
	}

	//Update:
	stress_s.set(sxx_s,sxy_s,syx_s,syy_s);
	stress_c.set(sxx_c,sxy_c,sxy_c,syy_c);
	//Overload division by double for Tensor2x2
	//Divide or not by volume here???
	stress_s = stress_s * (1. / cell_->getVolume());
	stress_c = stress_c * (1. / cell_->getVolume());
	
	//Total stress:
	stress_ = stress_s + stress_c;
}


void Interactions::debug(const int k) const{
	double dnaverage = 0. ;
	double dnmax = 0. ;
	cout<<"ncontacts = "<<clist_.size()<<endl;
	ofstream os("debugInteractions.txt",ios::app);
	for(vector<int>::const_iterator it = clist_.begin(); it != clist_.end(); it++){
		dnaverage += fabs(vlist_[*it].getdn());
		dnmax = max(fabs(vlist_[*it].getdn()),fabs(dnmax));
		//os<<vlist_[*it].getdn()<<" "<<vlist_[*it].geti()->getId()<<" "<<vlist_[*it].getj()->getId()<<" "<<vlist_[*it].getfn()<<" "<<vlist_[*it].getbranch().getNorme()<<endl;
	}
	if(clist_.size()!=0) dnaverage /= (double)clist_.size();
	os<<k<<" "<<dnaverage<<" "<<dnmax<<" "<<stress_s.getyy()<<" "<<stress_s.getxx()<<endl;
}
