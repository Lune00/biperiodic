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

	ofstream debug("internalstress.txt");
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
}


bool Interactions::checkDEMparameters() const{
	return (initkt_ && initkt_ && initgn_ && initgt_ && initmus_);
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

	svlist_.clear();

	//If less than two particles, no interaction possible
	if(spl_->getsize() < 2 ) return ;
	

	vector<Particle>* ps = spl_->getSample();
	//Get h:
	Tensor2x2 h = spl_->getCell()->geth();

	//Const iterator? But how, after we need these pointers to modify
	//Particles properties... 
	for(std::vector<Particle>::iterator iti = ps->begin(); iti!=ps->end();iti++){
		for(std::vector<Particle>::iterator itj = iti+1; itj!=ps->end();itj++){
			if( near( *(iti), *(itj), h , dsv_) ) {
				//cout<<"Near!"<<endl;
				particle_pair O__O = { &(*iti), &(*itj) };
				svlist_.push_back(O__O);
			}
		}
	}

	//cout<<"Super verlet list size: "<<svlist_.size()<<endl;
}

//True if distance between "surface" of particle i and j are lower than d
bool Interactions::near(const Particle& i, const Particle& j,const Tensor2x2& h,const double d) const{

	double sijx = j.getx() - i.getx();
	double sijy = j.gety() - i.gety();
	//Shortest branch through periodicity:
	sijx -= floor(sijx + 0.5);
	sijy -= floor(sijy + 0.5);

	Vecteur sij(sijx,sijy);
	//Absolute vector branch:
	sij = h * sij ;
	//Test distance compared to dsv
	if( sij.getNorme() - d < j.getRadius() + i.getRadius() ) return true;
	else return false;

}

void Interactions::updatevlist(){

	vlist_.clear();

	Tensor2x2 h = spl_->getCell()->geth();

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

//Build contact list (activated interactions)
void Interactions::detectContacts(){

	//cout<<"Contact detection..."<<endl;
	clist_.clear();

	for(vector<Contact>::iterator it = vlist_.begin(); it != vlist_.end(); it++){
		it->Frame();
		int k = distance( vlist_.begin(), it);

		//WIP
		//TEMPORAIRE!!! TEST SUR LISTEVERLET
		//clist_.push_back(k);

		if(it->isActif()){
			//cout<<"Le contact "<<k<<" est actif."<<endl;
			clist_.push_back(k);
		}
		//TEMP
		else{
			//cout<<"Le contact "<<k<<" est inactif."<<endl;
		}
	}


}

void Interactions::computeForces(const double dt){

	for(vector<int>::iterator it = clist_.begin(); it != clist_.end();it++){
		vlist_[*it].updateRelativeVelocity();
		vlist_[*it].computeForce(kn_,kt_,gn_,gt_,mus_,dt);
		vlist_[*it].updateAccelerations();
	}

}


void Interactions::askNumberOfContacts() const{
	cerr<<"Nombre de contacts: "<<clist_.size()<<endl;
}

void Interactions::writeContacts(int k) const {

	//if(clist_.size()==0) cerr<<"step "<<k<<" Il n'y a pas de contact"<<endl;
	string filename = formatfile(folder_, fInteractions_, k);
	ofstream file(filename.c_str());

	for(vector<int>::const_iterator it = clist_.begin(); it != clist_.end(); it++){
		vlist_[*it].write(file);
		vlist_[*it].print();
	}
}

double Interactions::getElasticEnergy() const {

	double E = 0.;
	for(vector<int>::const_iterator it = clist_.begin(); it != clist_.end(); it++){
		double dn = vlist_[*it].getdn();
		double mj = vlist_[*it].getj()->getMasse();
		double mi = vlist_[*it].geti()->getMasse();
		E += 0.5 * (dn) * (dn) * kn_  ;
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

	//Need to recaculculate the branch vector
	//Maybe store it in the contact

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
	//stress_s = stress_s * (1. / cell_->getVolume());
	//stress_c = stress_c * (1. / cell_->getVolume());
	
	//Total stress:
	stress_ = stress_s + stress_c;
}

void Interactions::debug(const int k) const{

	ofstream os("internalstress.txt",ios::app);
	//DEBUG
	os<<k<<" "<<stress_s.getxx()<<" "<<stress_s.getyy()<<endl;
	os.close();


}
