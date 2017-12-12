#include"Interactions.hpp"
#include"Sample.hpp"
#include"Cell.hpp"
#include"Particle.hpp"

using namespace std;

Interactions::Interactions(){
	//Temp
	//Distances a mettre en Rmax
	dv_ = 0. ;
	dsv_ = 0. ;
	nv_ = 0 ;
	nsv_ = 0;
	scale_=string();
	initScale_ = false;
	initdv_ = false;
	initdsv_ = false;
	checkInteractions_ = false;
	folder_ = string();
	fInteractions_ = "inter.txt";
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
		if(token=="}") break;
		is >> token;
	}


	//TMP:
	kn_ = 1.e7 ;
	kt_ = 1.e7 ;
	gn_ = 0. ;
	gt_ = 0. ;
	mus_ = 0. ;
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


bool Interactions::initcheck(){
	bool initNiter = false;
	if( nv_ != 0 && nsv_ != 0 ) initNiter = true; 
	checkInteractions_ = initNiter && initScale_ && initdv_ && initdsv_;


	cout<<"d verlet (scale)  = "<<dv_<<endl;
	cout<<"d sverlet (scale)  = "<<dsv_<<endl;
	return checkInteractions_;
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

	cout<<"Super verlet list size: "<<svlist_.size()<<endl;
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

	cout<<"Verlet list size: "<<vlist_.size()<<endl;
}

//Build contact list (activated interactions)
void Interactions::detectContacts(){

	cout<<"Contact detection..."<<endl;
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

void Interactions::computeForces(){

	for(vector<int>::iterator it = clist_.begin(); it != clist_.end();it++){
		vlist_[*it].updateRelativeVelocity();
		vlist_[*it].computeForce(kn_,kt_,gn_,gt_,mus_);
		vlist_[*it].updateAccelerations();
	}

}


void Interactions::writeContacts(int k) const {

	string filename = formatfile(folder_, fInteractions_, k);
	ofstream file(filename.c_str());

	for(vector<int>::const_iterator it = clist_.begin(); it != clist_.end(); it++){
		vlist_[*it].write(file);
	}
}
