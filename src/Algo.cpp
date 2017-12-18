#include"Algo.hpp"
#include"Cell.hpp"
#include"Sample.hpp"
#include"Interactions.hpp"
#include"Analyse.hpp"

using namespace std;

void Algo::init(ifstream& is){
	string token;
	is >> token;
	while(is){
		if(token=="dt") is >> dt_;
		if(token=="ns") is >> ns_;
		if(token=="nana") is >> nana_;
		if(token=="nrecord") is >> nrecord_;
		if(token=="}") break;
		is >> token;
	}
}

//A adapter comme on souhaite en terme de tests
bool Algo::initcheck(){
	computedtmax();
	compute_gnmax_restitution();
	if( t_ > 0. ) return false;
	if( !checkSimulationParameters() ) return false;
	if (dt_ < 1. && ns_ != 0 && nrecord_ != 0 ) return true;
	else return false;
}

void Algo::plug(Cell& cell, Sample& spl, Interactions& Int, Analyse& ana){
	Int_ = &Int;
	cell_ = &cell;
	spl_ = &spl;
	ana_ = &ana;
}

//Compute dtmax for the simulation (kn,m)
void Algo::computedtmax(){
	double const epsilon = 0.01;
	double const rho = spl_->getrho();
	double const kn = Int_->getkn();
	double m = rho * M_PI * spl_->getrmin() * spl_->getrmin(); 
	dtmax_ = epsilon * sqrt( m / kn ); 
}


//Compare dt to dtmax
bool Algo::checktimestep()const{
	cout<<"dtmax = "<<dtmax_<<endl;
	if(dt_>dtmax_) {
		//TMP
		cerr<<"Choose a lower dt."<<endl;
		return false;
	}
	else return true;
}

//Compute gnmax and restitution coeffcient
void Algo::compute_gnmax_restitution(){
	const double gn = Int_->getgn();
	double rho = spl_->getrho();
	double m = rho * M_PI * spl_->getrmin() * spl_->getrmin(); 
	double kn = Int_->getkn();
	gnmax_ = 2. * sqrt(kn * m);
	double ds = gn/gnmax_;
	//Compute resitution coefficient:
	e_ = exp(- M_PI * ds/(2.*sqrt(1.-ds*ds)) ) ;
}

bool Algo:: checkNormalViscosity()const{
	cout<<"Max normal viscosity gnmax = "<<gnmax_<<endl;
	cout<<"Restitution coefficient e = "<<e_<<endl;
	if(Int_->getgn()>gnmax_){
		cerr<<"Keeping all other parameters constants, choose a lower normal viscosity"<<endl;
		return false;
	}
	else return true;
}


//Check, according DEM parameters, if dt is set correctly
//Tangential viscosity check expression...???
bool Algo::checkSimulationParameters(){
	bool is_dt_ok = checktimestep();
	bool is_gn_ok = checkNormalViscosity();
	return (is_dt_ok && is_gn_ok);
}


void Algo::run(){

	double tfinal = ns_ * dt_ ;

	cout<<"Simulation:"<<endl;
	cout<<"dt = "<<dt_<<endl;
	cout<<"ns = "<<ns_<<endl;
	cout<<"tfinal = "<<tfinal<<endl;

	//Be precautious, restart clock for a new run:
	t_ = 0. ;
	tic_ = 0;
	ticw_ = 0 ;

	int tica = 0 ;
	int nprint = 5000;

	writesetup();

	//Tmp for debug:
	ofstream file("follow0.txt");
	ofstream file2("follow1.txt");

	while(t_<tfinal){
		//Update verlet list
		Int_->updateverlet(tic_);

		//Time step: integration & periodicity
		verletalgo2();

		if( tic_ % nrecord_ == 0){
			write();
			ticw_++;
		}

		if( tic_ % nana_ == 0 ) {
			ana_->analyse(tica,t_);
			tica++;
		}
		if( tic_ % nprint == 0) {
			std::streamsize ss = std::cout.precision();
			std::cout.precision(3);
			cout<<"t = "<<t_<<" - "<<t_/tfinal*100.<<"\% simulation"<<endl;
			std::cout.precision(ss);
		}
		//TMP
		if( tic_ % 5000 == 0 ){
			Int_->debug(tic_);
			cell_->debug(tic_);
		}

		t_+=dt_;
		tic_++;
		//TODO
		//Need to write a time.txt file which make the correspondance betwwen tics and time
	}
	file.close();
	file2.close();
}


//Write setup of the simulation
//Put here what matters to be remained later after forgeting things...
void Algo::writesetup() const{

	ofstream os(fsetup_.c_str());
	os <<"#Simulation set-up"<<endl;
	os <<"dt "<<dt_<<endl;
	os <<"tfinal "<<ns_ * dt_ <<endl;
	os <<"nrecord "<<nrecord_<<endl;
	os <<"nana "<<nana_<<endl;
	os <<"dtmax "<<dtmax_<<endl;
	os <<"gnmax "<<gnmax_<<endl;
	os <<"en "<<e_<<endl;
	os <<"kn "<<Int_->getkn()<<endl;
	os.close();
}


void Algo::write(){
	cout<<"Writing outputs..."<<endl;
	//spl_->writeAbsolute(ticw_);
	spl_->write(ticw_);
	Int_->writeContacts(ticw_);
	cell_->write(ticw_);
	ticw_++;
}

//Pour l'instant on l'implemente de maniere naive
//et non optimale (ecriture condensee a l'aide tenseurs/vecteurs)
//On verra apres comment rendre ca plus compacte
void Algo::verletalgo2(){

	double dt2_2 = 0.5 * dt_ * dt_ ;
	double dt_2 = 0.5 * dt_ ;

	// ------------- FIRST STEP VERLET ALGO STARTS HERE
	spl_->firstStepVerlet(dt_);

	//Periodicite en position des particules
	//Peut etre a bouger dans Sample plutot
	//Ca me parait plus etre un taff de sample de modifier les positions
	vector<Particle>* ps = spl_->getSample();

	//TODO: move this function to Sample
	cell_->PeriodicBoundaries2(ps);

	//Integrate cell motion
	cell_->firstStepVerlet(dt_);

	// ------------- FIRST STEP VERLET ALGO END HERE

	Int_->detectContacts();

	//Calcul des forces entre particules a la nouvelle position fin du pas de temps
	Int_->computeForces(dt_);

	//Calcul du tenseur de contraintes internes: stress_
	//Somme stress_c et stress_s
	Int_->computeInternalStress();

	//------------- SECOND STEP VERLET ALGO STARTS HERE
	//Calcul des vitesses a la fin du pas de temps:
	spl_->secondStepVerlet(dt_),

	//Apply stress_ext: si controle en force
	//stress ext est egal a celui impose
	//Sinon il est egal a -stressint et acc nulle
	cell_->computeExternalStress(Int_->stress());
	cell_->updatehdd(Int_->stress());
	cell_->updatehd(dt_);

	//Cell deformation
	cell_->CalculStrainTensor();
}


