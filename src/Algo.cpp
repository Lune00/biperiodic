#include"Algo.hpp"
#include"Cell.hpp"
#include"Sample.hpp"
#include"Interactions.hpp"
#include"Analyse.hpp"

using namespace std;

Algo::Algo(){
	dt_ = 1. ;
	ns_ = 0 ;
	nrecord_ = 0;
	nana_ = 0 ;
	nprint_ = 0 ;
	tic_= 0; 
	t_=0.; 
	ticw_ = 0; 
	fsetup_ = "simusetup.txt";
	damping_ = false;
	dampCoeff_ = 0. ;
	tf_ = 0. ;
	init_t_ = false;
}


//Maybe we can choose between ns and tfinal (to avoid large numbers)
void Algo::init(ifstream& is){
	string token;
	is >> token;
	while(is){
		if(token=="dt") is >> dt_;
		if(token=="ns") is >> ns_;
		if(token=="tf") {
			is >> tf_;
			init_t_ = true ;
			ns_ = round(tf_/dt_);
		}
		if(token=="nana") is >> nana_;
		if(token=="nrecord") is >> nrecord_;
		if(token=="nprint") is >> nprint_;
		if(token=="damping"){
			is >> dampCoeff_;
			damping_ = true ;
		}
		if(token=="}") break;
		is >> token;
	}
}

//A adapter comme on souhaite en terme de tests
bool Algo::initcheck(){

	computedtmax();
	compute_gmax();

	//Set initial values for tic (load vs build case)
	//If build automatically start at zero
	initTics();
	if( t_ > 0. ) return false;
	if( !checkSimulationParameters() ) return false;
	if (ns_ != 0 && nrecord_ != 0 ) return true;
	else return false;
}

void Algo::initTics(){
	//Calling load in config.cfg determines starting tic
	ticw_ = spl_->startingTic();
	tica_ = ticw_;
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
		cerr<<"Choose a lower dt."<<endl;
		return false;
	}
	else return true;
}

//Compute gnmax, gtmax  and restitution coeffcient
void Algo::compute_gmax(){

	double rho = spl_->getrho();
	double m = rho * M_PI * spl_->getrmin() * spl_->getrmin(); 
	double kn = Int_->getkn();
	double kt = Int_->getkt();

	//TMP
	gnmax_ = sqrt( 2. * kn * m);
	gtmax_ = sqrt( 2. * kt * m);

	//Should be sqrt(2knm)?

	if(Int_->setgnmax()){
		double gn = gnmax_ * 0.95 ;
		Int_->setgn(gn);
	}

	if(Int_->setgtmax()){
	  double gt = gtmax_ * 0.5 ;
	  Int_->setgt(gt);
	}

	//From Radjai Book (not sure)
	//But that's true we need a limit, maybe just a little above that
	//Compute resitution coefficient:
	double ds = Int_->getgn()/gnmax_;
	e_ = exp(- M_PI * ds/(2.*sqrt(1.-ds*ds)) ) ;
}

bool Algo::checkNormalViscosity()const{

	cout<<"max normal viscosity gnmax = "<<gnmax_<<endl;
	cout<<"restitution coefficient e = "<<e_<<endl;
	if(Int_->getgn()>gnmax_){
		cerr<<"Keeping all other parameters constants, choose a lower normal viscosity."<<endl;
		return false;
	}
	else return true;
}


//Check, according DEM parameters
bool Algo::checkSimulationParameters(){

	bool is_dt_ok = checktimestep();
	bool is_gn_ok = checkNormalViscosity();
	return (is_dt_ok && is_gn_ok);
}


void Algo::run(){

	double tfinal = (ns_) * dt_ ;

	cout<<"Simulation:"<<endl;
	cout<<"dt = "<<dt_<<endl;
	cout<<"ns = "<<ns_<<endl;
	cout<<"tfinal = "<<tfinal<<endl;

	//Be precautious, restart clock for a new run:
	t_ = 0. ;
	tic_ = 0;

	if(nprint_ == 0 ) nprint_ = 500 ;

	//ticw_ (writing tick cell/sample/network) and tica_ (analysis) should start at the same initial tic, set in initTics()

	writesetup();

	//return ;

	ofstream timefile("time.txt");

	//while(t_ <= tfinal){
	while(tic_ <= ns_){

		//Update verlet list
		Int_->updateverlet(tic_);

		//Time step: integration & periodicity
		verletalgo2();

		if(damping_) damping(dampCoeff_);

		if( tic_ % nrecord_ == 0) {
			write();
			timefile<<tic_<<" "<<t_<<endl;
		}

		if( tic_ % nana_ == 0 ) {
			ana_->analyse(tica_,t_,false);
			tica_++;
		}
		if( tic_ % nprint_ == 0) {
			std::streamsize ss = std::cout.precision();
			std::cout.precision(4);
			cout<<"step : "<<tic_<<" t = "<<t_<<" - "<<t_/tfinal*100.<<"\% simulation"<<endl;
			Int_->print();
			std::cout.precision(ss);
		}

		//TMP for debug
		if( tic_ % 8000 == 0 ){
			//Int_->debug(tic_);
			cell_->debug(tic_);
			//spl_->debug(tic_);
		}

		t_+=dt_;
		tic_++;
	}

	timefile.close();
}


//Write setup of the simulation
//Put here what matters to be remained later after forgeting things...
//Used for post-processing
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
	os <<"#Parameters:"<<endl;
	os <<"density "<<spl_->getrho()<<endl;
	os <<"kn "<<Int_->getkn()<<endl;
	os <<"kt "<<Int_->getkt()<<endl;
	os <<"gn "<<Int_->getgn()<<endl;
	os <<"gt "<<Int_->getgt()<<endl;
	os <<"mu "<<Int_->getmu()<<endl;
	os <<"mh "<<cell_->getMasse()<<endl;
	os <<"h0 ";
	cell_->geth().write(os);
	os <<endl;
	os.close();
}


void Algo::write(){
	spl_->write(ticw_);
	Int_->writeContacts(ticw_);
	cell_->write(ticw_);
	ticw_++;
}

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

	// ------------- FIRST STEP VERLET ENDS HERE
	Int_->detectContacts();

	//Calcul des forces entre particules a la nouvelle position fin du pas de temps
	Int_->computeForces(dt_);

	//Calcul du tenseur de contraintes internes:
	//Has been moved to computeForces(dt)
	//Int_->computeInternalStress();

	//------------- SECOND STEP VERLET STARTS HERE

	//Calcul des vitesses a la fin du pas de temps:
	spl_->secondStepVerlet(dt_);

	cell_->computeExternalStress(Int_->stress());
	cell_->updatehdd(Int_->stress());
	cell_->updatehd(dt_);

	//Cell deformation
	cell_->CalculStrainTensor();
}

void Algo::damping(const double e) {
  //Damp particle acceleration
  spl_->damp(e);
  //Damp cell acceleration
  cell_->damp(e);
}

