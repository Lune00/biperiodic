#include"Cell.hpp"
#include"Sample.hpp"

using namespace std;


Cell::Cell(){
	//Stress_ext et hd egal a 0 par construction par defaut
	Lx_ = 1. ;
	Ly_ = 1. ;
	xc_ = 0.5 ;
	yc_ = 0.5 ;
	h_.set(Lx_,0.,0.,Ly_);
	h0_ = h_;
	initCG_ = false ;
	initGeometry_ = false;
	initMass_ = false;
	mh_ = 1.;
	//Multiplier of sample mass
	mh_factor_ = 1. ;
	L_auto_ = false;
	imposeForce_ = false;
	stopForce_ = false;
	tstop_=0.;
	exy_stop_ = 0. ;
	stop_Shear_ = false;
	reverse_Shear_ = false;
	amp_shear_ = 0. ;

	folder_ = string();
	fcell_ = "cell.txt";

	//DEBUG
	ofstream debug("h.txt");
	debug.close();
	debug.open("hd.txt");
	debug.close();
	debug.open("hdd.txt");
	debug.close();
}

//Initialisation from configuration file
//Careful while loading if changing Boundary Conditions
void Cell::init(ifstream& is){

	//Check if every direction is initialised

	//The mode (v vitesse f stress) is stored in Control_ array
	//The value of the mode is stored in Control_values_Init array

	bool ixx = false;
	bool ixy = false;
	bool iyx = false;
	bool iyy = false;

	string token;
	is >> token;

	while(is){
		//Build manually cell geometry (assumed to be rectangular)
		if(token=="Lx") is >> Lx_;
		if(token=="Ly") is >> Ly_;
		//Auto: defined from bounding box around particles
		if(token=="L_auto") L_auto_ = true ;
		if(token=="mh_factor") is >> mh_factor_;

		//Boundary conditions:
		//direction v(speed)/f(stress) value
		//By default:
		if(token=="xx"){
			is >> Control_[0];
			is >> loadXX_; 
			ixx = true ;
		}
		if(token=="xy") {
			is >> Control_[1];
			is >> loadXY_; 
			ixy = true ;
		}
		if(token=="yx") {
			is >> Control_[2];
			is >> loadYX_; 
			iyx = true;
		}
		if(token=="yy") {
			is >> Control_[3];
			is >> loadYY_; 
			iyy = true;
		}

		if(token=="imposeForce"){
			is >> amplitude_;
			is >> mode_;
			imposeForce_ = true ;
		}

		//Driving variation
		if(token=="stopForce"){
			stopForce_ = true ;
			is >> tstop_;
		}
		if(token=="stopShear"){
			stop_Shear_ = true ;
			is >> exy_stop_;
		}
		if(token=="reverseShear"){
			reverse_Shear_ = true ;
			is >> amp_shear_;
		}

		//user : isocompression P
		if(token=="isocompression"){
			Control_[0] = 'f';
			Control_[3] = 'f';
			Control_[1] = 'v';
			Control_[2] = 'v';
			is >> loadXX_;
			loadYY_ = loadXX_;
			loadXY_ = 0. ;
			loadYX_ = 0. ;
		}
		//user: shear P dotgamma
		//user: forceDriven F mode P

		if(token=="}") break;

		is >> token ;
	}

	//Valid initialisation
	if(ixx && ixy && iyx && iyy){
		initCG_ = true ;
	}

	if(!L_auto_){
	h_.set(Lx_,0.,0.,Ly_);
	h0_ = h_ ;
	xc_ = 0.5 * Lx_;
	yc_ = 0.5 * Ly_;
	initGeometry_ = true;
	}

}


void Cell::talkinit(Sample& spl){

	//Load : if true continue a previous simulation from file
	//if fale start from a new sample (freshly created)

	//Mmax
	double rmax = spl.getrmax();
	double Mmax = M_PI * rmax * rmax * spl.getrho(); 

	if(spl.loaded()){
		//Ask which one.
		unsigned int filetoload = spl.filetoload();
		//Test if exists inside load
		//Ici on charge h, hd, hdd
		load(filetoload);
		//New reference geometry for strain computation
		h0_ = h_ ;
		//Assign mass
		mh_ = mh_factor_ * Mmax;

		initGeometry_ = true;
		initMass_ = true;
	}
	else{
		//If not, L_auto or not?
		//Initial Cell Geometry is rectangular
		//By def hd and stress_ext are set by user in config file
		if(L_auto_){
			Lx_ = spl.getxmax() - spl.getxmin();
			Ly_ = spl.getymax() - spl.getymin();
			//TMP
			h_.set(Lx_,0.0,0.,Ly_);
			//h_.set(Lx_,0.,0.,Ly_);
			xc_ = 0.5 * Lx_;
			yc_ = 0.5 * Ly_;
			h0_ = h_ ;
			initGeometry_ = true;
		}

		//Mass: sample mass for inertia
		mh_ = mh_factor_ * Mmax;
		initMass_ = true;
	}

	//Les nouvelles CLS seront appliquees au cours du premier
	//pas de temps, lors de l'appel de firstStepVerlet

	if(Control_[0] == 'v'){
		hd_.setxx(loadXX_);
	}
	else{
		stress_ext.setxx(loadXX_);
	}
	if(Control_[1] == 'v'){
		hd_.setxy(loadXY_);
	}
	else{
		stress_ext.setxy(loadXY_);
	}
	if(Control_[2] == 'v'){
		hd_.setyx( loadYX_);
	}
	else{
		stress_ext.setyx(loadYX_);
	}
	if(Control_[3] == 'v'){
		hd_.setyy(loadYY_);
	}
	else{
		stress_ext.setyy(loadYY_);
	}
}

bool Cell::initcheck(){
	return (initCG_ && initGeometry_ && initMass_);
}

//Volume of the cell is given by its determinant 
double Cell::getVolume() const{
	return h_.getDet();
}

void Cell::write(const int k) const{
	string filename = formatfile( folder_, fcell_, k );
	ofstream file(filename.c_str());
	file.precision(12);
	h_.write(file);
	file<<" ";
	hd_.write(file);
	file<<" ";
	hdd_.write(file);
	file<<endl;
	file.close();
}

void Cell::load(const int k) {
	string filename = formatfile( folder_, fcell_, k );
	ifstream is(filename.c_str());
	if(!is){
		cerr<<"Cell::load "<<filename<<" fail."<<endl;
	}
	else{
		h_.load(is);
		hd_.load(is);
		hdd_.load(is);
	}
	is.close();
}

void Cell::readh0(ifstream& is){
	h0_.load(is);
	return;
}

//Construit la cellule
//No need for now...
void Cell::writeGeometry(const int k) const{
	//string filename = formatfile( folder_, fcell_, k );
	//ofstream file(filename.c_str());
	//double ux=h_.getxx();
	//double uy=h_.getyx();
	//double vx=h_.getxy();
	//double vy=h_.getyy();
	//file<<k<<" "<<xc_- Lx_/2.<<" "<<yc_ - Ly_/2<<" "<<ux<<" "<<uy<<endl;
	//file<<k<<" "<<xc_- Lx_/2.+vx<<" "<<yc_- Ly_/2.+vy<<" "<<ux<<" "<<uy<<endl;
	//file<<k<<" "<<xc_- Lx_/2.<<" "<<yc_- Ly_/2.<<" "<<vx<<" "<<vy<<endl;
	//file<<k<<" "<<xc_- Lx_/2.+ux<<" "<<yc_- Ly_/2.+uy<<" "<<vx<<" "<<vy<<endl;

	//	file.close();
}

//Particle periodicity in position
//On travaille sur les coordonnees reduites
//Si elles sont plus petites que 0 ou plus grandes que 1 on periodise
//A bouger dans Sample peut etre
void Cell::PeriodicBoundaries2(std::vector<Particle>* sp){

	for(std::vector<Particle>::iterator it = sp->begin() ; it != sp->end(); it++){

		while(it->getx() > 1.){
			it->addrx(-1.);
		}
		while(it->gety() > 1.){
			it->addry(-1.);
		}
		while(it->gety() < 0.){
			it->addry(1.);
		}
		while(it->getx() < 0.){
			it->addrx(1.);
		}
	}
}

//Strain tensor of the cell
//Enginering strain
void Cell::CalculStrainTensor(){

	double Lx0= h0_.getxx();
	double Ly0= h0_.getyy();
	double Lxt= h_.getxx();
	double Lyt= h_.getyy();

	double sxx = (Lxt-Lx0)/Lx0;
	double syy = (Lyt-Ly0)/Ly0;
	double sxy = (h_.getxy() - h0_.getxy())/Ly0;
	double syx = (h_.getyx() - h0_.getyx())/Lx0; 

	s_.set(sxx,sxy,syx,syy);
	s_.eigenVectors();
}



double Cell::getxc() const{
	return ((h_.getxx() + h_.getxy()) * 0.5);
}

double Cell::getyc() const{
	return ((h_.getyx() + h_.getyy()) * 0.5);
}

//Original cell w/h ratio for printing sample
double Cell::get_width() const{
	return h0_.getxx();
}

double Cell::get_height() const{
	return h0_.getyy();
}

//Integrate h and hd at the begining of the step
void Cell::firstStepVerlet(const double dt) {

	const double dt_2 = dt * 0.5 ;
	const double dt2_2 = dt * dt * 0.5 ;

	if(getControlxx() == 'f' ) {
		h_.setxx( h_.getxx() + hd_.getxx() * dt + hdd_.getxx() * dt2_2);
		hd_.setxx( hd_.getxx() + hdd_.getxx() * dt_2);
	}
	else{
		h_.setxx( h_.getxx() + hd_.getxx() * dt) ; 
		hd_.setxx(loadXX_);
		hdd_.setxx(0.);
	}

	if(getControlxy() == 'f' ) {
		h_.setxy( h_.getxy() + hd_.getxy() * dt  + hdd_.getxy() * dt2_2 ); 
		hd_.setxy( hd_.getxy() + hdd_.getxy() * dt_2 ); 
	}
	else{
		h_.setxy(h_.getxy() + hd_.getxy() * dt ); 
		hd_.setxy(loadXY_);
		hdd_.setxy(0.);
	}

	if(getControlyx() == 'f' ) {
		h_.setyx(h_.getyx() + hd_.getyx() * dt  + hdd_.getyx() * dt2_2 );
		hd_.setyx(hd_.getyx() + hdd_.getyx() * dt_2); 
	}
	else{
		h_.setyx(h_.getyx() + hd_.getyx() * dt ); 
		hd_.setyx(loadYX_);
		hdd_.setyx(0.);
	}

	if(getControlyy() == 'f' ) {
		h_.setyy(h_.getyy() + hd_.getyy() * dt  + hdd_.getyy() * dt2_2 );
		hd_.setyy(hd_.getyy() + hdd_.getyy() * dt_2 );
	}
	else{
		h_.setyy(h_.getyy() + hd_.getyy() * dt ); 
		hd_.setyy(loadYY_);
		hdd_.setyy(0.);
	}
}

void Cell::computeExternalStress(const Tensor2x2& stress_int){
	//Si controle en vitesse alors hdd 0 ds cette direction
	//Composante xx:
	if(Control_[0] == 'v' ) {
		stress_ext.setxx(-stress_int.getxx());
	}
	//Composante xy:
	if(Control_[1] == 'v' ) {
		stress_ext.setxy(-stress_int.getxy());
	}
	//Composante yx:
	if(Control_[2] == 'v' ) {
		stress_ext.setyx(-stress_int.getyx());
	}
	//Composante yy:
	if(Control_[3] == 'v' ) {
		stress_ext.setyy(-stress_int.getyy());
	}
}

//Compute "acceleration" at the end of the time step
//using stress_ext and stress_int at the end of the time step
//If controled in velocity, hdd is zero in that direction
void Cell::updatehdd(const Tensor2x2 stress_int){

	Tensor2x2 TotalStress = stress_int + stress_ext;
	Tensor2x2 hinv = h_.getInverse();
	double Vmh = getVolume()/mh_;

	hdd_ = hinv * (TotalStress) * Vmh;

	if(Control_[0] == 'v' ) {
		hdd_.setxx(0.);
	}
	//Composante xy:
	if(Control_[1] == 'v' ) {
		hdd_.setxy(0.);
	}
	//Composante yx:
	if(Control_[2] == 'v' ) {
		hdd_.setyx(0.);
	}
	//Composante yy:
	if(Control_[3] == 'v' ) {
		hdd_.setyy(0.);
	}
}

//Second verlet step in velocity
//On a l'accleration en fin de pas, on peut integrer en vitesse a la fin du pas
void Cell::updatehd(const double dt){
	const double dt_2 = dt * 0.5 ;
	hd_ += hdd_ * dt_2;
}

void Cell::debug(const int k)const{

	ofstream debug("h.txt",ios::app);
	ofstream debug1("hd.txt",ios::app);
	ofstream debug2("hdd.txt",ios::app);

	debug<<k<<" "<<h_.getxx()<<" "<<h_.getxy()<<" "<<h_.getyx()<<" "<<h_.getyy()<<endl;

	debug1<<k<<" "<<hd_.getxx()<<" "<<hd_.getxy()<<" "<<hd_.getyx()<<" "<<hd_.getyy()<<endl;

	debug2<<k<<" "<<hdd_.getxx()<<" "<<hdd_.getxy()<<" "<<hdd_.getyx()<<" "<<hdd_.getyy()<<" "<<stress_ext.getxx()<<" "<<stress_ext.getyy()<<endl;

	debug.close();
	debug1.close();
	debug2.close();
}

void Cell::damp(const double e){


	if(getControlxx() == 'f' ) {

		if(hdd_.getxx() * hd_.getxx() >= 0.0){
			hdd_.setxx(hdd_.getxx() * ( 1. - e));
		}
		else  hdd_.setxx(hdd_.getxx() * ( 1. + e));

	}

	if(getControlyy() == 'f' ) {

		if(hdd_.getyy() * hd_.getyy() >= 0.0){
			hdd_.setyy(hdd_.getyy() * ( 1. - e));
		}
		else  hdd_.setyy(hdd_.getyy() * ( 1. + e));

	}

}

void Cell::reverseShear(){
	if( s_.getxy() > amp_shear_){
		loadXY_ = -fabs(loadXY_);
	}
	if( s_.getxy() < -amp_shear_){
		loadXY_ = fabs(loadXY_);
	}
}
