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
	//Multiplier of sample mass
	mh_ = 1.;
	mh_auto_ = false;
	L_auto_ = false;

	folder_ = string();
	fcell_ = "cell.txt";

	//DEBUG
	ofstream debug("stressext.txt");
	debug.close();
	ofstream debug2("hd.txt");
	debug2.close();
}
//Initialisation a partir du fichier de configuration
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
		//Manuellement on construit la geometrie de la cellule
		if(token=="Lx") is >> Lx_;
		if(token=="Ly") is >> Ly_;
		//Auto: defini a partir du sample initial
		if(token=="L_auto") {
			L_auto_ = true ;
		}
		//This option should be avoided
		if(token=="m"){
			is >> mh_;
			initMass_ = true;
		}
		if(token=="m_auto") mh_auto_ = true;

		if(token=="xx"){
			is >> Control_[0];
			is >> Control_values_Init[0]; 
			ixx = true ;
		}
		if(token=="xy") {
			is >> Control_[1];
			is >> Control_values_Init[1]; 
			ixy = true ;
		}
		if(token=="yx") {
			is >> Control_[2];
			is >> Control_values_Init[2]; 
			iyx = true;
		}
		if(token=="yy") {
			is >> Control_[3];
			is >> Control_values_Init[3]; 
			iyy = true;
		}
		if(token=="}") break;

		is >> token ;
	}

	//Valide initialisation
	if(ixx && ixy && iyx && iyy){
		initCG_ = true ;
	}

	if(!L_auto_){
	//Géometrie rectangle initiale definie manuellement
	h_.set(Lx_,0.,0.,Ly_);
	h0_ = h_ ;
	xc_ = 0.5 * Lx_;
	yc_ = 0.5 * Ly_;
	initGeometry_ = true;
	}

}


void Cell::talkinit(Sample& spl){

	//Load from file or no?
	if(spl.loaded()){
		//Ask which one.
		unsigned int filetoload = spl.filetoload();
		//Test if exists inside load
		load(filetoload);
		h0_ = h_ ;
		initGeometry_ = true;
		//Assign mass
		//the mutliplier should be defined somewhere...
		mh_ = 10. * spl.getMass();
		initMass_ = true;
		//APPLY NEW CL!!! car la on a loadé les anciennes et
	}
	else{

		//If not, L_auto or not?
		//Initial Cell Geometry is rectangular
		//By def hd and stress_ext are set by user in config file
		if(L_auto_){
			Lx_ = spl.getxmax() - spl.getxmin();
			Ly_ = spl.getymax() - spl.getymin();
			h_.set(Lx_,0.,0.,Ly_);
			xc_ = 0.5 * Lx_;
			yc_ = 0.5 * Ly_;
			h0_ = h_ ;
			initGeometry_ = true;
		}

		if(mh_auto_){
			//Mass: sample mass for inertia
			mh_ = 10. * spl.getMass();
			initMass_ = true;
		}
	}

	if(Control_[0] == 'v'){
		hd_.setxx(Control_values_Init[0]);
	}
	else{
		stress_ext.setxx(Control_values_Init[0]);
	}
	if(Control_[1] == 'v'){
		hd_.setxy(Control_values_Init[1]);
	}
	else{
		stress_ext.setxy(Control_values_Init[1]);
	}
	if(Control_[2] == 'v'){
		hd_.setyx( Control_values_Init[2]);
	}
	else{
		stress_ext.setyx(Control_values_Init[2]);
	}
	if(Control_[3] == 'v'){
		hd_.setyy(Control_values_Init[3]);
	}
	else{
		stress_ext.setyy(Control_values_Init[3]);
	}
}

//If load called by Sample, load cell geometry and dynamics
//Load h, hd and hdd
//check if problems with init, continuity with previous and new CL

bool Cell::initcheck(){
	return (initCG_ && initGeometry_ && initMass_);
}

//Le volume est donné par det(h) (deux vecteurs de base de la cellule)
double Cell::getVolume() const{
	return h_.getDet();
}

void Cell::write(const int k) const{
	string filename = formatfile( folder_, fcell_, k );
	cout<<"cell output: "<<filename<<endl;
	ofstream file(filename.c_str());
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

//On travaille sur les coordonnees reduites
//Si elles sont plus petites que 0 ou plus grandes que 1 on periodise
//A bouger dans Sample peut etre
void Cell::PeriodicBoundaries2(std::vector<Particle>* sp){

	for(std::vector<Particle>::iterator it = sp->begin() ; it != sp->end(); it++){
		double dx = 0.;
		double dy = 0.;

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

//Ici toujours un probelem pour le calcul
void Cell::CalculStrainTensor(){

	double Lx0= h0_.getxx();
	double Ly0= h0_.getyy();
	double Lxt= h_.getxx();
	double Lyt= h_.getyy();

	double sxx = (Lxt-Lx0)/Lx0;
	double syy = (Lyt-Ly0)/Ly0;
	double sxy = (h_.getxy() - h0_.getxy())/Ly0;
	double syx = (h_.getyx() - h0_.getyx())/Lx0; 

	//Si train doit etre calcule comme ca:
	//double syx = sxy ;

	s_.set(sxx,sxy,syx,syy);
	s_.eigenVectors();
	//Ecrire tenseur deformartion
	//s_.write(of2);
}

void Cell::writeStrainTensor(ofstream& os, double t){
	os<<t<<" ";
	s_.write(os);
	s_.eigenVectors();
}

//Ici on applique les BC definis par User: controle force/vitesse qui ensuite se repercute dans schema integration
//Suppose d'avoir stressInt au temps t
//Defaut: test a chaque fois ce qui est controle ou non, alors qu'on le sait depuis le début...
void Cell::computeExternalStress(const Tensor2x2& stress_int){

	//cout<<"Stress ext:"<<endl;
	//stress_ext.print();

	//cout<<"Stress int:"<<endl;
	//stress_int.print();

	//Si controle en vitesse alors hdd 0 ds cette direction
	double xx, xy, yx, yy;
	//Composante xx:

	//Voir s'il vaut mieux forcer a chaque fois
	//ou si les erreurs de troncature foutent le bordel

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

	//WIP WIP WIP
	if(getControlxx() == 'f' ) {
		h_.setxx( h_.getxx() + hd_.getxx() * dt + hdd_.getxx() * dt2_2);
		hd_.setxx( hd_.getxx() + hdd_.getxx() * dt_2);
	}
	else{
		//Vitesse imposee reste la meme (hdd, definit par Ld au debut)
		//hdd par definition nulle si control en vitesse sur hdd
		//Normalement hd_ dans cette direction ne doit jamais etre modifie si c bien fait
		h_.setxx( h_.getxx() + hd_.getxx() * dt) ; 
	}

	if(getControlxy() == 'f' ) {
		h_.setxy( h_.getxy() + hd_.getxy() * dt  + hdd_.getxy() * dt2_2 ); 
		hd_.setxy( hd_.getxy() + hdd_.getxy() * dt_2 ); 
	}
	else{
		h_.setxy(h_.getxy() + hd_.getxy() * dt ); 
	}

	if(getControlyx() == 'f' ) {
		h_.setyx(h_.getyx() + hd_.getyx() * dt  + hdd_.getyx() * dt2_2 );
		hd_.setyx(hd_.getyx() + hdd_.getyx() * dt_2); 
	}
	else{
		h_.setyx(h_.getyx() + hd_.getyx() * dt ); 
	}

	if(getControlyy() == 'f' ) {
		h_.setyy(h_.getyy() + hd_.getyy() * dt  + hdd_.getyy() * dt2_2 );
		hd_.setyy(hd_.getyy() + hdd_.getyy() * dt_2 );
	}
	else{
		h_.setyy(h_.getyy() + hd_.getyy() * dt ); 
	}
}

//Compute "acceleration" at the end of the time step
//using stress_ext and stress_int at the end of the time step
void Cell::updatehdd(const Tensor2x2 stress_int){
	Tensor2x2 TotalStress = stress_int + stress_ext;
	Tensor2x2 hinv = h_.getInverse();
	//cout<<" "<<endl;
	//cout<<"stress ext:"<<endl;
	//stress_ext.print();
	//cout<<"internal stress:"<<endl;
	//stress_int.print();
	//cout<<"total stress:"<<endl;
	//TotalStress.print();
	//cout<<" "<<endl;
	double V = getVolume();
	hdd_ = hinv * (V/mh_) * (TotalStress);
}

//Second verlet step in velocity
//On a l'accleration en fin de pas, on peut integrer l'espace en vitesse a la fin du pas
void Cell::updatehd(const double dt){
	double dt_2 = dt * 0.5 ;
	hd_ = hd_ + hdd_ * dt_2;
}

void Cell::debug(const int k)const{
	ofstream debug("stressext.txt",ios::app);
	ofstream debug2("hd.txt",ios::app);
	debug<<stress_ext.getxx()<<" "<<stress_ext.getyy()<<endl;
	debug2<<h_.getxx()<<" "<<hd_.getxx()<<" "<<h_.getyy()<<" "<<hd_.getyy()<<endl;
	debug.close();
	debug2.close();
}

