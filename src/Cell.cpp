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
	mh_auto_ = false;
	L_auto_ = false;

	folder_ = string();
	fcell_ = "cell.txt";
}
//Initialisation a partir du fichier de configuration
void Cell::init(ifstream& is){

	//Check si tout est initialise
	bool ixx,ixy,iyx,iyy;
	ixx=ixy=iyx=iyy=false;
	double xx,xy,yx,yy;
	xx=xy=yx=yy=0.;

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
		if(token=="m"){
			is >> mh_;
			initMass_ = true;
		}
		if(token=="m_auto") mh_auto_ = true;

		if(token=="xx"){
			is >> Control_[0];
			is >>xx; 
			ixx = true ;
		}
		if(token=="xy") {
			is >> Control_[1];
			is >>xy; 
			ixy = true ;
		}
		if(token=="yx") {
			is >> Control_[2];
			is >>yx; 
			iyx = true;
		}
		if(token=="yy") {
			is >> Control_[3];
			is >>yy; 
			iyy = true;
		}
		if(token=="}") break;

		is >> token ;
	}



	//Definir tenseur contraintes ext / hd (qu'on peut trans en L)
	if(Control_[0] == 'v'){
		hd_.setxx(xx);
	}
	else{
		stress_ext.setxx(xx);
	}
	if(Control_[1] == 'v'){
		hd_.setxy(xy);
	}
	else{
		stress_ext.setxy(xy);
	}
	if(Control_[2] == 'v'){
		hd_.setyx(yx);
	}
	else{
		stress_ext.setyx(yx);
	}
	if(Control_[3] == 'v'){
		hd_.setyy(yy);
	}
	else{
		stress_ext.setyy(yy);
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

void Cell::initFromSample(Sample& spl){
	//Initial Cell Geometry
	if(L_auto_){
	Lx_ = spl.getxmax() - spl.getxmin();
	Ly_ = spl.getymax() - spl.getymin();
	h_.set(Lx_,0.,0.,Ly_);
	h0_ = h_ ;
	xc_ = 0.5 * Lx_;
	yc_ = 0.5 * Ly_;
	initGeometry_ = true;
	}

	if(mh_auto_){
	//Mass: double sample mass for inertia
	mh_ = 2. * spl.getMass();
	initMass_ = true;
	}
}


bool Cell::initcheck(){
	return (initCG_ && initGeometry_ && initMass_);
}

//Renvoie vrai si on a besoin du sample pour initialiser cellule

bool Cell::needSample(){
	if( L_auto_ || mh_auto_) return true;
	else return false;
}

//Le volume est donné par det(h) (deux vecteurs de base de la cellule)
double Cell::getVolume() const{
	return h_.getDet();
}


//Construit la cellule
void Cell::write(ofstream& of,double t) const{

	double ux=h_.getxx();
	double uy=h_.getyx();
	double vx=h_.getxy();
	double vy=h_.getyy();
	of<<t<<" "<<xc_- Lx_/2.<<" "<<yc_ - Ly_/2<<" "<<ux<<" "<<uy<<endl;
	of<<t<<" "<<xc_- Lx_/2.+vx<<" "<<yc_- Ly_/2.+vy<<" "<<ux<<" "<<uy<<endl;
	of<<t<<" "<<xc_- Lx_/2.<<" "<<yc_- Ly_/2.<<" "<<vx<<" "<<vy<<endl;
	of<<t<<" "<<xc_- Lx_/2.+ux<<" "<<yc_- Ly_/2.+uy<<" "<<vx<<" "<<vy<<endl;

}

//On travaille sur les coordonnees reduites
//Si elles sont plus petites que 0 ou plus grandes que 1 on periodise
//A bouger dans Sample peut etre
void Cell::PeriodicBoundaries2(std::vector<Particle>* sp){
	//Vecteur a1(h_.getxx(),h_.getyx());
	//Vecteur a2(h_.getxy(),h_.getyy());

	//double Cx = a1.getNorme()  ;
	//double Cy = a2.getNorme()  ;
	double Lx2 = getLx() * 0.5 ;
	double Ly2 = getLy() * 0.5 ;

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
		//if(it->getR().getx() < 0.) dx +=  1.;
		//if(it->getR().gety() > 1.) dy += -1.;
		//if(it->getR().gety() < 0.) dy +=  1.;
		//it->Periodize(dx,dy);

	}
}
//A partir de la maj de Ld on peut facilement calculer le tenseur de deformations, oui mais ca marche pas... ????
void Cell::update(Tensor2x2 h, Tensor2x2 hd){

	h_ = h;
	hd_ = hd;
	//	Tensor2x2 hinv = h_.getInverse();
	//
	//	Ld_ = hd_*hinv;
	//	Tensor2x2 LdT = Ld_.getTranspose();
	//	Ld_.affiche();
	//Calcul du tenseur de deformation de maniere generale cumule
	//s_ = s_ +  (Ld_ + LdT) * dt;

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
void Cell::ApplyBC(const Tensor2x2& stress_int){

	//TEMPORAIRE:
	//	stress_int.set(0.,0.,0.,4.);

	double xx, xy, yx, yy;
	//Composante xx:
	if(Control_[0] == 'v' ) {
		xx = -stress_int.getxx();
	}
	//Sinon elle reste egale a elle meme (imposee depuis le debut)
	else{
		xx = stress_ext.getxx();
	}
	//Composante xy:
	if(Control_[1] == 'v' ) {
		xy = -stress_int.getxy();
	}
	//Sinon elle reste egale a elle meme (imposee depuis le debut)
	else{
		xy = stress_ext.getxy();
	}
	//Composante yx:
	if(Control_[2] == 'v' ) {
		yx = -stress_int.getyx();
	}
	//Sinon elle reste egale a elle meme (imposee depuis le debut)
	else{
		yx = stress_ext.getyx();
	}
	//Composante yy:
	if(Control_[3] == 'v' ) {
		yy = -stress_int.getyy();
	}
	//Sinon elle reste egale a elle meme (imposee depuis le debut)
	else{
		yy = stress_ext.getyy();
	}
	//Mise a jour tenseur contraintes:
	stress_ext.set(xx,xy,yx,yy);
}


double Cell::getxc() const{
	return ((h_.getxx() + h_.getxy()) * 0.5);
}

double Cell::getyc() const{
	return ((h_.getyx() + h_.getyy()) * 0.5);
}
