#include"Cell.hpp"

using namespace std;


Cell::Cell(){
	L_ = 1. ;
	xc_ = 0.5 ;
	yc_ = 0.5 ;
	h_.set(L_,0.,0.,L_);
	h0_ = h_;
	//Stress_ext et hd egal a 0 par construction par defaut
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
		//Automatise ensuite avec prop du sample
		if(token=="L") is >> L_;
		if(token=="m") is >> mh_;

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

	if(ixx && ixy && iyx && iyy){
		cout<<"Initialisation coordonnees generalisees...ok"<<endl;
	}

	//Géometrie initiale
	h_.set(L_,0.,0.,L_);
	h0_ = h_ ;
	//Definir masse a partir de l'échantillon

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
}


//Le volume est donné par det(h) (deux vecteurs de base de la cellule)
double Cell::getVolume(){
	return h_.getDet();
}

//Construit la cellule
void Cell::write(ofstream& of,ofstream& of2,double T){

	double ux=h_.getxx();
	double uy=h_.getyx();
	double vx=h_.getxy();
	double vy=h_.getyy();

	of<<T<<" "<<xc_- L_/2.<<" "<<yc_ - L_/2<<" "<<ux<<" "<<uy<<endl;
	of<<T<<" "<<xc_- L_/2.+vx<<" "<<yc_- L_/2.+vy<<" "<<ux<<" "<<uy<<endl;
	of<<T<<" "<<xc_- L_/2.<<" "<<yc_- L_/2.<<" "<<vx<<" "<<vy<<endl;
	of<<T<<" "<<xc_- L_/2.+ux<<" "<<yc_- L_/2.+uy<<" "<<vx<<" "<<vy<<endl;

}

//On travaille sur les coordonnees reduites
//Si elles sont plus petites que 0 ou plus grandes que 1 on periodise
void Cell::PeriodicBoundaries2(std::vector<Particle>& sp){
	//Vecteur a1(h_.getxx(),h_.getyx());
	//Vecteur a2(h_.getxy(),h_.getyy());

	//double Cx = a1.getNorme()  ;
	//double Cy = a2.getNorme()  ;

	for(std::vector<Particle>::iterator it = sp.begin() ; it != sp.end(); it++){
		double dx = 0.;
		double dy = 0.;
		if(it->getR().getx() > 1.) dx = -1.;
		if(it->getR().getx() < 0.) dx =  1.;
		if(it->getR().gety() > 1.) dy = -1.;
		if(it->getR().gety() < 0.) dy =  1.;
		it->Periodize(dx,dy);
		//if( it->getR().getx() > Cx) dx = - a1.getNorme();
		//if( it->getR().getx() < -Cx) dx =  a1.getNorme();
		//if( it->getR().gety() > Cy) dy = - a2.getNorme();
		//if( it->getR().gety() < -Cy) dy =  a2.getNorme();

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
	//Ecrire tenseur deformartion
	//s_.write(of2);
}

//Ici on applique les BC definis par User: controle force/vitesse qui ensuite se repercute dans schema integration
//Suppose d'avoir stressInt au temps t
//Defaut: test a chaque fois ce qui est controle ou non, alors qu'on le sait depuis le début...
void Cell::ApplyBC(){

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
