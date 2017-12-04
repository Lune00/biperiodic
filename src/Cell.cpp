#include"Cell.hpp"

using namespace std;

Cell::Cell(Config& config){
	L_ = 1. ;
	//xc_ = 0.5 * L_ ;
	//yc_ = 0.5 * L_ ;
	initCell(config);
	xc_ = 0.5;
	yc_ = 0.5;
}


//Init h et hdot
void Cell::initCell(Config& config){

	//Doit scaler avec masse echantillon
	mh_ = 2.;

	//Imposed BC by user:
	for(int i=0;i<4;i++){
		Control_[i] = config.getBCU(i);
	}

	//Geometrie initiale:
	h_.set(L_,0.,0.,L_);
	h0_ = h_;

	//Utilisateur
	//On initialise les tenseurs Ld et StressExt avec ceux de config
	//Cinematique:
	Ld_ = config.returnLd();
	//Dynamique:
	stress_ext = config.returnStress();
	hdd_ = stress_ext * (1./mh_);
	//On impose plutot un tenseur de gradient de vitesse (en eliminant la rotation) 
	//Ld_.set(0.,0.,shearrate,0.);
	//Cinematique:
	//Transformation pour le calcul des BC:
	//Vitesse de deformation de la cellue:
	//+tension,-compression
	hd_=Ld_*h_;
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
//PeriodicBoundary Conditions
//?????
void Cell::PeriodicBoundaries(std::vector<Particle>& sp){
	for(std::vector<Particle>::iterator it = sp.begin() ; it != sp.end(); it++){
		double dx = 0.;
		double dy = 0.;
		if( it->getR().getx() > getL2() )  dx = -1.;
		if( it->getR().getx() < -getL2() ) dx = 1.;
		if( it->getR().gety() < -getL2() ) dy = 1.;
		if( it->getR().gety() > getL2() )  dy = -1.;
		it->Periodize(dx,dy);

		//Other method without IF
		//Bug quand centre en 0,0 je comprends pas pourquoi
		//it->Periodize( - getL() *  floor( it->getR().getx()/ getL2()) ,0. );
		//it->Periodize(0., - getL() * floor( it->getR().gety()/ getL2()) );
	}
}

void Cell::PeriodicBoundaries2(std::vector<Particle>& sp){
	Vecteur a1(h_.getxx(),h_.getyx());
	Vecteur a2(h_.getxy(),h_.getyy());

	double Cx = a1.getNorme() * 0.5 ;
	double Cy = a2.getNorme() * 0.5 ;

	for(std::vector<Particle>::iterator it = sp.begin() ; it != sp.end(); it++){
		double dx = 0.;
		double dy = 0.;
		if( it->getR().getx() > Cx) dx = - a1.getNorme();
		if( it->getR().getx() < -Cx) dx =  a1.getNorme();
		it->Periodize(dx,dy);

	}
}
//Normalement si on change Ld ca ne changera pas le Ld initial
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
