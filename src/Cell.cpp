#include"Cell.hpp"

using namespace std;

Cell::Cell(Config& config){
	L_ = 1. ;
	xc_ = 0.5 * L_ ;
	yc_ = 0.5 * L_ ;
	initCell(config);
}


//Init h et hdot
void Cell::initCell(Config& config){

	//Imposed BC by user:
	for(int i=0;i<4;i++){
		Control_[i] = config.getBCU(i);
	}

	//Geometrie initiale:
	h_.set(L_,0.,0.,L_);
	h0_ = h_;

	//Utilisateur
	//On initialise les tenseurs Ld et StressExt avec ceux de config
	mh_ = 2.;
	//Cinematique:
	Ld_ = config.returnLd();
	//Dynamique:
	stress_ext = config.returnStress();

	//On impose plutot un tenseur de gradient de vitesse (en eliminant la rotation) 
	//Ld_.set(0.,0.,shearrate,0.);

	//Cinematique:

	//Transformation pour le calcul des BC:
	//Vitesse de deformation de la cellue:
	//+tension,-compression
	hd_=Ld_*h_;
	//On transforme ca en impose Ld et hd
	//hd_.set(0.,0.,1.,0.);
	hdd_ = stress_ext * (1./mh_);

	//Acceleration/Stress ext
	//stress_ext.set(0.,0.,0.,0.);
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

	of<<T<<" "<<xc_<<" "<<yc_<<" "<<ux<<" "<<uy<<endl;
	of<<T<<" "<<xc_+vx<<" "<<yc_+vy<<" "<<ux<<" "<<uy<<endl;
	of<<T<<" "<<xc_<<" "<<yc_<<" "<<vx<<" "<<vy<<endl;
	of<<T<<" "<<xc_+ux<<" "<<yc_+uy<<" "<<vx<<" "<<vy<<endl;

	//of<<T<<" 0 "<<" "<<"0"<<" "<<ux<<" "<<uy<<endl;
	//of<<T<<" "<<vx<<" "<<vy<<" "<<ux<<" "<<uy<<endl;
	//of<<T<<" 0 "<<" "<<"0"<<" "<<vx<<" "<<vy<<endl;
	//of<<T<<" "<<ux<<" "<<uy<<" "<<vx<<" "<<vy<<endl;
	s_.eigenVectors();
	s_.write(of2);
}
//PeriodicBoundary Conditions
void Cell::PeriodicBoundaries(Particle& p){

	//Periodic position
	//if( p.getR().getx() > getL2() ) p.getR().add(-getL(),0.);
	//if( p.getR().getx() < -getL2() ) p.getR().add(getL(),0.);
	//if( p.getR().gety() < -getL2() ) p.getR().add(0.,getL());
	//if( p.getR().gety() > getL2() ) p.getR().add(0.,-getL());

	//Other method without IF
	p.getR().add( - getL() *  floor( p.getR().getx()/ getL()) ,0. );
	p.getR().add(0., - getL() * floor( p.getR().gety()/ getL()) );
}

//Update h et hd de la cellule, ainsi que L et strain tensor
//On peut meme updater Ld, car par construction le control en vitesse se fait par le controle en force
//Normalement si on change Ld ca ne changera pas le Ld initial
//A partir de la maj de Ld on peut facilement calculer le tenseur de deformations, oui mais ca marche pas... ????
//UP
void Cell::update(Tensor2x2 h, Tensor2x2 hd, double dt){

	h_ = h;
	hd_ = hd;
	Tensor2x2 hinv = h_.getInverse();

	Ld_ = hd_*hinv;
	Tensor2x2 LdT = Ld_.getTranspose();
//	Ld_.affiche();
	//Calcul du tenseur de deformation de maniere generale cumule
	//s_ = s_ +  (Ld_ + LdT) * dt;

	//Provisoire:
	double Lx0= h0_.getxx();
	double Ly0= h0_.getyy();
	double Lxt= h_.getxx();
	double Lyt= h_.getyy();

	double sxx = (Lxt-Lx0)/Lx0;
	double syy = (Lyt-Ly0)/Ly0;
	double sxy = (h_.getxy() - h0_.getxy())/Ly0;
	double syx = (h_.getyx() - h0_.getyx())/Lx0; 

	s_.set(sxx,sxy,syx,syy);
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
