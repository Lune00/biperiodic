#include<iostream>
#include<string>
#include<cmath>
#include<fstream>
#include<vector>

class Tensor2x2;
class Vecteur;
class Cell;
class Particle;

using namespace std;



class Vecteur{
	private:
		double x_;
		double y_;
	public:
		Vecteur(): x_(0.), y_(0.){};
		Vecteur(double x, double y): x_(x), y_(y){};
		~Vecteur(){};

		Vecteur productTensor(Tensor2x2&);
		void print(){cout<<x_<<" "<<y_<<endl;}
		void update(double x,double y){ x_=x; y_=y;}
		void add(double dx, double dy) { x_ += dx ; y_ += dy;}

		//Accessors:
		double getNorme() const {return sqrt( x_*x_ + y_ * y_);}
		double getx() const {return x_;}
		double gety() const {return y_;}
		//Mutators:
		void setx(double x) {x_ = x;}
		void sety(double y) {y_ = y;}
};

class Tensor2x2{
	private:
		double xx_;
		double xy_;
		double yx_;
		double yy_;
		//Valeurs propres:
		double l1_;
		double l2_;
		//Vecteurs propres:
		Vecteur u1_;
		Vecteur u2_;
	public:
		Tensor2x2(): xx_(0.0), xy_(0.0),yx_(0.0), yy_(0.0){};
		Tensor2x2(double xx, double xy, double yx, double yy): xx_(xx), xy_(xy), yx_(yx), yy_(yy){};
		~Tensor2x2(){};
		
		void set(double xx, double xy, double yx, double yy);
		double getDet(){ return (xx_*yy_ - xy_*yx_) ;}
		Tensor2x2 getInverse();
		void affiche(){cout<<xx_<<" "<<xy_<<endl;
			cout<<yx_<<" "<<yy_<<endl;}

		//Accesseurs:
		double getxx() const { return xx_;} 
		double getxy() const { return xy_;} 
		double getyx() const { return yx_;} 
		double getyy() const { return yy_;} 

		void eigenValues();
		void eigenVectors();
		void write(ofstream&);

		Tensor2x2 operator * (Tensor2x2 a){
			Tensor2x2 p;
			p.xx_ = xx_*a.xx_ + xy_*a.yx_;
			p.xy_ = xx_*a.xy_ + xy_*a.yy_;
			p.yx_ = yx_*a.xx_ + yy_*a.yx_;
			p.yy_ = yx_*a.xy_ + yy_*a.yy_;
			return p;
		}

		Tensor2x2 operator * (double a){
			Tensor2x2 p;
			p.xx_ = xx_*a;
			p.xy_ = xy_*a;
			p.yx_ = yx_*a;
			p.yy_ = yy_*a;
			return p;
		}

		Tensor2x2 operator + (Tensor2x2 a){
			Tensor2x2 p;
			p.xx_ = xx_+a.xx_;
			p.xy_ = xy_+a.xy_;
			p.yx_ = yx_+a.yx_;
			p.yy_ = yy_+a.yy_;
			return p;
		}
};


void Tensor2x2::eigenValues(){
	double b= -yy_ - xx_;
	double c=-xy_*yx_+xx_*yy_;
	double d=b*b-4*c;
	l1_=(-b-sqrt(d))/2.;
	l2_=(-b+sqrt(d))/2.;
}

void Tensor2x2::eigenVectors(){

	eigenValues();

	double norme;
	double epsilon=0.0001;
	double vx1_,vy1_,vx2_,vy2_;

	if (abs(l1_)<epsilon && abs(l2_)<epsilon) 
	{
		return;
	}
	else
	{
	vx1_=xy_ / ( - xx_+ l1_);
	vy1_=1.;
	norme=sqrt(vx1_*vx1_ + vy1_*vy1_);
	vx1_/=norme;
	vy1_/=norme;
	
	vx2_=xy_/( - xx_ + l2_);
	vy2_=1.;
	norme=sqrt(vx2_*vx2_ + vy2_*vy2_);
	vx2_/=norme;
	vy2_/=norme;

	cout<<vx1_<<" "<<vy1_<<endl;
	u1_.setx(vx1_);
	u1_.sety(vy1_);
	u2_.setx(vx2_);
	u2_.sety(vy2_);
	}
}

void Tensor2x2::write(ofstream& of){
	of<<"0. 0. "<<u1_.getx()<<" "<<u1_.gety()<<" "<<u2_.getx()<<" "<<u2_.gety()<<endl;
}

void Tensor2x2::set(double xx, double xy, double yx, double yy){
	xx_ = xx;
	xy_ = xy;
	yx_ = yx;
	yy_ = yy;
}

Tensor2x2 Tensor2x2::getInverse(){
	double epsilon=0.00001;
	//Matrice nulle par defaut
	Tensor2x2 Inverse;
	double det = getDet();
	if(abs(det)<epsilon){
		cerr<<"@La matrice n'a pas d'inverse!"<<endl;
		return Inverse;
	}
	else{
		Inverse.set(yy_/det,-yx_/det,-xy_/det,xx_/det);
		return Inverse;
	}
}



Vecteur Vecteur::productTensor(Tensor2x2& T){
	Vecteur r;
	r.x_ = T.getxx() * this->x_ + T.getxy() * this->y_;
	r.y_ = T.getyx() * this->x_ + T.getyy() * this->y_;
	return r;
}

class Particle{
	private:
		//Position absolue
		Vecteur r_;
		//Reduced coordinate
		Vecteur s_;
		//Vitesse:
		Vecteur v_;
		//Acceleration:
		Vecteur a_;
	public:
		Particle(): r_(), v_() {}; 
		Particle(Vecteur r, double L, Vecteur v) {r_ = r; v_ = v; calculReduced(L);}
		~Particle(){};
		void updateSimple(double dt);
		Vecteur& getR() { return r_;}
		void write(ofstream&);
		void calculReduced(double L);
		void affiche();
};

//sx et sy doivent etre compris entre 0 et 1 (initialiement dans la cellule)
void Particle::calculReduced(double L){
	s_.update(r_.getx() / L,r_.gety() / L);
}

void Particle::write(ofstream& of){
	of<< r_.getx()<<" "<<r_.gety()<<endl;
}

void Particle::affiche(){
	cout<<"Coordonnees absolues: "<<r_.getx()<<" "<<r_.gety()<<endl;
	cout<<"Coordonnees reduites: "<<s_.getx()<<" "<<s_.gety()<<endl;
}	


//Euler methode: update acceleration, vitesse and position
//Kinematics: update juste position a partir de la vitesse
void Particle::updateSimple(double dt){
	double x = r_.getx() + v_.getx() * dt;
	double y = r_.gety() + v_.gety() * dt;
	r_.update(x,y);
}

//Permet de definir l'essai mecanique:Conditions aux limites (force/vitesse), mais aussi les parametres cellule etc...
class Config{
	private:
		//Tableau de 4 parametres de controle (chaque direction) "v" -> vitesse "f"-> force
		 char BCU[4];
		 //Controle en vitesse:
		 Tensor2x2 LdUser_;
		 //Controle en force:
		 Tensor2x2 StressUser_;
	public:
		 Config();
		 ~Config(){};
		 Tensor2x2 returnLd() const { return LdUser_;}
		 Tensor2x2 returnStress() const { return StressUser_;}
};
//Essai mecanique: 2 fichiers input parametres pour Config, echantillon initial (xmin,xmax...)
//Prendra input eventuellement dans un fichier
Config::Config(){
	//Pour le moment tout en vitesse (cinematique)
	for(int i=0;i<4;i++){
		BCU[i] = 'v';
	}
	LdUser_.set(0.,0.,1.,0.);
	StressUser_.set(0.,0.,0.,0.);
}

//Cellule periodique
class Cell{

	private:
		//Dimension:
		double L_;
		//Coordonnees centre:
		double xc_;
		double yc_;

		//Metrics: collective degrees of freedom
		Tensor2x2 h_;
		Tensor2x2 hd_;

		//Strain rate tensor (BC), on applique ca plutot que hd_ (utilisateur)
		Tensor2x2 sd_;
		//Strain tensor: cumulative (time measure of the essai)
		Tensor2x2 s_;
		//Velocity gradient tensor (a voir plus tard, mais BC plus generale possible)
		Tensor2x2 Ld_;
		//Imposed boundary conditions: garde en memoire ce qui est controle et libre
		//On impose soit une vitesse hdd, soit une contrainte stress_ext
		//Il faut que la partie antisymetrique de Ldot soit nulle, equiv a h(hdot) symetrique impose

		//Stress: int et ext
		Tensor2x2 stress_ext;
		Tensor2x2 stress_int;
		//Cell mass:(a ajuster plus tard)
		double mh_;
	public:
		Cell(Config&);
		//A l'avenir constructeur prend en argument un truc qui init/garde en memoire  les DOF controles
		Cell(double L,Config& config): L_(L),xc_(L/2.),yc_(L/2.){initAffine(config);}
		~Cell(){};

		//Init les parametres BC et de la grille par utilisateur
		void initAffine(Config&);
		void update(Tensor2x2,Tensor2x2);
		void computeStrainTensor();

		//Met a jours la periodicite des particules
		void PeriodicBoundaries(Particle&);
		double getVolume();
		void write(ofstream&,ofstream&,double);

		//Acces
		double getMasse() const { return mh_;}
		double getL() const { return L_;}
		double getL2() const { return L_/2.;}
		Tensor2x2 geth() const { return h_;}
		Tensor2x2 gethd() const { return hd_;}
		Tensor2x2 getStressInt() const { return stress_int;}
		Tensor2x2 getStressExt() const { return stress_ext;}

		//Debug & track:
		void affiche(){cout<<xc_<<" "<<yc_<<" "<<getVolume()<<endl;}

};

Cell::Cell(Config& config){
	L_ = 1. ;
	xc_ = 0.5 * L_ ;
	yc_ = 0.5 * L_ ;
	initAffine(config);
}

//Construit la cellule
void Cell::write(ofstream& of,ofstream& of2,double T){

	double ux=h_.getxx();
	double uy=h_.getxy();
	double vx=h_.getyx();
	double vy=h_.getyy();
	of<<T<<" "<<xc_<<" "<<yc_<<" "<<ux<<" "<<uy<<endl;
	of<<T<<" "<<xc_+vx<<" "<<yc_+vy<<" "<<ux<<" "<<uy<<endl;
	of<<T<<" "<<xc_<<" "<<yc_<<" "<<vx<<" "<<vy<<endl;
	of<<T<<" "<<xc_+ux<<" "<<yc_+uy<<" "<<vx<<" "<<vy<<endl;

	//of<<T<<" 0 "<<" "<<"0"<<" "<<ux<<" "<<uy<<endl;
	//of<<T<<" "<<vx<<" "<<vy<<" "<<ux<<" "<<uy<<endl;
	//of<<T<<" 0 "<<" "<<"0"<<" "<<vx<<" "<<vy<<endl;
	//of<<T<<" "<<ux<<" "<<uy<<" "<<vx<<" "<<vy<<endl;
	hd_.eigenVectors();
	hd_.write(of2);
}

//Init h et hdot
void Cell::initAffine(Config& config){

	//Geometrie initiale:
	h_.set(L_,0.,0.,L_);

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
	hd_.affiche();

	//On transforme ca en impose Ld et hd
	//hd_.set(0.,0.,1.,0.);

	//Acceleration/Stress ext
	//stress_ext.set(0.,0.,0.,0.);
}

//Le volume est donné par det(h) (deux vecteurs de base de la cellule)
double Cell::getVolume(){
	return h_.getDet();
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


void Cell::computeStrainTensor(){

}
	

//Update h et hd de la cellule
void Cell::update(Tensor2x2 h, Tensor2x2 hd){
	h_ = h;
	hd_ = hd;
}

//Le couplage passe par le tenseur de contrainte entre cellule et particules

//Update particules, calcul tenseur contrainte interne, update cellule

//Integre le mouvement des particules avec l etat cellule au temps t
void verletalgo_particles(Cell& cell,std::vector<Particle>& ps,double dt){

	//On recupere h,hd,hdd
	Tensor2x2 h = cell.geth();
	Tensor2x2 hd = cell.gethd();
	Tensor2x2 hdd = cell.gethd();
	//Calcul de la resultante des forces entre particules
	//Calcul du tenseur de contraintes interne
	//Calcul acceleration pour chaque particule (force et acc de l'espace via sigma_ext)
	//Integration par verletalgo

}

//On integre l'equation de la dynamique pour l'espace avec un algo de verlet
//Integre la dynamique de l'espace, de la cellule au temps t, a besoin tenseur contraintes
//On recupere le tenseur de contraintes internes calcule par verletalgo_particles
void verletalgo_space(Cell& cell,std::vector<Particle>& ps, double dt){

	//Les resultantddes du tenseur de contrainte ou de hd peuvent être imposees (a voir)
	Tensor2x2 h = cell.geth();
	Tensor2x2 hd = cell.gethd();

	Tensor2x2 hdd0;
	Tensor2x2 hdd1 = hdd0;

	Tensor2x2 stressint = cell.getStressInt();
	Tensor2x2 stressext = cell.getStressExt();
	double V0=cell.getVolume();
	double mh=cell.getMasse();

	//Calcul de hdd au debut du pas, hdd0
	//
	Tensor2x2 hinv = h.getInverse();
	//Eq (24)
	//Impose BC() en force: si une direction est controlee en vitesse, on impose hdd0 nulle dans cette dir
	//Et on calcule le tenseur de contrainte internes a partir de l'equation ac terme de droite nul
	//Cad que stressext=-stressint dans cette direction
	hdd0=hinv*(V0/mh)*(stressint+stressext);

	//Calcul position(h) avec hdd0 a la fin du pas de temps
	h = h + hd*dt + hdd0*0.5*dt*dt;
	//Volume fin du pas de temps
	double V1=h.getDet();
	//Calcul hdd fin du pas, hdd1
	hinv=h.getInverse();

	//Acceleration fin du pas: pareil que pour hdd0, on ompose hdd1 nul dans direction controle en vitesse
	hdd1=hinv*(V1/mh)*(stressint+stressext);

	//Calcul de hd fin du pas de temps
	//Alors la vitesse ds dir controle en vitesse est bien constante (hdd0+hdd1=0)
	hd=hd+(hdd0+hdd1)*dt*0.5;

	//Impose BC(): on maintient les BC imposees en appliquant les coordonees de L_ imposee par utilisateur
	//Du coup c'est redondant, je n'ai pas besoin de maintenir les BC en vitesse elles le seront automatiquement avec
	//le traitement en force du dessus
	//imposeBC(hd);

	//Update Cell:on fait remonter h,hd
	cell.update(h,hd);
//	cell.affiche();
}

int main(){

	//Parametres:
	double const L = 1.;
	double dt = .01 ;
	double T = 1.;

	vector<Particle> sample;
	Config config;
	Cell cell(L,config);

	//Initialise coordonnees reduites, absolues, vitesse fluctuante
	Vecteur r0(L/3,L/3);
	Vecteur v0(0.,0.);
	Particle p(r0,L,v0);
	sample.push_back(p);

	ofstream tmp("particle.txt");
	ofstream tmp2("cell.txt");
	ofstream tmp3("strainrateTensor.txt");

	double t=0.;
	int k=0;

	do{
		//Check for boundary:
		cell.PeriodicBoundaries(p);
		//Update position:
		verletalgo_particles(cell,sample,dt);
		verletalgo_space(cell,sample,dt);

		//outputs:
		if( k%1 == 0 ){
			cell.write(tmp2,tmp3,T);
			p.write(tmp);
		}

		t+=dt;
		k++;
	}while(t<T);

	tmp.close();
	tmp2.close();
	tmp3.close();

	return 0;
}

