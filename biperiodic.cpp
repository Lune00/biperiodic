#include<iostream>
#include<string>
#include<cmath>
#include<fstream>
#include<vector>

using namespace std;

class Tensor2x2{
	private:
		double xx_;
		double xy_;
		double yx_;
		double yy_;
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
		//Imposed boundary conditions: garde en memoire ce qui est controle et libre
		//On impose soit une vitesse hdd, soit une contrainte stress_ext
		//Il faut que la partie antisymetrique de Ldot soit nulle, equiv a h(hdot) symetrique impose

		//Stress: int et ext
		Tensor2x2 stress_ext;
		Tensor2x2 stress_int;
		//Cell mass:(a ajuster plus tard)
		double mh_;
	public:
		Cell();
		//A l'avenir constructeur prend en argument un truc qui init/garde en memoire  les DOF controles
		Cell(double L): L_(L),xc_(L/2.),yc_(L/2.){initAffine();}
		~Cell(){};

		void initAffine();
		void update(Tensor2x2,Tensor2x2);
		void affiche(){cout<<xc_<<" "<<yc_<<" "<<getVolume()<<endl;}
		//Met a jours la periodicite 
		void PeriodicBoundaries(Particle&);
		double getVolume();
		void write(ofstream&,double);
		//Acces
		double getMasse() const { return mh_;}
		double getL() const { return L_;}
		double getL2() const { return L_/2.;}
		Tensor2x2 geth() const { return h_;}
		Tensor2x2 gethd() const { return hd_;}
		Tensor2x2 getStressInt() const { return stress_int;}
		Tensor2x2 getStressExt() const { return stress_ext;}
};

Cell::Cell(){
	L_ = 1. ;
	xc_ = 0.5 * L_ ;
	yc_ = 0.5 * L_ ;
}

//Construit la cellule
void Cell::write(ofstream& of,double T){

	double ux=h_.getxx();
	double uy=h_.getxy();
	double vx=h_.getyx();
	double vy=h_.getyy();
	of<<T<<" "<<xc_<<" "<<yc_<<" "<<ux<<" "<<uy<<endl;
	of<<T<<" "<<xc_+vx<<" "<<yc_+vy<<" "<<ux<<" "<<uy<<endl;
	of<<T<<" "<<xc_<<" "<<yc_<<" "<<vx<<" "<<vy<<endl;
	of<<T<<" "<<xc_+ux<<" "<<yc_+uy<<" "<<vx<<" "<<vy<<endl;
}

//Init h et hdot
void Cell::initAffine(){
	mh_ = 2.;
	//Geometrie initiale:
	h_.set(L_,0.,0.,L_);
	//Vitesse de deformation de la cellue:
	hd_.set(0.,0.,1.,0.);
	//Acceleration/Stress ext
	stress_ext.set(0.,0.,0.,0.);
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
	Tensor2x2 hinv = h.getInverse();
	cout<<"V "<<V0<<endl;
	//Eq (24)
	hdd0=hinv*(V0/mh)*(stressint+stressext);
	//Calcul position(h) avec hdd0 a la fin du pas de temps
	h = h + hd*dt + hdd0*0.5*dt*dt;
	//Volume fin du pas de temps
	double V1=h.getDet();
	//Calcul hdd fin du pas, hdd1
	hinv=h.getInverse();
	hdd1=hinv*(V1/mh)*(stressint+stressext);
	//Calcul de hd fin du pas de temps
	hd=hd+(hdd0+hdd1)*dt*0.5;
	//Update Cell:on fait remonter h,hd
	cell.update(h,hd);
	cell.affiche();
}

int main(){

	//Initialisation:
	double const L = 1.;
	double dt = .05 ;
	double T = 0.;
	vector<Particle> sample;
	Cell cell(L);

	Vecteur r0(L/3,L/3);
	Vecteur v0(0.,0.);
	//Initialise coordonnees reduites, absolues, vitesse fluctuante
	Particle p(r0,L,v0);
	sample.push_back(p);

	ofstream tmp("particle.txt");
	ofstream tmp2("cell.txt");
	//Kinetics Time integration
	for(int t = 0; t < 10; t++){
		//Check for boundary:
		cell.PeriodicBoundaries(p);
		//Update position:
		verletalgo_particles(cell,sample,dt);
		verletalgo_space(cell,sample,dt);
		T+=dt;
		p.write(tmp);
		cell.write(tmp2,T);
	}
	tmp.close();
}

