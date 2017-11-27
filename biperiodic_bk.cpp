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
		void set(double x,double y){ x_=x; y_=y;}
		void add(double dx, double dy) { x_ += dx ; y_ += dy;}

		//Accessors:
		double getNorme() const {return sqrt( x_*x_ + y_ * y_);}
		double getx() const {return x_;}
		double gety() const {return y_;}
		//Mutators:
		void setx(double x) {x_ = x;}
		void sety(double y) {y_ = y;}

		//Surcharge operator:
		//Scalar product:
		double operator * (Vecteur a){
			return (x_*a.x_ + y_*a.y_);
		}

		Vecteur operator + (Vecteur a){
			Vecteur p;
			p.x_ = x_ + a.x_;
			p.y_ = y_ + a.y_;
			return p;
		}

		Vecteur operator * (double a){
			Vecteur p;
			p.x_ = a*x_ ;
			p.y_ = a*y_ ;
			return p;
		}
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
		Tensor2x2 getTranspose();
		void affiche(){cout<<endl;cout<<xx_<<" "<<xy_<<endl;
			cout<<yx_<<" "<<yy_<<endl;cout<<endl;}

		//Accesseurs:
		double getxx() const { return xx_;} 
		double getxy() const { return xy_;} 
		double getyx() const { return yx_;} 
		double getyy() const { return yy_;} 

		//Renvoie la composante i du tenseur:
		//0->xx 1->xy 2->yx 3->yy

		void eigenValues();
		void eigenVectors();
		void write(ofstream&);
		void writeEigenVectors(ofstream&);

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

		Vecteur operator * (Vecteur a){
			Vecteur p;
			double x = xx_*a.getx()+ xy_*a.gety();
			double y = yx_*a.getx() + yy_*a.gety();
			p.set(x,y);
			return p;
		}
};


void Tensor2x2::eigenValues(){
	double b= -yy_ - xx_;
	double c=-xy_*yx_+xx_*yy_;
	double d=b*b-4*c;
	l1_=(-b-sqrt(d))/2.;
	l2_=(-b+sqrt(d))/2.;
	double tmp=max(l1_,l2_);
	l2_=min(l1_,l2_);
	l1_=tmp;
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

	u1_.setx(vx1_);
	u1_.sety(vy1_);
	u2_.setx(vx2_);
	u2_.sety(vy2_);
	}
}

void Tensor2x2::writeEigenVectors(ofstream& of){
	of<<"0. 0. "<<u1_.getx()<<" "<<u1_.gety()<<" "<<u2_.getx()<<" "<<u2_.gety()<<endl;
}

void Tensor2x2::write(ofstream& of){
	eigenValues();
	double spherique=l1_+l2_;
	double deviatorique=l1_-l2_;
	of<<xx_<<" "<<xy_<<" "<<yx_<<" "<<yy_<<" "<<spherique<<" "<<deviatorique<<" "<<l1_<<" "<<l2_<<xx_+yy_<<endl;
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

Tensor2x2 Tensor2x2::getTranspose(){
	Tensor2x2 Transpose(xx_,yx_,xy_,yy_);
	return Transpose;
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
		//Acceleration debut pas de temps:
		Vecteur a0_;
		//Forces:
		Vecteur f_;
		//Torque:
		double t_;
		//Vrot:
		double vrot_;
		//Rot: angular position
		double rot_;
		//mass:
		double m_;
		//Rayon:
		double R_;
		
	public:
		Particle(): r_(), v_(), s_(), a0_(), f_(), t_(0.), vrot_(0.), rot_(0.), m_(1.), R_(1.) {}; 
		Particle(Vecteur r, double L, Vecteur v) {r_ = r; v_ = v; init(L);}
		~Particle(){};

		void updateSimple(double dt);
		void write(ofstream&);
		void init(double L);
		void affiche();

		double getMasse() const {return m_;}
		double getRadius() const { return R_;}
		double getTorque() const {return t_;}

		Vecteur getR() const { return r_;}
		Vecteur getV() const { return v_;}
		Vecteur getForce() const { return f_;}
		Vecteur geta0() const {return a0_;}
		void setV(Vecteur v) { v_=v;}
		void setR(Vecteur r, Vecteur a0) { r_ = r; a0_=a0;}
		
};

//sx et sy doivent etre compris entre 0 et 1 (initialiement dans la cellule)
void Particle::init(double L){
	//Les vecteurs appellent constructeur, tout egal a 0
	m_=1.;
	R_=0.1;
	vrot_=0.;
	rot_=0.;
	t_=0.;
	s_.set(r_.getx() / L,r_.gety() / L);
}

void Particle::write(ofstream& of){
	of<< r_.getx()<<" "<<r_.gety()<<" "<<v_.getx()<<" "<<v_.gety()<<endl;
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
	r_.set(x,y);
}

//Permet de definir l'essai mecanique:Conditions aux limites (force/vitesse), mais aussi les parametres cellule etc...
class Config{
	private:
		//Tableau de 4 parametres de controle (chaque direction) "v" -> vitesse "f"-> force
		//0 xx 1 xy 2 yx 3 yy
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

		 char getBCxx() const { return BCU[0];}
		 char getBCxy() const { return BCU[1];}
		 char getBCyx() const { return BCU[2];}
		 char getBCyy() const { return BCU[3];}
		 char getBCU(int i) const { return BCU[i];}
};
//Essai mecanique: 2 fichiers input parametres pour Config, echantillon initial (xmin,xmax...)
//Prendra input eventuellement dans un fichier
Config::Config(){
	//Pour le moment tout en vitesse (cinematique)
	for(int i=0;i<4;i++){
		BCU[i] = 'v';
	}
	BCU[3]='f';
	//Par souci de convention, les DL non controles sont init a 0
	//Par definition, le taux de cisaillement pur est egale a la moitie de LdUser.yx(ou xy)
	LdUser_.set(0.,0.,0.,0.);
	StressUser_.set(0.,0.,0.,4.);
}

//Cellule periodique
class Cell{

	private:
		//Dimension:
		double L_;
		//Coordonnees centre:
		double xc_;
		double yc_;
		//Cell mass:(a ajuster plus tard)
		double mh_;

		//Imposed boundary conditions: garde en memoire ce qui est controle et libre
		//On impose soit une vitesse hdd, soit une contrainte stress_ext
		//Il faut que la partie antisymetrique de Ldot soit nulle, equiv a h(hdot) symetrique impose
		char Control_[4];

		//Metrics: collective degrees of freedom
		Tensor2x2 h_;
		//Garde en memoire la forme originelle cellule pour calcul d'engeniring strain (L(t)-L0/LO))
		Tensor2x2 h0_;
		Tensor2x2 hd_;
		//Strain tensor: cumulative (time measure of the essai)
		Tensor2x2 s_;
		//Velocity gradient tensor (a voir plus tard, mais BC plus generale possible)
		Tensor2x2 Ld_;
		//Stress: int et ext
		Tensor2x2 stress_ext;
		//Maj par particules
		Tensor2x2 stress_int;
	public:
		Cell(Config&);
		//A l'avenir constructeur prend en argument un truc qui init/garde en memoire  les DOF controles
		Cell(double L,Config& config): L_(L),xc_(L/2.),yc_(L/2.){initCell(config);}
		~Cell(){};

		//Init les parametres BC et de la grille par utilisateur
		void initCell(Config&);
		//Controle force/vitesse
		void ApplyBC();
		//Maj de h,hdd,Ld,s apres a la fin du pas de temps
		void update(Tensor2x2,Tensor2x2,double);

		//Met a jours la periodicite des particules
		void PeriodicBoundaries(Particle&);

		void write(ofstream&,ofstream&,double);

		//Acces
		double getVolume();
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


//Le couplage passe par le tenseur de contrainte entre cellule et particules

//Update particules, calcul tenseur contrainte interne, update cellule

//On integre l'equation de la dynamique pour l'espace avec un algo de verlet
//Integre la dynamique de l'espace, de la cellule au temps t, a besoin tenseur contraintes
//Il y a couplage entre debut et fin du pas de temps
//Il faut integrer les deux en meme temps
void verletalgo(Cell& cell,std::vector<Particle>& ps, double dt){

	//Attention: variable globale algo hinv, stress ext et int
	//Init:
	Tensor2x2 h = cell.geth();
	Tensor2x2 hd = cell.gethd();
	Tensor2x2 hinv = h.getInverse();

	Tensor2x2 hdd0;
	Tensor2x2 hdd1 = hdd0;

	double mh=cell.getMasse();
	double V0=cell.getVolume();

	//Impose BC() en force: si une direction est controlee en vitesse, on impose hdd0 nulle dans cette dir
	//Et on calcule le tenseur de contrainte internes a partir de l'equation ac terme de droite nul
	//Cad que stressext=-stressint dans cette direction

	//Si dir controle en vitesse, stressext dans cette dir = stress int
	//On calcule les composantes de stressext dans cette direction
	//Apres on fait le produit matriciel
	//On recupere le stressint au debut du pas de temps (forces particules au debut du pas)
	// computeInternalStress()
	// Pour l'istant il est nul

	//On applique les BCU pour recuperer les contraintes ext
	cell.ApplyBC();
	Tensor2x2 TotalStress = cell.getStressInt() +cell.getStressExt();


	//Calcul de hdd au debut du pas, hdd0
	//Eq (24)
	//On a besoin de stressint au debut du pas de temps
	hdd0=hinv*(V0/mh)*(TotalStress);

	//Calcul position(h) avec hdd0 a la fin du pas de temps
	h = h + hd*dt + hdd0*0.5*dt*dt;

	//Calcul des positions particules fin du pas de temps

	for(std::vector<Particle>::iterator it = ps.begin(); it != ps.end(); it++){

		//Force pas de temps precedent:
		Vecteur F = it->getForce(); 
		Vecteur r = it->getR();
		Vecteur v = it->getV();
		double m = it->getMasse();

		//Acceleration debut pas de temps
		Vecteur tmp = hinv * r;
		tmp = TotalStress * tmp; 
		tmp = hinv * tmp;
		tmp = tmp * (V0/mh);
		//a0 a stocker
		Vecteur a0 = F * (1./m) + tmp; 
		r = r + v * dt + a0 * dt * dt * 0.5;
		//update position et acceleration debut pas temps
		it->setR(r,a0);
	}

	//On peut a present recalculer le tenseur des contraintes internes a partir des forces fin du pas de temps
	//Calculer les forces a partir nouvelle position
	//MD()
	//computeInternalStress()
	//Mettre a jout stressint de cell

	//Volume fin du pas de temps
	double V1=h.getDet();
	//Calcul hdd fin du pas, hdd1
	hinv=h.getInverse();

	//Ensuite on calcule l'accleration de chaque particule fin du pas de temps en utilisant le nouveau h,V
	//On a les vitesses des particules: position et vitesses particules fin pas de temps

	cell.ApplyBC(); // permet de definir le stressext fin pas de temps
	//TotalStress fin du pas de temps
	TotalStress = cell.getStressInt() +cell.getStressExt(); 

	TotalStress.affiche();
	//Derniere boucle pour calcul des vitesses:
	//Comment on applique les BCL en vitesse???
	for(std::vector<Particle>::iterator it = ps.begin(); it != ps.end(); it++){

		//Force recalculee fin du pas de temps:
		Vecteur F = it->getForce(); 
		//Position fin pas de temps
		Vecteur r = it->getR();
		//Acceleration debut pas de temps
		Vecteur a0 = it->geta0();
		double m = it->getMasse();

		//Acceleration fin pas de temps
		Vecteur tmp = hinv * r;
		tmp = TotalStress * tmp; 
		tmp = hinv * tmp;
		tmp = tmp * (V1/mh);
		Vecteur a1 = F * (1./m) + tmp; 
		a1.print();
		Vecteur v = it->getV() + (a1 + a0) * dt * 0.5 ;
		v.print();
		cout<<"vitesse:";
		//update vitesse
		it->setV(v);
	}
	//Fin integration particules


	//Acceleration fin du pas: pareil que pour hdd0, on ompose hdd1 nul dans direction controle en vitesse
	//Il faut prendre le nouveau tenseur contrainte interne et externe en appliqand les BCU fin pas de temps

	hdd1=hinv*(V1/mh)*(TotalStress);

	//Calcul de hd fin du pas de temps
	//Alors la vitesse ds dir controle en vitesse est bien constante (hdd0+hdd1=0)
	hd=hd+(hdd0+hdd1)*dt*0.5;

	//Update Cell:on fait remonter h,hd
	cell.update(h,hd,dt);
}

void write(std::vector<Particle> sp,ofstream& of,double t){
	of<<t<<" ";
	for(std::vector<Particle>::iterator it = sp.begin(); it != sp.end(); it++){
		it->write(of);
	}
}

int main(){

	//Parametres:
	double const L = 1.;
	double dt = .001 ;
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
	ofstream tmp3("strain.txt");

	double t=0.;
	int k=0;

	do{
		//Check for boundary:
		cell.PeriodicBoundaries(p);
		//Update :
		verletalgo(cell,sample,dt);

		//outputs:
		if( k% 1 == 0 ){
			cell.write(tmp2,tmp3,t);
			write(sample,tmp,t);
			//cout<<t<<" "<<cell.getVolume()<<endl;
		}

		t+=dt;
		k++;
	}while(t<T);

	tmp.close();
	tmp2.close();
	tmp3.close();

	return 0;
}

