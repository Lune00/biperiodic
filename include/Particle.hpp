#ifndef Particle_hpp
#define Particle_hpp

#include"Vecteur.hpp" 
#include<vector>
#include<fstream>

class Particle{
	private:
		//Position absolue
		Vecteur r_;
		//Vitesse:
		Vecteur v_;

		//TEMPORAIRE
		//Acceleration debut pas de temps:
		Vecteur a0_;
		//Forces:
		Vecteur f_;
		//END TEMMPORAIRE

		//Vecteur acc:
		Vecteur a_;
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
		Particle(): r_(), v_(), a_(), a0_(), f_(), t_(0.), vrot_(0.), rot_(0.), m_(1.), R_(1.) {}; 
		Particle(Vecteur r, Vecteur v) {r_ = r; v_ = v; init();}
		~Particle(){};

		void write(std::ofstream&);
		void init();
		void affiche();

		double getMasse() const {return m_;}
		double getRadius() const { return R_;}
		double getTorque() const {return t_;}

		Vecteur getR() const { return r_;}
		Vecteur getV() const { return v_;}
		Vecteur getForce()  { return a_*(1./m_);}
		Vecteur geta() const {return a_;}
		void setV(Vecteur v) { v_=v;}
		void setA(Vecteur a) { a_ = a;}
		void setR(Vecteur r, Vecteur a0) { r_ = r; a0_=a0;}
		void setRV(Vecteur r, Vecteur v) { r_ = r ; v_ = v; }
		void Periodize( double lx, double ly) { r_.add(lx,ly)  ;}
		
};



#endif
