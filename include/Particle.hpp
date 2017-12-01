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
		//Reduced coordinate
		Vecteur s_;
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
		Particle(): r_(), v_(), s_(),a_(), a0_(), f_(), t_(0.), vrot_(0.), rot_(0.), m_(1.), R_(1.) {}; 
		Particle(Vecteur r, double L, Vecteur v) {r_ = r; v_ = v; init(L);}
		~Particle(){};

		void write(std::ofstream&);
		void init(double L);
		void affiche();

		double getMasse() const {return m_;}
		double getRadius() const { return R_;}
		double getTorque() const {return t_;}

		Vecteur getR() const { return r_;}
		Vecteur getV() const { return v_;}
		Vecteur getForce() const { return f_;}
		Vecteur geta0() const {return a0_;}
		Vecteur geta() const {return a_;}
		void setV(Vecteur v) { v_=v;}
		void setR(Vecteur r, Vecteur a0) { r_ = r; a0_=a0;}
		void setRV(Vecteur r, Vecteur v) { r_ = r ; v_ = v }
		
};



#endif
