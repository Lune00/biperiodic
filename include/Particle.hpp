#ifndef Particle_hpp
#define Particle_hpp

#include"Vecteur.hpp" 
#include<vector>
#include<fstream>

class Particle{
	private:
		//identificateur:
		int id_;
		//mass:
		double m_;
		//Moment of inertia around particle center
		double I_;
		//Rayon:
		double R_;

		//Coordonnees reduites:
		Vecteur r_;
		//Vitesses:
		Vecteur v_;
		//Vecteur acc:
		Vecteur a_;
		//Torque:
		double t_;
		//Vrot:
		double vrot_;
		//Rot: angular position
		double rot_;

		
	public:
		Particle(): r_(), v_(), a_(), t_(0.), vrot_(0.), rot_(0.), m_(1.), R_(1.) {}; 
		Particle(Vecteur r, Vecteur v) {r_ = r; v_ = v; init();}
		~Particle(){};

		Particle(std::ifstream&);

		void write(std::ofstream&) const;
		void init();
		void affiche() const;

		double getMasse() const {return m_;}
		double getRadius() const { return R_;}
		double getTorque() const {return t_;}

		Vecteur getR() const { return r_;}
		Vecteur getV() const { return v_;}
		Vecteur getForce()  { return a_*(1./m_);}
		Vecteur geta() const {return a_;}
		void setV(Vecteur v) { v_=v;}
		void setA(Vecteur a) { a_ = a;}
		void setRV(Vecteur r, Vecteur v) { r_ = r ; v_ = v; }
		void Periodize( double lx, double ly) { r_.add(lx,ly)  ;}
		void setMasse(double m) { m_ = m ;}
		void setInertia(double I) {I_ = I ;}
		
};



#endif
