#ifndef Particle_hpp
#define Particle_hpp

#include"Vecteur.hpp" 
#include"Tenseur.hpp"
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
		void write(std::ofstream&,Tensor2x2&) const;
		void init();
		void affiche() const;
		void Periodize( double lx, double ly) { r_.add(lx,ly)  ;}

		void setV(Vecteur v) { v_=v;}
		void setA(Vecteur a) { a_ = a;}
		void setRV(Vecteur r, Vecteur v) { r_ = r ; v_ = v; }
		void setr(double sx, double sy) { r_.set(sx,sy) ;}
		void setMasse(double m) { m_ = m ;}
		void setInertia(double I) {I_ = I ;}

		void setrx(double x) { r_.setx(x);}
		void setry(double y) { r_.sety(y);}
		void addrx(double dx) {r_.addx(dx);}
		void addry(double dy) {r_.addy(dy);}


		double getx() const { return r_.getx() ; }
		double gety() const { return r_.gety() ; }
		double getMasse() const {return m_;}
		double getRadius() const { return R_;}
		double getTorque() const {return t_;}
		int getId() const { return id_;}

		Vecteur getR() const { return r_;}
		Vecteur getV() const { return v_;}
		Vecteur getForce()  { return a_*(1./m_);}
		Vecteur geta() const {return a_;}
		
};



#endif
