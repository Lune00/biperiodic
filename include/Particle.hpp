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
		//radius:
		double R_;

		//I should have done 3d vector to store rotation too
		//But not consistent with Cell tensor operations...

		//Reduced coordinates:
		Vecteur r_;
		//Vitesses:
		Vecteur v_;
		//Vecteur acc:
		Vecteur a_;

		//Rotation:
		//Time rate of change of rotational velocity
		double arot_;
		//Vrot:
		double vrot_;
		//Rot: angular position
		double rot_;

		
	public:
		Particle(): r_(), v_(), a_(), arot_(0.), vrot_(0.), rot_(0.), m_(1.), R_(1.), I_(0.5*m_*R_*R_) {}; 
		Particle(Vecteur r, Vecteur v) {r_ = r; v_ = v; init();}
		~Particle(){};

		Particle(std::ifstream&);

		void write(std::ofstream&) const;
		void write(std::ofstream&,const Tensor2x2&,const Tensor2x2&) const;
		void init();
		void load(std::ifstream&);
		void print() const;
		void Periodize( double lx, double ly) { r_.add(lx,ly)  ;}

		void setr(double sx, double sy) { r_.set(sx,sy) ;}
		void setInertia(const double m);

		void setrx(double x) { r_.setx(x);}
		void setry(double y) { r_.sety(y);}
		void addrx(double dx) {r_.addx(dx);}
		void addry(double dy) {r_.addy(dy);}

		//NOT USED YET
		void removevmean(const Vecteur& vmean) { v_ = v_ - vmean;}

		//Integration methods:
		void updateA(Vecteur f) { a_ = a_ + f / m_;}
		//void updateA(Vecteur a) { a_ = a_ + a ;}
		void updateArot(double Torque) { arot_ = arot_ + Torque/I_;}
		void updateR(const double dt);
		void updateV(const double dt);
		void updateVrot(const double dt);
		void updateRot(const double dt);
		void resetA(); 
		void removeHpart(const Tensor2x2&);

		double getx() const { return r_.getx() ; }
		double gety() const { return r_.gety() ; }
		double getMasse() const {return m_;}
		double getI() const { return I_;}
		double getRadius() const { return R_;}
		//Rotation accessors
		double getArot() const {return arot_;}
		double getRot() const { return rot_;}
		double getVrot() const { return vrot_;}
		int getId() const { return id_;}

		Vecteur getR() const { return r_;}
		Vecteur getV() const { return v_;}
		Vecteur getA() const {return a_;}
};



#endif
