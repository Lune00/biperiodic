#ifndef hpp_CONTACT_hpp
#define hpp_CONTACT_hpp

#include<iostream>
#include"Vecteur.hpp"

class Particle;
class Tensor2x2;
class Cell;

class Contact{

	private:
		//Should be Particle * const i_/j_
		Particle * i_;
		Particle * j_;
		//Contact position
		Vecteur r_;
		//Relative velocity at contact (normal,tangent) in the
		//contact frame
		Vecteur v_;

		//Contact frame:
		//Normal to the contact
		Vecteur n_;
		//Tangent to the contact
		Vecteur t_;

		//Force at the contact
		Vecteur f_;

		//Interpenetration:
		double dn_;
		double dt_;

		//Relative rate of rotation at contact
		double rvrot_;

		//Pointer on the cell (to access h and hd)
		Cell * cell_;

		bool isActif_;
		//Avoid truncature errors
		//Used by frame for threshold for interpenetration
		static const double tolerance_;

	public:
		Contact(){ i_ = NULL; j_ =NULL, cell_=NULL ; isActif_ = false; }
		Contact(Particle* i, Particle* j, Cell&);
		~Contact(){};

		bool isActif() const { return isActif_;}
		void activate() { isActif_ = true;}
		void Frame();
		void write(std::ofstream&) const;
		void updateRelativeVelocity();
		void computeForce(const double,const double,const double,const double, const double, const double);
		void updateAccelerations();
		double sign(double x){if(x<0.) return -1.;else return 1.;}
		Vecteur getfxy() const { return (n_ *f_.getx() + t_ * f_.gety());} 
		double getdn() const { return dn_;}
		//Debug:
		void print() const;
		const Particle* getj() const { return j_;}
		const Particle* geti() const { return i_;}
};



#endif
