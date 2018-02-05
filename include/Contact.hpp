#ifndef hpp_CONTACT_hpp
#define hpp_CONTACT_hpp

#include<iostream>
#include"Vecteur.hpp"
#include<algorithm>
#include<utility>
#include<iterator>

class Particle;
class Tensor2x2;
class Cell;

class Contact{

	private:

		Particle * i_;
		Particle * j_;

		//Contact position
		Vecteur r_;
		//Smalest branch vector through periodicity bet i and j
		Vecteur branch_;
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

		double rvrot_;

		//Pointer on the cell (to access h and hd)
		Cell * cell_;

		bool isActif_;

		//For a contact between an original and image particle
		//give indexes i in the x direction, and j in y direction
		std::pair<int,int> indexes_;

		//Avoid truncature errors
		//Used by frame for threshold for interpenetration
		static const double tolerance_;

	public:
		Contact(){ i_ = NULL; j_ =NULL, cell_=NULL ; isActif_ = false; }
		Contact(Particle* i, Particle* j, Cell*);
		Contact(std::ifstream&);
		//Contact(std::ifstream&,std::vector<Contact>&);
		~Contact(){};

		bool isActif() const { return isActif_;}
		void activate() { isActif_ = true;}
		void Frame();
		void write(std::ofstream&) const;
		void updateRelativeVelocity();
		void computeForce(const double,const double,const double,const double, const double, const double);
		void updateAccelerations();
		double sign(double x){if(x<0.) return -1.;else return 1.;}
		double getdn() const { return dn_;}
		double getdt() const { return dt_;}
		Vecteur getfxy() const { return (n_ *f_.getx() + t_ * f_.gety());} 

		//Debug:
		void print() const;
		const Particle* getj() const { return j_;}
		const Particle* geti() const { return i_;}

		void computeShortestBranch() ;
		Vecteur getbranch() const { return branch_;}
		double getfn() const { return f_.getx();}
		double getft() const { return f_.gety();}
		Vecteur getrv() const { return v_;}
		void setdt(const double dt) { dt_ = dt ;}

		//For post-processing:
		void setrv(const Vecteur rv) { v_ = rv ;}
		void setf(const Vecteur f) { f_ = f ;}

};



#endif
