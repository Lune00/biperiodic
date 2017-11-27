#ifndef Vecteur_hpp
#define Vecteur_hpp

#include<iostream>
#include<cmath>

class Vecteur{
	private:
		double x_;
		double y_;
	public:
		Vecteur(): x_(0.), y_(0.){};
		Vecteur(double x, double y): x_(x), y_(y){};
		~Vecteur(){};

		void print(){std::cout<<x_<<" "<<y_<<std::endl;}
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

#endif
