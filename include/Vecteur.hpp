#ifndef Vecteur_hpp
#define Vecteur_hpp

#include<iostream>
#include<fstream>
#include<cmath>

class Vecteur{
	private:
		double x_;
		double y_;
	public:
		Vecteur(): x_(0.), y_(0.){};
		Vecteur(double x, double y): x_(x), y_(y){};
		~Vecteur(){};

		void print() const{std::cerr<<x_<<" "<<y_<<" "<<getNorme()<<std::endl;}
		void add(double dx, double dy) { x_ += dx ; y_ += dy;}
		void addx(double dx) { x_ += dx;}
		void addy(double dy) { y_ += dy;}

		//Accessors:
		double getNorme() const {return sqrt( x_*x_ + y_ * y_);}
		double getNorme2() const {return  x_*x_ + y_ * y_;}

		double getx() const {return x_;}
		double gety() const {return y_;}
		//Mutators:
		void setx(double x) {x_ = x;}
		void sety(double y) {y_ = y;}
		void set(double x,double y){ x_=x; y_=y;}
		//Read from ifstream
		void load(std::ifstream& is);
		//Write in ofstream
		void write(std::ofstream& os) const;

		//Surcharge operator:

		//Scalar product:
		//double operator * (Vecteur a){
		//	return (x_*a.x_ + y_*a.y_);
		//}

		friend double operator * (const Vecteur& a, const Vecteur& b){
		  return (a.getx() * b.getx() + a.gety() * b.gety());
		}

		//Vecteur operator + (Vecteur a){
		//	Vecteur p;
		//	p.x_ = x_ + a.x_;
		//	p.y_ = y_ + a.y_;
		//	return p;
		//}

		friend Vecteur operator + (const Vecteur& a, const Vecteur& b)
		{
		  return Vecteur(a.x_ + b.x_, a.y_ + b.y_ );
		}
		friend Vecteur operator - (const Vecteur& a, const Vecteur& b)
		{
		  return Vecteur(a.x_ - b.x_, a.y_ - b.y_ );
		}

		//Vecteur operator - (Vecteur a){
		//	Vecteur p;
		//	p.x_ = x_ - a.x_;
		//	p.y_ = y_ - a.y_;
		//	return p;
		//}

		//Vecteur operator * (double a){
		//	Vecteur p;
		//	p.x_ = a*x_ ;
		//	p.y_ = a*y_ ;
		//	return p;
		//}

		//Vecteur operator * (double a) const{
		//	Vecteur p;
		//	p.x_ = a*x_ ;
		//	p.y_ = a*y_ ;
		//	return p;
		//}

		friend Vecteur operator * (const Vecteur& a, const double k)
		{
		  return Vecteur(a.getx()*k, a.gety()*k );
		}

		friend Vecteur operator / (const Vecteur& a, const double k)
		{
		  double invk = 1. / k ;
		  return Vecteur(a.getx()*invk, a.gety()*invk );
		}

		friend Vecteur operator * (const Vecteur& a, const int k)
		{
		  return Vecteur(a.getx()*(double)k, a.gety()*(double)k );
		}

		//Vecteur operator * (int a) {
		//	Vecteur p;
		//	p.x_ = (double)a * x_ ;
		//	p.y_ = (double)a * y_ ;
		//	return p;
		//}

		//Vecteur operator / (double a){
		//	Vecteur p;
		//	p.x_ = x_ / a ;
		//	p.y_ = y_ / a ;
		//	return p;
		//}

		Vecteur operator - () {
			Vecteur p;
			p.x_ = - x_;
			p.y_ = - y_;
			return p;
		}
		Vecteur& operator += (const Vecteur& a){
		  x_ += a.x_;
		  y_ += a.y_;
		  return *this;
		}
		Vecteur& operator *= (const double e){
		  x_ *= e;
		  y_ *= e;
		  return *this;
		}
};

#endif
