#ifndef Tensor_hpp
#define Tensor_hpp

#include<fstream>
#include<iostream>

#include"Vecteur.hpp"

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

		Tensor2x2 getInverse() const;
		Tensor2x2 getTranspose() const;

		void print() const{std::cout<<std::endl;std::cout<<xx_<<" "<<xy_<<std::endl;
			std::cout<<yx_<<" "<<yy_<<std::endl;std::cout<<std::endl;}

		//Accesseurs:
		double getxx() const { return xx_;} 
		double getxy() const { return xy_;} 
		double getyx() const { return yx_;} 
		double getyy() const { return yy_;} 

		double getl1() const { return l1_;}
		double getl2() const { return l2_;}
		double getMajorDirection() const;

		double getDet() const{ return (xx_*yy_ - xy_*yx_) ;}

		void setxx(double xx) { xx_ = xx;}
		void setxy(double xy) { xy_ = xy;}
		void setyx(double yx) { yx_ = yx;}
		void setyy(double yy) { yy_ = yy;}

		void eigenValues();
		void eigenVectors();
		void write(std::ofstream&) const;
		void load(std::ifstream&);
		void writeEigenVectors(std::ofstream&);

		friend Tensor2x2 operator * (const Tensor2x2& a, const Tensor2x2& b){

		  return Tensor2x2(a.getxx() * b.getxx() +a.getxy()*b.getxy(), a.getxx() * b.getxy() + a.getxy()*b.getyy(), a.getyx()*b.getxx() + a.getyy() * b.getyx(), a.getyx()*b.getxy() + a.getyy()*b.getyy());
		}

		friend Tensor2x2 operator * (const Tensor2x2& a, const double k){
		  return Tensor2x2(a.getxx() * k, a.getxy()* k , a.getyx() * k , a.getyy() * k);
		}

		friend Tensor2x2 operator + (const Tensor2x2& a, const Tensor2x2& b){
		  return Tensor2x2(a.getxx() + b.getxx(), a.getxy()+b.getxy(), a.getyx() + b.getyx() , a.getyy() + b.getyy());
		}

		friend Vecteur operator * (const Tensor2x2& t, const Vecteur& v){
		  return Vecteur(t.getxx() * v.getx() + t.getxy()*v.gety(), t.getyx() * v.getx() + t.getyy() * v.gety());
		}

		Tensor2x2& operator += (const Tensor2x2& a){
		  xx_ += a.getxx();
		  xy_ += a.getxy();
		  yx_ += a.getyx();
		  yy_ += a.getyy();
		  return *this;
		}

};

#endif
