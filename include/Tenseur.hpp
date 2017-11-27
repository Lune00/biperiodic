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
		double getDet(){ return (xx_*yy_ - xy_*yx_) ;}
		Tensor2x2 getInverse();
		Tensor2x2 getTranspose();
		void affiche(){std::cout<<std::endl;std::cout<<xx_<<" "<<xy_<<std::endl;
			std::cout<<yx_<<" "<<yy_<<std::endl;std::cout<<std::endl;}

		//Accesseurs:
		double getxx() const { return xx_;} 
		double getxy() const { return xy_;} 
		double getyx() const { return yx_;} 
		double getyy() const { return yy_;} 

		//Renvoie la composante i du tenseur:
		//0->xx 1->xy 2->yx 3->yy

		void eigenValues();
		void eigenVectors();
		void write(std::ofstream&);
		void writeEigenVectors(std::ofstream&);

		Tensor2x2 operator * (Tensor2x2 a){
			Tensor2x2 p;
			p.xx_ = xx_*a.xx_ + xy_*a.yx_;
			p.xy_ = xx_*a.xy_ + xy_*a.yy_;
			p.yx_ = yx_*a.xx_ + yy_*a.yx_;
			p.yy_ = yx_*a.xy_ + yy_*a.yy_;
			return p;
		}

		Tensor2x2 operator * (double a){
			Tensor2x2 p;
			p.xx_ = xx_*a;
			p.xy_ = xy_*a;
			p.yx_ = yx_*a;
			p.yy_ = yy_*a;
			return p;
		}

		Tensor2x2 operator + (Tensor2x2 a){
			Tensor2x2 p;
			p.xx_ = xx_+a.xx_;
			p.xy_ = xy_+a.xy_;
			p.yx_ = yx_+a.yx_;
			p.yy_ = yy_+a.yy_;
			return p;
		}

		Vecteur operator * (Vecteur a){
			Vecteur p;
			double x = xx_*a.getx()+ xy_*a.gety();
			double y = yx_*a.getx() + yy_*a.gety();
			p.set(x,y);
			return p;
		}
};

#endif
