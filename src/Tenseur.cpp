#include"Tenseur.hpp"

using namespace std;

void Tensor2x2::eigenValues(){
	double b= -yy_ - xx_;
	double c=-xy_*yx_+xx_*yy_;
	double d=b*b-4*c;
	l1_=(-b-sqrt(d))/2.;
	l2_=(-b+sqrt(d))/2.;
	double tmp=max(l1_,l2_);
	l2_=min(l1_,l2_);
	l1_=tmp;
}

//Suppose eigenVectors have been called before
double Tensor2x2::getMajorDirection() const{
	return (  u1_.gety() < 0. ? 2.*M_PI-acos(u1_.getx()): acos(u1_.getx()));
}

void Tensor2x2::eigenVectors(){

	eigenValues();

	double norme;
	double epsilon=0.0001;
	double vx1_,vy1_,vx2_,vy2_;

	if (abs(l1_)<epsilon && abs(l2_)<epsilon) 
	{
		return;
	}
	else
	{
	vx1_=xy_ / ( - xx_+ l1_);
	vy1_=1.;
	norme=sqrt(vx1_*vx1_ + vy1_*vy1_);
	vx1_/=norme;
	vy1_/=norme;
	
	vx2_=xy_/( - xx_ + l2_);
	vy2_=1.;
	norme=sqrt(vx2_*vx2_ + vy2_*vy2_);
	vx2_/=norme;
	vy2_/=norme;

	u1_.setx(vx1_);
	u1_.sety(vy1_);
	u2_.setx(vx2_);
	u2_.sety(vy2_);
	}
}

void Tensor2x2::writeEigenVectors(ofstream& of){
	of<<"0. 0. "<<u1_.getx()<<" "<<u1_.gety()<<" "<<u2_.getx()<<" "<<u2_.gety()<<endl;
}

//All information can be extracted from comp if needed later...
void Tensor2x2::write(ofstream& of) const{
	of<<xx_<<" "<<xy_<<" "<<yx_<<" "<<yy_;
}

void Tensor2x2::load(ifstream& is){
	is >> xx_ >> xy_ >> yx_ >> yy_;
}

void Tensor2x2::set(double xx, double xy, double yx, double yy){
	xx_ = xx;
	xy_ = xy;
	yx_ = yx;
	yy_ = yy;
}

Tensor2x2 Tensor2x2::getInverse() const{
	double epsilon=0.00001;
	//Matrice nulle par defaut
	Tensor2x2 Inverse;
	double det = getDet();
	if(abs(det)<epsilon){
		cerr<<"@La matrice n'a pas d'inverse!"<<endl;
		return Inverse;
	}
	else{
		Inverse.set(yy_/det,-xy_/det,-yx_/det,xx_/det);
		return Inverse;
	}
}

Tensor2x2 Tensor2x2::getTranspose() const{
	Tensor2x2 Transpose(xx_,yx_,xy_,yy_);
	return Transpose;
}
