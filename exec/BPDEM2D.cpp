#include<iostream>
#include<string>
#include<cmath>
#include<fstream>
#include<vector>

#include"Particle.hpp"
#include"Cell.hpp"
#include"Config.hpp"
#include"Algo.hpp"

using namespace std;

//Appartient a Sample
void write(std::vector<Particle> sp,ofstream& of,double t){
	of<<t<<" ";
	for(std::vector<Particle>::iterator it = sp.begin(); it != sp.end(); it++){
		it->write(of);
	}
}
//Temporaire
//Ecrit coordonnees absolues:
void write(std::vector<Particle> sp,Cell& cell,ofstream& of,double t){
	Tensor2x2 h = cell.geth();
	of<<t<<" ";
	for(std::vector<Particle>::iterator it = sp.begin(); it != sp.end(); it++){
		//it->affiche();
		h.affiche();
		Vecteur r = it->getR();
		r = h * r ;
		of<<r.getx()<<" "<<r.gety()<<endl;
	}

}

int main(){

	cout<<"BPDEM2D"<<endl;

	//Parametres:
	double const L = 1.;
	double dt = .01 ;
	double T = 1.;

	vector<Particle> sample;
	Config config;
	Cell cell(L,config);
	Algo algo(dt);

	//Initialise coordonnees(pos & vit) reduites
	Vecteur r0(0.5,0.5);
	Vecteur v0(0.5,0.);

	Particle p(r0,v0);
	sample.push_back(p);

	ofstream tmp("particle.txt");
	ofstream tmp2("cell.txt");
	ofstream tmp3("strain.txt");

	double t=0.;
	int k=0;

	do{
		//Check for boundary:
		cell.PeriodicBoundaries2(sample);
		//Update :
		algo.verletalgo2(cell,sample);

		//outputs:
		if( k% 1 == 0 ){
			cell.write(tmp2,tmp3,t);
			//write(sample,tmp,t);
			write(sample,cell,tmp,t);
			//cout<<t<<" "<<cell.getVolume()<<endl;
		}

		t+=dt;
		k++;
	}while(t<T);

	tmp.close();
	tmp2.close();
	tmp3.close();

	return 0;
}

