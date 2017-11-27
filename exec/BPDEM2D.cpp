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

int main(){

	//Parametres:
	double const L = 1.;
	double dt = .001 ;
	double T = 1.;

	vector<Particle> sample;
	Config config;
	Cell cell(L,config);
	Algo algo(dt);

	//Initialise coordonnees reduites, absolues, vitesse fluctuante
	Vecteur r0(L/3,L/3);
	Vecteur v0(0.,0.);
	Particle p(r0,L,v0);
	sample.push_back(p);

	ofstream tmp("particle.txt");
	ofstream tmp2("cell.txt");
	ofstream tmp3("strain.txt");

	double t=0.;
	int k=0;

	do{
		//Check for boundary:
		cell.PeriodicBoundaries(p);
		//Update :
		algo.verletalgo(cell,sample);

		//outputs:
		if( k% 1 == 0 ){
			cell.write(tmp2,tmp3,t);
			write(sample,tmp,t);
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

