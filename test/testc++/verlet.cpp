#include<iostream>
#include<fstream>
#include<cmath>

using namespace std;

//Integration du mouvement du pendule simple

int main(){

	double m=1.;
	double theta;
	double thetadot;
	double a1, a0;
	double dt=0.01;
	double t=0.;
	double T=100.;

	//Condition initiale:
	theta=3.;
	thetadot=0.;

	ofstream pendule("pendule.txt");
	pendule<<t<<" "<<theta<<" "<<thetadot<<endl;
	do{
		//Calcul de l acceleration au debut du pas de temps
		a0=-sin(theta)/m;
		//Update position:
		theta += dt*thetadot + 0.5*dt*dt*a0;
		//Calcul de l acceleration a fin du pas de temps:
		a1= -sin(theta)/m;
		//Update vitesse:
		thetadot += 0.5*dt*(a0+a1);
		t+=dt;
		pendule<<t<<" "<<theta<<" "<<thetadot<<endl;
	}while(t<T);
	pendule.close();
	return 0;
}
