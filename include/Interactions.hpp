#ifndef ListeInteraction_hpp
#define ListeInteraction_hpp

//Manage verlet list, contact list


//Peut etre algo de calcul des forces egalement
//Car c'est different de l'algorithme (schema d'integration de Algo)

//A besoin de communiquer avec Sample
//On va le plug a Sample


//liste interactions potentielles

class Interactions{

	private:

		double dverlet_;
		double Rmax_;

	public:
		Interactions();
		~Interactions();



};



#endif
