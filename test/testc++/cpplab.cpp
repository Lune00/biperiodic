#include<vector>
#include<iostream>
#include<string>

using namespace std;


class A{

	private:
	vector<int> v_;
	public:
	A(){v_.push_back(1);}
	~A(){}
	vector<int>* getvector() { return &v_;}
	int getelement(int i) { return v_[i];}
	int getsize() { return v_.size();}

};

int main(){

	A a;
	vector<int> * u = a.getvector();
	u->push_back(2);
	cout<<a.getsize()<<endl;
	cout<<a.getelement(1)<<endl;
	cout<<u->at(1)<<endl;


	//Double acces:
	vector<int> s;
	//s.push_back(1);
//	s.push_back(2);
//	s.push_back(3);
	for(vector<int>::const_iterator iti = s.begin(); iti != s.end(); iti++){
		for(vector<int>::const_iterator itj = iti + 1; itj != s.end(); itj++){
			cout<<*iti<<" "<<*itj<<endl;
		}
	}


	int k = 2 ;
	string file = "file.txt";
	std::string kk = to_string(k);
	file = kk + file ;
	cout<<file<<endl;

	return 0;
}
