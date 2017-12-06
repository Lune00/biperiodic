#include<vector>
#include<iostream>

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


	return 0;
}
