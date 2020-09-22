#include "help.h"
#include<iostream>

#define READ 1
#define WRITE 0

using namespace std;

int main(){
//I'm just calculating everything by hand.
//Initializing========================================================================

	Vector *class1 = new Vector[100000];
	Vector *class2 = new Vector[100000];

	Vector *class1_2 = new Vector[100000];
	Vector *class2_2 = new Vector[100000];

	if(!READ){
		for(int i=0; i<100000; i++){
			class1[i].v[0][0] = box_muller(1,1);
			class1[i].v[1][0] = box_muller(1,1);	
			class2[i].v[0][0] = box_muller(4,1);
			class2[i].v[1][0] = box_muller(4,1);

			class1_2[i].v[0][0] = box_muller(1,1);
			class1_2[i].v[1][0] = box_muller(1,1);
			class2_2[i].v[0][0] = box_muller(4,2);
			class2_2[i].v[1][0] = box_muller(4,sqrt(8));
		}
	}

	if(READ){
		ifstream source("data.txt");
		if(source.is_open()){
			for(int i=0; i<100000; i++){
				source >> class1[i].v[0][0];
				source >> class1[i].v[1][0];
				source >> class2[i].v[0][0];
				source >> class2[i].v[1][0];

				source >> class1_2[i].v[0][0];
				source >> class1_2[i].v[1][0];
				source >> class2_2[i].v[0][0];
				source >> class2_2[i].v[1][0];
			}
			source.close();
		}
		else{
			cout << "File not found."<<endl;
		}
	}

	if(WRITE){
		ofstream dest("data.txt");
		if(dest.is_open()){
			for(int i=0; i<100000; i++){
				dest << class1[i].v[0][0] << ' ' << class1[i].v[1][0] << endl;
				dest << class2[i].v[0][0] << ' ' << class2[i].v[1][0] << endl;
				dest << class1_2[i].v[0][0] << ' ' << class1_2[i].v[1][0] << endl;
				dest << class2_2[i].v[0][0] << ' ' << class2_2[i].v[1][0] << endl;
			}
		}
		dest.close();
	}

//1========================================================================
//Assuming P(w1) = P(w2):
//Decision boundary:
	int c1_error=0, c2_error=0;

	float mu1[2][1], mu2[2][1], Sigma[2][2] ={0};
	mu1[0][0] = 1;
	mu1[1][0] = 1;
	mu2[0][0] = 4;
	mu2[1][0] = 4;
	Sigma[0][0] = 1;
	Sigma[1][1] = 1;

	for(int i=0; i<100000; i++){
				c1_error += determinant(class1[i], mu1, Sigma, 0.5) < determinant(class1[i], mu2, Sigma, 0.5);
				c2_error += determinant(class2[i], mu2, Sigma, 0.5) < determinant(class2[i], mu1, Sigma, 0.5);
	}

	cout << "Problem 1:"<<endl;
	cout << "P(w1) = P(w2)"<<endl;
	cout << "Class 1 errors: "<< c1_error << endl << "Class 2 errors: " << c2_error << endl << "Total errors: "<<c1_error+c2_error<<endl<<endl;

	int c1_error_b=0, c2_error_b=0;

	for(int i=0; i<100000; i++){
		c1_error_b += determinant(class1[i], mu1, Sigma, 0.2) < determinant(class1[i], mu2, Sigma, 0.8);
		c2_error_b += determinant(class2[i], mu2, Sigma, 0.8) < determinant(class2[i], mu1, Sigma, 0.2);
	}

	cout << "P(w1) = 0.2, P(w2) = 0.8"<<endl;
	cout << "Class 1 errors: "<< c1_error_b << endl << "Class 2 errors: " << c2_error_b << endl << "Total errors: "<<c1_error_b+c2_error_b<<endl<<endl;

	cout<<"#1 Bhattacharryya = "<<Bhatt(mu1, mu2, Sigma, Sigma)<<endl<<endl;

//2========================================================================
//P(w1) = P(w2)

	int c1_2_error=0, c2_2_error=0;
	float Sigma2[2][2] = {0};
	Sigma2[0][0] = 4;
	Sigma2[1][1] = 8;

	for(int i=0; i<100000; i++){
		c1_2_error += determinant(class1_2[i], mu1, Sigma, 0.5) < determinant(class1_2[i], mu2, Sigma2, 0.5);
		c2_2_error += determinant(class2_2[i], mu2, Sigma2, 0.5) < determinant(class2_2[i], mu1, Sigma, 0.5);
	}

	cout << "Problem 2:"<<endl;
	cout << "P(w1) = P(w2)"<<endl;
	cout << "Class 1 errors: "<< c1_2_error << endl << "Class 2 errors: " << c2_2_error << endl << "Total errors: "<<c1_2_error+c2_2_error<<endl<<endl;

//P(w1) = 0.2, P(w2) = 0.8

	int c1_2_b_error=0, c2_2_b_error=0;

	for(int i=0; i<100000; i++){
		c1_2_b_error += determinant(class1_2[i], mu1, Sigma, 0.2) < determinant(class1_2[i], mu2, Sigma2, 0.8);
		c2_2_b_error += determinant(class2_2[i], mu2, Sigma2, 0.8) < determinant(class2_2[i], mu1, Sigma, 0.2);		
	}

	cout << "P(w1) = 0.2, P(w2) = 0.8"<<endl;
	cout << "Class 1 errors: "<< c1_2_b_error << endl << "Class 2 errors: " << c2_2_b_error << endl << "Total errors: "<<c1_2_b_error+c2_2_b_error<<endl<<endl;

	cout<< "#2 Bhattacharyya = "<<Bhatt(mu1, mu2, Sigma, Sigma2)<<endl<<endl;

	int c1_3_error=0, c2_3_error=0;

	for(int i=0; i<100000; i++){
		c1_3_error += Euclidean_distance(class1_2[i], mu1) < Euclidean_distance(class1_2[i], mu2);
		c2_3_error += Euclidean_distance(class2_2[i], mu2) < Euclidean_distance(class2_2[i], mu1);
	}

	cout<<"Problem 3:"<<endl;
	cout<<"Class 1 errors: "<<c1_3_error<<endl<<"Class 2 errors: "<< c2_3_error<<endl<<"Total errors: "<<c1_3_error+c2_3_error<<endl;


	delete class1;
	delete class2;
	delete class1_2;
	delete class2_2;

	return 0;
}
