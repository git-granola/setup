#ifndef HELP_H
#define HELP_H

#include <math.h>
#include <iostream>
#include <fstream>
using namespace std;

//help.h implementing box muller

/* boxmuller.c           Implements the Polar form of the Box-Muller
                         Transformation

                      (c) Copyright 1994, Everett F. Carter Jr.
                          Permission is granted by the author to use
			  this software for any application provided this
			  copyright notice is preserved.

*/

struct Vector{//As a compact way of storing coordinate vectors
	float v[2][1];
};

double ranf(double m){
	return (m*rand())/(double)RAND_MAX;
};

float box_muller(float m, float s)	/* normal random variate generator */
{				        /* mean m, standard deviation s */
	float x1, x2, w, y1;
	static float y2;
	static int use_last = 0;

	if (use_last)		        /* use value from previous call */
	{
		y1 = y2;
		use_last = 0;
	}
	else
	{
		do {
			x1 = 2.0 * ranf(m) - 1.0;
			x2 = 2.0 * ranf(m) - 1.0;
			w = x1 * x1 + x2 * x2;
		} while ( w >= 1.0 );

		w = sqrt( (-2.0 * log( w ) ) / w );
		y1 = x1 * w;
		y2 = x2 * w;
		use_last = 1;
	}

	return( m + y1 * s );
}

float determinant(Vector coord, float mu[2][1], float Sigma[2][2], float Pw){//calculates determinant

	float x_minus_mu_T[1][2];
	float x_minus_mu[2][1];

	x_minus_mu_T[0][0] = coord.v[0][0] - mu[0][0];
	x_minus_mu_T[0][1] = coord.v[1][0] - mu[1][0];

	x_minus_mu[0][0] = coord.v[0][0] - mu[0][0];
	x_minus_mu[1][0] = coord.v[1][0] - mu[1][0];

	float det_Sigma = 1.0/(Sigma[0][0]*Sigma[1][1] - Sigma[1][0]*Sigma[0][1]);

	float Sigma_inverse[2][2];
	Sigma_inverse[0][0] = Sigma[1][1];
	Sigma_inverse[1][1] = Sigma[0][0];
	Sigma_inverse[0][1] = -1*Sigma[1][0];
	Sigma_inverse[1][0] = -1*Sigma[0][1];

	float B[1][2];
	B[0][0] = x_minus_mu_T[0][0]*Sigma_inverse[0][0];
	B[0][1] = x_minus_mu_T[0][1]*Sigma_inverse[1][1];

	float allTogether = B[0][0]*x_minus_mu[0][0] + B[0][1]*x_minus_mu[1][0];//(x-mu)^T * Sigma^-1 * (x-mu), 1x2 * 2x2 * 2x1 vectors are transformed into a scalar

	return -0.5*det_Sigma*allTogether - 0.5 * log(det_Sigma) + log(Pw);
}

float det(float arr[2][2]){//just returns the determinant of an array, assuming that it is diagonal
	return 1.0/(arr[0][0]*arr[1][1] - arr[1][0]*arr[0][1]);
}

float Bhatt(float mu1[2][1], float mu2[2][1], float Sigma1[2][2], float Sigma2[2][2]){//assumes Sigma is diagonal
	//assuming beta = 0.5
	float mu[2][1], mu_T[1][2];
	mu[0][0] = mu1[0][0] - mu2[0][0];
	mu[1][0] = mu1[1][0] - mu2[1][0];

	mu_T[0][0] = mu[0][0];
	mu_T[0][1] = mu[1][0];

	float bSigma1[2][2]={0}, bSigma2[2][2] = {0}, bSigmaSum[2][2] = {0}, bSigmaSumI[2][2] = {0};
	
	bSigma1[0][0] = 0.5 * Sigma1[0][0];
	bSigma1[1][1] = 0.5 * Sigma1[1][1];

	bSigma2[0][0] = 0.5 * Sigma2[0][0];
	bSigma2[1][1] = 0.5 * Sigma2[1][1];

	bSigmaSum[0][0] = bSigma1[0][0] + bSigma2[0][0];
	bSigmaSum[1][1] = bSigma1[1][1] + bSigma2[1][1];

	bSigmaSumI[0][0] = det(bSigmaSum)*bSigmaSum[1][1];
	bSigmaSumI[1][1] = det(bSigmaSum)*bSigmaSum[0][0];

	float mu_SigmaI[1][2];
	mu_SigmaI[0][0] = mu_T[0][0]*bSigmaSumI[0][0]+mu_T[0][1]*bSigmaSumI[1][0];
	mu_SigmaI[0][1] = mu_T[0][0]*bSigmaSumI[0][1]+mu_T[0][1]*bSigmaSumI[1][1];

	float allTogether = mu_SigmaI[0][0]*mu[0][0] + mu_SigmaI[0][1]*mu[1][0];

	return exp(-(1.0/8)*allTogether + 0.5 * log(det(bSigmaSum) / ( sqrt(det(Sigma1)) * sqrt(det(Sigma2)))));//returns e^-k(.5), the Bhatt... bound
}

float Euclidean_distance(Vector coord, float mu[2][1]){//assumes Pw1 = Pw2, minimum distance classifier
	float x_minus_mu[2][1], x_minus_mu_T[1][2];
	x_minus_mu[0][0] = coord.v[0][0] - mu[0][0];
	x_minus_mu[1][0] = coord.v[1][0] - mu[1][0];

	x_minus_mu_T[0][0] = x_minus_mu[0][0];
	x_minus_mu_T[0][1] = x_minus_mu[1][0];

	return -1.0 * (x_minus_mu_T[0][0]*x_minus_mu[0][0] + x_minus_mu_T[0][1]*x_minus_mu[1][0]);
}


#endif

