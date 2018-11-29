#include <math.h>
#include "navier.h"
#include <iostream>
#include <armadillo>

void initalize(arma::vec &h)
{
	int N = h.n_rows;
	int index = 0;
	int index_2 = N-1;
	double shift = 7.8;
	arma::vec half_grid = arma::linspace<arma::vec>(-10,10,(int) N/2);
	for (int i = 0; i < (int)N/2; i++)
	{
		h(index++) = (-tanh(half_grid(i) - shift)+1.0)/2.0;
		h(index_2--) = (-tanh(half_grid(i) - shift)+1.0)/2.0;
	}
	//h((int)N/2) = (-tanh(half_grid((int)N/2 - shift))+1.0)/2.0;
}

int main()
{
	const unsigned int size = 500;
	arma::vec inital = arma::zeros<arma::vec>(size);
	initalize(inital);
	navier simulation(-10,10,inital,2);
	simulation.dt() = .001;
	simulation.sigma() = 1.0;
	simulation.mu() = 1.0;
	simulation.gravity() = 0;
	simulation.rho() = 1;
	simulation.samples() = 10;
	simulation.lambda_1() = 20.0;
	simulation.h_star() = .01;
	simulation.n() = 3;
	simulation.m() = 2;
	simulation.contact_angle() = 0.0;
	simulation.set_tolerence() = .001;
	simulation.set_max_iterations() = 10;
	simulation.del_sig() = 0.0;
	simulation.run();
}
