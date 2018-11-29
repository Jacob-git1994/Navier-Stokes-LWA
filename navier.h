#pragma once
#include "coe.h"
#include "newton_roots.h"
#include "stats.h"
#include <stdio.h>
#include <armadillo>
#include "profiles.h"

class navier
{
private:
	/*Parameters*/
	/*When to sample*/
	int samp;
	/*Size of the system*/
	int size;
	/*When to stop time*/
	double end_time;
	/*Tolerence and max iteration count for the root finding method*/
	double tol;
	int max_iterations;
	/*Store the data of the paremeters*/
	coe data;
	/*Initalize the grid and vectors*/
	arma::vec grid;
	arma::vec inital;
public:
	/*Constructors*/
	navier(void);
	navier(double,double,double (*func)(double),int,double);
	navier(double,double,arma::vec&,double);
	/*Return reference to the variables*/
	int & samples(void);
	double & set_tolerence(void);
	int & set_max_iterations(void);
	double & sigma(void) {return data.get_sigma();}
	double & mu(void) {return data.get_mu();}
	double & theta(void) {return data.get_theta();}
	double & gravity(void) {return data.get_gravity();}
	double & rho(void) {return data.get_rho();}
	double & h_star(void) {return data.get_h_star();}
	double & contact_angle(void) {return data.get_contact_angle();}
	double & n(void) {return data.get_n();}
	double & m(void) {return data.get_m();}
	double & lambda_1(void) {return data.get_lambda_1();}
	double & del_sig(void) {return data.get_del_sig();}
	double & kappa(void) {return data.get_kappa();}
	double & alpha(void) {return data.get_alpha();}
	double & dt(void) {return data.get_dt();}
	/*Start solving*/
	void run(void);
};

navier::navier(void)
{
	size = 0;
}

navier::navier(double a,double b,double (*func)(double),int N,double time_)
{
	size = N;
	tol = .0001;
	max_iterations = 10;
	end_time = time_;
	grid = arma::linspace<arma::vec>(a,b,N);
	data.get_dx() = grid(1) - grid(0); 
	inital = arma::zeros<arma::vec>(arma::size(grid));
#pragma omp parallel num_threads(8)
	{
#pragma omp for simd
		for (int i = 0; i < grid.n_rows; i++)
		{
			inital(i) = func(grid(i));
		}
	}
}

navier::navier(double a,double b,arma::vec &intial_condition,double time_)
{
	size = intial_condition.n_rows;
	tol = .0001;
	max_iterations = 10;
	end_time = time_;
	grid = arma::linspace<arma::vec>(a,b,size);
	data.get_dx() = grid(1)- grid(0);
	inital = intial_condition;
}

int & navier::samples(void)
{
	return samp;
}

double & navier::set_tolerence(void)
{
	return tol;
}

int & navier::set_max_iterations(void)
{
	return max_iterations;
}

void navier::run(void)
{
	int count = 0;
	int count_2 = (int) (end_time/data.get_dt())/samp;
	FILE *f = fopen("output.txt","w");
	newton_roots sim(inital);
	stats stat(inital);
	profiles pro;
	sim.tolerance() = tol;
	sim.steps() = max_iterations;
	while ((count++)*data.get_dt() <= end_time)
	{
		stat.record((count-1)*data.get_dt(),sim.profile());
		if (count_2 == (int) (end_time/data.get_dt())/samp)
		{
			pro.project(grid,sim.profile(),(count-1)*data.get_dt());
			std::cout << "Running\t" << ((count-1)*data.get_dt()/end_time)*100 << "% Done\n";
			for (int i = 0; i < size; i++)
			{ 
				fprintf(f,"%.14f\t%.14f\n\n",grid(i),(sim.profile())(i));
				count_2 = 0;
			}
		}
		else 
		{
			count_2++;
		}

		sim.step(data);
		sim.update();
	}	
	/*Finish the Program*/
	std::cout << "Finished! " << 100.0 << "% Done\n";
	fclose(f);
}
	
	
	
	
	
	
	
	
