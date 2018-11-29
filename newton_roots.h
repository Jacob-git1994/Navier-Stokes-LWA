#pragma once
#include "coe.h"
#include "functions.h"
#include <armadillo>

class newton_roots
{
private:
	double tol;
	int max_iterations;
	functions func;
	arma::vec init;
	arma::vec guess;
	arma::vec sol;
	//arma::sp_mat jacobian;
	arma::mat jacobian;
	arma::vec function_vector;
	void update_jacobian(arma::vec&,coe&);
	void update_function_vector(arma::vec&,coe&);
public:
	newton_roots(void);
	newton_roots(arma::vec&);
	double & tolerance(void);
	int & steps(void);
	void step(coe&);
	void update(void);
	arma::vec & profile(void);
};

newton_roots::newton_roots(void)
{
	tol = .00001;
	max_iterations = 20;
}

newton_roots::newton_roots(arma::vec &inital)
{
	tol = .0001;
	max_iterations = 20;
	init = inital;
    sol = arma::zeros<arma::vec>(size(init));
	guess = inital;
	jacobian = arma::zeros<arma::sp_mat>(init.n_rows,init.n_rows);
	//jacobian = arma::zeros<arma::mat>(init.n_rows,init.n_rows);
	function_vector = arma::zeros<arma::vec>(init.n_rows);
}

double & newton_roots::tolerance(void)
{
	return tol;
}

int & newton_roots::steps(void)
{
	return max_iterations;
}

void newton_roots::update_function_vector(arma::vec &g,coe &data)
{
	int end = init.n_rows - 1;
	/*Boundary i = 0*/
	function_vector(0) = func.eval_functions(g(end-2),g(end-1),g(0),g(1),g(2),init(end-2),init(end-1),init(0),init(1),init(2),data);
	/*Boundary i = 1*/
	function_vector(1) = func.eval_functions(g(end-1),g(0),g(1),g(2),g(3),init(end-1),init(0),init(1),init(2),init(3),data);
	/*Boundary i = N-1*/
	function_vector(end-1) = func.eval_functions(g(end-3),g(end-2),g(end-1),g(end),g(1),init(end-3),init(end-2),init(end-1),init(end),init(1),data);
	/*Boundary i = N*/
	function_vector(end) = func.eval_functions(g(end-2),g(end-1),g(end),g(1),g(2),init(end-2),init(end-1),init(end),init(1),init(2),data);
	
	/*Iterate Interior Grid*/
#pragma omp parallel num_threads(8)
	{
#pragma omp for simd schedule(static)
		for(int i = 2; i < end-1; i++)
		{
			function_vector(i) = func.eval_functions(g(i-2),g(i-1),g(i),g(i+1),g(i+2),init(i-2),init(i-1),init(i),init(i+1),init(i+2),data);
		}
	}	
}

void newton_roots::update_jacobian(arma::vec &g,coe &data)
{
	int end = init.n_rows-1;
	
	/*Boundaryi = 0*/
	jacobian(0,0) = func.eval_f(g(end-2),g(end-1),g(0),g(1),g(2),init(end-2),init(end-1),init(0),init(1),init(2),data);
	jacobian(0,1) = func.eval_fp1(g(end-2),g(end-1),g(0),g(1),g(2),init(end-2),init(end-1),init(0),init(1),init(2),data);
	jacobian(0,2) = func.eval_fp2(g(end-2),g(end-1),g(0),g(1),g(2),init(end-2),init(end-1),init(0),init(1),init(2),data);
	jacobian(0,end-1) = func.eval_fm1(g(end-2),g(end-1),g(0),g(1),g(2),init(end-2),init(end-1),init(0),init(1),init(2),data);
	jacobian(0,end-2) = func.eval_fm2(g(end-2),g(end-1),g(0),g(1),g(2),init(end-2),init(end-1),init(0),init(1),init(2),data);
	
	/*Boundary i = 1*/
	jacobian(1,0) = func.eval_fm1(g(end-1),g(0),g(1),g(2),g(3),init(end-1),init(0),init(1),init(2),init(3),data);
	jacobian(1,1) = func.eval_f(g(end-1),g(0),g(1),g(2),g(3),init(end-1),init(0),init(1),init(2),init(3),data);
	jacobian(1,2) = func.eval_fp1(g(end-1),g(0),g(1),g(2),g(3),init(end-1),init(0),init(1),init(2),init(3),data);
	jacobian(1,3) = func.eval_fp2(g(end-1),g(0),g(1),g(2),g(3),init(end-1),init(0),init(1),init(2),init(3),data);
	jacobian(1,end-1) = func.eval_fm2(g(end-1),g(0),g(1),g(2),g(3),init(end-1),init(0),init(1),init(2),init(3),data);
	
	/*Boundary i = N-1*/
	jacobian(end-1,end) = func.eval_fp1(g(end-3),g(end-2),g(end-1),g(end),g(1),init(end-3),init(end-2),init(end-1),init(end),init(1),data);
	jacobian(end-1,end-1) = func.eval_f(g(end-3),g(end-2),g(end-1),g(end),g(1),init(end-3),init(end-2),init(end-1),init(end),init(1),data);
	jacobian(end-1,end-2) = func.eval_fm1(g(end-3),g(end-2),g(end-1),g(end),g(1),init(end-3),init(end-2),init(end-1),init(end),init(1),data);
	jacobian(end-1,end-3) = func.eval_fm2(g(end-3),g(end-2),g(end-1),g(end),g(1),init(end-3),init(end-2),init(end-1),init(end),init(1),data);
	jacobian(end-1,1) = func.eval_fp2(g(end-3),g(end-2),g(end-1),g(end),g(1),init(end-3),init(end-2),init(end-1),init(end),init(1),data);
			
	/*Boundary i = N*/
	jacobian(end,end) = func.eval_f(g(end-2),g(end-1),g(end),g(1),g(2),init(end-2),init(end-1),init(end),init(1),init(2),data);
	jacobian(end,end-1) = func.eval_fm1(g(end-2),g(end-1),g(end),g(1),g(2),init(end-2),init(end-1),init(end),init(1),init(2),data);
	jacobian(end,end-2) = func.eval_fm2(g(end-2),g(end-1),g(end),g(1),g(2),init(end-2),init(end-1),init(end),init(1),init(2),data);
	jacobian(end,1) = func.eval_fp1(g(end-2),g(end-1),g(end),g(1),g(2),init(end-2),init(end-1),init(end),init(1),init(2),data);
    jacobian(end,2) = func.eval_fp2(g(end-2),g(end-1),g(end),g(1),g(2),init(end-2),init(end-1),init(end),init(1),init(2),data);

	/*Iterate over internal grid*/
#pragma omp parallel num_threads(8)
	{
#pragma omp for simd schedule(static)
		for (int i = 2; i < end-1; i++)
		{
			jacobian(i,i-2) = func.eval_fm2(g(i-2),g(i-1),g(i),g(i+1),g(i+2),init(i-2),init(i-1),init(i),init(i+1),init(i+2),data);
			jacobian(i,i-1) = func.eval_fm1(g(i-2),g(i-1),g(i),g(i+1),g(i+2),init(i-2),init(i-1),init(i),init(i+1),init(i+2),data);
			jacobian(i,i) = func.eval_f(g(i-2),g(i-1),g(i),g(i+1),g(i+2),init(i-2),init(i-1),init(i),init(i+1),init(i+2),data);
			jacobian(i,i+1) = func.eval_fp1(g(i-2),g(i-1),g(i),g(i+1),g(i+2),init(i-2),init(i-1),init(i),init(i+1),init(i+2),data);
			jacobian(i,i+2) = func.eval_fp2(g(i-2),g(i-1),g(i),g(i+1),g(i+2),init(i-2),init(i-1),init(i),init(i+1),init(i+2),data);
		}
	}
	//std::cout << jacobian(0,0);
}

void newton_roots::step(coe &data)
{
	int count = 0;
	arma::vec sol_2 = arma::zeros<arma::vec>(size(guess));
	while (count++ < max_iterations && arma::norm(sol_2 - guess) > tol)
	{
		/*Update the function vector and jacobian*/
		this->update_function_vector(guess,data);
		this->update_jacobian(guess,data);
		
		/*(Not working) Checks to see if the matrix is singular*/
		//if (rcond(jacobian) < 1e-8) {data.get_dt() /= 10.0;}
		//else {}
		
		sol = guess - arma::solve(jacobian,function_vector,arma::solve_opts::fast);//arma::solve_opts::fast//arma::spsolve(jacobian,function_vector);
		sol_2 = guess;
		guess = sol;
	}
	guess = sol;
}

void newton_roots::update(void)
{
	init = guess;
	guess = init;
}

arma::vec & newton_roots::profile(void)
{
	return init;
}

