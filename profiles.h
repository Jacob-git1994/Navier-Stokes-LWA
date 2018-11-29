#pragma once
#include <armadillo>
#include <stdio.h>
#include <math.h>

class profiles
{
private:
	FILE *fluid;
	FILE *profile;
public:
	profiles(void);
	~profiles(void);
	void project(arma::vec&,arma::vec&,double);
};

profiles::profiles(void)
{
	fluid = fopen("fluid.txt","w");
	profile = fopen("profiles.txt","w");
}

profiles::~profiles(void)
{
	fclose(fluid);
	fclose(profile);
}

void profiles::project(arma::vec &x,arma::vec &data,double t)
{
	const int size = data.n_rows;
	arma::vec copy_x = x;
	arma::vec copy_data = data;
	copy_x /= pow(arma::min(data),2.0);
	copy_data /= arma::min(data);
	for (int i = 0; i < size; i++)
	{
		fprintf(fluid,"%.14f\t%.14f\t%.14f\n",x(i),data(i),t);
		fprintf(profile,"%.14f\t%.14f\t%.14f\n",copy_x(i),copy_data(i),t);
	}
	
}