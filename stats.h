#pragma once
#include "coe.h"
#include <armadillo>
#include <stdio.h>

class stats
{
private:
	FILE *min;
	FILE *max;
	FILE *middle;
	FILE *mean;
	int index_min;
	int index_max;
	int index_middle;
public:
	stats(arma::vec&);
	~stats(void);
	void record(double,arma::vec&);
};

stats::stats(arma::vec &data)
{
	min = fopen("min_v_time.txt","w");
	max = fopen("max_v_time.txt","w");
	middle = fopen("middle_v_time.txt","w");
	mean = fopen("mean_v_time.txt","w");
	
	index_min = (int) data.index_min();
	index_max = (int) data.index_max();
	index_middle = (int) data.n_rows/2;
}

stats::~stats(void)
{
	fclose(min);
	fclose(max);
	fclose(middle);
}

void stats::record(double t,arma::vec &data)
{
	fprintf(min,"%.14f\t%.14f\n",t,data(index_min));
	fprintf(max,"%.14f\t%.14f\n",t,data(index_max));
	fprintf(middle,"%.14f\t%.14f\n",t,data(index_middle));
	fprintf(mean,"%.14f\t%.14f\n",t,arma::mean(data));
}