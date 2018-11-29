#pragma once

class coe
{
private:
	double sigma;
	double mu;
	double theta;
	double gravity;
	double rho;
	double h_star;
	double contact_angle;
	double n;
	double m;
	double lambda_1;
	double del_sig;
	double kappa;
	double alpha;
	double dt;
	double dx;
public:
	coe(void);
	double & get_sigma(void) {return sigma;}
	double & get_mu(void) {return mu;}
	double & get_theta(void) {return theta;}
	double & get_gravity(void) {return gravity;}
	double & get_rho(void) {return rho;}
	double & get_h_star(void) {return h_star;}
	double & get_contact_angle(void) {return contact_angle;}
	double & get_n(void) {return n;}
	double & get_m(void) {return m;}
	double & get_lambda_1(void) {return lambda_1;}
	double & get_del_sig(void) {return del_sig;}
	double & get_kappa(void) {return kappa;}
	double & get_alpha(void) {return alpha;}
	double & get_dt(void) {return dt;}
	double & get_dx(void) {return dx;}
};

coe::coe(void)
{
	sigma = 1.0;
	mu = 1.0;
	theta = 0.0;
	gravity = 1.0;
	rho = 1.0;
	h_star = .0001;
	contact_angle = 0.0;
	n = 3.0;
	m = 2.0;
	lambda_1 = 0.0;
	del_sig = 0.0;
	kappa = 1.0;
	alpha = 1.0;
	dt = .001;
	dx = .001;
}

