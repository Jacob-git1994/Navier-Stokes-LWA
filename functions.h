#pragma once
#include "coe.h"
#include <cmath>
#include <math.h>

class functions
{
private:
	double main_function(double hm2, double hm1, double h, double hp1, double hp2, double hIm2, double hIm1, double hI, double hIp1, double hIp2, double mew,double sigma, double n, double m, double h_star, double theta,double thetaI, double lam1, double dx, double dt, double gr, double p, double delsig, double alpha, double kapa);
	double fm2(double hm2, double hm1, double h, double hp1, double hp2, double hIm2, double hIm1, double hI, double hIp1, double hIp2, double mew,double sigma, double n, double m, double h_star, double theta,double thetaI, double lam1, double dx, double dt, double gr, double p, double delsig, double alpha, double kapa);
	double fm1(double hm2, double hm1, double h, double hp1, double hp2, double hIm2, double hIm1, double hI, double hIp1, double hIp2, double mew,double sigma, double n, double m, double h_star, double theta,double thetaI, double lam1, double dx, double dt, double gr, double p, double delsig, double alpha, double kapa);
	double f(double hm2, double hm1, double h, double hp1, double hp2, double hIm2, double hIm1, double hI, double hIp1, double hIp2, double mew,double sigma, double n, double m, double h_star, double theta,double thetaI, double lam1, double dx, double dt, double gr, double p, double delsig, double alpha, double kapa);
	double fp1(double hm2, double hm1, double h, double hp1, double hp2, double hIm2, double hIm1, double hI, double hIp1, double hIp2, double mew,double sigma, double n, double m, double h_star, double theta,double thetaI, double lam1, double dx, double dt, double gr, double p, double delsig, double alpha, double kapa);
	double fp2(double hm2, double hm1, double h, double hp1, double hp2, double hIm2, double hIm1, double hI, double hIp1, double hIp2, double mew,double sigma, double n, double m, double h_star, double theta,double thetaI, double lam1, double dx, double dt, double gr, double p, double delsig, double alpha, double kapa);
public:
	double eval_functions(double,double,double,double,double,double,double,double,double,double,coe&);
	double eval_fm2(double,double,double,double,double,double,double,double,double,double,coe&);
	double eval_fm1(double,double,double,double,double,double,double,double,double,double,coe&);
	double eval_f(double,double,double,double,double,double,double,double,double,double,coe&);
	double eval_fp1(double,double,double,double,double,double,double,double,double,double,coe&);
	double eval_fp2(double,double,double,double,double,double,double,double,double,double,coe&);
};

double functions::main_function(double hm2, double hm1, double h, double hp1, double hp2, double hIm2, double hIm1, double hI, double hIp1, double hIp2, double mew,double sigma, double n, double m, double h_star, double theta,double thetaI, double lam1, double dx, double dt, double gr, double p, double delsig, double alpha, double kapa)
{
    double a;
     double b_a;
     double c_a;
     double d_a;
     double e_a;
     double f_a;
     double g_a;
     double h_a;
     double i_a;
     double j_a;
     double k_a;
     double l_a;
     double m_a;
     double n_a;
     double o_a;
     double p_a;
     double q_a;
     double r_a;
     double s_a;
     double t_a;
     double u_a;
     double v_a;
     double w_a;
     double x_a;
     double y_a;
     double ab_a;
     double bb_a;
     double cb_a;
     double db_a;
     double eb_a;
     a = h / 2.0 + hm1 / 2.0;
     b_a = h / 2.0 + hm1 / 2.0;
     c_a = h / 2.0 + hp1 / 2.0;
     d_a = h / 2.0 + hp1 / 2.0;
     e_a = hI / 2.0 + hIm1 / 2.0;
     f_a = hI / 2.0 + hIm1 / 2.0;
     g_a = hI / 2.0 + hIp1 / 2.0;
     h_a = hI / 2.0 + hIp1 / 2.0;
     i_a = h / 2.0 + hm1 / 2.0;
     j_a = h / 2.0 + hm1 / 2.0;
     k_a = h / 2.0 + hp1 / 2.0;
     l_a = h / 2.0 + hp1 / 2.0;
     m_a = hI / 2.0 + hIm1 / 2.0;
     n_a = hI / 2.0 + hIm1 / 2.0;
     o_a = hI / 2.0 + hIp1 / 2.0;
     p_a = hI / 2.0 + hIp1 / 2.0;
     q_a = h + hp1;
     r_a = h / 2.0 + hp1 / 2.0;
     s_a = h / 2.0 + hp1 / 2.0;
     t_a = h + hm1;
     u_a = h / 2.0 + hm1 / 2.0;
     v_a = h / 2.0 + hm1 / 2.0;
     w_a = h + hm1;
     x_a = h + hp1;
     y_a = hI + hIm1;
     ab_a = hI + hIp1;
     bb_a = h + hm1;
     cb_a = h + hp1;
     db_a = h + hm1;
     eb_a = h + hp1;
     return (((((((((mew * (h - hI) / dt - (((sigma * pow(h + hm1, 3.0) * (std::cos
       (thetaI) - 1.0) * (h - hm1) * (m - 1.0) * (n - 1.0) * (h_star * m * pow
       (h_star / (h / 2.0 + hm1 / 2.0), m - 1.0) / (a * a) - h_star * n * pow
       (h_star / (h / 2.0 + hm1 / 2.0), n - 1.0) / (b_a * b_a)) / (h_star * (m - n))
       + sigma * pow(h + hp1, 3.0) * (std::cos(thetaI) - 1.0) * (h - hp1) * (m -
       1.0) * (n - 1.0) * (h_star * m * pow(h_star / (h / 2.0 + hp1 / 2.0), m - 1.0)
                           / (c_a * c_a) - h_star * n * pow(h_star / (h / 2.0 + hp1
       / 2.0), n - 1.0) / (d_a * d_a)) / (h_star * (m - n))) + sigma * pow(hI +
       hIm1, 3.0) * (std::cos(thetaI) - 1.0) * (hI - hIm1) * (m - 1.0) * (n - 1.0) *
       (h_star * m * pow(h_star / (hI / 2.0 + hIm1 / 2.0), m - 1.0) / (e_a * e_a) -
        h_star * n * pow(h_star / (hI / 2.0 + hIm1 / 2.0), n - 1.0) / (f_a * f_a)) /
       (h_star * (m - n))) + sigma * pow(hI + hIp1, 3.0) * (std::cos(thetaI) - 1.0)
       * (hI - hIp1) * (m - 1.0) * (n - 1.0) * (h_star * m * pow(h_star / (hI / 2.0
       + hIp1 / 2.0), m - 1.0) / (g_a * g_a) - h_star * n * pow(h_star / (hI / 2.0
       + hIp1 / 2.0), n - 1.0) / (h_a * h_a)) / (h_star * (m - n))) / (48.0 * (dx *
       dx))) + sigma * (((pow(h + hm1, 3.0) * (((3.0 * h - 3.0 * hm1) + hm2) - hp1)
                          + pow(hI + hIm1, 3.0) * (((3.0 * hI - 3.0 * hIm1) + hIm2)
       - hIp1)) + pow(hI + hIp1, 3.0) * (((3.0 * hI - hIm1) - 3.0 * hIp1) + hIp2))
                        + pow(h + hp1, 3.0) * (((3.0 * h - hm1) - 3.0 * hp1) + hp2))
                    / (48.0 * pow(dx, 4.0))) - gr * p * (((pow(h + hm1, 3.0) * (std::
       sin(theta) - std::cos(theta) * (h - hm1) / dx) - pow(h + hp1, 3.0) * (std::
       sin(theta) + std::cos(theta) * (h - hp1) / dx)) + pow(hI + hIm1, 3.0) * (std::
       sin(theta) - std::cos(theta) * (hI - hIm1) / dx)) - pow(hI + hIp1, 3.0) *
       (std::sin(theta) + std::cos(theta) * (hI - hIp1) / dx)) / (48.0 * dx)) -
                  lam1 * (((sigma * pow(h + hm1, 3.0) * (std::cos(thetaI) - 1.0) *
       (h - hm1) * (m - 1.0) * (n - 1.0) * (h_star * m * pow(h_star / (h / 2.0 +
       hm1 / 2.0), m - 1.0) / (i_a * i_a) - h_star * n * pow(h_star / (h / 2.0 +
       hm1 / 2.0), n - 1.0) / (j_a * j_a)) / (h_star * (m - n)) + sigma * pow(h +
       hp1, 3.0) * (std::cos(thetaI) - 1.0) * (h - hp1) * (m - 1.0) * (n - 1.0) *
       (h_star * m * pow(h_star / (h / 2.0 + hp1 / 2.0), m - 1.0) / (k_a * k_a) -
        h_star * n * pow(h_star / (h / 2.0 + hp1 / 2.0), n - 1.0) / (l_a * l_a)) /
       (h_star * (m - n))) - sigma * pow(hI + hIm1, 3.0) * (std::cos(thetaI) - 1.0)
                           * (hI - hIm1) * (m - 1.0) * (n - 1.0) * (h_star * m *
       pow(h_star / (hI / 2.0 + hIm1 / 2.0), m - 1.0) / (m_a * m_a) - h_star * n *
       pow(h_star / (hI / 2.0 + hIm1 / 2.0), n - 1.0) / (n_a * n_a)) / (h_star * (m
       - n))) - sigma * pow(hI + hIp1, 3.0) * (std::cos(thetaI) - 1.0) * (hI - hIp1)
                          * (m - 1.0) * (n - 1.0) * (h_star * m * pow(h_star / (hI /
       2.0 + hIp1 / 2.0), m - 1.0) / (o_a * o_a) - h_star * n * pow(h_star / (hI /
       2.0 + hIp1 / 2.0), n - 1.0) / (p_a * p_a)) / (h_star * (m - n))) / (24.0 *
       dt * (dx * dx))) + lam1 * (sigma * (q_a * q_a) * (std::cos(thetaI) - 1.0) *
       (h - hp1) * (m - 1.0) * (n - 1.0) * (h_star * m * pow(h_star / (h / 2.0 +
       hp1 / 2.0), m - 1.0) / (r_a * r_a) - h_star * n * pow(h_star / (h / 2.0 +
       hp1 / 2.0), n - 1.0) / (s_a * s_a)) * (((h - hI) - hIp1) + hp1) / (h_star *
       (m - n)) + sigma * (t_a * t_a) * (std::cos(thetaI) - 1.0) * (h - hI) * (h -
       hm1) * (m - 1.0) * (n - 1.0) * (h_star * m * pow(h_star / (h / 2.0 + hm1 /
       2.0), m - 1.0) / (u_a * u_a) - h_star * n * pow(h_star / (h / 2.0 + hm1 /
       2.0), n - 1.0) / (v_a * v_a)) / (h_star * (m - n))) / (16.0 * dt * (dx * dx)))
                - alpha * delsig * (((w_a * w_a * (h - hm1) + x_a * x_a * (h - hp1))
       + y_a * y_a * (hI - hIm1)) + ab_a * ab_a * (hI - hIp1)) / (16.0 * (dx * dx) *
                 kapa)) + lam1 * sigma * (((pow(h + hm1, 3.0) * (((3.0 * h - 3.0 *
       hm1) + hm2) - hp1) - pow(hI + hIm1, 3.0) * (((3.0 * hI - 3.0 * hIm1) + hIm2)
       - hIp1)) - pow(hI + hIp1, 3.0) * (((3.0 * hI - hIm1) - 3.0 * hIp1) + hIp2))
                + pow(h + hp1, 3.0) * (((3.0 * h - hm1) - 3.0 * hp1) + hp2)) /
               (24.0 * dt * pow(dx, 4.0))) - lam1 * sigma * (bb_a * bb_a * (((3.0 *
       h - 3.0 * hm1) + hm2) - hp1) * (((h - hI) - hIm1) + hm1) + cb_a * cb_a *
               (((3.0 * h - hm1) - 3.0 * hp1) + hp2) * (((h - hI) - hIp1) + hp1)) /
              (16.0 * dt * pow(dx, 4.0))) + gr * lam1 * p * (db_a * db_a * (std::
               sin(theta) - std::cos(theta) * (h - hm1) / dx) * (((h - hI) - hIm1)
               + hm1) - eb_a * eb_a * (std::sin(theta) + std::cos(theta) * (h - hp1)
               / dx) * (((h - hI) - hIp1) + hp1)) / (16.0 * dt * dx)) - gr * lam1 *
       p * (((pow(h + hm1, 3.0) * (std::sin(theta) - std::cos(theta) * (h - hm1) /
               dx) - pow(h + hp1, 3.0) * (std::sin(theta) + std::cos(theta) * (h -
                hp1) / dx)) - pow(hI + hIm1, 3.0) * (std::sin(theta) - std::cos
              (theta) * (hI - hIm1) / dx)) + pow(hI + hIp1, 3.0) * (std::sin(theta)
             + std::cos(theta) * (hI - hIp1) / dx)) / (24.0 * dt * dx);
}

double functions::fm2(double hm2, double hm1, double h, double hp1, double hp2, double hIm2, double hIm1, double hI, double hIp1, double hIp2, double mew,double sigma, double n, double m, double h_star, double theta,double thetaI, double lam1, double dx, double dt, double gr, double p, double delsig, double alpha, double kapa)
{
    double a;
    a = h + hm1;
    return (sigma * pow(h + hm1, 3.0) / (48.0 * pow(dx, 4.0)) + lam1 * sigma * pow
            (h + hm1, 3.0) / (24.0 * dt * pow(dx, 4.0))) - lam1 * sigma * (a * a) *
      (((h - hI) - hIm1) + hm1) / (16.0 * dt * pow(dx, 4.0));
}

double functions::fm1(double hm2, double hm1, double h, double hp1, double hp2, double hIm2, double hIm1, double hI, double hIp1, double hIp2, double mew,double sigma, double n, double m, double h_star, double theta,double thetaI, double lam1, double dx, double dt, double gr, double p, double delsig, double alpha, double kapa)
{
	double a;
	  double b_a;
	  double c_a;
	  double d_a;
	  double e_a;
	  double f_a;
	  double g_a;
	  double h_a;
	  double i_a;
	  double j_a;
	  double k_a;
	  double l_a;
	  double m_a;
	  double n_a;
	  double o_a;
	  double p_a;
	  double q_a;
	  double r_a;
	  double s_a;
	  double t_a;
	  double u_a;
	  double v_a;
	  double w_a;
	  double x_a;
	  double y_a;
	  double ab_a;
	  a = h / 2.0 + hm1 / 2.0;
	  b_a = h / 2.0 + hm1 / 2.0;
	  c_a = h + hm1;
	  d_a = h / 2.0 + hm1 / 2.0;
	  e_a = h / 2.0 + hm1 / 2.0;
	  f_a = h + hm1;
	  g_a = h / 2.0 + hm1 / 2.0;
	  h_a = h / 2.0 + hm1 / 2.0;
	  i_a = h + hm1;
	  j_a = h / 2.0 + hm1 / 2.0;
	  k_a = h / 2.0 + hm1 / 2.0;
	  l_a = h + hm1;
	  m_a = h / 2.0 + hm1 / 2.0;
	  n_a = h / 2.0 + hm1 / 2.0;
	  o_a = h / 2.0 + hm1 / 2.0;
	  p_a = h / 2.0 + hm1 / 2.0;
	  q_a = h + hm1;
	  r_a = h + hm1;
	  s_a = h + hm1;
	  t_a = h + hm1;
	  u_a = h + hm1;
	  v_a = h + hm1;
	  w_a = h + hp1;
	  x_a = h + hm1;
	  y_a = h + hm1;
	  ab_a = h + hm1;
	  return ((((((((((sigma * pow(h + hm1, 3.0) * (std::cos(thetaI) - 1.0) * (m -
	    1.0) * (n - 1.0) * (h_star * m * pow(h_star / (h / 2.0 + hm1 / 2.0), m - 1.0)
	                        / (a * a) - h_star * n * pow(h_star / (h / 2.0 + hm1 /
	    2.0), n - 1.0) / (b_a * b_a)) / (h_star * (m - n)) + sigma * pow(h + hm1,
	    3.0) * (std::cos(thetaI) - 1.0) * (h - hm1) * (m - 1.0) * (n - 1.0) *
	                   (((h_star * m * pow(h_star / (h / 2.0 + hm1 / 2.0), m - 1.0) /
	                      pow(h / 2.0 + hm1 / 2.0, 3.0) - h_star * n * pow(h_star /
	    (h / 2.0 + hm1 / 2.0), n - 1.0) / pow(h / 2.0 + hm1 / 2.0, 3.0)) + h_star *
	                     h_star * m * pow(h_star / (h / 2.0 + hm1 / 2.0), m - 2.0) *
	                     (m - 1.0) / (2.0 * pow(h / 2.0 + hm1 / 2.0, 4.0))) - h_star
	                    * h_star * n * pow(h_star / (h / 2.0 + hm1 / 2.0), n - 2.0) *
	                    (n - 1.0) / (2.0 * pow(h / 2.0 + hm1 / 2.0, 4.0))) / (h_star
	    * (m - n))) - 3.0 * sigma * (c_a * c_a) * (std::cos(thetaI) - 1.0) * (h -
	    hm1) * (m - 1.0) * (n - 1.0) * (h_star * m * pow(h_star / (h / 2.0 + hm1 /
	    2.0), m - 1.0) / (d_a * d_a) - h_star * n * pow(h_star / (h / 2.0 + hm1 /
	    2.0), n - 1.0) / (e_a * e_a)) / (h_star * (m - n))) / (48.0 * (dx * dx)) -
	                 sigma * ((3.0 * pow(h + hm1, 3.0) - 3.0 * (f_a * f_a) * (((3.0 *
	    h - 3.0 * hm1) + hm2) - hp1)) + pow(h + hp1, 3.0)) / (48.0 * pow(dx, 4.0)))
	                + lam1 * ((sigma * pow(h + hm1, 3.0) * (std::cos(thetaI) - 1.0) *
	    (m - 1.0) * (n - 1.0) * (h_star * m * pow(h_star / (h / 2.0 + hm1 / 2.0), m
	    - 1.0) / (g_a * g_a) - h_star * n * pow(h_star / (h / 2.0 + hm1 / 2.0), n -
	    1.0) / (h_a * h_a)) / (h_star * (m - n)) + sigma * pow(h + hm1, 3.0) * (std::
	    cos(thetaI) - 1.0) * (h - hm1) * (m - 1.0) * (n - 1.0) * (((h_star * m * pow
	    (h_star / (h / 2.0 + hm1 / 2.0), m - 1.0) / pow(h / 2.0 + hm1 / 2.0, 3.0) -
	    h_star * n * pow(h_star / (h / 2.0 + hm1 / 2.0), n - 1.0) / pow(h / 2.0 +
	    hm1 / 2.0, 3.0)) + h_star * h_star * m * pow(h_star / (h / 2.0 + hm1 / 2.0),
	    m - 2.0) * (m - 1.0) / (2.0 * pow(h / 2.0 + hm1 / 2.0, 4.0))) - h_star *
	    h_star * n * pow(h_star / (h / 2.0 + hm1 / 2.0), n - 2.0) * (n - 1.0) / (2.0
	    * pow(h / 2.0 + hm1 / 2.0, 4.0))) / (h_star * (m - n))) - 3.0 * sigma * (i_a
	    * i_a) * (std::cos(thetaI) - 1.0) * (h - hm1) * (m - 1.0) * (n - 1.0) *
	    (h_star * m * pow(h_star / (h / 2.0 + hm1 / 2.0), m - 1.0) / (j_a * j_a) -
	     h_star * n * pow(h_star / (h / 2.0 + hm1 / 2.0), n - 1.0) / (k_a * k_a)) /
	    (h_star * (m - n))) / (24.0 * dt * (dx * dx))) - lam1 * ((sigma * (l_a * l_a)
	    * (std::cos(thetaI) - 1.0) * (h - hI) * (m - 1.0) * (n - 1.0) * (h_star * m *
	    pow(h_star / (h / 2.0 + hm1 / 2.0), m - 1.0) / (m_a * m_a) - h_star * n *
	    pow(h_star / (h / 2.0 + hm1 / 2.0), n - 1.0) / (n_a * n_a)) / (h_star * (m -
	    n)) - sigma * (std::cos(thetaI) - 1.0) * (h - hI) * (h - hm1) * (m - 1.0) *
	    (n - 1.0) * (h_star * m * pow(h_star / (h / 2.0 + hm1 / 2.0), m - 1.0) /
	                 (o_a * o_a) - h_star * n * pow(h_star / (h / 2.0 + hm1 / 2.0),
	    n - 1.0) / (p_a * p_a)) * (2.0 * h + 2.0 * hm1) / (h_star * (m - n))) +
	    sigma * (q_a * q_a) * (std::cos(thetaI) - 1.0) * (h - hI) * (h - hm1) * (m -
	    1.0) * (n - 1.0) * (((h_star * m * pow(h_star / (h / 2.0 + hm1 / 2.0), m -
	    1.0) / pow(h / 2.0 + hm1 / 2.0, 3.0) - h_star * n * pow(h_star / (h / 2.0 +
	    hm1 / 2.0), n - 1.0) / pow(h / 2.0 + hm1 / 2.0, 3.0)) + h_star * h_star * m *
	    pow(h_star / (h / 2.0 + hm1 / 2.0), m - 2.0) * (m - 1.0) / (2.0 * pow(h /
	    2.0 + hm1 / 2.0, 4.0))) - h_star * h_star * n * pow(h_star / (h / 2.0 + hm1 /
	    2.0), n - 2.0) * (n - 1.0) / (2.0 * pow(h / 2.0 + hm1 / 2.0, 4.0))) /
	    (h_star * (m - n))) / (16.0 * dt * (dx * dx))) - gr * p * (3.0 * (r_a * r_a)
	    * (std::sin(theta) - std::cos(theta) * (h - hm1) / dx) + std::cos(theta) *
	    pow(h + hm1, 3.0) / dx) / (48.0 * dx)) + alpha * delsig * (s_a * s_a - (h -
	    hm1) * (2.0 * h + 2.0 * hm1)) / (16.0 * (dx * dx) * kapa)) - lam1 * sigma *
	            ((3.0 * pow(h + hm1, 3.0) - 3.0 * (t_a * t_a) * (((3.0 * h - 3.0 *
	    hm1) + hm2) - hp1)) + pow(h + hp1, 3.0)) / (24.0 * dt * pow(dx, 4.0))) -
	           lam1 * sigma * (((u_a * u_a * (((3.0 * h - 3.0 * hm1) + hm2) - hp1) -
	              3.0 * (v_a * v_a) * (((h - hI) - hIm1) + hm1)) - w_a * w_a * (((h
	    - hI) - hIp1) + hp1)) + (2.0 * h + 2.0 * hm1) * (((3.0 * h - 3.0 * hm1) +
	              hm2) - hp1) * (((h - hI) - hIm1) + hm1)) / (16.0 * dt * pow(dx,
	             4.0))) + gr * lam1 * p * ((x_a * x_a * (std::sin(theta) - std::cos
	             (theta) * (h - hm1) / dx) + (std::sin(theta) - std::cos(theta) * (h
	              - hm1) / dx) * (2.0 * h + 2.0 * hm1) * (((h - hI) - hIm1) + hm1))
	           + std::cos(theta) * (y_a * y_a) * (((h - hI) - hIm1) + hm1) / dx) /
	          (16.0 * dt * dx)) - gr * lam1 * p * (3.0 * (ab_a * ab_a) * (std::sin
	    (theta) - std::cos(theta) * (h - hm1) / dx) + std::cos(theta) * pow(h + hm1,
	    3.0) / dx) / (24.0 * dt * dx);
}

double functions::f(double hm2, double hm1, double h, double hp1, double hp2, double hIm2, double hIm1, double hI, double hIp1, double hIp2, double mew,double sigma, double n, double m, double h_star, double theta,double thetaI, double lam1, double dx, double dt, double gr, double p, double delsig, double alpha, double kapa)
{
    double a;
     double b_a;
     double c_a;
     double d_a;
     double e_a;
     double f_a;
     double g_a;
     double h_a;
     double i_a;
     double j_a;
     double k_a;
     double l_a;
     double m_a;
     double n_a;
     double o_a;
     double p_a;
     double q_a;
     double r_a;
     double s_a;
     double t_a;
     double u_a;
     double v_a;
     double w_a;
     double x_a;
     double y_a;
     double ab_a;
     double bb_a;
     double cb_a;
     double db_a;
     double eb_a;
     double fb_a;
     double gb_a;
     double hb_a;
     double ib_a;
     double jb_a;
     double kb_a;
     double lb_a;
     double mb_a;
     double nb_a;
     double ob_a;
     double pb_a;
     double qb_a;
     double rb_a;
     double sb_a;
     double tb_a;
     double ub_a;
     double vb_a;
     double wb_a;
     double xb_a;
     double yb_a;
     double ac_a;
     double bc_a;
     double cc_a;
     double dc_a;
     double ec_a;
     double fc_a;
     a = h / 2.0 + hm1 / 2.0;
     b_a = h / 2.0 + hm1 / 2.0;
     c_a = h / 2.0 + hp1 / 2.0;
     d_a = h / 2.0 + hp1 / 2.0;
     e_a = h + hm1;
     f_a = h / 2.0 + hm1 / 2.0;
     g_a = h / 2.0 + hm1 / 2.0;
     h_a = h + hp1;
     i_a = h / 2.0 + hp1 / 2.0;
     j_a = h / 2.0 + hp1 / 2.0;
     k_a = h + hm1;
     l_a = h + hp1;
     m_a = h + hm1;
     n_a = h / 2.0 + hm1 / 2.0;
     o_a = h / 2.0 + hm1 / 2.0;
     p_a = h + hm1;
     q_a = h / 2.0 + hm1 / 2.0;
     r_a = h / 2.0 + hm1 / 2.0;
     s_a = h + hp1;
     t_a = h / 2.0 + hp1 / 2.0;
     u_a = h / 2.0 + hp1 / 2.0;
     v_a = h + hp1;
     w_a = h / 2.0 + hp1 / 2.0;
     x_a = h / 2.0 + hp1 / 2.0;
     y_a = h + hp1;
     ab_a = h / 2.0 + hm1 / 2.0;
     bb_a = h / 2.0 + hm1 / 2.0;
     cb_a = h / 2.0 + hp1 / 2.0;
     db_a = h / 2.0 + hp1 / 2.0;
     eb_a = h + hm1;
     fb_a = h + hp1;
     gb_a = h + hm1;
     hb_a = h / 2.0 + hm1 / 2.0;
     ib_a = h / 2.0 + hm1 / 2.0;
     jb_a = h / 2.0 + hp1 / 2.0;
     kb_a = h / 2.0 + hp1 / 2.0;
     lb_a = h + hm1;
     mb_a = h / 2.0 + hm1 / 2.0;
     nb_a = h / 2.0 + hm1 / 2.0;
     ob_a = h + hp1;
     pb_a = h / 2.0 + hp1 / 2.0;
     qb_a = h / 2.0 + hp1 / 2.0;
     rb_a = h + hm1;
     sb_a = h + hp1;
     tb_a = h + hm1;
     ub_a = h + hp1;
     vb_a = h + hm1;
     wb_a = h + hp1;
     xb_a = h + hm1;
     yb_a = h + hp1;
     ac_a = h + hp1;
     bc_a = h + hm1;
     cc_a = h + hp1;
     dc_a = h + hm1;
     ec_a = h + hm1;
     fc_a = h + hp1;
     return (((((((((mew / dt - (((((sigma * pow(h + hm1,m) * (std::cos(thetaI) -
       1.0) * (m - 1.0) * (n - 1.0) * (h_star * m * pow(h_star / (h / 2.0 + hm1 /
       2.0), m - 1.0) / (a * a) - h_star * n * pow(h_star / (h / 2.0 + hm1 / 2.0),
       n - 1.0) / (b_a * b_a)) / (h_star * (m - n)) + sigma * pow(h + hp1,m) *
       (std::cos(thetaI) - 1.0) * (m - 1.0) * (n - 1.0) * (h_star * m * pow(h_star /
       (h / 2.0 + hp1 / 2.0), m - 1.0) / (c_a * c_a) - h_star * n * pow(h_star / (h
       / 2.0 + hp1 / 2.0), n - 1.0) / (d_a * d_a)) / (h_star * (m - n))) - sigma *
       pow(h + hm1, 3.0) * (std::cos(thetaI) - 1.0) * (h - hm1) * (m - 1.0) * (n -
       1.0) * (((h_star * m * pow(h_star / (h / 2.0 + hm1 / 2.0), m - 1.0) / pow(h /
       2.0 + hm1 / 2.0, 3.0) - h_star * n * pow(h_star / (h / 2.0 + hm1 / 2.0), n -
       1.0) / pow(h / 2.0 + hm1 / 2.0, 3.0)) + h_star * h_star * m * pow(h_star /
       (h / 2.0 + hm1 / 2.0), m - 2.0) * (m - 1.0) / (2.0 * pow(h / 2.0 + hm1 / 2.0,
       4.0))) - h_star * h_star * n * pow(h_star / (h / 2.0 + hm1 / 2.0), n - 2.0) *
               (n - 1.0) / (2.0 * pow(h / 2.0 + hm1 / 2.0, 4.0))) / (h_star * (m -
       n))) - sigma * pow(h + hp1, 3.0) * (std::cos(thetaI) - 1.0) * (h - hp1) * (m
       - 1.0) * (n - 1.0) * (((h_star * m * pow(h_star / (h / 2.0 + hp1 / 2.0), m -
       1.0) / pow(h / 2.0 + hp1 / 2.0, 3.0) - h_star * n * pow(h_star / (h / 2.0 +
       hp1 / 2.0), n - 1.0) / pow(h / 2.0 + hp1 / 2.0, 3.0)) + h_star * h_star * m *
       pow(h_star / (h / 2.0 + hp1 / 2.0), m - 2.0) * (m - 1.0) / (2.0 * pow(h /
       2.0 + hp1 / 2.0, 4.0))) - h_star * h_star * n * pow(h_star / (h / 2.0 + hp1 /
       2.0), n - 2.0) * (n - 1.0) / (2.0 * pow(h / 2.0 + hp1 / 2.0, 4.0))) /
       (h_star * (m - n))) + 3.0 * sigma * (e_a * e_a) * (std::cos(thetaI) - 1.0) *
       (h - hm1) * (m - 1.0) * (n - 1.0) * (h_star * m * pow(h_star / (h / 2.0 +
       hm1 / 2.0), m - 1.0) / (f_a * f_a) - h_star * n * pow(h_star / (h / 2.0 +
       hm1 / 2.0), n - 1.0) / (g_a * g_a)) / (h_star * (m - n))) + 3.0 * sigma *
       (h_a * h_a) * (std::cos(thetaI) - 1.0) * (h - hp1) * (m - 1.0) * (n - 1.0) *
       (h_star * m * pow(h_star / (h / 2.0 + hp1 / 2.0), m - 1.0) / (i_a * i_a) -
        h_star * n * pow(h_star / (h / 2.0 + hp1 / 2.0), n - 1.0) / (j_a * j_a)) /
       (h_star * (m - n))) / (48.0 * (dx * dx))) + sigma * (((3.0 * (k_a * k_a) *
       (((3.0 * h - 3.0 * hm1) + hm2) - hp1) + 3.0 * (l_a * l_a) * (((3.0 * h - hm1)
       - 3.0 * hp1) + hp2)) + 3.0 * pow(h + hm1, 3.0)) + 3.0 * pow(h + hp1, 3.0)) /
                    (48.0 * pow(dx, 4.0))) + lam1 * (((((((sigma * (m_a * m_a) *
       (std::cos(thetaI) - 1.0) * (h - hI) * (m - 1.0) * (n - 1.0) * (h_star * m *
       pow(h_star / (h / 2.0 + hm1 / 2.0), m - 1.0) / (n_a * n_a) - h_star * n *
       pow(h_star / (h / 2.0 + hm1 / 2.0), n - 1.0) / (o_a * o_a)) / (h_star * (m -
       n)) + sigma * (p_a * p_a) * (std::cos(thetaI) - 1.0) * (h - hm1) * (m - 1.0)
       * (n - 1.0) * (h_star * m * pow(h_star / (h / 2.0 + hm1 / 2.0), m - 1.0) /
                      (q_a * q_a) - h_star * n * pow(h_star / (h / 2.0 + hm1 / 2.0),
       n - 1.0) / (r_a * r_a)) / (h_star * (m - n))) + sigma * (s_a * s_a) * (std::
       cos(thetaI) - 1.0) * (h - hp1) * (m - 1.0) * (n - 1.0) * (h_star * m * pow
       (h_star / (h / 2.0 + hp1 / 2.0), m - 1.0) / (t_a * t_a) - h_star * n * pow
       (h_star / (h / 2.0 + hp1 / 2.0), n - 1.0) / (u_a * u_a)) / (h_star * (m - n)))
       + sigma * (v_a * v_a) * (std::cos(thetaI) - 1.0) * (m - 1.0) * (n - 1.0) *
       (h_star * m * pow(h_star / (h / 2.0 + hp1 / 2.0), m - 1.0) / (w_a * w_a) -
        h_star * n * pow(h_star / (h / 2.0 + hp1 / 2.0), n - 1.0) / (x_a * x_a)) *
       (((h - hI) - hIp1) + hp1) / (h_star * (m - n))) - sigma * (y_a * y_a) * (std::
       cos(thetaI) - 1.0) * (h - hp1) * (m - 1.0) * (n - 1.0) * (((h - hI) - hIp1)
       + hp1) * (((h_star * m * pow(h_star / (h / 2.0 + hp1 / 2.0), m - 1.0) / pow
                   (h / 2.0 + hp1 / 2.0, 3.0) - h_star * n * pow(h_star / (h / 2.0
       + hp1 / 2.0), n - 1.0) / pow(h / 2.0 + hp1 / 2.0, 3.0)) + h_star * h_star *
                  m * pow(h_star / (h / 2.0 + hp1 / 2.0), m - 2.0) * (m - 1.0) /
                  (2.0 * pow(h / 2.0 + hp1 / 2.0, 4.0))) - h_star * h_star * n *
                 pow(h_star / (h / 2.0 + hp1 / 2.0), n - 2.0) * (n - 1.0) / (2.0 *
       pow(h / 2.0 + hp1 / 2.0, 4.0))) / (h_star * (m - n))) + sigma * (std::cos
       (thetaI) - 1.0) * (h - hI) * (h - hm1) * (m - 1.0) * (n - 1.0) * (h_star * m
       * pow(h_star / (h / 2.0 + hm1 / 2.0), m - 1.0) / (ab_a * ab_a) - h_star * n *
       pow(h_star / (h / 2.0 + hm1 / 2.0), n - 1.0) / (bb_a * bb_a)) * (2.0 * h +
       2.0 * hm1) / (h_star * (m - n))) + sigma * (std::cos(thetaI) - 1.0) * (h -
       hp1) * (m - 1.0) * (n - 1.0) * (h_star * m * pow(h_star / (h / 2.0 + hp1 /
       2.0), m - 1.0) / (cb_a * cb_a) - h_star * n * pow(h_star / (h / 2.0 + hp1 /
       2.0), n - 1.0) / (db_a * db_a)) * (2.0 * h + 2.0 * hp1) * (((h - hI) - hIp1)
       + hp1) / (h_star * (m - n))) - sigma * (eb_a * eb_a) * (std::cos(thetaI) -
       1.0) * (h - hI) * (h - hm1) * (m - 1.0) * (n - 1.0) * (((h_star * m * pow
       (h_star / (h / 2.0 + hm1 / 2.0), m - 1.0) / pow(h / 2.0 + hm1 / 2.0, 3.0) -
       h_star * n * pow(h_star / (h / 2.0 + hm1 / 2.0), n - 1.0) / pow(h / 2.0 +
       hm1 / 2.0, 3.0)) + h_star * h_star * m * pow(h_star / (h / 2.0 + hm1 / 2.0),
       m - 2.0) * (m - 1.0) / (2.0 * pow(h / 2.0 + hm1 / 2.0, 4.0))) - h_star *
       h_star * n * pow(h_star / (h / 2.0 + hm1 / 2.0), n - 2.0) * (n - 1.0) / (2.0
       * pow(h / 2.0 + hm1 / 2.0, 4.0))) / (h_star * (m - n))) / (16.0 * dt * (dx *
       dx))) + gr * p * (((3.0 * (fb_a * fb_a) * (std::sin(theta) + std::cos(theta)
       * (h - hp1) / dx) - 3.0 * (gb_a * gb_a) * (std::sin(theta) - std::cos(theta)
       * (h - hm1) / dx)) + std::cos(theta) * pow(h + hm1, 3.0) / dx) + std::cos
                         (theta) * pow(h + hp1, 3.0) / dx) / (48.0 * dx)) - lam1 *
                 (((((sigma * pow(h + hm1, 3.0) * (std::cos(thetaI) - 1.0) * (m -
       1.0) * (n - 1.0) * (h_star * m * pow(h_star / (h / 2.0 + hm1 / 2.0), m - 1.0)
                           / (hb_a * hb_a) - h_star * n * pow(h_star / (h / 2.0 +
       hm1 / 2.0), n - 1.0) / (ib_a * ib_a)) / (h_star * (m - n)) + sigma * pow(h +
       hp1, 3.0) * (std::cos(thetaI) - 1.0) * (m - 1.0) * (n - 1.0) * (h_star * m *
       pow(h_star / (h / 2.0 + hp1 / 2.0), m - 1.0) / (jb_a * jb_a) - h_star * n *
       pow(h_star / (h / 2.0 + hp1 / 2.0), n - 1.0) / (kb_a * kb_a)) / (h_star * (m
       - n))) - sigma * pow(h + hm1, 3.0) * (std::cos(thetaI) - 1.0) * (h - hm1) *
                     (m - 1.0) * (n - 1.0) * (((h_star * m * pow(h_star / (h / 2.0
       + hm1 / 2.0), m - 1.0) / pow(h / 2.0 + hm1 / 2.0, 3.0) - h_star * n * pow
       (h_star / (h / 2.0 + hm1 / 2.0), n - 1.0) / pow(h / 2.0 + hm1 / 2.0, 3.0)) +
       h_star * h_star * m * pow(h_star / (h / 2.0 + hm1 / 2.0), m - 2.0) * (m -
       1.0) / (2.0 * pow(h / 2.0 + hm1 / 2.0, 4.0))) - h_star * h_star * n * pow
       (h_star / (h / 2.0 + hm1 / 2.0), n - 2.0) * (n - 1.0) / (2.0 * pow(h / 2.0 +
       hm1 / 2.0, 4.0))) / (h_star * (m - n))) - sigma * pow(h + hp1, 3.0) * (std::
       cos(thetaI) - 1.0) * (h - hp1) * (m - 1.0) * (n - 1.0) * (((h_star * m * pow
       (h_star / (h / 2.0 + hp1 / 2.0), m - 1.0) / pow(h / 2.0 + hp1 / 2.0, 3.0) -
       h_star * n * pow(h_star / (h / 2.0 + hp1 / 2.0), n - 1.0) / pow(h / 2.0 +
       hp1 / 2.0, 3.0)) + h_star * h_star * m * pow(h_star / (h / 2.0 + hp1 / 2.0),
       m - 2.0) * (m - 1.0) / (2.0 * pow(h / 2.0 + hp1 / 2.0, 4.0))) - h_star *
       h_star * n * pow(h_star / (h / 2.0 + hp1 / 2.0), n - 2.0) * (n - 1.0) / (2.0
       * pow(h / 2.0 + hp1 / 2.0, 4.0))) / (h_star * (m - n))) + 3.0 * sigma *
                   (lb_a * lb_a) * (std::cos(thetaI) - 1.0) * (h - hm1) * (m - 1.0)
                   * (n - 1.0) * (h_star * m * pow(h_star / (h / 2.0 + hm1 / 2.0),
       m - 1.0) / (mb_a * mb_a) - h_star * n * pow(h_star / (h / 2.0 + hm1 / 2.0),
       n - 1.0) / (nb_a * nb_a)) / (h_star * (m - n))) + 3.0 * sigma * (ob_a * ob_a)
                  * (std::cos(thetaI) - 1.0) * (h - hp1) * (m - 1.0) * (n - 1.0) *
                  (h_star * m * pow(h_star / (h / 2.0 + hp1 / 2.0), m - 1.0) /
                   (pb_a * pb_a) - h_star * n * pow(h_star / (h / 2.0 + hp1 / 2.0),
       n - 1.0) / (qb_a * qb_a)) / (h_star * (m - n))) / (24.0 * dt * (dx * dx))) +
                lam1 * sigma * (((3.0 * (rb_a * rb_a) * (((3.0 * h - 3.0 * hm1) +
       hm2) - hp1) + 3.0 * (sb_a * sb_a) * (((3.0 * h - hm1) - 3.0 * hp1) + hp2)) +
       3.0 * pow(h + hm1, 3.0)) + 3.0 * pow(h + hp1, 3.0)) / (24.0 * dt * pow(dx,
       4.0))) - alpha * delsig * (((tb_a * tb_a + ub_a * ub_a) + (h - hm1) * (2.0 *
       h + 2.0 * hm1)) + (h - hp1) * (2.0 * h + 2.0 * hp1)) / (16.0 * (dx * dx) *
                kapa)) - lam1 * sigma * (((((vb_a * vb_a * (((3.0 * h - 3.0 * hm1)
       + hm2) - hp1) + wb_a * wb_a * (((3.0 * h - hm1) - 3.0 * hp1) + hp2)) + 3.0 *
       (xb_a * xb_a) * (((h - hI) - hIm1) + hm1)) + 3.0 * (yb_a * yb_a) * (((h - hI)
       - hIp1) + hp1)) + (2.0 * h + 2.0 * hm1) * (((3.0 * h - 3.0 * hm1) + hm2) -
                 hp1) * (((h - hI) - hIm1) + hm1)) + (2.0 * h + 2.0 * hp1) * (((3.0
       * h - hm1) - 3.0 * hp1) + hp2) * (((h - hI) - hIp1) + hp1)) / (16.0 * dt *
               pow(dx, 4.0))) + gr * lam1 * p * (((3.0 * (ac_a * ac_a) * (std::sin
                 (theta) + std::cos(theta) * (h - hp1) / dx) - 3.0 * (bc_a * bc_a) *
                (std::sin(theta) - std::cos(theta) * (h - hm1) / dx)) + std::cos
               (theta) * pow(h + hm1, 3.0) / dx) + std::cos(theta) * pow(h + hp1,
               3.0) / dx) / (24.0 * dt * dx)) - gr * lam1 * p * (((((cc_a * cc_a *
       (std::sin(theta) + std::cos(theta) * (h - hp1) / dx) - dc_a * dc_a * (std::
       sin(theta) - std::cos(theta) * (h - hm1) / dx)) - (std::sin(theta) - std::
       cos(theta) * (h - hm1) / dx) * (2.0 * h + 2.0 * hm1) * (((h - hI) - hIm1) +
       hm1)) + (std::sin(theta) + std::cos(theta) * (h - hp1) / dx) * (2.0 * h +
       2.0 * hp1) * (((h - hI) - hIp1) + hp1)) + std::cos(theta) * (ec_a * ec_a) *
       (((h - hI) - hIm1) + hm1) / dx) + std::cos(theta) * (fc_a * fc_a) * (((h -
       hI) - hIp1) + hp1) / dx) / (16.0 * dt * dx);
}

double functions::fp1(double hm2, double hm1, double h, double hp1, double hp2, double hIm2, double hIm1, double hI, double hIp1, double hIp2, double mew,double sigma, double n, double m, double h_star, double theta,double thetaI, double lam1, double dx, double dt, double gr, double p, double delsig, double alpha, double kapa)
{
	double a;
	  double b_a;
	  double c_a;
	  double d_a;
	  double e_a;
	  double f_a;
	  double g_a;
	  double h_a;
	  double i_a;
	  double j_a;
	  double k_a;
	  double l_a;
	  double m_a;
	  double n_a;
	  double o_a;
	  double p_a;
	  double q_a;
	  double r_a;
	  double s_a;
	  double t_a;
	  double u_a;
	  double v_a;
	  double w_a;
	  double x_a;
	  double y_a;
	  double ab_a;
	  double bb_a;
	  double cb_a;
	  double db_a;
	  a = h / 2.0 + hp1 / 2.0;
	  b_a = h / 2.0 + hp1 / 2.0;
	  c_a = h + hp1;
	  d_a = h / 2.0 + hp1 / 2.0;
	  e_a = h / 2.0 + hp1 / 2.0;
	  f_a = h + hp1;
	  g_a = h / 2.0 + hp1 / 2.0;
	  h_a = h / 2.0 + hp1 / 2.0;
	  i_a = h + hp1;
	  j_a = h / 2.0 + hp1 / 2.0;
	  k_a = h / 2.0 + hp1 / 2.0;
	  l_a = h + hp1;
	  m_a = h / 2.0 + hp1 / 2.0;
	  n_a = h / 2.0 + hp1 / 2.0;
	  o_a = h + hp1;
	  p_a = h / 2.0 + hp1 / 2.0;
	  q_a = h / 2.0 + hp1 / 2.0;
	  r_a = h + hp1;
	  s_a = h / 2.0 + hp1 / 2.0;
	  t_a = h / 2.0 + hp1 / 2.0;
	  u_a = h + hp1;
	  v_a = h + hp1;
	  w_a = h + hp1;
	  x_a = h + hp1;
	  y_a = h + hm1;
	  ab_a = h + hp1;
	  bb_a = h + hp1;
	  cb_a = h + hp1;
	  db_a = h + hp1;
	  return ((((((((((sigma * pow(h + hp1, 3.0) * (std::cos(thetaI) - 1.0) * (m -
	    1.0) * (n - 1.0) * (h_star * m * pow(h_star / (h / 2.0 + hp1 / 2.0), m - 1.0)
	                        / (a * a) - h_star * n * pow(h_star / (h / 2.0 + hp1 /
	    2.0), n - 1.0) / (b_a * b_a)) / (h_star * (m - n)) + sigma * pow(h + hp1,
	    3.0) * (std::cos(thetaI) - 1.0) * (h - hp1) * (m - 1.0) * (n - 1.0) *
	                   (((h_star * m * pow(h_star / (h / 2.0 + hp1 / 2.0), m - 1.0) /
	                      pow(h / 2.0 + hp1 / 2.0, 3.0) - h_star * n * pow(h_star /
	    (h / 2.0 + hp1 / 2.0), n - 1.0) / pow(h / 2.0 + hp1 / 2.0, 3.0)) + h_star *
	                     h_star * m * pow(h_star / (h / 2.0 + hp1 / 2.0), m - 2.0) *
	                     (m - 1.0) / (2.0 * pow(h / 2.0 + hp1 / 2.0, 4.0))) - h_star
	                    * h_star * n * pow(h_star / (h / 2.0 + hp1 / 2.0), n - 2.0) *
	                    (n - 1.0) / (2.0 * pow(h / 2.0 + hp1 / 2.0, 4.0))) / (h_star
	    * (m - n))) - 3.0 * sigma * (c_a * c_a) * (std::cos(thetaI) - 1.0) * (h -
	    hp1) * (m - 1.0) * (n - 1.0) * (h_star * m * pow(h_star / (h / 2.0 + hp1 /
	    2.0), m - 1.0) / (d_a * d_a) - h_star * n * pow(h_star / (h / 2.0 + hp1 /
	    2.0), n - 1.0) / (e_a * e_a)) / (h_star * (m - n))) / (48.0 * (dx * dx)) -
	                 sigma * ((pow(h + hm1, 3.0) - 3.0 * (f_a * f_a) * (((3.0 * h -
	    hm1) - 3.0 * hp1) + hp2)) + 3.0 * pow(h + hp1, 3.0)) / (48.0 * pow(dx, 4.0)))
	                + lam1 * ((sigma * pow(h + hp1, 3.0) * (std::cos(thetaI) - 1.0) *
	    (m - 1.0) * (n - 1.0) * (h_star * m * pow(h_star / (h / 2.0 + hp1 / 2.0), m
	    - 1.0) / (g_a * g_a) - h_star * n * pow(h_star / (h / 2.0 + hp1 / 2.0), n -
	    1.0) / (h_a * h_a)) / (h_star * (m - n)) + sigma * pow(h + hp1, 3.0) * (std::
	    cos(thetaI) - 1.0) * (h - hp1) * (m - 1.0) * (n - 1.0) * (((h_star * m * pow
	    (h_star / (h / 2.0 + hp1 / 2.0), m - 1.0) / pow(h / 2.0 + hp1 / 2.0, 3.0) -
	    h_star * n * pow(h_star / (h / 2.0 + hp1 / 2.0), n - 1.0) / pow(h / 2.0 +
	    hp1 / 2.0, 3.0)) + h_star * h_star * m * pow(h_star / (h / 2.0 + hp1 / 2.0),
	    m - 2.0) * (m - 1.0) / (2.0 * pow(h / 2.0 + hp1 / 2.0, 4.0))) - h_star *
	    h_star * n * pow(h_star / (h / 2.0 + hp1 / 2.0), n - 2.0) * (n - 1.0) / (2.0
	    * pow(h / 2.0 + hp1 / 2.0, 4.0))) / (h_star * (m - n))) - 3.0 * sigma * (i_a
	    * i_a) * (std::cos(thetaI) - 1.0) * (h - hp1) * (m - 1.0) * (n - 1.0) *
	    (h_star * m * pow(h_star / (h / 2.0 + hp1 / 2.0), m - 1.0) / (j_a * j_a) -
	     h_star * n * pow(h_star / (h / 2.0 + hp1 / 2.0), n - 1.0) / (k_a * k_a)) /
	    (h_star * (m - n))) / (24.0 * dt * (dx * dx))) + lam1 * (((sigma * (l_a *
	    l_a) * (std::cos(thetaI) - 1.0) * (h - hp1) * (m - 1.0) * (n - 1.0) *
	    (h_star * m * pow(h_star / (h / 2.0 + hp1 / 2.0), m - 1.0) / (m_a * m_a) -
	     h_star * n * pow(h_star / (h / 2.0 + hp1 / 2.0), n - 1.0) / (n_a * n_a)) /
	    (h_star * (m - n)) - sigma * (o_a * o_a) * (std::cos(thetaI) - 1.0) * (m -
	    1.0) * (n - 1.0) * (h_star * m * pow(h_star / (h / 2.0 + hp1 / 2.0), m - 1.0)
	                        / (p_a * p_a) - h_star * n * pow(h_star / (h / 2.0 + hp1
	    / 2.0), n - 1.0) / (q_a * q_a)) * (((h - hI) - hIp1) + hp1) / (h_star * (m -
	    n))) - sigma * (r_a * r_a) * (std::cos(thetaI) - 1.0) * (h - hp1) * (m - 1.0)
	    * (n - 1.0) * (((h - hI) - hIp1) + hp1) * (((h_star * m * pow(h_star / (h /
	    2.0 + hp1 / 2.0), m - 1.0) / pow(h / 2.0 + hp1 / 2.0, 3.0) - h_star * n *
	    pow(h_star / (h / 2.0 + hp1 / 2.0), n - 1.0) / pow(h / 2.0 + hp1 / 2.0, 3.0))
	    + h_star * h_star * m * pow(h_star / (h / 2.0 + hp1 / 2.0), m - 2.0) * (m -
	    1.0) / (2.0 * pow(h / 2.0 + hp1 / 2.0, 4.0))) - h_star * h_star * n * pow
	    (h_star / (h / 2.0 + hp1 / 2.0), n - 2.0) * (n - 1.0) / (2.0 * pow(h / 2.0 +
	    hp1 / 2.0, 4.0))) / (h_star * (m - n))) + sigma * (std::cos(thetaI) - 1.0) *
	    (h - hp1) * (m - 1.0) * (n - 1.0) * (h_star * m * pow(h_star / (h / 2.0 +
	    hp1 / 2.0), m - 1.0) / (s_a * s_a) - h_star * n * pow(h_star / (h / 2.0 +
	    hp1 / 2.0), n - 1.0) / (t_a * t_a)) * (2.0 * h + 2.0 * hp1) * (((h - hI) -
	    hIp1) + hp1) / (h_star * (m - n))) / (16.0 * dt * (dx * dx))) + gr * p *
	              (3.0 * (u_a * u_a) * (std::sin(theta) + std::cos(theta) * (h - hp1)
	    / dx) - std::cos(theta) * pow(h + hp1, 3.0) / dx) / (48.0 * dx)) + alpha *
	             delsig * (v_a * v_a - (h - hp1) * (2.0 * h + 2.0 * hp1)) / (16.0 *
	              (dx * dx) * kapa)) - lam1 * sigma * ((pow(h + hm1, 3.0) - 3.0 *
	              (w_a * w_a) * (((3.0 * h - hm1) - 3.0 * hp1) + hp2)) + 3.0 * pow(h
	              + hp1, 3.0)) / (24.0 * dt * pow(dx, 4.0))) - lam1 * sigma * (((x_a
	              * x_a * (((3.0 * h - hm1) - 3.0 * hp1) + hp2) - y_a * y_a * (((h -
	    hI) - hIm1) + hm1)) - 3.0 * (ab_a * ab_a) * (((h - hI) - hIp1) + hp1)) +
	            (2.0 * h + 2.0 * hp1) * (((3.0 * h - hm1) - 3.0 * hp1) + hp2) * (((h
	    - hI) - hIp1) + hp1)) / (16.0 * dt * pow(dx, 4.0))) - gr * lam1 * p * ((bb_a
	            * bb_a * (std::sin(theta) + std::cos(theta) * (h - hp1) / dx) + (std::
	             sin(theta) + std::cos(theta) * (h - hp1) / dx) * (2.0 * h + 2.0 *
	             hp1) * (((h - hI) - hIp1) + hp1)) - std::cos(theta) * (cb_a * cb_a)
	           * (((h - hI) - hIp1) + hp1) / dx) / (16.0 * dt * dx)) + gr * lam1 * p
	    * (3.0 * (db_a * db_a) * (std::sin(theta) + std::cos(theta) * (h - hp1) / dx)
	       - std::cos(theta) * pow(h + hp1, 3.0) / dx) / (24.0 * dt * dx);
}

double functions::fp2(double hm2, double hm1, double h, double hp1, double hp2, double hIm2, double hIm1, double hI, double hIp1, double hIp2, double mew,double sigma, double n, double m, double h_star, double theta,double thetaI, double lam1, double dx, double dt, double gr, double p, double delsig, double alpha, double kapa)
{
	double a;
	 a = h + hp1;
	 return (sigma * pow(h + hp1, 3.0) / (48.0 * pow(dx, 4.0)) + lam1 * sigma * pow
	          (h + hp1, 3.0) / (24.0 * dt * pow(dx, 4.0))) - lam1 * sigma * (a * a) * (((h - hI) - hIp1) + hp1) / (16.0 * dt * pow(dx, 4.0));
}

double functions::eval_functions(double hm2, double hm1, double h, double hp1, double hp2, double hIm2, double hIm1, double hI, double hIp1, double hIp2,coe &data)
{
	return this->main_function(hm2,hm1,h,hp1,hp2,hIm2,hIm1,hI,hIp1,hIp2,data.get_mu(),data.get_sigma(),data.get_n(),data.get_m(), data.get_h_star(), data.get_theta(),data.get_contact_angle(),data.get_lambda_1(),data.get_dx(),data.get_dt(),data.get_gravity(),data.get_rho(),data.get_del_sig(),data.get_alpha(),data.get_kappa());
}

double functions::eval_fm2(double hm2, double hm1, double h, double hp1, double hp2, double hIm2, double hIm1, double hI, double hIp1, double hIp2,coe &data)
{
	return this->fm2(hm2,hm1,h,hp1,hp2,hIm2,hIm1,hI,hIp1,hIp2,data.get_mu(),data.get_sigma(),data.get_n(),data.get_m(), data.get_h_star(), data.get_theta(),data.get_contact_angle(),data.get_lambda_1(),data.get_dx(),data.get_dt(),data.get_gravity(),data.get_rho(),data.get_del_sig(),data.get_alpha(),data.get_kappa());
}

double functions::eval_fm1(double hm2, double hm1, double h, double hp1, double hp2, double hIm2, double hIm1, double hI, double hIp1, double hIp2,coe &data)
{
	return this->fm1(hm2,hm1,h,hp1,hp2,hIm2,hIm1,hI,hIp1,hIp2,data.get_mu(),data.get_sigma(),data.get_n(),data.get_m(), data.get_h_star(), data.get_theta(),data.get_contact_angle(),data.get_lambda_1(),data.get_dx(),data.get_dt(),data.get_gravity(),data.get_rho(),data.get_del_sig(),data.get_alpha(),data.get_kappa());
}

double functions::eval_f(double hm2, double hm1, double h, double hp1, double hp2, double hIm2, double hIm1, double hI, double hIp1, double hIp2,coe &data)
{
	return this->f(hm2,hm1,h,hp1,hp2,hIm2,hIm1,hI,hIp1,hIp2,data.get_mu(),data.get_sigma(),data.get_n(),data.get_m(), data.get_h_star(), data.get_theta(),data.get_contact_angle(),data.get_lambda_1(),data.get_dx(),data.get_dt(),data.get_gravity(),data.get_rho(),data.get_del_sig(),data.get_alpha(),data.get_kappa());
}

double functions::eval_fp1(double hm2, double hm1, double h, double hp1, double hp2, double hIm2, double hIm1, double hI, double hIp1, double hIp2,coe &data)
{
	return this->fp1(hm2,hm1,h,hp1,hp2,hIm2,hIm1,hI,hIp1,hIp2,data.get_mu(),data.get_sigma(),data.get_n(),data.get_m(), data.get_h_star(), data.get_theta(),data.get_contact_angle(),data.get_lambda_1(),data.get_dx(),data.get_dt(),data.get_gravity(),data.get_rho(),data.get_del_sig(),data.get_alpha(),data.get_kappa());
}

double functions::eval_fp2(double hm2, double hm1, double h, double hp1, double hp2, double hIm2, double hIm1, double hI, double hIp1, double hIp2,coe &data)
{
	return this->fp2(hm2,hm1,h,hp1,hp2,hIm2,hIm1,hI,hIp1,hIp2,data.get_mu(),data.get_sigma(),data.get_n(),data.get_m(), data.get_h_star(), data.get_theta(),data.get_contact_angle(),data.get_lambda_1(),data.get_dx(),data.get_dt(),data.get_gravity(),data.get_rho(),data.get_del_sig(),data.get_alpha(),data.get_kappa());
}


