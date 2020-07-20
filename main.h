#ifndef MAIN_H_
#define MAIN_H_

#define PROGRAM "twoPoisModel"
#define VERSION "v1.0"
#define AUTHORS "yerui"
#define CONTACT "yerui@genomics.cn"
#define REMARKS "(two-Poisson Model for LFMD)"
#define DEBUG false

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <getopt.h>
#include <map>

#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "optimization.h"
#include "ap.h"

using namespace std;
using namespace alglib;

void usage();

struct fn_data {
	double N;
	double Z;
	double precision;
	vector<double> y;
};

double llh(
		const double &N,
		const double &Z,
		const double &precision,
		const vector<double> &y,
		const double &m,
		const double &p
		);

double llh2(
		const double &N,
		const double &Z,
		const double &precision,
		const vector<double> &y,
		const double &m,
		const double &p
		);

double llh_3var(
		const double &N,
		const double &precision,
		const vector<double> &y,
		const double &lambda,
		const double &m,
		const double &p
		);

double llh2_3var(
		const double &N,
		const double &precision,
		const vector<double> &y,
		const double &lambda,
		const double &m,
		const double &p
		);

double llh_2lambda (
		const double &N,
		const double &precision,
		const vector<double> &y,
		const double &lambda,
		const double &lam2
	);

double prob_y(
		const double &lambda,
		const double &m,
		const double &p,
		const double &y,
		const double &precision
		);

double prob_y2(
		const double &lambda,
		const double &m,
		const double &p,
		const double &y,
		const double &precision
		);

double prob_y_2lambda(
		const double &lambda,
		const double &lam2,
		const double &y,
		const double &precision
		);


// given readN y in a bin, return prob that has k fragments
double prob_k_given_y (
		const double &lambda,
		const double &m,
		const double &p,
		const double &y,
		const double &precision,
		const double &k
		);

double prob_k_given_y2 (
		const double &lambda,
		const double &m,
		const double &p,
		const double &y,
		const double &precision,
		const double &k
		);

// given readN y in a bin, return prob that has k fragments
double prob_k_given_y_2lambda (
		const double &lambda,
		const double &lam2,
		const double &y,
		const double &precision,
		const double &k
		);

double p_k_given_y_2lambda (
		const double &lambda,
		const double &lam2,
		const double &y,
		const double &precision,
		const double &k
		);


void LogLikelihoodFunc (
		const real_1d_array &x,
		double &func,
		void   *opt_data
		);

void LogLikelihoodFunc2 (
		const real_1d_array &x,
		double &func,
		void   *opt_data
		);

void LogLikelihoodFunc_3var (
		const real_1d_array &x,
		double &func,
		void   *opt_data
		);

void LogLikelihoodFunc2_3var (
		const real_1d_array &x,
		double &func,
		void   *opt_data
		);

void LogLikelihoodFunc_2lambda (
		const real_1d_array &x,
		double &func,
		void   *opt_data
		);

void diff_prob_y0_2lambda (
		const real_1d_array &x,
		double &func,
		void   *opt_data
		);

// probability that a mutation with frequency 'f', can not be detected in a read family with size 's'
//
double notDetect(
		const double &f,
		const double &s
		);

#endif // MAIN_H_
