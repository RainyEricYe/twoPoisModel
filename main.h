/*
 *
 */

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
//#include "ap.h"

using namespace std;
using namespace alglib;


void usage() {
    cout << "Program: " << PROGRAM << "  " << REMARKS << "\n"
            "Version: " << VERSION << "\n"
            "Authors: " << AUTHORS << "\n"
            "Contact: " << CONTACT << "\n\n"

            "Options: " << PROGRAM << "\n"
            "\n"
            "\n";

    exit(1);
}

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

double prob_y(
		const double &lambda,
		const double &m,
		const double &p,
		const double &y,
		const double &precision
		);

void LogLikelihoodFunc (
		const real_1d_array &x,
		double &func,
		void   *opt_data
		)
{
	fn_data* obj = reinterpret_cast<fn_data*>(opt_data);

	// minimize func == maximize llh
	func = -llh(obj->N, obj->Z, obj->precision, obj->y, x[0], x[1]);
}

// log_likelihood() of mu and phi given Nb, Z, and vector of y>0
double llh(
		const double &N,
		const double &Z,
		const double &precision,
		const vector<double> &y,
		const double &m,
		const double &p
		)
{
    double mp = m * p;
	double lambda = Z/m/N;
    double llh(0);
    double k_max = 100; // 10 * lambda

	if (DEBUG)
		cout << "~~" << m << ' ' << p << endl;

    for ( size_t i(0); i != y.size(); i++ ) { // elements in y are > 0

        double sum(0.0);
        for ( double k(1.0); k < k_max; k++ ) {
			double tmp = -k*log(1+mp)/p - y[i]*log(1+1/mp) + k*log(lambda) - lambda
				+ lgamma(y[i]+k/p) - lgamma(k/p) - lgamma(k+1) -lgamma(y[i]+1);

			if (DEBUG)
				cout << "1~" << i << ' ' << y[i] << ' ' << k << ' ' << exp(tmp) << ' ' << sum << endl;

            sum += exp(tmp);
            if ( exp(tmp) < precision ) break;
        }

        llh += log(sum);
    }

    double y0sum = exp(-lambda);
    for ( double k(1.0); k < k_max; k++ ) {
        double tmp = -k*log(1+mp)/p + k*log(lambda) - lambda - lgamma(k+1);

		if (DEBUG)
			cout << "2~" << k << ' ' << exp(tmp) << ' ' << y0sum << endl;

        y0sum += exp(tmp);
        if ( exp(tmp) < precision ) break;
    }

    llh += y0sum * ( N - y.size() ); // sum * Num_of_empty_bins

	return llh;
}

double prob_y(
		const double &lambda,
		const double &m,
		const double &p,
		const double &y,
		const double &precision
		)
{
	double pr(0.0);

	if ( y == 0.0 ) pr += exp(-lambda);

	double k_max(200.00);
	double mp = m * p;

	for ( double k(1.0); k < k_max; k++ ) {
		double tmp = -k*log(1+mp)/p - y*log(1+1/mp) + k*log(lambda) - lambda
			+ lgamma(y+k/p) - lgamma(k/p) - lgamma(k+1) -lgamma(y+1);

		pr += exp(tmp);
		if ( exp(tmp) < precision ) break;
	}

	return pr;
}

#endif // MAIN_H_
