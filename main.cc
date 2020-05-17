#include "main.h"

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
	double lambda = Z/m/N;
    double llh(0);

	if (DEBUG)
		cout << "~~" << m << ' ' << p << endl;

    for ( size_t i(0); i != y.size(); i++ ) { // elements in y are > 0
        llh += log( prob_y(lambda, m, p, y[i], precision) );
    }

	llh += ( N - y.size() ) * prob_y(lambda, m, p, 0.0, precision);
	return llh;
}

double llh_2lambda (
		const double &N,
		const double &precision,
		const vector<double> &y,
		const double &lambda,
		const double &lam2
	)
{
    double llh(0);

    for ( size_t i(0); i != y.size(); i++ ) { // elements in y are > 0
        llh += log( prob_y_2lambda(lambda, lam2, y[i], precision) );
    }

	llh += ( N - y.size() ) * prob_y_2lambda(lambda, lam2, 0.0, precision);
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
	double k_max(300.00);
	double mp = m * p;

	if ( y == 0.0 ) pr += exp(-lambda);

	for ( double k(1.0); k < k_max; k++ ) {
		double tmp = -k*log(1+mp)/p + k*log(lambda) - lambda - lgamma(k+1);

		if ( y != 0.0 )
			tmp += lgamma(y+k/p) - lgamma(k/p) -lgamma(y+1) - y*log(1+1/mp);

		pr += exp(tmp);
		if ( exp(tmp) < precision ) break;
	}

	return pr;
}

double prob_y_2lambda(
		const double &lambda,
		const double &lam2,
		const double &y,
		const double &precision
		)
{
	double pr(0.0);
	double k_max(300.00);

	if ( y == 0.0 ) pr += exp(-lambda);
	
	for ( double k(1.0); k < k_max; k++ ) {
		double tmp = -k*lam2 + k*log(lambda) - lambda - lgamma(k+1);

		if ( y != 0.0 )
			tmp += y*log(k*lam2) -lgamma(y+1);

		pr += exp(tmp);
		if ( exp(tmp) < precision ) break;
	}

	return pr;
}

double prob_k_given_y (
		const double &lambda,
		const double &m,
		const double &p,
		const double &y,
		const double &precision,
		const double &k
		)
{
		double mp = m * p;

		if ( k == 0 ) {
			return ( y == 0 ? exp(-lambda)/prob_y(lambda, m, p, y, precision) : 0.0 );
		}
		else if ( k > 0 ) {
			double tmp = -k*log(1+mp)/p + k*log(lambda) - lambda - lgamma(k+1);

			if ( y != 0.0 )
				tmp += lgamma(y+k/p) - lgamma(k/p) -lgamma(y+1) - y*log(1+1/mp);

			return exp(tmp) / prob_y(lambda, m, p, y, precision);
		}
		else {
			cerr << "k can not < 0" << endl;
			exit(1);
		}
}

// given readN y in a bin, return prob that has k fragments
double prob_k_given_y_2lambda (
		const double &lambda,
		const double &lam2,
		const double &y,
		const double &precision,
		const double &k
		)
{
	if ( k == 0 ) {
		return ( y == 0 ? exp(-lambda)/prob_y_2lambda(lambda, lam2, y, precision) : 0.0 );
	}
	else if ( k > 0 ) {
		double tmp = -k*lam2 + k*log(lambda) - lambda - lgamma(k+1);

		if ( y != 0.0 )
			tmp += y*log(k*lam2) -lgamma(y+1);

		return exp(tmp) / prob_y_2lambda(lambda, lam2, y, precision);
	}
	else {
		cerr << "k can not < 0" << endl;
		exit(1);
	}
}

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

void LogLikelihoodFunc_2lambda (
		const real_1d_array &x,
		double &func,
		void   *opt_data
		)
{
	fn_data* obj = reinterpret_cast<fn_data*>(opt_data);
	// minimize func == maximize llh
	func = -llh_2lambda(obj->N, obj->precision, obj->y, x[0], x[1]);
}

// diff = P(Y=0) = Nempty/Nb, given lambda2 = Z/Nb/lambda
// return diff^2, to find lambda to minimize func
//
void diff_prob_y0_2lambda (
		const real_1d_array &x,
		double &func,
		void   *opt_data
		)
{
	fn_data* obj = reinterpret_cast<fn_data*>(opt_data);
	double diff = prob_y_2lambda(x[0], obj->Z/obj->N/x[0], 0.0, obj->precision)
		-( obj->N - obj->y.size() )/obj->N;

	func = pow(diff, 2);
}
