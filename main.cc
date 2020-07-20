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

// pois + pois_dis  2var
// log_likelihood() of mu and phi given Nb, Z, and vector of y>0
//
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
	return llh_3var(N, precision, y, lambda, m, p);
}

// pois_dis + pois  2var
// log_likelihood() of mu and phi given Nb, Z, and vector of y>0 ; prob_y2
//
double llh2(
		const double &N,
		const double &Z,
		const double &precision,
		const vector<double> &y,
		const double &m,
		const double &p
		)
{
	double lambda = Z/m/N;
	return llh2_3var(N, precision, y, lambda, m, p);
}

// pois + pois_dis  3var
// log_likelihood() of lambda mu and phi
//
double llh_3var(
                const double &N,
                const double &precision,
                const vector<double> &y,
                const double &lambda,
                const double &m,
                const double &p
                )
{
    double llh(0);

	if (DEBUG)
		cout << "~~" << m << ' ' << p << endl;

    for ( size_t i(0); i != y.size(); i++ ) { // elements in y are > 0
        llh += log( prob_y(lambda, m, p, y[i], precision) );
    }

	llh += ( N - y.size() ) * prob_y(lambda, m, p, 0.0, precision);
	return llh;
}

// pois_dis + pois  3var
// log_likelihood() of lambda mu and phi
//
double llh2_3var(
                const double &N,
                const double &precision,
                const vector<double> &y,
                const double &lambda,
                const double &m,
                const double &p
                )
{
    double llh(0);

	if (DEBUG)
		cout << "~~" << m << ' ' << p << endl;

    for ( size_t i(0); i != y.size(); i++ ) { // elements in y are > 0
        llh += log( prob_y2(lambda, m, p, y[i], precision) );
    }

	llh += ( N - y.size() ) * prob_y2(lambda, m, p, 0.0, precision);
	return llh;
}

// pois + pois
// log_likelihood() of lambda mu and phi
//
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

// pois + pois_dis
// prob of observe y reads in the family, given 3var: lambda mu phi
//
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

// pois_dis + pois
// prob of observe y reads in the family, given 3var: lambda mu phi
//
double prob_y2(
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

	if ( y == 0.0 ) pr += pow( 1/(1+mp), 1/p );

	for ( double k(1.0); k < k_max; k++ ) {
		double tmp = -k*lambda + y*log(k*lambda) - lgamma(y+1)
			+ lgamma(k+1/p) - lgamma(1/p) - lgamma(k+1) - log(1+mp)/p - k*log(1+1/mp);

		pr += exp(tmp);
		if ( exp(tmp) < precision ) break;
	}

	return pr;
}

// pois + pois
// prob of observe y reads in the family, given 2var: lambda lam2
//
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

// pois + pois_dis
// given readN y in a bin, return prob that has k fragments
//
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

			return exp(tmp)/prob_y(lambda, m, p, y, precision);
		}
		else {
			cerr << "k can not < 0" << endl;
			exit(1);
		}
}

// pois_dis + pois
// given readN y in a bin, return prob that has k fragments
//
double prob_k_given_y2 (
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
			return ( y == 0 ? pow(1/(1+mp), 1/p) /prob_y2(lambda, m, p, y, precision) : 0.0 );
		}
		else if ( k > 0 ) {
			double tmp = -k*lambda + y*log(k*lambda) - lgamma(y+1)
				+ lgamma(k+1/p) - lgamma(1/p) - lgamma(k+1) - log(1+mp)/p - k*log(1+1/mp);

			return exp(tmp)/prob_y2(lambda, m, p, y, precision);
		}
		else {
			cerr << "k can not < 0" << endl;
			exit(1);
		}
}

// pois + pois
// given readN y in a bin, return prob that has k fragments
//
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

		return exp(tmp)/prob_y_2lambda(lambda, lam2, y, precision);
	}
	else {
		cerr << "k can not < 0" << endl;
		exit(1);
	}
}

// pois + pois_dis  2var
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

// pois_dis + pois  2var
void LogLikelihoodFunc2 (
		const real_1d_array &x,
		double &func,
		void   *opt_data
		)
{
	fn_data* obj = reinterpret_cast<fn_data*>(opt_data);
	// minimize func == maximize llh
	func = -llh2(obj->N, obj->Z, obj->precision, obj->y, x[0], x[1]);
}

// pois + pois_dis  3var
void LogLikelihoodFunc_3var (
                const real_1d_array &x,
                double &func,
                void   *opt_data
                )
{
	fn_data* obj = reinterpret_cast<fn_data*>(opt_data);
	func = -llh_3var(obj->N, obj->precision, obj->y, x[0], x[1], x[2]);
}

// pois_dis + pois  3var
void LogLikelihoodFunc2_3var (
                const real_1d_array &x,
                double &func,
                void   *opt_data
                )
{
	fn_data* obj = reinterpret_cast<fn_data*>(opt_data);
	func = -llh2_3var(obj->N, obj->precision, obj->y, x[0], x[1], x[2]);
}

// pois + pois  2var
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

// pois + pois  1var
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

// probability that a mutation with frequency 'f', can not be detected in a read family with size 's'
//
double notDetect(
		const double &f,
		const double &s
		)
{
	if ( f == 1.0 ) return 0.0;

	double lnF = log(f);
	double lnC = log(1-f);

	return (
			exp( s*lnC )
			+ exp( log(s) + lnF + (s-1)*lnC )
			+ exp( log(0.5) + log(s) + log(s-1) + 2*lnF + (s-2)*lnC )
			+ exp( log(0.005) + log(s) + log(s-1) + log(s-2) + 3*lnF + (s-3)*lnC )
		   );
}


