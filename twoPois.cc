
#include "main.h"
#include "optimization.h"

using namespace std;
using namespace alglib;

int main( int argc, char **argv )
{
	if ( argc < 2 ) usage(), exit(1);

	ifstream inf(argv[1]);

	string line;
	while ( getline(inf, line) ) {
		istringstream stm(line);
		double isert, t;
		stm >> isert;

		fn_data data;
		data.precision = 1.0e-12;
		data.N = 146.0;
		data.Z = 0.0;

		double epsg = 1.0e-6;
		double epsf = 0;
		double epsx = 0;
		ae_int_t maxits = 0;
		double diffstep = 1.0e-6;

		while (stm >> t) {
			data.y.push_back(t);
		}

		double max_y(0.0);

		for ( size_t i(0); i != data.y.size(); i++) {
			data.Z += data.y[i];
			if (data.y[i] > max_y) max_y = data.y[i];
		}


		map<double, vector<double> > mLP;
		cerr << endl;
/*
		// optimizaion log-likelihood function: 3var: lambda, mu and phi
		//
		cerr << isert << ' ' << data.y.size() << ' ' << data.Z << ' ';
		try {
			real_1d_array x = "[0.1, 6.0, 0.1]";
			real_1d_array bndl = "[0.0, 0.0, 0.0]";
			real_1d_array bndu = "[50.0, 50.0, 50.0]";
			minbleicstate state;
			minbleicreport rep;

			minbleiccreatef(x, diffstep, state);
			minbleicsetbc(state, bndl, bndu);
			minbleicsetcond(state, epsg, epsf, epsx, maxits);
			alglib::minbleicoptimize(state, LogLikelihoodFunc_3var, NULL, &data);
			minbleicresults(state, x, rep);

			cerr << "0~ " << x[0] << ' ' << x[1] << ' ' << x[2] << ' ' << rep.terminationtype << ' ';

			if ( rep.terminationtype != -8 ) {
				double llh = llh_3var( data.N, data.precision, data.y, x[0], x[1], x[2]);
				cerr << llh;

				for (size_t i(0); i != 3; i++ )
					mLP[llh].push_back(x[i]);

				mLP[llh].push_back(-1.0); // as a seperateor
			}
			cerr << endl;
		}
		catch (alglib::ap_error &e) {
			cerr << "catch error: " << e.msg << " at insertSize=" << isert << endl;
		}


		// optimizaion log-likelihood function: 2var: mu and phi
		//
		cerr << isert << ' ' << data.y.size() << ' ' << data.Z << ' ';
		try {
			real_1d_array x = "[6.0, 0.1]";
			real_1d_array bndl = "[0.0, 0.0]";
			real_1d_array bndu = "[50.0, 50.0]";
			minbleicstate state;
			minbleicreport rep;

			minbleiccreatef(x, diffstep, state);
			minbleicsetbc(state, bndl, bndu);
			minbleicsetcond(state, epsg, epsf, epsx, maxits);
			alglib::minbleicoptimize(state, LogLikelihoodFunc, NULL, &data);
			minbleicresults(state, x, rep);

			double lambda = data.Z/data.N/x[0];
			cerr << "1~ " << lambda << ' ' << x[0] << ' ' << x[1] << ' ' << rep.terminationtype << ' ';

			if ( rep.terminationtype != -8 ) {
				double llh = llh_3var( data.N, data.precision, data.y, lambda, x[0], x[1]);
				cerr << llh;

				mLP[llh].push_back(lambda);
				for (size_t i(0); i != 2; i++ )
					mLP[llh].push_back(x[i]);

				mLP[llh].push_back(-1.0); // as a seperateor
			}
			cerr << endl;
		}
		catch (alglib::ap_error &e) {
			cerr << "catch error: " << e.msg << " at insertSize=" << isert << endl;
		}

		// maximize the 2lambda model; 2var: lam and lam2
		//
		cerr << isert << ' ' << data.y.size() << ' ' << data.Z << ' ';
		cout << "~ " << isert << ' ' << data.y.size() << ' ' << data.Z << ' ';
		try {
			real_1d_array x = "[0.1, 6.0]";
			real_1d_array bndl = "[0.0, 0.0]";
			real_1d_array bndu = "[50.0, 50.0]";
			minbleicstate state;
			minbleicreport rep;

			minbleiccreatef(x, diffstep, state);
			minbleicsetbc(state, bndl, bndu);
			minbleicsetcond(state, epsg, epsf, epsx, maxits);
			alglib::minbleicoptimize(state, LogLikelihoodFunc_2lambda, NULL, &data);
			minbleicresults(state, x, rep);

			cerr << "2~ " << x[0] << ' ' << x[1] << ' ' << rep.terminationtype << ' ';

			if ( rep.terminationtype != -8 ) {
				double llh = llh_2lambda(data.N, data.precision, data.y, x[0], x[1]);
				cerr << llh;
				cout << x[0] << ' ' << x[1] << ' ' << llh;

				for (size_t i(0); i != 2; i++ )
					mLP[llh].push_back(x[i]);

				mLP[llh].push_back(-1.0); // as a seperateor
			}
			cerr << endl;
			cout << endl;
		}
		catch (alglib::ap_error &e) {
			cerr << "catch error: " << e.msg << " at insertSize=" << isert << endl;
		}

		// solve diff_prob_y0_2lamda to get lambda
		cerr << isert << ' ' << data.y.size() << ' ' << data.Z << ' ';
		cout << "~ " << isert << ' ' << data.y.size() << ' ' << data.Z << ' ';
		try {
			real_1d_array x = "[1.0]";
			real_1d_array bndl = "[0.0]";
			real_1d_array bndu = "[50.0]";
			minbleicstate state;
			minbleicreport rep;

			minbleiccreatef(x, diffstep, state);
			minbleicsetbc(state, bndl, bndu);
			minbleicsetcond(state, epsg, epsf, epsx, maxits);
			alglib::minbleicoptimize(state, diff_prob_y0_2lambda, NULL, &data);
			minbleicresults(state, x, rep);

			double expect_lam2 = data.Z/data.N/x[0];
			cerr << "3~ " << x[0] << ' ' << expect_lam2 << ' ' << rep.terminationtype << ' ';

			if ( rep.terminationtype != -8 ) {
				double llh = llh_2lambda(data.N, data.precision, data.y, x[0], expect_lam2);
				cerr << llh;
				cout << x[0] << ' ' << expect_lam2 << ' ' << llh;

				mLP[llh].push_back(x[0]);
				mLP[llh].push_back(expect_lam2);

				mLP[llh].push_back(-1.0); // as a seperateor
			}
			cerr << endl;
			cout << endl;

		}
		catch (alglib::ap_error &e) {
			cerr << "catch error: " << e.msg << " at insertSize=" << isert << endl;
		}
		*/
		//
		//
		// optimizaion log-likelihood function: 3var: lambda, mu and phi   Pois_dis + Pois model
		//
		cerr << isert << ' ' << data.y.size() << ' ' << data.Z << ' ';
		try {
			real_1d_array x = "[ 6.0, 0.1, 0.1]";
			real_1d_array bndl = "[0.0, 0.0, 0.0]";
			real_1d_array bndu = "[50.0, 50.0, 50.0]";
			minbleicstate state;
			minbleicreport rep;

			minbleiccreatef(x, diffstep, state);
			minbleicsetbc(state, bndl, bndu);
			minbleicsetcond(state, epsg, epsf, epsx, maxits);
			alglib::minbleicoptimize(state, LogLikelihoodFunc2_3var, NULL, &data);
			minbleicresults(state, x, rep);

			cerr << "*3~ " << x[0] << ' ' << x[1] << ' ' << x[2] << ' ' << rep.terminationtype << ' ';

			if ( rep.terminationtype != -8 ) {
				double llh = llh2_3var( data.N, data.precision, data.y, x[0], x[1], x[2]);
				cerr << llh;

				for (size_t i(0); i != 3; i++ )
					mLP[llh].push_back(x[i]);

				mLP[llh].push_back(-1.0); // as a seperateor
			}
			cerr << endl;
		}
		catch (alglib::ap_error &e) {
			cerr << "catch error: " << e.msg << " at insertSize=" << isert << endl;
		}


		// optimizaion log-likelihood function2: 2var: mu and phi  for Pois_dis + Pois model
		//
		cerr << isert << ' ' << data.y.size() << ' ' << data.Z << ' ';
		try {
			real_1d_array x = "[0.1, 0.1]";
			real_1d_array bndl = "[0.0, 0.0]";
			real_1d_array bndu = "[50.0, 50.0]";
			minbleicstate state;
			minbleicreport rep;

			minbleiccreatef(x, diffstep, state);
			minbleicsetbc(state, bndl, bndu);
			minbleicsetcond(state, epsg, epsf, epsx, maxits);
			alglib::minbleicoptimize(state, LogLikelihoodFunc2, NULL, &data);
			minbleicresults(state, x, rep);

			double lambda = data.Z/data.N/x[0];
			cerr << "*2~ " << lambda << ' ' << x[0] << ' ' << x[1] << ' ' << rep.terminationtype << ' ';

			if ( rep.terminationtype != -8 ) {
				double llh = llh2_3var( data.N, data.precision, data.y, lambda, x[0], x[1]);
				cerr << llh;

				mLP[llh].push_back(lambda);
				for (size_t i(0); i != 2; i++ )
					mLP[llh].push_back(x[i]);

				mLP[llh].push_back(-1.0); // as a seperateor
			}
			cerr << endl;
		}
		catch (alglib::ap_error &e) {
			cerr << "catch error: " << e.msg << " at insertSize=" << isert << endl;
		}


/*
		if ( mLP.size() > 0 ) {
			map<double, vector<double> >::reverse_iterator it = mLP.rbegin();
			vector<double> &v = it->second;

			cout << "~~~ " << isert << ' ' << data.y.size() << ' ' << data.Z
				<< ' ' << it->first << ' ' << v[0] << ' ' << v[1];

			if ( v[2] != -1.0 )
				cout << ' ' << v[2] << ' ' << "3var" << endl;
			else
				cout << ' ' << "2var" << endl;

			for ( double y(1.0); y <= max_y; y++ ) {
				cout << y;

				double pr;
				for ( double k(1.0); k < max_y; k++ ) {
					if ( v[2] != -1.0 ) {
						pr = prob_k_given_y(v[0], v[1], v[2], y, data.precision, k);
					}
					else {
						pr = prob_k_given_y_2lambda(v[0], v[1], y, data.precision, k);
					}

					if ( pr > 1e-3 ) {
						cout << ' ' << k << ':' << pr;
					}
				}
				cout << endl;
			}
		}
		*/

	}
	inf.close();

	return 0;
}


