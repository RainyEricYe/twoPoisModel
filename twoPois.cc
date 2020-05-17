
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

		while (stm >> t) {
			data.y.push_back(t);
		}

		for ( size_t i(0); i != data.y.size(); i++) {
			data.Z += data.y[i];
		}
/*
		cerr << isert << ' ' << data.y.size() << ' ' << data.Z << ' ';

		// optimizaion log-likelihood function
		//
		try {
			real_1d_array x = "[7.0, 0.5]";
			real_1d_array bndl = "[0.0, 0.0]";
			real_1d_array bndu = "[300.0, 300.0]";
			minbleicstate state;
			minbleicreport rep;

			double epsg = 0.000001;
			double epsf = 0;
			double epsx = 0;
			ae_int_t maxits = 0;

			// This variable contains differentiation step
			double diffstep = 1.0e-6;

			// Now we are ready to actually optimize something:
			// * first we create optimizer
			// * we add boundary constraints
			// * we tune stopping conditions
			// * and, finally, optimize and obtain results...
			//
			minbleiccreatef(x, diffstep, state);
			minbleicsetbc(state, bndl, bndu);
			minbleicsetcond(state, epsg, epsf, epsx, maxits);
			alglib::minbleicoptimize(state, LogLikelihoodFunc, NULL, &data);
			minbleicresults(state, x, rep);

			printf("\t%d\t", int(rep.terminationtype)); // EXPECTED: 4
			//		printf("%s\n", x.tostring(4).c_str()); // EXPECTED: [-1,1]

			double lambda = data.Z/data.N/x[0];
			cerr << "1~ " << lambda << ' ' << x[0] << ' ' << x[1] << ' ' << rep.terminationtype << ' ';
			if ( rep.terminationtype != 8 ) {
				cerr << llh(lambda, x[0], x[1], data.y, data.precision, data.N);
			}
			cerr << endl;
		}
		catch (alglib::ap_error &e) {
			cerr << "catch error: " << e.msg << " at insertSize=" << isert << endl;
		}
*/
		cerr << isert << ' ' << data.y.size() << ' ' << data.Z << ' ';

		double f_lam(0.0), f_lam2(0.0);
		// maximize the 2lambda model
		//
		try {
			real_1d_array x = "[1.0, 6.0]";
			real_1d_array bndl = "[0.0, 0.0]";
			real_1d_array bndu = "[300.0, 300.0]";
			minbleicstate state;
			minbleicreport rep;

			double epsg = 0.000001;
			double epsf = 0;
			double epsx = 0;
			ae_int_t maxits = 0;

			// This variable contains differentiation step
			double diffstep = 1.0e-6;

			minbleiccreatef(x, diffstep, state);
			minbleicsetbc(state, bndl, bndu);
			minbleicsetcond(state, epsg, epsf, epsx, maxits);
			alglib::minbleicoptimize(state, LogLikelihoodFunc_2lambda, NULL, &data);
			minbleicresults(state, x, rep);

//			printf("\t%d\t", int(rep.terminationtype)); // EXPECTED: 4
			//		printf("%s\n", x.tostring(4).c_str()); // EXPECTED: [-1,1]

			double expect_lam2 = data.Z/data.N/x[0];
			cerr << "2~ " << expect_lam2 << ' ' << x[0] << ' ' << x[1] << ' ' << rep.terminationtype << ' ';
			if ( rep.terminationtype != 8 ) {
				cerr << llh_2lambda(data.N, data.precision, data.y, x[0], x[1]);
			}
			cerr << endl;
		}
		catch (alglib::ap_error &e) {
			cerr << "catch error: " << e.msg << " at insertSize=" << isert << endl;
		}

		// solve diff_prob_y0_2lamda to get lambda
		cerr << isert << ' ' << data.y.size() << ' ' << data.Z << ' ';

		try {
			real_1d_array x = "[1.0]";
			real_1d_array bndl = "[0.0]";
			real_1d_array bndu = "[300.0]";
			minbleicstate state;
			minbleicreport rep;

			double epsg = 0.000001;
			double epsf = 0;
			double epsx = 0;
			ae_int_t maxits = 0;

			// This variable contains differentiation step
			double diffstep = 1.0e-6;

			minbleiccreatef(x, diffstep, state);
			minbleicsetbc(state, bndl, bndu);
			minbleicsetcond(state, epsg, epsf, epsx, maxits);
			alglib::minbleicoptimize(state, diff_prob_y0_2lambda, NULL, &data);
			minbleicresults(state, x, rep);

//			printf("\t%d\t", int(rep.terminationtype)); // EXPECTED: 4
			//		printf("%s\n", x.tostring(4).c_str()); // EXPECTED: [-1,1]

			double expect_lam2 = data.Z/data.N/x[0];
			cerr << "3~ " << expect_lam2 << ' ' << x[0] << ' ' << expect_lam2 << ' ' << rep.terminationtype << ' ';
			if ( rep.terminationtype != 8 ) {
				cerr << llh_2lambda(data.N, data.precision, data.y, x[0], expect_lam2);
			}
			cerr << endl;

		}
		catch (alglib::ap_error &e) {
			cerr << "catch error: " << e.msg << " at insertSize=" << isert << endl;
		}

	}
	inf.close();

	return 0;

}


