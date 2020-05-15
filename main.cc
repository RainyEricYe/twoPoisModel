
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

		//cerr << "is: " << isert << " Z: " << data.Z << " fBin: " << data.y.size() << ' ';
		cerr << isert << ' ' << data.y.size() << ' ' << data.Z << ' ';

		// optimizaion log-likelihood function
		//
		real_1d_array x = "[7.0, 0.5]";
		real_1d_array bndl = "[0.0, 0.0]";
		real_1d_array bndu = "[300.0, 300.0]";
		minbleicstate state;
		minbleicreport rep;

		double epsg = 0.000001;
		double epsf = 0;
		double epsx = 0;
		ae_int_t maxits = 0;

		//
		// This variable contains differentiation step
		//
		double diffstep = 1.0e-6;

		//
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

		//
		// ...and evaluate these results
		//
		printf("\t%d\t", int(rep.terminationtype)); // EXPECTED: 4
//		printf("%s\n", x.tostring(4).c_str()); // EXPECTED: [-1,1]
		double lambda = data.Z/data.N/x[0];
		cerr << lambda << ' ' << x[0] << ' ' << x[1] <<  endl;
/*
		map<double, double> y_pr, y_raw_pr;
		for ( size_t i(0); i != data.y.size(); i++ ) {
			y_pr[ data.y[i] ] = prob_y( lambda, x[0], x[1], data.y[i], data.precision );
			y_raw_pr[ data.y[i] ] += 1.0;
		}

		for ( size_t i(0); i != data.y.size(); i++ ) {
			cerr << data.y[i] << ' ' << y_pr[ data.y[i] ] << ' ' << y_raw_pr[ data.y[i] ] / data.N << endl;
		}
		*/
	}
	inf.close();

    return 0;

}


