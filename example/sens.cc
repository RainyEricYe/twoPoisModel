/*
 * Estimate sensitivity for each point mutation
 *
 * considering depth, read distribution, read family, strand bias, mutation frequency
 * ignoring indel yet
 *
 *  v1.0    20200720    yerui@genomics.cn
 *
 */

#define VERSION "v1.0"

#include "../main.h"
#include "Option.h"

#include "SeqLib/BamReader.h"
#include "SeqLib/BamWriter.h"
#include <boost/tokenizer.hpp>

using namespace std;
using namespace SeqLib;

void usage(char *prog) {
	cout << "Program: \n"
		"Version: " << VERSION << "\n"
		"Contact: yerui <yerui@genomics.cn>\n\n"
		"Usage: " << prog << " -o outf mut.f dcs.sort.bam raw.sort.bam\n\n"
		"       -o [s]    outfile\n"
		"\n" << endl;
}

int main( int argc, char **argv )
{
	Option opt;

	int c;
	while ( (c=getopt(argc,argv,"o:h")) != -1 ) {
		switch (c) {
			case 'o': opt.outFile = optarg; break;
			case 'h':
			default:  usage(argv[0]);   exit(1);
		}
	}

	if ( argc < optind + 3 ) usage(argv[0]), exit(1);

	string mutFile( argv[optind]   );
	string dcsFile( argv[optind+1] );
	string rawFile( argv[optind+2] );

	ifstream mutf( mutFile.c_str() );
	ofstream outf( opt.outFile.c_str() );

	SeqLib::BamReader dcsBam, rawBam;

	if ( !mutf.is_open() ) cerr << "open error: " << mutFile << endl, exit(1);
	if ( !outf.is_open() ) cerr << "open error: " << opt.outFile << endl, exit(1);

	if ( !dcsBam.Open(dcsFile) ) cerr << "open error: " << dcsFile << endl, exit(1);
	if ( !rawBam.Open(rawFile) ) cerr << "open error: " << rawFile << endl, exit(1);

	// read each line of mutFile
	string mutLine;
	getline(mutf, mutLine); // skip head line

	while ( getline(mutf, mutLine) ) {
		istringstream ist(mutLine);
		string chr, ref;
		double pos, depth, mutN;

		ist >> chr >> ref >> pos >> depth >> mutN;
		if ( mutN == 0 ) continue;

		double Nb(0);
		map<double, map<double, double> > wBin, cBin, tBin;  // watson, crick, total, keys: insertSize, startPos
		map<double, vector<double> > para; // parameters of model. key: insertSize
		map<double, double> wSize, cSize, tSize;

		rawBam.Reset();
		GenomicRegion gr(chr, to_string(pos-1), to_string(pos), rawBam.Header() ); // region is 0-based
		rawBam.SetRegion(gr);

		BamRecord br;
		while ( rawBam.GetNextRecord(br)) {
			if ( !br.ProperPair()
					or br.QCFailFlag()
					or br.Interchromosomal()
					or br.SecondaryFlag()
			   ) continue;
			/*
			   if ( opt.filtSoftClip
			   && ( br.Position() != br.PositionWithSClips() || br.PositionEnd() != br.PositionWithSClips() )
			   ) continue;
			   */
			double insert( br.InsertSize() ) ;
			double startPos( br.Position() ) ;
			double endPos( br.PositionEnd() );

			if ( startPos + opt.softEndTrim > pos
					or endPos - opt.softEndTrim <= pos
			   ) continue;

			if ( insert < 0 ) {
				startPos = br.MatePosition();
				insert *= -1;
			}

			if ( Nb == 0 )
				Nb = br.Sequence().size() * 2 - opt.softEndTrim * 4;

			if ( ( br.FirstFlag() && !br.ReverseFlag() )
					or ( !br.FirstFlag() && br.ReverseFlag() )
			   ) {
				wBin[insert][startPos]++;
				wSize[insert]++;
			}
			else {
				cBin[insert][startPos]++;
				cSize[insert]++;
			}

			tBin[insert][startPos]++;
			tSize[insert]++;
		}

		for ( auto it : tBin ) {

			double isert(it.first);

			fn_data data;
			data.precision = 1.0e-12;
			data.N = 146.0;
			data.Z = 0.0;

			double epsg = 1.0e-6;
			double epsf = 0;
			double epsx = 0;
			ae_int_t maxits = 0;
			double diffstep = 1.0e-6;

			for ( auto t : it.second ) {
				data.y.push_back( t.second );
			}

			double max_y(0.0);

			for ( size_t i(0); i != data.y.size(); i++) {
				data.Z += data.y[i];
				if (data.y[i] > max_y) max_y = data.y[i];
			}

			map<double, vector<double> > mLP;

			// solve diff_prob_y0_2lamda to get lambda; 1var ; lam
			// pois + pois (pp)
			//
			cerr << isert << ' ' << data.y.size() << ' ' << data.Z << ' ';
			//cout << "~ " << isert << ' ' << data.y.size() << ' ' << data.Z << ' ';
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
				cerr << "pp1~ " << x[0] << ' ' << expect_lam2 << ' ' << rep.terminationtype << ' ';
				para[isert].push_back(x[0]);
				para[isert].push_back(expect_lam2);


				if ( rep.terminationtype != -8 ) {
					double llh = llh_2lambda(data.N, data.precision, data.y, x[0], expect_lam2);
					cerr << llh;
					//cout << x[0] << ' ' << expect_lam2 << ' ' << llh;

					mLP[llh].push_back(x[0]);
					mLP[llh].push_back(expect_lam2);

					mLP[llh].push_back(-1.0); // as a seperateor
				}
				cerr << endl;
				//cout << endl;

			}
			catch (alglib::ap_error &e) {
				cerr << "catch error: " << e.msg << " at insertSize=" << isert << endl;
			}
		}

		double pNotDetect(0.0);

		dcsBam.Reset();
		dcsBam.SetRegion(gr);

		while ( dcsBam.GetNextRecord(br) ) {
			double insert( br.InsertSize() ) ;
			if ( insert < 0 ) insert *= -1;

			map<double, vector<double> >::const_iterator it = para.find(insert);
			if ( it == para.end() ) continue;

			string spType("sp"), spStr("");
			if ( br.GetZTag(spType, spStr) ) {
				if ( spStr != "1" ) continue;
			}

			double pNotDetectInFam(0.0);

			string fsTyps("fs"), fsStr("");
			if ( br.GetZTag(fsTyps, fsStr) ) {
				typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
				std::string x = ",";
				boost::char_separator<char> sep{ x.c_str() };
				tokenizer tok{fsStr, sep};
				cout << "fsStr " << fsStr << endl;
				for (const auto &t : tok) std::cout << t << '\n';
			}
		}
	}

	dcsBam.Close();
	rawBam.Close();

	mutf.close();
	outf.close();

	exit(0);
}
