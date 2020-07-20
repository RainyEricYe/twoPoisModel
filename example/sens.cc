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

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>

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
	while ( (c=getopt(argc,argv,"o:m:dh")) != -1 ) {
		switch (c) {
			case 'o': opt.outFile = optarg; 		break;
			case 'm': opt.mutFreq = atof(optarg); 	break;
			case 'd': opt.debug = true; 			break;
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

	string mutLine;
	getline(mutf, mutLine); // skip head line

	while ( getline(mutf, mutLine) ) {
		istringstream ist(mutLine);
		string chr, ref;
		double pos, depth, mutN;

		ist >> chr >> ref >> pos >> depth >> mutN;
		if ( mutN == 0 ) continue;

//		if (opt.debug) cout << chr << ' ' << ref << ' ' << pos << ' ' << depth << ' ' << mutN << endl;

		double Nb(0);
		map<double, map<double, double> > wBin, cBin, tBin;
		// watson, crick, total --  keys: insertSize, startPos

		map<double, vector<double> > para;
		// parameters of model --  key: insertSize

		map<double, double> wSize, cSize, tSize;

		// read raw.sort.bam to collect original reads distribution
		//
		rawBam.Reset();
		GenomicRegion gr(chr, to_string(pos-1), to_string(pos), rawBam.Header() ); // region is 0-based
		rawBam.SetRegion(gr);

		BamRecord br;
		while ( rawBam.GetNextRecord(br)) {
			uint32_t flag = br.AlignmentFlag();

			if ( br.ChrID() != br.MateChrID()
					or !br.ProperPair()
					or flag & 0x4 or flag & 0x8
					or flag & 0x100 or flag & 0x200 or flag & 0x400 or flag & 0x800
			   ) continue;
			/*
			   if ( !br.ProperPair()
			   or br.QCFailFlag()
			   or br.Interchromosomal()
			   or br.SecondaryFlag()
			   ) continue;

			   if ( opt.filtSoftClip
			   && ( br.Position() != br.PositionWithSClips() || br.PositionEnd() != br.PositionWithSClips() )
			   ) continue;
			   */

			double insert( br.InsertSize() ) ;
			double startPos( br.Position() ) ;
			double endPos( br.PositionEnd() );

//			if (opt.debug) cout << insert << ' ' << startPos+1 << ' ' << endPos+1 << endl;

			if ( startPos + opt.softEndTrim + 1 > pos
					or endPos - opt.softEndTrim + 1 <= pos
			   ) continue;

			if ( insert < 0 ) {
				startPos = br.MatePosition();
				insert *= -1;
			}

//			if (opt.debug) cout << br;

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

		// estimate parameters based on collected data
		//
		for ( auto tb : tBin ) {

			double isert(tb.first);
			const map<double, double> &mPosY = tb.second;

			fn_data data;
			data.precision = opt.dataPrecision;
			data.N = 146.0;
			data.Z = 0.0;

			double epsg = 1.0e-6;
			double epsf = 0;
			double epsx = 0;
			ae_int_t maxits = 0;
			double diffstep = 1.0e-6;

			for ( auto t : mPosY ) {
				data.y.push_back( t.second );
			}

			if (opt.debug) {
				cout << isert;
				for ( auto t : mPosY )
					cout << ' ' << t.second;
				cout << endl;
			}

			double max_y(0.0);
			for ( size_t i(0); i != data.y.size(); i++) {
				data.Z += data.y[i];
				if (data.y[i] > max_y) max_y = data.y[i];
			}

			// solve diff_prob_y0_2lamda to get lambda; 1var ; lam
			// pois + pois (pp)
			//
			//cerr << isert << ' ' << data.y.size() << ' ' << data.Z << ' ';
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
			//	cerr << "pp1~ " << x[0] << ' ' << expect_lam2 << ' ' << rep.terminationtype << ' ';
				para[isert].push_back(x[0]);
				para[isert].push_back(expect_lam2);
//				if (opt.debug) cout << "isert " << isert << " lambda " << x[0] << " lambda2 " << expect_lam2 << endl;

				if ( rep.terminationtype != -8 ) {
					//double llh = llh_2lambda(data.N, data.precision, data.y, x[0], expect_lam2);
					//cerr << llh;
				}
			//	cerr << endl;
			}
			catch (alglib::ap_error &e) {
				cerr << "catch error: " << e.msg << " at insertSize=" << isert << endl;
			}
		}


		// read dcs.bam file
		double pNotDetect(0.0);

		dcsBam.Reset();
		dcsBam.SetRegion(gr);

		while ( dcsBam.GetNextRecord(br) ) {
			double insert( br.InsertSize() ) ;
			if ( insert < 0 ) insert *= -1;

			if (opt.debug) cout << "insertSize " << insert << endl;

			map<double, vector<double> >::const_iterator it = para.find(insert);
			if ( it == para.end() ) continue;
			const vector<double> &v = it->second;

			string spType("sp"), spStr("");
			if ( br.GetZTag(spType, spStr) ) {
				if ( spStr != "1" ) continue;
			}

			double pNotDetectInFam(0.0);
			string fsTyps("fs"), fsStr("");
			if ( br.GetZTag(fsTyps, fsStr) ) {
				vector<string> fsVec;
				boost::split(fsVec, fsStr, boost::is_any_of(","), boost::token_compress_on );

				double wS = stod(fsVec[0]);
				double cS = stod(fsVec[1]);
				double tS = wS + cS;

				if (opt.debug) cout << "w_s " << wS << " c_s " << cS << " s " << tS << endl;

				double sumP = prob_y_2lambda(v[0], v[1], tS, opt.dataPrecision);

				double target_k = int(tS / v[1]);
				double start_k = ( target_k - 30 <= 1.0 ? 1.0 : target_k - 30 );

				for ( double k(start_k); k < target_k + 30; k++ ) {

					double pr = p_k_given_y_2lambda(v[0], v[1], tS, opt.dataPrecision, k) / sumP;
					if ( pr < opt.dataPrecision ) continue;

					if (opt.debug) cout << "~" << k << ' ' << tS << ' ' << v[0] << ' ' << v[1] << ' ' << pr << endl;

					double wNotDet = notDetect(1/k, wS);
					double cNotDet = notDetect(1/k, cS);

					pNotDetectInFam += pr * ( wNotDet + cNotDet - wNotDet * cNotDet );
					if (opt.debug) cout << "~~w_n " << wNotDet << " c_n " << cNotDet << " p_fam " << pNotDetectInFam << endl;
				}
			}
			else {
				cerr << "NO fs:Z tag in the bam file: " << dcsFile << endl;
				exit(1);
			}

			pNotDetect += log( 1 - opt.mutFreq + opt.mutFreq * pNotDetectInFam );
		}

		double sensitivity = 1 - exp(pNotDetect);
		cout << "freq:sensitivity => " << opt.mutFreq << '\t' << sensitivity << endl;
	}

	dcsBam.Close();
	rawBam.Close();

	mutf.close();
	outf.close();

	exit(0);
}
