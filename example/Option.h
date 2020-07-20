#ifndef OPTION_H_
#define OPTION_H_

#include <cstdio>
#include <string>

using namespace std;

class Option {
    public:
        Option():
            baseQuaCutoff(20),
            mapQuaCutoff(30),
            minSupOnEachStrand(3),
            maxSupOnEachStrand(3000),
            Ncutoff(0.1),
            minFractionInFam(0.002),
            freqPrecision(0.00001),
            lhrGapCutoff(2.0),
            phredOffset(33),
            minSupOnHaplotype(3),
            filtSoftClip(false),
            outFile("out"),
            sscsOut(false),
            singleOut(false),
            debug(false),
            pvalue(0.001),
            pcrError(1.0e-5),
            softEndTrim(5),
            randNread(0) {}

        ~Option(){}

        int baseQuaCutoff;
        int mapQuaCutoff;
        long minSupOnEachStrand;
        long maxSupOnEachStrand;

        double Ncutoff;
        double minFractionInFam;
        double freqPrecision;
        double lhrGapCutoff;

        int phredOffset;
        long minSupOnHaplotype;
        bool filtSoftClip;
        string outFile;
        bool sscsOut;
        bool singleOut;
        bool debug;
        double pvalue;
        double pcrError;
        long softEndTrim;
        size_t randNread;
};

#endif // OPTION_H_
