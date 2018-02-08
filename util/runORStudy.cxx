//SusyNtuple
#include "ORStudy/ORStudy.h"
#include "SusyNtuple/ChainHelper.h"
#include "SusyNtuple/string_utils.h"

//std/stl
#include <iostream>
#include <cstdlib>
#include <string>
using namespace std;

//ROOT
#include "TChain.h"

//////////////////////////////////////////////////////
//
// runORStudy 
// Executable auto-generated with SusyNtuple/make_susy_skeleton on 2018-02-01 00:26
//
//
//////////////////////////////////////////////////////


void help()
{
    cout << "----------------------------------------------------------" << endl;
    cout << " runORStudy" << endl;
    cout << endl;
    cout << "  Options:" << endl;
    cout << "   -n          number of events to process (default: all)" << endl;
    cout << "   -d          debug level (integer) (default: 0)" << endl;
    cout << "   -i          input file (ROOT file, *.txt file, or directory)" << endl;
    cout << "   --tagger    b-tagging algorithm to use (mv2 or dl1) [default: mv2]" << endl;
    cout << "   --release   release 20 or 21 (sets which version of mv2c10 to use) [default: none]" << endl;
    cout << "   --lower-mu  lower requirement on avgMu [default: 0]" << endl;
    cout << "   --upper-mu  higher bound on avgMu [default: 100]" << endl;
    cout << "   -h          print this help message" << endl;
    cout << endl;
    cout << "  Example Usage:" << endl;
    cout << "   runORStudy -i susyNt.root -n 500" << endl;
    cout << "----------------------------------------------------------" << endl;
}

int main(int argc, char** argv)
{

    /////////////////////////
    // cmd line options
    /////////////////////////

    int n_events = -1;
    int dbg = 0;
    string input = "";
    string tagger_name = "mv2";
    string release = "";
    int lower_mu = 0;
    int higher_mu = 100;

    for(int i = 1; i < argc; i++) {
        if      (strcmp(argv[i], "-n") == 0) n_events = atoi(argv[++i]);
        else if (strcmp(argv[i], "-d") == 0) dbg = atoi(argv[++i]);
        else if (strcmp(argv[i], "-i") == 0) input = argv[++i];
        else if (strcmp(argv[i], "--tagger") == 0) tagger_name = argv[++i];
        else if (strcmp(argv[i], "--release") == 0) release = argv[++i];
        else if (strcmp(argv[i], "--lower-mu") == 0) lower_mu = atoi(argv[++i]);
        else if (strcmp(argv[i], "--upper-mu") == 0) higher_mu = atoi(argv[++i]);
        else if (strcmp(argv[i], "-h") == 0) { help(); return 0; }
        else {
            cout << "runORStudy    Unknown command line argument '" << argv[i] << "', exiting" << endl;
            help();
            return 1;
        }
    } // i

    if(input.empty()) {
        cout << "runORStudy    You must specify an input" << endl;
        return 1;
    }

    if(! (tagger_name=="mv2" || tagger_name=="dl1")) {
        cout << "runORStudy    ERROR Invalid tagger name (=" << tagger_name << ") provided" << endl;
        return 1;
    }

    if(release.empty()) {
        cout << "runORStudy    You must specify a release (20 or 21)" << endl;
        return 1;
    }

    if(! (release == "20" || release == "21")) {
        cout << "runORStudy   ERROR Invalid release (=" << release << ") provided" << endl;
        return 1;
    }


    /////////////////////////////////////////////////////////
    // Build the TChain object
    // For SusyNtuple analysis, the chain name is susyNt
    /////////////////////////////////////////////////////////
    TChain* chain = new TChain("susyNt");

    // use ChainHelper to infer the input type (ROOT file, *.txt, or dir/)
    // and build the full chain of the input files
    // (c.f. SusyNtuple/ChainHelper.h)
    ChainHelper::addInput(chain, input, dbg>0);
    Long64_t n_entries_in_chain = chain->GetEntries();
    // let's see what it looks like
    chain->ls();

    /////////////////////////////////////////////////////////
    // Build the TSelector object
    // SusyNt analyses inheriting from SusyNtAna must
    // build their own TSelector looper
    /////////////////////////////////////////////////////////
    ORStudy* analysis = new ORStudy();

    // set to do the 2 lepton analysis object selection (c.f. SusyNtuple/AnalysisType.h)
    // the AnalysisType configures all of the selector tools (c.f. SusyNtuple/SusyNtTools.h)
    analysis->setAnaType(AnalysisType::Ana_Stop2L);

    analysis->set_debug(dbg);
    analysis->setSampleName(ChainHelper::sampleName(input, dbg>0)); // SusyNtAna setSampleName (c.f. SusyNtuple/SusyNtAna.h)
    analysis->set_chain(chain); // propagate the TChain to the analysis
    analysis->set_tagger(tagger_name);
    analysis->set_release(release);
    analysis->set_mu_bounds(lower_mu, higher_mu);

    // for using the TriggerTools (c.f. SusyNtuple/TriggerTools.h) we
    // must provide the first file in our chain to initialize the
    // underlying TriggerTool object (to inspect the "trig" histogram
    // stored in the susyNt file)

    // we use the inherited SusyNtTools object from SusyNtAna base class
    // to initialize the underlying TriggerTool object
    analysis->nttools().initTriggerTool(ChainHelper::firstFile(input, dbg>0));

    if(n_events < 0) n_events = n_entries_in_chain;

    cout << "---------------------------------------------------------" << endl;
    cout << " Total entries in input chain          : " << n_entries_in_chain << endl;
    cout << " Total entries to process for analysis : " << n_events << endl;
    cout << "---------------------------------------------------------" << endl;

    
    // call TChain Process to star the TSelector looper over the input TChain
    if(n_events > 0) chain->Process(analysis, input.c_str(), n_events);

    cout << endl;
    cout << "runORStudy    Analysis loop done" << endl;

    delete chain;
    return 0;
}
