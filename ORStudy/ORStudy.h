#ifndef SusyNtuple_ORStudy_h
#define SusyNtuple_ORStudy_h

//ROOT
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"

//SusyNtuple
#include "SusyNtuple/SusyNtAna.h"
#include "SusyNtuple/SusyNtTools.h"

//std/stl
#include <string>
#include <fstream>

/////////////////////////////////////////////////////////////
//
// ORStudy
// Class auto-generated with SusyNtuple/make_susy_skeleton on 2018-02-01 00:26
//
/////////////////////////////////////////////////////////////

// for TSelector analysis loopers processing susyNt you MUST inherit from SusyNtAna
// in order to pick up the susyNt class objects
class ORStudy : public SusyNtAna
{

    public :
        ORStudy();
        virtual ~ORStudy() {};

        void set_debug(int dbg) { m_dbg = dbg; }
        int dbg() { return m_dbg; }

        void set_mu_bounds(float lower, float higher) { m_lower_mu = lower; m_higher_mu = higher; }

        void set_chain(TChain* chain) { m_input_chain = chain; }
        TChain* chain() { return m_input_chain; }

        void set_tagger(std::string tagger_name) { m_tagger = tagger_name; }
        void set_release(std::string release) { m_release = release; }

        ////////////////////////////////////////////
        // analysis methods
        ////////////////////////////////////////////

        // histograms
        void initialize_histograms();
        void write_histograms();

        // fill truth containers
        void selectTruthObjects();
        void clearTruthObjects();
        void performOverlap();
        void j_e_overlap(ElectronVector& electrons, JetVector& jets);
        void e_j_overlap(ElectronVector& electrons, JetVector& jets);
        void look_at_jet(const Jet* jet);
        void look_at_high_pt_ratio_electron(const Electron* el);
        void look_at_overlapping_electron(const Electron* el);
        void look_at_photon_electron_isolation();
        float get_photon_sum_et(const Susy::Particle* particle, float dRcone, bool require_fsr_photon = false);
        Susy::Jet* electron_matches_bjet(const Susy::Electron* ele, bool& is_matched);

        // standard ATLAS event cleaning
        bool passEventCleaning(const MuonVector& preMuons, const MuonVector& baseMuons,
                const JetVector& baseJets);

        ////////////////////////////////////////////
        // TSelector methods override
        ////////////////////////////////////////////
        virtual void Begin(TTree* tree); // Begin is called before looping on entries
        virtual Bool_t Process(Long64_t entry); // Main event loop function called on each event
        virtual void Terminate(); // Terminate is called after looping has finished


    private :
        int m_dbg;
        float m_lower_mu;
        float m_higher_mu;
        TChain* m_input_chain; // the TChain object we are processing
        std::string m_tagger;
        std::string m_release;
        float m_mc_weight;

        // Truth objects
        TruthParticleVector m_truth_particles;
        TruthJetVector m_truth_jets;
        TruthMet m_truth_met;

        TruthParticleVector m_truth_electrons;
        TruthParticleVector m_truth_muons;
        TruthParticleVector m_truth_leptons;

        TruthParticleVector m_truth_photons;

        ElectronVector m_original_base_electrons;
        JetVector m_original_base_jets;

        // Histograms and the like
        TH1F* h_l0pt;
        TH1F* h_l0pt_lowmu;
        TH1F* h_l0pt_highmu;
        TH1F* h_l0eta;
        TH1F* h_l1pt;
        TH1F* h_l1eta;

        TH1F* h_l0pt300;
        TH1F* h_l1pt300;
        TH1F* h_m0pt;
        TH1F* h_m1pt;
        TH1F* h_dphill;
        TH1F* h_drll;

        // ele and jet
        TH1F* h_min_drll_ele_jet;
        TH1F* h_pt_ratio_ele_jet;
        TH2F* h2_pt_ratio_min_drll_ele_jet;
        TH2F* h2_pt_ratio_ntracks;
        TH2F* h2_mv2_min_drll_ele_jet;
        TH2F* h2_dl1_min_drll_ele_jet;
        TH1F* h_closest_truth_obj_jet;
        TH1F* h_type_closest_ele;
        TH2F* h2_type_closest_ele_min_drjet;
        TH1F* h_origin_closest_ele;
        TH2F* h2_origin_closest_ele_min_drjet;
        TH1F* h_flavor_closest_jet;
        TH2F* h2_flavor_closest_jet_min_drjet;
        TH2F* h2_flavor_closest_jet_ntracks;

        TH2F* h2_mv2_ntracks;
        TH2F* h2_dl1_ntracks;
        TH2F* h2_mindr_ntracks;
        TH2F* h2_mindr_ntracks_bjet;

        int n_kin;
        int n_kin_b;

        TH1F* h_nbjets_survive;
        TH1F* h_nbjets_survive_matched;

        TH1F* h_hptrele_min_dr_truth;
        TH1F* h_hptrele_nearest_truth_ele_type;
        TH1F* h_hptrele_nearest_truth_ele_origin;
        TH1F* h_hptrele_nearest_truth_ele_status;
        TH1F* h_hptrele_relptcone20;
        TH1F* h_hptrele_relptcone30;
        TH1F* h_hptrele_d0;
        TH1F* h_hptrele_d0sig;

        TH1F* h_overele_min_dr_truth;
        TH1F* h_overele_nearest_truth_ele_type;
        TH1F* h_overele_nearest_truth_ele_origin;
        TH1F* h_overele_nearest_truth_ele_status;
        TH1F* h_overele_relptcone20;
        TH1F* h_overele_relptcone30;
        TH1F* h_overele_d0;
        TH1F* h_overele_d0sig;

        TH1F* h_dr_eljet_photon;
        TH1F* h_dr_eljet_fsr_brem_photon;
        TH1F* h_dr_eljet_fsr_brem_photon_pt0_100;
        TH1F* h_dr_eljet_fsr_brem_photon_pt100_200;
        TH1F* h_dr_eljet_fsr_brem_photon_pt200_300;
        TH1F* h_dr_eljet_fsr_brem_photon_pt300_400;

        TH1F* h_dr_eljet_fsr_brem_conv_photon;
        TH1F* h_dr_eljet_brem_photon;
        TH1F* h_dr_eljet_brem_photon_pt0_100;
        TH1F* h_dr_eljet_brem_photon_pt100_200;
        TH1F* h_dr_eljet_brem_photon_pt200_300;
        TH1F* h_dr_eljet_brem_photon_pt300_400;
        TH1F* h_dr_eljet_brem_conv_photon;

        TH2F* h_dr_eljet_photon_vs_lpt;
        TH2F* h_dr_eljet_fsr_brem_photon_vs_lpt;

        TH1F* h_etcone_gamma_electron_bjet;
        TH1F* h_etcone_gamma_electron;
        TH1F* h_etcone_gamma_electron_bjet_fsr;
        TH1F* h_etcone_gamma_electron_fsr;
        TH1F* h_etcone_gamma_electron_pt200_bjet;
        TH1F* h_etcone_gamma_electron_pt200;
        TH1F* h_etcone_gamma_electron_pt200_bjet_fsr;
        TH1F* h_etcone_gamma_electron_pt200_fsr;
        TH1F* h_etcone_gamma_electron_pt300_bjet;
        TH1F* h_etcone_gamma_electron_pt300;
        TH1F* h_etcone_gamma_electron_pt300_bjet_fsr;
        TH1F* h_etcone_gamma_electron_pt300_fsr;

}; //class


#endif
