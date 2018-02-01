#include "ORStudy/ORStudy.h"

// SusyNtuple
#include "SusyNtuple/KinematicTools.h"
#include "SusyNtuple/SusyDefs.h"
using namespace Susy; // everything in SusyNtuple is in this namespace

//ROOT
#include "TFile.h"

// std/stl
#include <iomanip> // setw
#include <iostream>
#include <string>
#include <vector>
#include <sstream> // stringstream, ostringstream
using namespace std;

//////////////////////////////////////////////////////////////////////////////
ORStudy::ORStudy() :
    m_dbg(0),
    m_input_chain(nullptr),
    m_tagger("mv2"),
    m_release(""),
    m_mc_weight(1.0)
{
}
//////////////////////////////////////////////////////////////////////////////
void ORStudy::Begin(TTree* /*tree*/)
{
    n_kin = 0;
    n_kin_b = 0;
    // call base class' Begin method
    SusyNtAna::Begin(0);
    if(dbg()) cout << "ORStudy::Begin" << endl;

    initialize_histograms();

    return;
}
//////////////////////////////////////////////////////////////////////////////
void ORStudy::initialize_histograms()
{
    // lead lepton pT
    h_l0pt = new TH1F("h_l0pt", "Lead Electron p_{T};p_{T} [GeV]; Entries", 40, 0, 600);
    h_l1pt = new TH1F("h_l1pt", "Sub-lead Electron p_{T};p_{T} [GeV]; Entries", 40, 0, 600);
    h_m0pt = new TH1F("h_m0pt", "Lead #mu p_{T};p_{T} [GeV]; Entries", 40, 0, 600);
    h_m1pt = new TH1F("h_m1pt", "Sub-lead #mu p_{T};p_{T} [GeV]; Entries", 40, 0, 600);
    h_dphill = new TH1F("h_dphill", "#Delta#phi_{ll};#Delta#phi_{ll};Entries", 64, -3.2, 3.2);
    h_drll = new TH1F("h_drll", "#DeltaR_{ll};#DeltaR_{ll};Entries", 60, 0, 6);
    h_min_drll_ele_jet = new TH1F("h_min_drll_ele_jet", "min #DeltaR_{e,b-jet};min #DeltaR_{l,j};Entries", 60, 0, 6);
    h_pt_ratio_ele_jet = new TH1F("h_pt_ratio_ele_jet", "p_{T}^{l}/p_{T}^{j};Ratio;Entries", 40, 0, 2);
    h2_pt_ratio_min_drll_ele_jet = new TH2F("h2_pt_ratio_min_drll_ele_jet",
        "min #DeltaR_{e,b-jet} vs p_{T}^{e}/p_{T}^{j};min #DeltaR_{e,b-jet};p_{T}^{l}/p_{T}^{j}", 60, 0, 6, 40, 0, 2);
    h2_mv2_min_drll_ele_jet = new TH2F("h2_mv2_min_drll_ele_jet", "mv2c10 vs min #DeltaR_{e,b-jet};min#DeltaR_{e,b-jet};mv2c10", 30, 0, 3, 40, 0, 1);
    h2_dl1_min_drll_ele_jet = new TH2F("h2_dl1_min_drll_ele_jet", "DL1 vs min #DeltaR_{e,b-jet};min#DeltaR_{e,b-jet};DL1", 30, 0, 3, 60, -0.5,1.2);

    h_closest_truth_obj_jet = new TH1F("h_closest_truth_obj_jet", ";Object (0=truth jet, 1=truth ele);Entries", 3, 0, 3);
    h_type_closest_ele = new TH1F("h_type_closest_ele", ";Truth Type Closest Truth Electron to b-Jet;Entries", 40, 0, 40);
    h_origin_closest_ele = new TH1F("h_origin_closest_ele", ";Truth Origin Closest Truth Electron to b-Jet;Entries", 50, 0, 50);
    h_flavor_closest_jet = new TH1F("h_flavor_closest_jet", ";Flavor of Closest Truth Jet to b-Jet;Entries", 10, 0, 10);

    h2_type_closest_ele_min_drjet = new TH2F("h2_type_closest_ele_min_drjet",
        ";min #DeltaR_{e,b-jet};MC Truth Type Closest Truth Electron", 30, 0, 3, 6, 0, 6);
    h2_origin_closest_ele_min_drjet = new TH2F("h2_origin_closest_ele_min_drjet",
        ";min #DeltaR_{e,b-jet};MC Truth Origin Closest Truth Electron", 30, 0, 3, 40, 0, 40);
    h2_flavor_closest_jet_min_drjet = new TH2F("h2_flavor_closest_jet_min_drjet",
        ";min #DeltaR_{e,b-jet};Closest Truth Jet Flavor", 30, 0, 30, 6,0,6);

    h2_pt_ratio_ntracks = new TH2F("h2_pt_ratio_ntracks",
        ";p_{T}^{l}/p_{T}^{j};n_{trk}^{j}", 30, 0, 3, 40, 0, 40);
    h2_flavor_closest_jet_ntracks = new TH2F("h2_flavor_closest_jet_ntracks",
        ";Flavor of Closest Truth Jet to Reco-bJet;Reco-bJet nTracks", 6, 0, 6, 30, 0, 30);

    h_nbjets_survive = new TH1F("h_nbjets_survive", ";n-bjets;Entries",6,0,6);
    h_nbjets_survive_matched = new TH1F("h_nbjets_survive_matched", ";n-bjets;Entries", 6, 0, 6);

    h2_mv2_ntracks = new TH2F("h2_mv2_ntracks", ";mv2c10 score;n-tracks", 40, 0, 1, 30, 0, 30);
    h2_dl1_ntracks = new TH2F("h2_dl1_ntracks", ";DL1 score;n-tracks", 60, -4, 2, 30, 0, 30);
}
//////////////////////////////////////////////////////////////////////////////
Bool_t ORStudy::Process(Long64_t entry)
{

    // calling "GetEntry" loads into memory the susyNt class objects for this event
    GetEntry(entry);
    SusyNtAna::clearObjects(); // clear the previous event's objects
    m_original_base_electrons.clear();
    m_original_base_jets.clear();

    // increment the chain entry (c.f. SusyNtuple/SusyNtAna.h)
    m_chainEntry++;

    // evt() provides pointer to the SusyNt::Event class object for this event
    int run_number = nt.evt()->run;
    int event_number = nt.evt()->eventNumber;

    if(dbg() || m_chainEntry%1000==0) {
        cout << "ORStudy::Process    **** Processing entry " << setw(6) << m_chainEntry
                << "  run " << run_number << "  event " << event_number << " **** " << endl;
    }

    // SusyNtAna::selectObject fills the baseline and signal objects
    // for the given AnalysisType
    // m_preX    = objects before any selection (as they are in susyNt)
    // m_baseX   = objects with the Analysis' baseline selection AND overlap removal applied
    // m_signalX = objects with the Analysis' signal selection applied (and baseline AND overlap removal)
    //SusyNtAna::selectObjects();

    // get the objects at baseline level prior to OR
    nttools().getPreObjects(&nt, NtSys::NOM, m_preElectrons, m_preMuons, m_preJets,
            m_preTaus, m_prePhotons);
    nttools().getBaselineObjects(m_preElectrons, m_preMuons, m_preJets, m_preTaus,
                m_prePhotons,
                m_baseElectrons, m_baseMuons, m_baseJets, m_baseTaus, m_basePhotons);

    m_original_base_electrons = m_baseElectrons;
    m_original_base_jets = m_baseJets;

    // get the truth level objects
    ORStudy::selectTruthObjects();

    ORStudy::performOverlap();

    int n_bjets = 0;
    int n_bjets_matched = 0;
    for(auto & x : m_baseJets) {
        bool pass_pt = x->Pt() > 20.0;
        bool pass_eta = (fabs(x->Eta()) < 2.5);
        bool pass_mv2 = true;
        if(m_release == "20") {
            pass_mv2 = (x->mv2c10 > JetSelector::mv2c10_77efficiency_rel20());
        }
        else if(m_release == "21") {
            pass_mv2 = (x->mv2c10 > JetSelector::mv2c10_77efficiency());
        }
        bool pass_dl1 = (x->dl1 > JetSelector::dl1_77efficiency());
        bool pass_tagger = true;
        if(m_tagger=="mv2") {
            pass_tagger = pass_mv2;
        }
        else if(m_tagger=="dl1") {
            pass_tagger = pass_dl1;
        }
        if(pass_tagger) {
            n_bjets++;
            if(x->truthLabel==5)
                n_bjets_matched++;
        }
    }
    h_nbjets_survive->Fill(n_bjets);
    h_nbjets_survive_matched->Fill(n_bjets_matched);

    // get signal objects after overlap
    nttools().getSignalObjects(m_baseElectrons, m_baseMuons, m_baseJets, m_baseTaus,
        m_basePhotons,
        m_signalElectrons, m_signalMuons, m_signalJets, m_signalTaus, m_signalPhotons);

    // get the MC weight using the inherited MCWeighter object
    // (c.f. SusyNtuple/MCWeighter.h)
    if(nt.evt()->isMC) {
        float lumi = 100000; // normalize the MC to 100 fb-1
        m_mc_weight = SusyNtAna::mcWeighter().getMCWeight(nt.evt(), lumi, NtSys::NOM);
    }
    else {
        m_mc_weight = 1.; // don't re-weight data
    }

    // check that the event passes the standard ATLAS event cleaning cuts
    if(!passEventCleaning(m_preMuons, m_baseMuons, m_baseJets)) return false;


    // for R20 vs R21 comparison
    if(m_signalElectrons.size() >=2) {
        h_l0pt->Fill(m_signalElectrons.at(0)->Pt(), m_mc_weight);
        h_l1pt->Fill(m_signalElectrons.at(1)->Pt(), m_mc_weight);
        float dphi = m_signalElectrons.at(0)->DeltaPhi(*m_signalElectrons.at(1));
        h_dphill->Fill(dphi, m_mc_weight);
        float drll = m_signalElectrons.at(0)->DeltaR(*m_signalElectrons.at(1));
        h_drll->Fill(drll, m_mc_weight);
    }


    return kTRUE;
}
//////////////////////////////////////////////////////////////////////////////
void ORStudy::selectTruthObjects()
{
    clearTruthObjects();

    // all particles
    nttools().getTruthParticles(&nt, m_truth_particles);

    // electrons
    nttools().getTruthParticles(&nt, m_truth_electrons, 11);

    // muons
    nttools().getTruthParticles(&nt, m_truth_muons, 13);

    // jets
    nttools().getTruthJets(&nt, m_truth_jets);

    // met
    nttools().getTruthMet(&nt, m_truth_met);

    // leptons
    for(auto x : m_truth_electrons)
        m_truth_leptons.push_back(x);

    for(auto x : m_truth_muons)
        m_truth_leptons.push_back(x);

    std::sort(m_truth_leptons.begin(), m_truth_leptons.end(), comparePt);    

    //cout << "n ele = " << m_truth_electrons.size() << ", n muo = " << m_truth_muons.size() << ", n lep = " << m_truth_leptons.size() << " :: ";
    //for(auto x : m_truth_leptons) cout << " " << x->Pt();
    //cout << endl;
}
//////////////////////////////////////////////////////////////////////////////
void ORStudy::performOverlap()
{
    nttools().overlapTool().m_e_overlap(m_baseMuons, m_baseElectrons);
    nttools().overlapTool().e_m_overlap(m_baseElectrons, m_baseMuons);
    //nttools().overlapTool().j_e_overlap(m_baseElectrons, m_baseJets, 0.2, true);
    //nttools().overlapTool().e_j_overlap(m_baseElectrons, m_baseJets, 0.4, false, true);
    ORStudy::j_e_overlap(m_baseElectrons, m_baseJets);
    ORStudy::e_j_overlap(m_baseElectrons, m_baseJets);
    nttools().overlapTool().j_m_overlap(m_baseJets, m_baseMuons, 0.2, true, true, false);
    nttools().overlapTool().m_j_overlap(m_baseMuons, m_baseJets, 0.4, false, true);
}
//////////////////////////////////////////////////////////////////////////////
void ORStudy::j_e_overlap(ElectronVector& electrons, JetVector& jets)
{
    float dRcone = 0.2;
    bool doBJetOR = true;

    if(electrons.size()==0 || jets.size()==0) return;
    for(int ij = jets.size()-1; ij >= 0; ij--) {
        const Jet* j = jets.at(ij);
        if(doBJetOR) {
            bool pass_pt = j->Pt() > 20.0;
            bool pass_eta = (fabs(j->Eta()) < 2.5);
            bool pass_mv2 = true;//(j->mv2c10 > JetSelector::mv2c10_77efficiency());
            if(m_release == "20") {
                pass_mv2 = (j->mv2c10 > JetSelector::mv2c10_85efficiency_rel20());
            }
            else if(m_release == "21") {
                pass_mv2 = (j->mv2c10 > JetSelector::mv2c10_85efficiency());
            }
            else {
                cout << "ORStudy::j_e_overlap    ERROR Invalid release encountered" << endl;
                exit(1);
            }
            bool pass_dl1 = (j->dl1 > JetSelector::dl1_85efficiency());
            bool pass_tagger = true;
            if(m_tagger=="mv2") {
                pass_tagger = pass_mv2;
            }
            else if(m_tagger=="dl1") {
                pass_tagger = pass_dl1;
            }
            else {
                cout << "ORStudy::j_e_overlap   ERROR Invalid tagger name encountered" << endl;
                exit(1);
            }
            bool pass_jvt = JetSelector::passJvt(j);

            if(pass_pt && pass_eta && pass_jvt) {
                h2_mv2_ntracks->Fill(j->mv2c10, j->nTracks);
                h2_dl1_ntracks->Fill(j->dl1, j->nTracks);

                n_kin++;
            }

            //look_at_jet(j);
            if(pass_pt && pass_eta && pass_tagger && pass_jvt) {
                n_kin_b++;
                //int min_ele_idx = -1;
                //float min_dr = 9999.;
                //for(size_t iel = 0; iel < electrons.size(); iel++) {
                //    const Electron* el = electrons.at(iel);
                //    float dr = j->DeltaR(*el);
                //    if(dr<min_dr) { min_dr = dr; min_ele_idx = iel; }
                //}

                //float pt_j = j->Pt();
                //float pt_e = electrons.at(min_ele_idx)->Pt();
                //float ratio = pt_e / pt_j;

                //if(ratio < 0.5)

                look_at_jet(j);

                    continue;
            }
        }
        for(int ie = 0; ie < (int)electrons.size(); ie++) {
            const Electron* el = electrons.at(ie);
            if(el->DeltaRy(*j) < dRcone) {
                jets.erase(jets.begin() + ij);
                break;
            }
        } // ie
    } // ij

}
//////////////////////////////////////////////////////////////////////////////
void ORStudy::e_j_overlap(ElectronVector& electrons, JetVector& jets)
{
    float dRcone = 0.4;
    bool boosted = false;
    bool applyJvt = true;

    if(electrons.size()==0 || jets.size()==0) return;
    for(int ie = electrons.size()-1; ie>=0; ie--) {
        const Electron* el = electrons.at(ie);
        for(int ij = 0; ij < (int)jets.size(); ij++) {
            const Jet* j = jets.at(ij);
            if(applyJvt && !JetSelector::passJvt(j)) continue;
            if(boosted) dRcone = nttools().overlapTool().getSlidingDRCone(el->Pt());
            if(el->DeltaRy(*j) < dRcone) {

                //look_at_jet_removing_electron(j, el);

                electrons.erase(electrons.begin() + ie);
                break;
            }
        } // ij
    } // ie
    
}
//////////////////////////////////////////////////////////////////////////////
void ORStudy::look_at_jet(const Jet* jet)
{
    if(!jet) return;

    int min_ele_idx = -1;
    float min_dr = 9999.;
    for(size_t iel = 0; iel < m_original_base_electrons.size(); iel++) {
        const Electron* el = m_original_base_electrons.at(iel);
        float dr = jet->DeltaR(*el);
        if(dr<min_dr) { min_dr = dr; min_ele_idx = iel; }
    }

    float pt_j = jet->Pt();
    float pt_e = m_original_base_electrons.at(min_ele_idx)->Pt();
    float ratio = pt_e / pt_j;

    h_min_drll_ele_jet->Fill(min_dr, m_mc_weight);
    h_pt_ratio_ele_jet->Fill(ratio, m_mc_weight);
    h2_pt_ratio_min_drll_ele_jet->Fill(min_dr, ratio);
    h2_mv2_min_drll_ele_jet->Fill(min_dr, jet->mv2c10);
    h2_dl1_min_drll_ele_jet->Fill(min_dr, jet->dl1);
    h2_pt_ratio_ntracks->Fill(ratio, jet->nTracks);


   // int min_truth_ele_idx = -1;
   // float min_te_dr = 9999.;
   // for(size_t ite = 0; ite < m_truth_electrons.size(); ite++) {
   //     const TruthParticle* te = m_truth_electrons.at(ite);
   //     float dr = jet->DeltaR(*te);
   //     if(dr < min_te_dr) { min_te_dr = dr; min_truth_ele_idx = ite; }
   // } // te

   // int type_closest_ele = m_truth_electrons.at(min_truth_ele_idx)->type;
   // int origin_closest_ele = m_truth_electrons.at(min_truth_ele_idx)->origin;
   // h_type_closest_ele->Fill(type_closest_ele);
   // h_origin_closest_ele->Fill(origin_closest_ele);
   // h2_type_closest_ele_min_drjet->Fill(min_te_dr, type_closest_ele);
   // h2_origin_closest_ele_min_drjet->Fill(min_te_dr, origin_closest_ele);

   // int min_truth_jet_idx = -1;
   // float min_tj_dr = 9999.;
   // for(size_t itj = 0; itj < m_truth_jets.size(); itj++) {
   //     const TruthJet* tj = m_truth_jets.at(itj);
   //     float dr = jet->DeltaR(*tj);
   //     if(dr < min_tj_dr) { min_tj_dr = dr; min_truth_jet_idx = itj; }
   // } // itj

   // int flavor = m_truth_jets.at(min_truth_jet_idx)->flavor;
   // h_flavor_closest_jet->Fill(flavor);
   // h2_flavor_closest_jet_min_drjet->Fill(min_tj_dr, flavor);
   // h2_flavor_closest_jet_ntracks->Fill(flavor, jet->nTracks);

   // bool closest_object_is_jet = true;
   // if(min_te_dr < min_tj_dr) closest_object_is_jet = false;
   // if(closest_object_is_jet) h_closest_truth_obj_jet->Fill(0.5, m_mc_weight);
   // else { h_closest_truth_obj_jet->Fill(1.5, m_mc_weight); }
    
}
//////////////////////////////////////////////////////////////////////////////
void ORStudy::clearTruthObjects()
{
    m_truth_particles.clear();
    m_truth_jets.clear();

    m_truth_electrons.clear();
    m_truth_muons.clear();
    m_truth_leptons.clear();
}
//////////////////////////////////////////////////////////////////////////////
bool ORStudy::passEventCleaning(const MuonVector& preMuons, const MuonVector& baseMuons,
            const JetVector& baseJets)
{
    int flags = nt.evt()->cutFlags[NtSys::NOM];

    if(!nttools().passGRL(flags))           return false;

    if(!nttools().passLarErr(flags))        return false;

    if(!nttools().passTileErr(flags))       return false;

    if(!nttools().passTTC(flags))           return false;

    if(!nttools().passSCTErr(flags))        return false;

    if(!nttools().passGoodVtx(flags))       return false;


    ///////////////////////////////////////////////////////
    // for bad muon, cosmic moun, and jet cleaning the
    // cuts depend on the baseline object defintion
    // (and in thec ase of the cosmic muon cut, it also
    // depends on the analysis' overlap removal
    // procedure) -- so we do not use the cutFlags but
    // rather use the objects that have passed the various
    // analysis selections to do the checks
    ///////////////////////////////////////////////////////
    if(!nttools().passBadMuon(preMuons))    return false;

    if(!nttools().passCosmicMuon(baseMuons)) return false;

    if(!nttools().passJetCleaning(baseJets)) return false;

    return true;
}
//////////////////////////////////////////////////////////////////////////////
void ORStudy::Terminate()
{
    // close SusyNtAna and print timers
    SusyNtAna::Terminate();

    write_histograms();

    cout << "----------------------" << endl;
    cout << " n_kin_b / n_kin = " << n_kin_b << " / " << n_kin << endl;

    return;
}
//////////////////////////////////////////////////////////////////////////////
void ORStudy::write_histograms()
{
    stringstream outname;
    outname << "or_study_rel" << m_release << "_" << m_tagger << ".root";
    TFile* out_file = new TFile(outname.str().c_str(), "RECREATE"); 
    if(out_file->IsZombie()) {
        cout << "ERROR Failed to open output file" << endl;
        return;
    }
    out_file->cd();
    h_l0pt->Write();
    h_l1pt->Write();
    h_dphill->Write();
    h_drll->Write();

    h_min_drll_ele_jet->Write();
    h_pt_ratio_ele_jet->Write();
    h2_pt_ratio_min_drll_ele_jet->Write();
    h2_mv2_min_drll_ele_jet->Write();
    h2_dl1_min_drll_ele_jet->Write();

    h_closest_truth_obj_jet->Write();
    h_type_closest_ele->Write();
    h_origin_closest_ele->Write();
    h_flavor_closest_jet->Write();

    h2_type_closest_ele_min_drjet->Write();
    h2_origin_closest_ele_min_drjet->Write();
    h2_flavor_closest_jet_min_drjet->Write();
    h2_pt_ratio_ntracks->Write();
    h2_flavor_closest_jet_ntracks->Write();

    h_nbjets_survive->Write();
    h_nbjets_survive_matched->Write();

    h2_mv2_ntracks->Write();
    h2_dl1_ntracks->Write();

    out_file->Write();
    out_file->Close();

}
//////////////////////////////////////////////////////////////////////////////
