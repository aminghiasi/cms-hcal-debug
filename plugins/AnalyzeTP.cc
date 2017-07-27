// -*- C++ -*-
//
// Package:    HcalDebug
// Class:      AnalyzeTP
// 
/**\class AnalyzeTP AnalyzeTP.cc HcalDebug/CompareChans/src/AnalyzeTP.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  matthias wolf
//         Created:  Fri Nov 27 11:21:58 CET 2015
// $Id$
//
//


// system include files
#include <memory>
#include <unordered_map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CalibFormats/CaloTPG/interface/CaloTPGTranscoder.h"
#include "CalibFormats/CaloTPG/interface/CaloTPGRecord.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"

#include "CondFormats/DataRecord/interface/HcalChannelQualityRcd.h"
#include "CondFormats/DataRecord/interface/L1CaloGeometryRecord.h"
#include "CondFormats/HcalObjects/interface/HcalChannelQuality.h"
#include "CondFormats/L1TObjects/interface/L1CaloGeometry.h"
#include "CondFormats/L1TObjects/interface/L1RCTParameters.h"
#include "CondFormats/DataRecord/interface/L1RCTParametersRcd.h"
#include "CondFormats/L1TObjects/interface/L1CaloHcalScale.h"
#include "CondFormats/DataRecord/interface/L1CaloHcalScaleRcd.h"

#include "DataFormats/Common/interface/SortedCollection.h"
#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalTriggerPrimitiveDigi.h"
#include "DataFormats/HcalDetId/interface/HcalTrigTowerDetId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/TriggerResults.h"


#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/HcalTowerAlgo/interface/HcalGeometry.h"
#include "Geometry/HcalTowerAlgo/interface/HcalTrigTowerGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TTree.h"
//
// class declaration
//

class AnalyzeTP : public edm::EDAnalyzer {
   public:
      explicit AnalyzeTP(const edm::ParameterSet&);
      ~AnalyzeTP();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void analyze(const edm::Event&, const edm::EventSetup&);

      // ----------member data ---------------------------
      edm::InputTag digis_;
      edm::InputTag offlinePrimaryVertices_;
      edm::EDGetTokenT<edm::TriggerResults>TrigTagTok_;
      edm::InputTag trigTagSrc_;


      double threshold_;

      int event_, ls_, bx_;
      int npv_;
      TTree *match_;
      int m_ieta_;
      int m_iphi_;
      double old_et_;
      double new_et_;
      int new_count_;
      int old_fg_;
      int new_fg_;

      TTree *tps_;

      int tp_ieta_;
      int tp_iphi_;
      int tp_depth_;
      int tp_version_;
      int tp_soi_;
      int tp_fg_;
      double tp_et_;
      bool HLT_Random_v2;
      bool HLT_ZeroBias_v5;

      TTree *ev_;
      double ev_tp_v0_et_;
      double ev_tp_v1_et_, sum_tp_et_hf, sum_tp_et_he, sum_tp_et_hb, sum_tp_et_all;
      int no_tps_g_0_0 = 0, no_tps_g_0_5 = 0, no_tps_g_1 = 0, no_tps_g_2 = 0, no_tps_g_5 = 0, no_tps_g_10 = 0;
      int no_hf_tps_g_0_0 = 0, no_hf_tps_g_0_5 = 0, no_hf_tps_g_1 = 0, no_hf_tps_g_2 = 0, no_hf_tps_g_5 = 0, no_hf_tps_g_10 = 0;
      int no_he_tps_g_0_0 = 0, no_he_tps_g_0_5 = 0, no_he_tps_g_1 = 0, no_he_tps_g_2 = 0, no_he_tps_g_5 = 0, no_he_tps_g_10 = 0;
      int no_hb_tps_g_0_0 = 0, no_hb_tps_g_0_5 = 0, no_hb_tps_g_1 = 0, no_hb_tps_g_2 = 0, no_hb_tps_g_5 = 0, no_hb_tps_g_10 = 0;



};

AnalyzeTP::AnalyzeTP(const edm::ParameterSet& config) :
   edm::EDAnalyzer(),
   digis_(config.getParameter<edm::InputTag>("triggerPrimitives")),
   offlinePrimaryVertices_(config.getParameter<edm::InputTag>("offlinePrimaryVertices")),
   threshold_(config.getUntrackedParameter<double>("threshold", 0.5))
{
   edm::Service<TFileService> fs;

   consumes<HcalTrigPrimDigiCollection>(digis_);
   consumes<std::vector<reco::Vertex>> (offlinePrimaryVertices_);
   trigTagSrc_ = config.getParameter<edm::InputTag> ("trigTagSrc");
   TrigTagTok_ = consumes<edm::TriggerResults>(trigTagSrc_);


   tps_ = fs->make<TTree>("tps", "Trigger primitives");
   tps_->Branch("event", &event_);
   tps_->Branch("BX", &bx_);
   tps_->Branch("ls", &ls_);
   tps_->Branch("ieta", &tp_ieta_);
   tps_->Branch("iphi", &tp_iphi_);
   tps_->Branch("depth", &tp_depth_);
   tps_->Branch("version", &tp_version_);
   tps_->Branch("soi", &tp_soi_);
   tps_->Branch("et", &tp_et_);
   tps_->Branch("fg", &tp_fg_);
   tps_->Branch("npv",&npv_);
   tps_->Branch("HLT_Random_v2"  ,&HLT_Random_v2);
   tps_->Branch("HLT_ZeroBias_v5",&HLT_ZeroBias_v5);

   ev_ = fs->make<TTree>("evs", "Event quantities");
   ev_->Branch("ls", &ls_);
   ev_->Branch("BX", &bx_);
   ev_->Branch("event", &event_);
   ev_->Branch("tp_v0_et", &ev_tp_v0_et_);
   ev_->Branch("tp_v1_et"   , &ev_tp_v1_et_);
   ev_->Branch("sum_tp_et_hf", &sum_tp_et_hf);
   ev_->Branch("sum_tp_et_he", &sum_tp_et_he);
   ev_->Branch("sum_tp_et_hb", &sum_tp_et_hb);
   ev_->Branch("sum_tp_et_all", &sum_tp_et_all);
   ev_->Branch("npv",&npv_);
   ev_->Branch("HLT_Random_v2"             ,&HLT_Random_v2);
   ev_->Branch("HLT_ZeroBias_v5"           ,&HLT_ZeroBias_v5);

   ev_->Branch("no_tps_g_0_0"           ,&no_tps_g_0_5);
   ev_->Branch("no_tps_g_0_5"           ,&no_tps_g_0_5);
   ev_->Branch("no_tps_g_1"           ,&no_tps_g_1);
   ev_->Branch("no_tps_g_2"           ,&no_tps_g_2);
   ev_->Branch("no_tps_g_5"           ,&no_tps_g_5);
   ev_->Branch("no_tps_g_10"           ,&no_tps_g_10);

   ev_->Branch("no_hf_tps_g_0_0"         ,&no_hf_tps_g_0_5);
   ev_->Branch("no_hf_tps_g_0_5"         ,&no_hf_tps_g_0_5);
   ev_->Branch("no_hf_tps_g_1"           ,&no_hf_tps_g_1);
   ev_->Branch("no_hf_tps_g_2"           ,&no_hf_tps_g_2);
   ev_->Branch("no_hf_tps_g_5"           ,&no_hf_tps_g_5);
   ev_->Branch("no_hf_tps_g_10"          ,&no_hf_tps_g_10);

   ev_->Branch("no_he_tps_g_0_0"         ,&no_he_tps_g_0_5);
   ev_->Branch("no_he_tps_g_0_5"         ,&no_he_tps_g_0_5);
   ev_->Branch("no_he_tps_g_1"           ,&no_he_tps_g_1);
   ev_->Branch("no_he_tps_g_2"           ,&no_he_tps_g_2);
   ev_->Branch("no_he_tps_g_5"           ,&no_he_tps_g_5);
   ev_->Branch("no_he_tps_g_10"          ,&no_he_tps_g_10);

   ev_->Branch("no_hb_tps_g_0_0"         ,&no_hb_tps_g_0_5);
   ev_->Branch("no_hb_tps_g_0_5"         ,&no_hb_tps_g_0_5);
   ev_->Branch("no_hb_tps_g_1"           ,&no_hb_tps_g_1);
   ev_->Branch("no_hb_tps_g_2"           ,&no_hb_tps_g_2);
   ev_->Branch("no_hb_tps_g_5"           ,&no_hb_tps_g_5);
   ev_->Branch("no_hb_tps_g_10"          ,&no_hb_tps_g_10);



   match_ = fs->make<TTree>("ms", "TP matches");
   match_->Branch("event", &event_);
   match_->Branch("BX", &bx_);
   match_->Branch("ls", &ls_);
   match_->Branch("ieta", &m_ieta_);
   match_->Branch("iphi", &m_iphi_);
   match_->Branch("et1x1", &new_et_);
   match_->Branch("et2x3", &old_et_);
   match_->Branch("n1x1", &new_count_);
   match_->Branch("fg1x1", &new_fg_);
   match_->Branch("fg2x3", &old_fg_);
   match_->Branch("npv",&npv_);

}

AnalyzeTP::~AnalyzeTP() {}

void
AnalyzeTP::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
   using namespace edm;
   using namespace std;

   event_ = event.id().event();
   ls_ = event.id().luminosityBlock();
   bx_ = event.bunchCrossing();

   Handle<HcalTrigPrimDigiCollection> digis;
   if (!event.getByLabel(digis_, digis)) {
      LogError("AnalyzeTP") <<
         "Can't find hcal trigger primitive digi collection with tag '" <<
         digis_ << "'" << std::endl;
      return;
   }

   int npv=0;

        edm::Handle < vector < reco::Vertex > > pvHandle;
   try {
     event.getByLabel( offlinePrimaryVertices_, pvHandle );
   } catch ( cms::Exception & e ) {
   }

   npv = pvHandle->size();
   npv_ = npv;


   edm::Handle<edm::TriggerResults> trigResults;
   event.getByToken(TrigTagTok_,trigResults);

   const edm::TriggerNames& trigNames = event.triggerNames(*trigResults);
   size_t pathIndex = trigNames.triggerIndex("HLT_Random_v2");
   HLT_Random_v2 = false;
   if ( pathIndex < trigResults->size() && trigResults->accept(pathIndex))
      HLT_Random_v2 = true;

   size_t pathIndex2 = trigNames.triggerIndex("HLT_ZeroBias_v5");
   HLT_ZeroBias_v5 = false;
   if ( pathIndex2 < trigResults->size() && trigResults->accept(pathIndex2))
      HLT_ZeroBias_v5 = true;



   ESHandle<CaloTPGTranscoder> decoder;

   setup.get<CaloTPGRecord>().get(decoder);

   std::unordered_map<int, std::unordered_map<int, double>> old_ets;
   std::unordered_map<int, std::unordered_map<int, double>> new_ets;
   std::unordered_map<int, std::unordered_map<int, int>> new_counts;

   ev_tp_v0_et_ = 0.;
   ev_tp_v1_et_ = 0.;

   sum_tp_et_hf = 0; sum_tp_et_he = 0; sum_tp_et_hb = 0; sum_tp_et_all = 0;

   no_tps_g_0_0 = 0; no_tps_g_0_5 = 0; no_tps_g_1 = 0; no_tps_g_2 = 0; no_tps_g_5 = 0; no_tps_g_10 = 0;
   no_hf_tps_g_0_0 = 0; no_hf_tps_g_0_5 = 0; no_hf_tps_g_1 = 0; no_hf_tps_g_2 = 0; no_hf_tps_g_5 = 0; no_hf_tps_g_10 = 0;
   no_he_tps_g_0_0 = 0; no_he_tps_g_0_5 = 0; no_he_tps_g_1 = 0; no_he_tps_g_2 = 0; no_he_tps_g_5 = 0; no_he_tps_g_10 = 0;
   no_hb_tps_g_0_0 = 0; no_hb_tps_g_0_5 = 0; no_hb_tps_g_1 = 0; no_hb_tps_g_2 = 0; no_hb_tps_g_5 = 0; no_hb_tps_g_10 = 0;


   ESHandle<HcalTrigTowerGeometry> tpd_geo;
   setup.get<CaloGeometryRecord>().get(tpd_geo);

   std::map<HcalTrigTowerDetId, HcalTriggerPrimitiveDigi> ttids;
   for (const auto& digi: *digis)
   {
      if (digi.id().version() == 1)
      {
         ttids[digi.id()] = digi;
      }

   }

   for (const auto& digi: *digis) {
      HcalTrigTowerDetId id = digi.id();

      if (id.version() == 1 and abs(id.ieta()) >= 40 and id.iphi() % 4 == 1)
         continue;

      tp_ieta_ = id.ieta();
      tp_iphi_ = id.iphi();
      tp_depth_ = id.depth();
      tp_version_ = id.version();
      tp_soi_ = digi.SOI_compressedEt();
      tp_et_ = decoder->hcaletValue(id, digi.t0());
      tp_fg_ = digi.SOI_fineGrain();


      if ( abs(tp_ieta_) == 28 || abs(tp_ieta_) == 16 ) if ( tp_et_ < 2 ) std::cout << "tp_et: " << tp_et_ << std::endl;
      if (tp_version_ == 0 && abs(tp_ieta_) >= 29) {
         ev_tp_v0_et_ += tp_et_;
      } else if (tp_version_ == 1) {
         ev_tp_v1_et_ += tp_et_;
      }
      if (  tp_et_ > 0.5 && HLT_ZeroBias_v5 == false ) 
      {
         if ( abs(tp_ieta_) >= 29 && tp_version_ == 1   ) sum_tp_et_hf += tp_et_;
         if ( abs(tp_ieta_) <  29 && abs(tp_ieta_) > 16 ) sum_tp_et_he += tp_et_;
         if ( abs(tp_ieta_) <=  16                      ) sum_tp_et_hb += tp_et_;

      }// if




      if ( tp_et_ > 0.0 && HLT_ZeroBias_v5 == false ) 
      {
         if ( abs(tp_ieta_) >= 29 && tp_version_ == 1   ) no_hf_tps_g_0_0++;
         if ( abs(tp_ieta_) <  29 && abs(tp_ieta_) > 16 ) no_he_tps_g_0_0++;
         if ( abs(tp_ieta_) <=  16                      ) no_hb_tps_g_0_0++;

         no_tps_g_0_0++;
      }

      if ( tp_et_ > 0.5 && HLT_ZeroBias_v5 == false ) 
      {
         if ( abs(tp_ieta_) >= 29 && tp_version_ == 1   ) no_hf_tps_g_0_5++;
         if ( abs(tp_ieta_) <  29 && abs(tp_ieta_) > 16 ) no_he_tps_g_0_5++;
         if ( abs(tp_ieta_) <=  16                      ) no_hb_tps_g_0_5++;
         no_tps_g_0_5++;
      }

      if ( tp_et_ > 1 && HLT_ZeroBias_v5 == false ) 
      {
         if ( abs(tp_ieta_) >= 29 && tp_version_ == 1   ) no_hf_tps_g_1++;
         if ( abs(tp_ieta_) <  29 && abs(tp_ieta_) > 16 ) no_he_tps_g_1++;
         if ( abs(tp_ieta_) <=  16                      ) no_hb_tps_g_1++;
         no_tps_g_1++;
      }

      if ( tp_et_ > 2 && HLT_ZeroBias_v5 == false ) 
      {
         if ( abs(tp_ieta_) >= 29 && tp_version_ == 1   ) no_hf_tps_g_2++;
         if ( abs(tp_ieta_) <  29 && abs(tp_ieta_) > 16 ) no_he_tps_g_2++;
         if ( abs(tp_ieta_) <=  16                      ) no_hb_tps_g_2++;
         no_tps_g_2++;
      }

      if ( tp_et_ > 5 && HLT_ZeroBias_v5 == false )  
      {
         if ( abs(tp_ieta_) >= 29 && tp_version_ == 1   ) no_hf_tps_g_5++;
         if ( abs(tp_ieta_) <  29 && abs(tp_ieta_) > 16 ) no_he_tps_g_5++;
         if ( abs(tp_ieta_) <=  16                      ) no_hb_tps_g_5++;
         no_tps_g_5++;
      }

      if ( tp_et_ > 10 && HLT_ZeroBias_v5 == false )  
      {
         if ( abs(tp_ieta_) >= 29 && tp_version_ == 1   ) no_hf_tps_g_10++;
         if ( abs(tp_ieta_) <  29 && abs(tp_ieta_) > 16 ) no_he_tps_g_10++;
         if ( abs(tp_ieta_) <=  16                      ) no_hb_tps_g_10++;
         no_tps_g_10++;
      }
      if (tp_et_ < threshold_)
         continue;

      tps_->Fill();



      if (abs(tp_ieta_) >= 29 and tp_version_ == 0) {
         std::set<HcalTrigTowerDetId> matches;
         for (const auto& detid: tpd_geo->detIds(id)) {
            for (const auto& ttid: tpd_geo->towerIds(detid)) {
               if (ttid.version() == 1)
                  matches.insert(ttid);
            }
         }

         m_ieta_ = tp_ieta_;
         m_iphi_ = tp_iphi_;
         new_et_ = 0;
         new_count_ = 0;
         old_et_ = tp_et_;
         old_fg_ = tp_fg_;
         new_fg_ = 0;
         for (const auto& m: matches) {
            if (m.version() == 1 and abs(m.ieta()) >= 40 and m.iphi() % 4 == 1)
               continue;

            new_et_ += decoder->hcaletValue(m, ttids[m].t0());
            ++new_count_;
            new_fg_ = new_fg_ || ttids[m].t0().fineGrain();
         }
         match_->Fill();
      }
   }

   ev_->Fill();
}

void
AnalyzeTP::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(AnalyzeTP);
