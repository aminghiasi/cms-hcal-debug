// -*- C++ -*-
//
// Package:    HcalCompareLegacyChains
// Class:      HcalCompareLegacyChains
// 
/**\class HcalCompareLegacyChains HcalCompareLegacyChains.cc HcalDebug/CompareChans/src/HcalCompareLegacyChains.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  matthias wolf
//         Created:  Fri Aug 26 11:37:21 CDT 2013
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CalibFormats/CaloTPG/interface/CaloTPGTranscoder.h"
#include "CalibFormats/CaloTPG/interface/CaloTPGRecord.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"

#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalTriggerPrimitiveDigi.h"
#include "DataFormats/HcalDetId/interface/HcalTrigTowerDetId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "DataFormats/HcalRecHit/interface/HFRecHit.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/HcalTowerAlgo/interface/HcalGeometry.h"
#include "Geometry/HcalTowerAlgo/interface/HcalTrigTowerGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TTree.h"
//
// class declaration
//

class HcalCompareLegacyChains : public edm::EDAnalyzer {
   public:
      explicit HcalCompareLegacyChains(const edm::ParameterSet&);
      ~HcalCompareLegacyChains();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void analyze(const edm::Event&, const edm::EventSetup&);

      double get_cosh(const HcalDetId&);

      // ----------member data ---------------------------
      bool first_;
      edm::EDGetTokenT<edm::TriggerResults>TrigTagTok_;
      edm::InputTag trigTagSrc_;
      edm::InputTag offlinePrimaryVertices_;

      std::vector<edm::InputTag> frames_;
      edm::InputTag digis_;
      std::vector<edm::InputTag> rechits_;

      edm::ESHandle<CaloGeometry> gen_geo_;

      bool swap_iphi_;

      TH2D *df_multiplicity_;
      TH2D *tp_multiplicity_;

      TTree *tps_;
      TTree *tpsplit_;
      TTree *events_;
      TTree *matches_;

      bool HLT_ZeroBias_v5_;

      double tp_energy_;
      int tp_ieta_;
      int tp_iphi_;
      int tp_soi_;

      double tpsplit_energy_;
      double tpsplit_oot_;
      int tpsplit_ieta_;
      int tpsplit_iphi_;
      int tpsplit_depth_;
      double tpsplit_ettot_;
      double tpsplit_rise_avg_;
      double tpsplit_rise_rms_;
      double tpsplit_fall_avg_;
      double tpsplit_fall_rms_;

      double ev_rh_energy0_;
      double ev_rh_energy2_;
      double ev_rh_energy3_;
      double ev_tp_energy_;
      int ev_rh_unmatched_;
      int ev_tp_unmatched_;

      double mt_rh_energy0_;
      double mt_rh_energy2_;
      double mt_rh_energy3_;
      double mt_tp_energy_;

      int mt_ieta_;
      int mt_iphi_;
      int mt_version_;
      int mt_tp_soi_;
};

HcalCompareLegacyChains::HcalCompareLegacyChains(const edm::ParameterSet& config) :
   edm::EDAnalyzer(),
   first_(true),
   offlinePrimaryVertices_(config.getParameter<edm::InputTag>("offlinePrimaryVertices")),
   frames_(config.getParameter<std::vector<edm::InputTag>>("dataFrames")),
   digis_(config.getParameter<edm::InputTag>("triggerPrimitives")),
   rechits_(config.getParameter<std::vector<edm::InputTag>>("recHits")),
   swap_iphi_(config.getParameter<bool>("swapIphi"))
{
   consumes<HcalTrigPrimDigiCollection>(digis_);
   consumes<HBHEDigiCollection>(frames_[0]);
   consumes<HFDigiCollection>(frames_[1]);
   consumes<edm::SortedCollection<HBHERecHit>>(rechits_[0]);
   consumes<edm::SortedCollection<HFRecHit>>(rechits_[1]);
   consumes<std::vector<reco::Vertex>> (offlinePrimaryVertices_);


   trigTagSrc_ = config.getParameter<edm::InputTag> ("trigTagSrc");
   TrigTagTok_ = consumes<edm::TriggerResults>(trigTagSrc_);

   edm::Service<TFileService> fs;

   df_multiplicity_ = fs->make<TH2D>("df_multiplicity", "DataFrame multiplicity;ieta;iphi", 65, -32.5, 32.5, 72, 0.5, 72.5);
   tp_multiplicity_ = fs->make<TH2D>("tp_multiplicity", "TrigPrim multiplicity;ieta;iphi", 65, -32.5, 32.5, 72, 0.5, 72.5);

   tps_ = fs->make<TTree>("tps", "Trigger primitives");
   tps_->Branch("et", &tp_energy_);
   tps_->Branch("ieta", &tp_ieta_);
   tps_->Branch("iphi", &tp_iphi_);
   tps_->Branch("soi", &tp_soi_);

   tpsplit_ = fs->make<TTree>("tpsplit", "Trigger primitives");
   tpsplit_->Branch("et", &tpsplit_energy_);
   tpsplit_->Branch("oot", &tpsplit_oot_);
   tpsplit_->Branch("ieta", &tpsplit_ieta_);
   tpsplit_->Branch("iphi", &tpsplit_iphi_);
   tpsplit_->Branch("depth", &tpsplit_depth_);
   tpsplit_->Branch("etsum", &tpsplit_ettot_);
   tpsplit_->Branch("rise_avg", &tpsplit_rise_avg_);
   tpsplit_->Branch("rise_rms", &tpsplit_rise_rms_);
   tpsplit_->Branch("fall_avg", &tpsplit_fall_avg_);
   tpsplit_->Branch("fall_rms", &tpsplit_fall_rms_);

   events_ = fs->make<TTree>("events", "Event quantities");
   events_->Branch("RH_energyM0", &ev_rh_energy0_);
   events_->Branch("RH_energyM2", &ev_rh_energy2_);
   events_->Branch("RH_energyM3", &ev_rh_energy3_);
   events_->Branch("TP_energy", &ev_tp_energy_);
   events_->Branch("RH_unmatched", &ev_rh_unmatched_);
   events_->Branch("TP_unmatched", &ev_tp_unmatched_);
   events_->Branch("HLT_ZeroBias_v5", &HLT_ZeroBias_v5_);



   matches_ = fs->make<TTree>("matches", "Matched RH and TP");
   matches_->Branch("RH_energyM0", &mt_rh_energy0_);
   matches_->Branch("RH_energyM2", &mt_rh_energy2_);
   matches_->Branch("RH_energyM3", &mt_rh_energy3_);
   matches_->Branch("TP_energy", &mt_tp_energy_);
   matches_->Branch("ieta", &mt_ieta_);
   matches_->Branch("iphi", &mt_iphi_);
   matches_->Branch("tp_version", &mt_version_);
   matches_->Branch("tp_soi", &mt_tp_soi_);
   matches_->Branch("HLT_ZeroBias_v5", &HLT_ZeroBias_v5_);

}

HcalCompareLegacyChains::~HcalCompareLegacyChains() {}

double
HcalCompareLegacyChains::get_cosh(const HcalDetId& id)
{
   const auto *sub_geo = dynamic_cast<const HcalGeometry*>(gen_geo_->getSubdetectorGeometry(id));
   auto eta = sub_geo->getPosition(id).eta();
   return cosh(eta);
}

void
HcalCompareLegacyChains::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
   using namespace edm;

   edm::ESHandle<HcalTrigTowerGeometry> tpd_geo_h;
   setup.get<CaloGeometryRecord>().get(tpd_geo_h);
   edm::ESHandle<HcalDbService> conditions;
   setup.get<HcalDbRecord>().get(conditions);
   const HcalTrigTowerGeometry& tpd_geo = *tpd_geo_h;

   // ==========
   // Dataframes
   // ==========

   HLT_ZeroBias_v5_ = false;
   edm::Handle<edm::TriggerResults> trigResults;
   event.getByToken(TrigTagTok_,trigResults);
   const edm::TriggerNames& trigNames = event.triggerNames(*trigResults);
   size_t pathIndex = trigNames.triggerIndex("HLT_ZeroBias_v5");
   if ( pathIndex < trigResults->size() && trigResults->accept(pathIndex))
      HLT_ZeroBias_v5_ = true;

   Handle<HBHEDigiCollection> frames;
   Handle<HFDigiCollection> hfframes;
   if (frames_.size() == 2) {
      if (first_ && event.getByLabel(frames_[0], frames) && event.getByLabel(frames_[1], hfframes)) {
         std::set<HcalTrigTowerDetId> ids;

         for (const auto& frame: *(frames.product())) {
            auto mapped = tpd_geo_h->towerIds(frame.id());

            for (const auto& id: mapped) {
               df_multiplicity_->Fill(id.ieta(), id.iphi());
               ids.insert(id);
            }
         }

         for (const auto& frame: *(hfframes.product())) {
            auto mapped = tpd_geo_h->towerIds(frame.id());

            for (const auto& id: mapped) {
               df_multiplicity_->Fill(id.ieta(), id.iphi());
               ids.insert(id);
            }
         }

         for (const auto& id: ids) {
            tp_multiplicity_->Fill(id.ieta(), id.iphi());
         }

         first_ = false;
      }
   }


   // ==============
   // Matching stuff
   // ==============

   ev_rh_energy0_ = 0.;
   ev_rh_energy2_ = 0.;
   ev_rh_energy3_ = 0.;
   ev_rh_unmatched_ = 0.;
   ev_tp_energy_ = 0.;
   ev_tp_unmatched_ = 0.;

   std::map<HcalTrigTowerDetId, std::vector<HBHERecHit>> rhits;
   std::map<HcalTrigTowerDetId, std::vector<HFRecHit>> fhits;
   std::map<HcalTrigTowerDetId, std::vector<HcalTriggerPrimitiveDigi>> tpdigis;

   Handle<HcalTrigPrimDigiCollection> digis;
   if (!event.getByLabel(digis_, digis)) {
      LogError("HcalTrigPrimDigiCleaner") <<
         "Can't find hcal trigger primitive digi collection with tag '" <<
         digis_ << "'" << std::endl;
      return;
   }

   edm::Handle< edm::SortedCollection<HBHERecHit> > hits;
   if (!event.getByLabel(rechits_[0], hits)) {
      edm::LogError("HcalCompareLegacyChains") <<
         "Can't find rec hit collection with tag '" << rechits_[0] << "'" << std::endl;
      /* return; */
   }

   edm::Handle< edm::SortedCollection<HFRecHit> > hfhits;
   if (!event.getByLabel(rechits_[1], hfhits)) {
      edm::LogError("HcalCompareLegacyChains") <<
         "Can't find rec hit collection with tag '" << rechits_[1] << "'" << std::endl;
      /* return; */
   }

   setup.get<CaloGeometryRecord>().get(gen_geo_);

   if (hits.isValid()) {
      for (auto& hit: *(hits.product())) {
         HcalDetId id(hit.id());
         ev_rh_energy0_ += hit.eraw() / get_cosh(id);
         ev_rh_energy2_ += hit.energy() / get_cosh(id);
         ev_rh_energy3_ += hit.eaux() / get_cosh(id);

         auto tower_ids = tpd_geo.towerIds(id);
         for (auto& tower_id: tower_ids) {
            tower_id = HcalTrigTowerDetId(tower_id.ieta(), tower_id.iphi(), 1);
            rhits[tower_id].push_back(hit);
         }
      }
   }

   if (hfhits.isValid()) {
      for (auto& hit: *(hfhits.product())) {
         HcalDetId id(hit.id());
         ev_rh_energy0_ += hit.energy() / get_cosh(id);
         ev_rh_energy2_ += hit.energy() / get_cosh(id);
         ev_rh_energy3_ += hit.energy() / get_cosh(id);

         auto tower_ids = tpd_geo.towerIds(id);
         for (auto& tower_id: tower_ids) {
            tower_id = HcalTrigTowerDetId(tower_id.ieta(), tower_id.iphi(), 1, tower_id.version());
            fhits[tower_id].push_back(hit);
         }
      }
   }

   // ESHandle<L1CaloHcalScale> hcal_scale;
   // setup.get<L1CaloHcalScaleRcd>().get(hcal_scale);
   // const L1CaloHcalScale* h = hcal_scale.product();

   ESHandle<CaloTPGTranscoder> decoder;
   setup.get<CaloTPGRecord>().get(decoder);

   for (const auto& digi: *digis) {
      HcalTrigTowerDetId id = digi.id();
      id = HcalTrigTowerDetId(id.ieta(), id.iphi(), 1, id.version());
      ev_tp_energy_ += decoder->hcaletValue(id, digi.t0());
      tpdigis[id].push_back(digi);

      tp_energy_ = decoder->hcaletValue(id, digi.t0());
      tp_ieta_ = id.ieta();
      tp_iphi_ = id.iphi();
      tp_soi_ = digi.SOI_compressedEt();

      tps_->Fill();
   }

   for (const auto& pair: tpdigis) {
      auto id = pair.first;

      auto new_id(id);
      if (swap_iphi_ and id.version() == 1 and id.ieta() > 28 and id.ieta() < 40) {
         if (id.iphi() % 4 == 1)
            new_id = HcalTrigTowerDetId(id.ieta(), (id.iphi() + 70) % 72, id.depth(), id.version());
         else
            new_id = HcalTrigTowerDetId(id.ieta(), (id.iphi() + 2) % 72 , id.depth(), id.version());
      }

      auto rh = rhits.find(new_id);
      auto fh = fhits.find(new_id);

      if (rh != rhits.end() and fh != fhits.end()) {
         assert(0);
      }

      mt_ieta_ = new_id.ieta();
      mt_iphi_ = new_id.iphi();
      mt_version_ = new_id.version();
      mt_tp_energy_ = 0;
      mt_tp_soi_ = 0;

      for (const auto& tp: pair.second) {
         mt_tp_energy_ += decoder->hcaletValue(new_id, tp.t0());
         mt_tp_soi_ = tp.SOI_compressedEt();
      }
      mt_rh_energy0_ = 0.;
      mt_rh_energy2_ = 0.;
      mt_rh_energy3_ = 0.;

      if (rh != rhits.end()) {
         for (const auto& hit: rh->second) {
            HcalDetId id(hit.id());
            auto tower_ids = tpd_geo.towerIds(id);
            auto count = std::count_if(std::begin(tower_ids), std::end(tower_ids),
                  [&](const auto& t) { return t.version() == new_id.version(); });
            mt_rh_energy0_ += hit.eraw() / get_cosh(id) / count;
            mt_rh_energy2_ += hit.energy() / get_cosh(id) / count;
            mt_rh_energy3_ += hit.eaux() / get_cosh(id) / count;
         }
         matches_->Fill();
         rhits.erase(rh);
      } else if (fh != fhits.end()) {
         for (const auto& hit: fh->second) {
            HcalDetId id(hit.id());
            auto tower_ids = tpd_geo.towerIds(id);
            auto count = std::count_if(std::begin(tower_ids), std::end(tower_ids),
                  [&](const auto& t) { return t.version() == new_id.version(); });
            mt_rh_energy0_ += hit.energy() / get_cosh(id) / count;
            mt_rh_energy2_ += hit.energy() / get_cosh(id) / count;
            mt_rh_energy3_ += hit.energy() / get_cosh(id) / count;
         }
         matches_->Fill();
         fhits.erase(fh);
      } else {
         ++ev_tp_unmatched_;
      }
   }

   for (const auto& pair: rhits) {
      ev_rh_unmatched_ += pair.second.size();
   }

   events_->Fill();
}

void
HcalCompareLegacyChains::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HcalCompareLegacyChains);
