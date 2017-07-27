import FWCore.ParameterSet.Config as cms

from Configuration.AlCa.GlobalTag import GlobalTag
from Configuration.StandardSequences.Eras import eras

process = cms.Process('PLOT', eras.Run2_2017)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag = GlobalTag(process.GlobalTag, '92X_dataRun2_HLT_v3', '')

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring())
process.source.fileNames.extend([

'/store/data/Run2017B/ZeroBias/RAW/v1/000/297/469/00000/02E754DD-2459-E711-BDA8-02163E01A583.root'
])

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.RawToDigi_Data_cff")
process.load("SimCalorimetry.Configuration.hcalDigiSequence_cff")
process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
process.load("EventFilter.L1TRawToDigi.caloStage2Digis_cfi")
process.load("EventFilter.L1TXRawToDigi.caloLayer1Stage2Digis_cfi")

process.TFileService = cms.Service("TFileService",
                                   closeFileFasO=cms.untracked.bool(True),
                                   fileName=cms.string('analyze_tps.root'))

process.emulTP = process.simHcalTriggerPrimitiveDigis.clone()
process.emulTP.upgradeHF = cms.bool(True)
process.emulTP.upgradeHE = cms.bool(True)
process.emulTP.inputLabel = cms.VInputTag("hcalDigis", "hcalDigis")
process.emulTP.inputUpgradeLabel = cms.VInputTag("hcalDigis", "hcalDigis")

process.emulTP.numberOfSamples = cms.int32(4)
process.emulTP.numberOfPresamples = cms.int32(2)
process.emulTP.numberOfSamplesHF = cms.int32(2)
process.emulTP.numberOfPresamplesHF = cms.int32(1)

process.GlobalTag.toGet = cms.VPSet(
  cms.PSet(record = cms.string("HcalElectronicsMapRcd"),
           tag = cms.string("HcalElectronicsMap_2017plan1_v3.0_data"),
           connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS")
          )
)

# process.hcalDigis.InputLabel = cms.InputTag("source")
process.analyzeRAW = cms.EDAnalyzer("AnalyzeTP",
                                    triggerPrimitives=cms.InputTag("hcalDigis", "", ""))
process.analyzeSIM = cms.EDAnalyzer("AnalyzeTP",
                                    triggerPrimitives=cms.InputTag("emulTP", "", ""))
process.compare = cms.EDAnalyzer("CompareTP",
                                 triggerPrimitives=cms.InputTag("hcalDigis"),
                                 emulTriggerPrimitives=cms.InputTag("emulTP"),
                                 swapIphi=cms.bool(False))

process.analyzeL1T = cms.EDAnalyzer("AnalyzeTP",
                                    triggerPrimitives = cms.InputTag("l1tCaloLayer1Digis", "" , "")
)


process.emulTP2016 = process.simHcalTriggerPrimitiveDigis.clone()
process.emulTP2016.upgradeHF = cms.bool(False)
process.emulTP2016.upgradeHE = cms.bool(True)
process.emulTP2016.inputLabel = cms.VInputTag("hcalDigis", "hcalDigis")
process.emulTP2016.inputUpgradeLabel = cms.VInputTag("hcalDigis", "hcalDigis")
process.emulTP2016.numberOfSamples = cms.int32(3)
process.emulTP2016.numberOfPresamples = cms.int32(1)

process.compare2016 = cms.EDAnalyzer("CompareTP",
                                     triggerPrimitives=cms.InputTag("hcalDigis"),
                                     emulTriggerPrimitives=cms.InputTag("emulTP2016"),
                                     swapIphi=cms.bool(False))

process.dump = cms.EDAnalyzer("EventContentAnalyzer")

process.p = cms.Path(
    process.hcalDigis *
    # process.dump *
    process.l1tCaloLayer1Digis *
    process.emulTP *
    process.analyzeRAW *
    process.analyzeSIM *
    process.compare 
#    process.emulTP2016 *
#    process.compare2016
)

# print process.dumpPython()
