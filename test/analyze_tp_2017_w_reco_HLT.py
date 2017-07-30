import FWCore.ParameterSet.Config as cms

from Configuration.AlCa.GlobalTag import GlobalTag
from Configuration.StandardSequences.Eras import eras

# Import of standard configurations
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100


process = cms.Process('PLOT', eras.Run2_2017)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag = GlobalTag(process.GlobalTag, '92X_dataRun2_HLT_v4', '')

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(100))

process.source = cms.Source("PoolSource",
                            skipBadFiles=cms.untracked.bool(True),
                            fileNames=cms.untracked.vstring(),
                            secondaryFileNames=cms.untracked.vstring(),
                            firstRun=cms.untracked.uint32(260627)
                            )


process.source.fileNames.extend([

'/store/data/Run2017B/ZeroBias/RECO/PromptReco-v1/000/297/469/00000/1C06F39D-B85A-E711-B8B3-02163E01A5DF.root'

])

process.source.secondaryFileNames.extend([

'/store/data/Run2017B/ZeroBias/RAW/v1/000/297/469/00000/F04F22DC-2459-E711-B1C7-02163E019D7B.root'

])

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.RawToDigi_Data_cff")
process.load("SimCalorimetry.Configuration.hcalDigiSequence_cff")
process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')

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


import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing ('standard')
options.register('hltName', 'HLT', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "HLT menu to use for trigger matching")

# process.hcalDigis.InputLabel = cms.InputTag("source")
process.analyzeRAW = cms.EDAnalyzer("AnalyzeTP",
                                    triggerPrimitives=cms.InputTag("hcalDigis", "", ""),

)
process.analyzeSIM = cms.EDAnalyzer("AnalyzeTP",
                                    triggerPrimitives=cms.InputTag("emulTP", "", ""),

)
process.compare = cms.EDAnalyzer("CompareTP",
                                    triggerPrimitives=cms.InputTag("hcalDigis"),
                                    emulTriggerPrimitives=cms.InputTag("emulTP"),
                                    swapIphi=cms.bool(False),

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

process.chainplotter = cms.EDAnalyzer("HcalCompareLegacyChains",
                                      triggerPrimitives=cms.InputTag('hcalDigis', '', 'PLOT'),
                                      recHits=cms.VInputTag('hbhereco', 'hfreco'),
                                      dataFrames=cms.VInputTag(
                                          cms.InputTag("hcalDigis", "", "PLOT"),
                                          cms.InputTag("hcalDigis", "", "PLOT")
                                      ),
                                      swapIphi=cms.bool(False),

                                      offlinePrimaryVertices = cms.InputTag("offlinePrimaryVertices","","RECO"),
                                      trigTagSrc = cms.InputTag("TriggerResults","",options.hltName)

                                      )

process.dump = cms.EDAnalyzer("EventContentAnalyzer")

process.p = cms.Path(
    process.hcalDigis *
    # process.dump *
    process.emulTP *
#    process.analyzeRAW *
#    process.analyzeSIM *
#    process.compare *
    process.chainplotter
#    process.emulTP2016 *
#    process.compare2016
)

# print process.dumpPython()
