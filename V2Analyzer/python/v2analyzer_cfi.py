import FWCore.ParameterSet.Config as cms

v2analyzer = cms.EDAnalyzer("V2Analyzer",
                            dzerr_ = cms.untracked.double(10.),
                            chi2_ = cms.untracked.double(40.),
                            jetAnal_ = cms.untracked.bool(True)
 )
