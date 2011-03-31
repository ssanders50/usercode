import FWCore.ParameterSet.Config as cms

v2analyzer = cms.EDAnalyzer("V2Analyzer",
                            dzerr_ = cms.untracked.double(14.),
                            chi2_ = cms.untracked.double(80.)
 )
