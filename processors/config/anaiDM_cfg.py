import HpstrConf
import sys
import os
import baseConfig as base

base.parser.add_argument("-V", "--splitVolume", type=int, dest="splitVolume",
                         help="Require positron in Top and Bottom", metavar="splitVolume", default=0)

options = base.parser.parse_args()

# Use the input file to set the output file name
infile = options.inFilename
outfile = options.outFilename

print('Input file: %s' % infile)
print('Output file: %s' % outfile)

p = HpstrConf.Process()

p.run_mode = 1
#p.max_events = 1000

# Library containing processors
p.add_library("libprocessors")

###############################
#          Processors         #
###############################

recoana = HpstrConf.Processor('idmpre', 'iDMPreSelection')

###############################
#   Processor Configuration   #
###############################
#RecoHitAna
recoana.parameters["anaName"] = "idmpre"
recoana.parameters["trkColl"] = "KalmanFullTracks"
recoana.parameters["tsColl"] = "TSData"
recoana.parameters["vtxColl"] = "UnconstrainedV0Vertices_KF"
recoana.parameters["mcColl"] = "MCParticle"
recoana.parameters["hitColl"] = "SiClusters"
recoana.parameters["ecalColl"] = "RecoEcalClusters"
recoana.parameters["vtxSelectionjson"] = os.environ['HPSTR_BASE']+"/analysis/selections/idm/vertexSelection_2016.json"
#####
recoana.parameters["beamE"] = base.beamE[str(options.year)]
recoana.parameters["isData"] = options.isData
recoana.parameters["analysis"] = options.analysis
recoana.parameters["debug"] = 0
CalTimeOffset = -999

if (options.isData == 1):
    CalTimeOffset = 56.
    print("Running on data file: Setting CalTimeOffset %d" % CalTimeOffset)
elif (options.isData == 0):
    CalTimeOffset = 43.
    print("Running on MC file: Setting CalTimeOffset %d" % CalTimeOffset)
else:
    print("Specify which type of ntuple you are running on: -t 1 [for Data] / -t 0 [for MC]")

recoana.parameters["CalTimeOffset"] = CalTimeOffset

p.sequence = [recoana]
p.skip_events = options.skip_events
p.max_events = options.nevents
p.input_files = infile
p.output_files = [outfile]

p.printProcess()
