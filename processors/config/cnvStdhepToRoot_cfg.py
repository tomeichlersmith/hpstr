import HpstrConf
import baseConfig as base
import os


options = base.parser.parse_args()

# Use the input file to set the output file name
stdhep_file = options.inFilename[0]
root_file = options.outFilename


p = HpstrConf.Process()

p.run_mode = 2
p.skip_events = options.skip_events
if(options.nevents>-1):
    p.max_events = options.skip_events+options.nevents
else:
    p.max_events = -1

# Library containing processors
p.libraries.append("libprocessors.so")

###############################
#          Processors         #
###############################

cnvStd = HpstrConf.Processor('cnvStd', 'StdhepMCParticleProcessor')

###############################
#   Processor Configuration   #
###############################
#MCParticles
cnvStd.parameters["mcPartCollStdhep"] = 'MCParticle'
cnvStd.parameters["mcPartCollRoot"] = 'MCParticle'
cnvStd.parameters["maxEvent"] = options.skip_events+options.nevents

# Sequence which the processors will run.
p.sequence = [cnvStd]

p.input_files=[stdhep_file]
#p.input_files=[]
p.output_files = [root_file]

p.printProcess()
