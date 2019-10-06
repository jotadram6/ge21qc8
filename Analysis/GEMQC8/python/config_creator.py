import csv
import os, sys, io

def configMaker(run_number,masks):

	minClusterSize = 1
	maxClusterSize = 10
	maxResidual = 5.0
	trackChi2 = 3
	trackResX = 0.3
	trackResY = 0.3697
	MulSigmaOnWindow = 5
	minRecHitsPerTrack = 4

	configTablesPath = os.path.abspath("config_creator.py").split('QC8Test')[0] + 'QC8Test/src/Analysis/GEMQC8/data/StandConfigurationTables/'

	runPath = os.path.abspath("config_creator.py").split('QC8Test')[0] + 'QC8Test/src/Analysis/GEMQC8/test/'

	infileName = configTablesPath + "StandGeometryConfiguration_run" + str(run_number) + ".csv"

	with open(infileName) as infile:
		for line in infile:
			line = line.split('\n')[0]
			SCtype = line.split(',')[0]
			if (SCtype!='CH_SERIAL_NUMBER'):
				if (int(line.split(',')[8])!=int(run_number)):
					sys.exit('StandGeometryConfiguration file has something wrong: run rumber not matching...')

	out_name = 'out_run_'
	for i in range(6-len(str(run_number))):
	    out_name = out_name + '0'
	out_name = out_name + str(run_number) + '.root'

	outfileName = runPath + "configureRun_cfi.py"

	outfile = open(outfileName,"w")

	outfile.write('RunNumber = ' + str(run_number) + '\n\n')

	outfile.write('# Output file name definition\n')
	outfile.write('OutputFileName = \'' + out_name + '\'\n\n')

	outfile.write('# Parameters definition\n')
	outfile.write('minClusterSize = {}\n'.format(minClusterSize))
	outfile.write('maxClusterSize = {}\n'.format(maxClusterSize))
	outfile.write('maxResidual = {} # cm\n'.format(maxResidual))
	outfile.write('trackChi2 = {}\n'.format(trackChi2))
	outfile.write('trackResX = {}\n'.format(trackResX))
	outfile.write('trackResY = {}\n'.format(trackResY))
	outfile.write('MulSigmaOnWindow = {}\n'.format(MulSigmaOnWindow))
	outfile.write('minRecHitsPerTrack = {}\n'.format(minRecHitsPerTrack))
	if (masks == "yesMasks") outfile.write('applyMasks = True\n')
	if (masks == "noMasks") outfile.write('applyMasks = False\n')

	outfile.write('# Stand configuration definition\n')
	StandConfiguration = ['0','0','0','0','0',\
	                      '0','0','0','0','0',\
	                      '0','0','0','0','0']

	with open(infileName) as infile:
		for line in infile:
			line = line.split('\n')[0]
			SCtype = line.split(',')[0]
			if (SCtype!='CH_SERIAL_NUMBER'):
				position = line.split(',')[2]
				row = int(position.split('/')[0])
				column = int(position.split('/')[1])
				SCnumber = (5 * (column - 1)) + (row - 1)
				StandConfiguration[SCnumber] = (SCtype)[8]

	outfile.write('StandConfiguration = [\\\n')
	for entry in range(15):
		if (entry==4 or entry==9):
			outfile.write('\'' + StandConfiguration[entry] + '\',\\\n')
		elif (entry==14):
			outfile.write('\'' + StandConfiguration[entry] + '\']')
		else:
			outfile.write('\'' + StandConfiguration[entry] + '\',')

	outfile.close()

	print("\n")
	print("Success: configuration file created for run " + str(run_number))
	print("\n")

if __name__ == '__main__':
	run_num = sys.argv[1]
	configMaker(run_num)
