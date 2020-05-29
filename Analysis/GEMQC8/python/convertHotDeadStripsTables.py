import csv
import os, sys, io

def convertHotDead(run_number,typeOfTable):

	if (typeOfTable == "hot"):
		hotTablesPath = os.path.abspath("convertHotDeadStripsTables.py").split('QC8Test')[0] + 'QC8Test/src/Analysis/GEMQC8/data/HotStripsTables/'
		infileName = hotTablesPath + "HotStrips_run" + str(run_number) + "_ToDB.csv"
	elif (typeOfTable == "dead"):
		deadTablesPath = os.path.abspath("convertHotDeadStripsTables.py").split('QC8Test')[0] + 'QC8Test/src/Analysis/GEMQC8/data/DeadStripsTables/'
		infileName = deadTablesPath + "DeadStrips_run" + str(run_number) + "_ToDB.csv"
	else: print("This mask type does not exist!")

	with open(infileName) as infile:
		for line in infile:
			line = line.split('\n')[0]
			SCtype = line.split(',')[0]
			if (SCtype=='RUN_NUMBER'):
				if (int(line.split(',')[1])!=int(run_number)):
					sys.exit('StandGeometryConfiguration file has something wrong: run rumber not matching...')

	outfileName = infileName[:-9] + '.csv'

	outfile = open(outfileName,"w")

	with open(infileName) as infile:
		for line in infile:
			line = line.split('\n')[0]
			FirstCol = line.split(',')[0]
			if (FirstCol!='RUN_NUMBER'):
				if (FirstCol=='CH_SERIAL_NUMBER'):
					line = line + ',RUN_NUMBER\n'
				else:
					line = line + ',' + str(run_number) + '\n'
				outfile.write(line)

	outfile.close()

	print("\n")
	print("Success: converted " + typeOfTable + " strips table for run " + str(run_number))
	print("\n")

if __name__ == '__main__':
	run_num = sys.argv[1]
	tableType = sys.argv[2]
	convertHotDead(run_num,tableType)
