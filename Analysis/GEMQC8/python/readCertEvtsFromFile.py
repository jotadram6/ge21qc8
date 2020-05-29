import os, sys, io

def GetCertifiedEvents(run_number):

    inputPath = os.path.abspath("readCertEvtsFromFile.py").split('QC8Test')[0]+'QC8Test/src/Analysis/GEMQC8/data/CertifiedEvents/'
    in_name = inputPath + "CertifiedEvents_run" + str(run_number) + ".csv"

    if (os.path.exists(in_name)):
        output = []
        with open(in_name) as infile:
            for line in infile:
                if ('\n' in line):
                    line = line.split('\n')[0]
                if ('\r' in line):
                    line = line.split('\r')[0]
                if (line.split(',')[0]!='RunNumber'):
                    output.append(line)

		print("\n")
		print("Success: the certified events have been read!")
		print("\n")

        return(output)

    else:
		print("\n")
		print("The input file does not exist!")
		print("\n")

if __name__ == '__main__':
    run_num = sys.argv[1]
    outList = GetCertifiedEvents(run_num)
