import xlrd, openpyxl
import csv
import os, sys, io

def conversion(in_name):

	if in_name.endswith('.xlsx'):
		infile = openpyxl.load_workbook(in_name)
		sh = infile.active
		out_name = in_name[:-5] + '.csv'
		if sys.version_info[0] == 2: # python2
			with open(out_name, 'wb') as outfile:
				c = csv.writer(outfile)
				for r in sh.rows:
					c.writerow([cell.value for cell in r])
		if sys.version_info[0] == 3: # python3
			with open(out_name, 'w', newline="") as outfile:
				c = csv.writer(outfile)
				for r in sh.rows:
					c.writerow([cell.value for cell in r])
		print("\n")
		print("Success: " + in_name + " converted to " + out_name)
		print("\n")

	else:
		print("\n")
		print("The input file does not have the correct extension!")
		print("\n")

if __name__ == '__main__':
    input_file = sys.argv[1]
    conversion(input_file)
