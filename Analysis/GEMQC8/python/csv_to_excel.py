import xlrd, openpyxl
import csv
import os, sys, io

def conversion(in_name):

	if in_name.endswith('.csv'):
		wb = openpyxl.Workbook()
		ws = wb.active

		with open(in_name) as infile:
		    reader = csv.reader(infile, delimiter=',')
		    for row in reader:
		        ws.append(row)

		out_name = in_name[:-4] + '.xlsx'
		wb.save(out_name)

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
