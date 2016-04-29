import sys
import os
import csv
import treeSearch


cells = ["P1", "P2", "P3", "P4", "AB", "EMS", "C", "D"]

id_cell_dict = {}

def initIDsToCellList():
	print "\nmatching IDs to cell names"

	tableS1_filename = "Tables1.csv"

	with open(tableS1_filename, 'rU') as csvfile:
		reader = csv.reader(csvfile, delimiter= ',')
		for row in reader:
			sample = str(row[0])

			try:
				idx = sample.index("_")
				ID = sample[0:idx]
				cell = sample[idx+1:sample.index("_", idx+1)]
					
				# only add cells of usable quality
				if cell != "tossed":
					id_cell_dict[ID] = cell

			except:
				pass

def readRPMVals(root):
	print "\nadding RPKM values for each sample to respective cells"

	dirname = "Tintori_2016"
	ext = ".txt"

	for filename in os.listdir(dirname):
		# check if correct extension
		if filename.endswith(ext):
			# extract the sample name
			idx1 = filename.index("_")
			idx2 = filename.index("_", idx1+1)

			sample = filename[idx1+1: idx2]

			try:
				targetCellName = id_cell_dict[sample]

				# see if cell isn't completely identified i.e. ends in 'x' or in #
				# figure out if ends with a,p or l,r
				char = targetCellName[len(targetCellName)-2]
				if targetCellName[len(targetCellName)-1].lower() == 'x':
					if char.lower() == 'a' or char.lower() == 'ap':
						cell1 = treeSearch.findCellInOrder(root, targetCellName[:len(targetCellName)-1] + "l")
						if cell1 != None:
							cell1.addSample(collectRPKMValues(dirname, filename))
						else:
							print "error 1"

						cell2 = treeSearch.findCellInOrder(root, targetCellName[:len(targetCellName)-1] + "r")
						if cell2 != None:
							cell2.addSample(collectRPKMValues(dirname, filename))
						else:
							print "error 2"

					elif char.lower() == 'l' or char.lower() == 'r':
						cell1 = treeSearch.findCellInOrder(root, targetCellName[:len(targetCellName)-1] + "a")
						if cell1 != None:
								cell1.addSample(collectRPKMValues(dirname, filename))
						else:
							print "error 3"

						cell2 = treeSearch.findCellInOrder(root, targetCellName[:len(targetCellName)-1] + "p")
						if cell2 != None:
							cell2.addSample(collectRPKMValues(dirname, filename))
						else:
							print "error 4"

				elif targetCellName[len(targetCellName)-1].lower().isdigit() and char.lower() == 'x':
					cell1 = treeSearch.findCellInOrder(root, targetCellName[:len(targetCellName)-2] + "a")
					if cell1 != None:
						cell1.addSample(collectRPKMValues(dirname, filename))
					else:
						print "error 5"

					cell2 = treeSearch.findCellInOrder(root, targetCellName[:len(targetCellName)-2] + "p")
					if cell2 != None:
						cell2.addSample(collectRPKMValues(dirname, filename))
					else:
						print "error 6"
				else:
					# find the cell in the lineage tree
					cell = treeSearch.findCellInOrder(root, targetCellName)

					# add the rpkms to the cell
					if cell != None:
						cell.addSample(collectRPKMValues(dirname, filename))
					else:
						print "error 7"

				


			except:
				pass

def collectRPKMValues(dirname, filename):
	# get the rpkm values
	gene_RPKM_dict = {}
	with open(dirname + "/" + filename) as ins:
		for line in ins:
			tokens = line.split("	")
	
			gene = tokens[0]
			rpkm = float(tokens[3])
											
			gene_RPKM_dict[gene] = rpkm						
					
	return gene_RPKM_dict				