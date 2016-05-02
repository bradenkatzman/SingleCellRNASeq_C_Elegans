import sys
import os
import time
import treeSearch
import initRPKMs
import cell
import mRNAEvol
import vennDiagram

preOrderTree = ["P0", "AB", "ABa", "ABal", "ABala", "ABalp", "ABar", "ABara", "ABarp", "ABp", "ABpl", "ABpla", "ABplp",
"ABpr", "ABpra", "ABprp", "P", "EMS", "MS", "MSa", "MSp", "E", "Ea", "Ep", "P2", "C", "Ca", "Cp", "P3", "D", "P4"]

headerLine = "Parent Cell - Child Cell, Gene Exp. Increased in Child, Gene Exp. Decreased in Child, Difference (Incr - Decr)"
headerLine2 = "Total Num Increases, Total Num Decreases"

newline = "\n"
dash = "-"

data = {}

# parameters
p = 0.0100000000
medRPKMThreshold = 25
logFCThreshold = 0.5
CPM = -1 # not currently in use


def buildTree():
	print "\nconstructing tree rooted at P0"

	# 16 cell stage --> leaves 
	# (left side of lineage)
	ABala = cell.Cell("ABala", None, None)
	ABalp = cell.Cell("ABalp", None, None)

	ABara = cell.Cell("ABara", None, None)
	ABarp = cell.Cell("ABarp", None, None)

	ABpla = cell.Cell("ABpla", None, None)
	ABplp = cell.Cell("ABplp", None, None)

	ABpra = cell.Cell("ABpra", None, None)
	ABprp = cell.Cell("ABprp", None, None)

	# (right side of lineage)
	MSa = cell.Cell("MSa", None, None)
	MSp = cell.Cell("MSp", None, None)
	
	Ea = cell.Cell("Ea", None, None)
	Ep = cell.Cell("Ep", None, None)

	Ca = cell.Cell("Ca", None, None)
	Cp = cell.Cell("Cp", None, None)

	D = cell.Cell("D", None, None)
	P4 = cell.Cell("P4", None, None)

	# 8 cell stage
	# (left side of lineage)
	ABal = cell.Cell("ABal", ABala, ABalp)
	ABar = cell.Cell("ABar", ABara, ABarp)

	ABpl = cell.Cell("ABpl", ABpla, ABplp)
	ABpr = cell.Cell("ABpr", ABpra, ABprp)

	# (right side of lineage)
	MS = cell.Cell("MS", MSa, MSp)
	E = cell.Cell("E", Ea, Ep)

	C = cell.Cell("C", Ca, Cp)
	P3 = cell.Cell("P3", D, P4)

	# 4 cell stage
	# (left side of lineage)
	ABa = cell.Cell("ABa", ABal, ABar)
	ABp = cell.Cell("ABp", ABpl, ABpr)

	# (right side of lineage)
	EMS = cell.Cell("EMS", MS, E)
	P2 = cell.Cell("P2", C, P3)

	# 2 cell stage
	AB = cell.Cell("AB", ABa, ABp) # left side
	P1 = cell.Cell("P1", EMS, P2) # right side

	# root
	P0 = cell.Cell("P0", AB, P1)

	return P0

def gatherRPKMData(root):
	initRPKMs.initIDsToCellList()

	initRPKMs.readRPMVals(root)


def addmRNADataToTree(root):
	print "\nadding mRNA data to tree"

	results = []

	dirName = "export_data"

	underscore = "_"
	ext = ".csv"


	totalNumIncreases = 0
	totalNumDecreases = 0

	for filename in os.listdir(dirName):
		# extract the cells name from the filename
		delineator = filename.index(underscore)
		end = filename.index(ext)

		parentCell = filename[0:delineator]
		childCell = filename[delineator+1:end]

		geneIncreaseCounter = 0
		geneDecreaseCounter = 0
		with open(dirName + "/" + filename) as ins:
			dataEntry = {}
			
			for line in ins:
				geneDataTokens = line.split(",")

				gene = geneDataTokens[0]
				p_value = float(geneDataTokens[1])
				logCPM = float(geneDataTokens[2])
				logFC = float(geneDataTokens[3])

				pc = treeSearch.findCellInOrder(root, parentCell)
				cc = treeSearch.findCellInOrder(root, childCell)

				if pc != None:
					if p_value <= p: # manually set p-value
						if logFC > logFCThreshold: # the parent cell expresses this gene at a higher level than the child cell
							geneDecreaseCounter += 1

							# make sure the gene meets the RPKM threshold in at least one of the cells
							if pc.getMedGeneRPKM(gene) > medRPKMThreshold or cc.getMedGeneRPKM(gene) > medRPKMThreshold:
								dataEntry[gene] = 1 # use 1 to denote increase on this gene from parent to child




						elif logFC < logFCThreshold: # the parent cell expresses this gene at a lower level than the child cell
							geneIncreaseCounter += 1

							if pc.getMedGeneRPKM(gene) > medRPKMThreshold or cc.getMedGeneRPKM(gene) > medRPKMThreshold:
								dataEntry[gene] = -1 # use -1 to denote decrease on this gene from parent to child


		data[parentCell + dash + childCell] = dataEntry

		# find the difference
		diff = geneIncreaseCounter - geneDecreaseCounter
		if diff < 0:
			totalNumDecreases += 1
		else:
			totalNumIncreases += 1

		# print the results
		print "{pc} - {cc}	{geneIncr}	{geneDecr}	{diff}".format(pc=parentCell, cc=childCell, 
			geneIncr=geneIncreaseCounter, geneDecr=geneDecreaseCounter, diff=diff)

		result = [parentCell + " - " + childCell, geneIncreaseCounter, geneDecreaseCounter, diff]
		results.append(result)
		

	# sort the results
	results.sort(key=lambda x: x[0])

	# write to file
	writeToFile(results, totalNumIncreases, totalNumDecreases)

def writeToFile(results, totalNumIncreases, totalNumDecreases):
	print "\nwriting results to file"

	resultsDirName = "Results"

	filename = "mRNA_Evolution_Data.csv"

	if not os.path.exists(resultsDirName):
		os.makedirs(resultsDirName)

	file = open(resultsDirName + "/" + filename, "w+")

	file.write(headerLine)

	size = len(results)
	iterator = 0

	while iterator < size:
		if len(results[iterator]) == 4:
			file.write(newline)
			file.write(str(results[iterator][0]) + ", " + str(results[iterator][1]) + ", " + str(results[iterator][2]) + ", " + str(results[iterator][3]))
			
		iterator += 1

	# write the global values to file
	file.write(newline)
	file.write(newline)

	file.write(headerLine2)
	file.write(newline)
	file.write(str(totalNumIncreases) + ", " + str(totalNumDecreases))

	file.close()



if __name__ == '__main__':
	t0 = time.clock()
	print "start"

	P0_root = buildTree()

	gatherRPKMData(P0_root)

	# also writes the data to file
	addmRNADataToTree(P0_root)

	vennDiagram.makeVennDiagram("P0", "AB", "P1", data["P0-P1"], data["P0-AB"], p, medRPKMThreshold, logFCThreshold)

	print "\nprogram execution: {t} seconds".format(t=time.clock()-t0)
	print "exiting"