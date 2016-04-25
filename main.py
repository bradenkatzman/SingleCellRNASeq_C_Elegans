import sys
import os
import time
import cell
import mRNAEvol

preOrderTree = ["P0", "AB", "ABa", "ABal", "ABala", "ABalp", "ABar", "ABara", "ABarp", "ABp", "ABpl", "ABpla", "ABplp",
"ABpr", "ABpra", "ABprp", "P", "EMS", "MS", "MSa", "MSp", "E", "Ea", "Ep", "P2", "C", "Ca", "Cp", "P3", "D", "P4"]

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

def addmRNADataToTree(root):
	print "\n adding mRNA data to tree"

	dirName = "export_data"

	underscore = "_"
	ext = ".csv"

	for filename in os.listdir(dirName):
		# extract the cells name from the filename
		delineator = filename.index(underscore)
		end = filename.index(ext)

		parentCell = filename[0:delineator]
		childCell = filename[delineator+1:end]

		with open(dirName + "/" + filename) as ins:
			for line in ins:
				geneDataTokens = line.split(",")

				p_value = float(geneDataTokens[1])
				logCPM = float(geneDataTokens[2])
				logFC = float(geneDataTokens[3])

				c = findCell(root, parentCell)

				if logFC > 0.: # the parent cell expresses this gene at a higher level than the child cell
					y = 1
				elif logFC < 0.: # the parent cell expresses this gene at a lower level than the child cell
					y = 0


def findCell(root, targetCellLineageName):
	if root == None:
		return 0

	if root.getCellLineageName() == targetCellLineageName:
		return root

	c = findCell(root.getCellL(), targetCellLineageName)
	if c != 0:
		return c
	else:
		c = findCell(root.getCellR(), targetCellLineageName)


if __name__ == '__main__':
	t0 = time.clock()
	print "start"

	P0_root = buildTree()

	addmRNADataToTree(P0_root)

	print "\nprogram execution: {t} seconds".format(t=time.clock()-t0)
	print "exiting"