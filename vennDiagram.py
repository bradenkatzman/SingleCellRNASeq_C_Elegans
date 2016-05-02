import os

comma = ","
newline = "\n"
headerLineIncr = "Left Cell Increase" + comma + "Shared Increase" + comma + "Right Cell Increase"
headerLineDecr = "Left Cell Decrease" + comma + "Shared Decrease" + comma + "Right Cell Decrease"

# Build Venn Diagrams for Pairwise Cell Comparison
#
#	Current Interests:
#		- Increase and descrease of genes between germline cells P1-P4 and somatic siblings from a germline parent
#


# parameters:
#	- the 3 cells: the parent, the left somatic child, and the right germline child
#	- the dictionary of gene increases and decreases from the germline parent to the germ line child
#	- the dictionary of gene increases and decreases from the germline parent to the somatic child
def makeVennDiagram(germlineParent, somaticChild, germlineChild, glpToglc_dict, glpTosc_dict, p_value, medRPKMThreshold, logFCThreshold):
	print "\nBuilding Venn Diagrams for siblings: germline cell {glc} and somatic cell {sc} from germline parent {glp} ".format(
		glc=germlineChild, sc=somaticChild, glp=germlineParent)

	title = germlineChild + "_" + somaticChild
	criteria = ", ,Criteria:, p = {p}, medRPKM = {medRPKMThreshold}, logFC Threshold = {logFCThreshold}".format(p=p_value, medRPKMThreshold=medRPKMThreshold, logFCThreshold=logFCThreshold)

	# initialize lists
	leftCellIncr = []
	sharedIncr = []
	rightCellIncr = []

	leftCellDecr = []
	sharedDecr = []
	rightCellDecr = []

	for gene in glpToglc_dict:
		valR = glpToglc_dict[gene]

		if gene in glpTosc_dict.keys(): # the gene was statistically significant in the somatic cell too
			valL = glpTosc_dict[gene]

			if valR == valL: # both genes increased or decreased
				# check if increase or decrease
				if valR == 1:
					sharedIncr.append(gene)
				elif valR == -1:
					sharedDecr.append(gene)
				
			elif valR > valL:
				rightCellIncr.append(gene)
				leftCellDecr.append(gene)
				
			elif valR < valL:
				rightCellDecr.append(gene)
				leftCellIncr.append(gene)

		else: # the gene wasn't statistically significant in the somatic cell
			if valR == 1:
				rightCellIncr.append(gene)
				
			elif valR == -1:
				rightCellDecr.append(gene)


	for gene in glpTosc_dict:
		valL = glpTosc_dict[gene]

		if gene not in sharedIncr or gene not in sharedDecr:
			if valL == 1:
				leftCellIncr.append(gene)
			elif valL == -1:
				leftCellDecr.append(gene)


		# sort the lists
	leftCellIncr.sort()
	sharedIncr.sort()
	rightCellIncr.sort()
	leftCellDecr.sort()
	sharedDecr.sort()
	rightCellDecr.sort()

	writeVennDiagramToFile(title, criteria, leftCellIncr, sharedIncr, rightCellIncr, leftCellDecr, sharedDecr, rightCellDecr)

def writeVennDiagramToFile(title, criteria, leftCellIncr, sharedIncr, rightCellIncr, leftCellDecr, sharedDecr, rightCellDecr):
	print "writing venn diagram to file"

	resultsDirName = "Results"

	if not os.path.exists(resultsDirName):
		os.makedirs(resultsDirName)


	filename = title + "_vd.csv"

	file = open(resultsDirName + "/" + filename, "w+")

	file.write(headerLineIncr)
	file.write(criteria)

	iterator = 0

	# the increase diagram
	while iterator < len(leftCellIncr) or iterator < len(sharedIncr) or iterator < len(rightCellIncr):
		line = ""

		# left
		if iterator < len(leftCellIncr):
			line += (leftCellIncr[iterator] + comma)
		else:
			line += (" " + comma)

		# shared
		if iterator < len(sharedIncr):
			line += (sharedIncr[iterator] + comma)
		else:
			line += (" " + comma)

		# right
		if iterator < len(rightCellIncr):
			line += (rightCellIncr[iterator])
		else:
			line += " "

		file.write(newline)
		file.write(line)
		iterator += 1


	file.write(newline)
	file.write(newline)
	file.write(headerLineDecr)
	
	iterator = 0

	# the decrease diagram
	while iterator < len(leftCellDecr) or iterator < len(sharedDecr) or iterator < len(rightCellDecr):
		line = ""

		# left
		if iterator < len(leftCellDecr):
			line += (leftCellDecr[iterator] + comma)
		else:
			line += (" " + comma)

		# shared
		if iterator < len(sharedDecr):
			line += (sharedDecr[iterator] + comma)
		else:
			line += (" " + comma)

		# right
		if iterator < len(rightCellDecr):
			line += (rightCellDecr[iterator])
		else:
			line += " "

		file.write(newline)
		file.write(line)
		iterator += 1


	file.close()