import os

comma = ","
newline = "\n"
headerLineIncr = "Left Cell Increase; logFC" + comma + "Shared Increase; logFC (L) - logFC (R)" + comma + "Right Cell Increase; logFC"
headerLineDecr = "Left Cell Decrease; logFC" + comma + "Shared Decrease; logFC (L) - logFC (R)" + comma + "Right Cell Decrease; logFC"

headerLine = "Left Cell Higher; logFC" + comma + "~Same Level; logFC" + comma + "Right Cell Higher; logFC"

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
	

	if glpTosc_dict == None: # just a single comparison case --> i.e. between 2 cells
		print "\nBuilding Venn Diagrams for siblings: germline cell {glc} and somatic cell {sc}".format(
			glc=germlineChild, sc=somaticChild)


		title = germlineChild + "_" + somaticChild
		criteria = ", ,Criteria:, p = {p}, medRPKM = {medRPKMThreshold}, logFC Threshold = {logFCThreshold}".format(p=p_value, medRPKMThreshold=medRPKMThreshold, logFCThreshold=logFCThreshold)
		
		# initialize lists
		leftCellHigher = []
		rightCellHigher = []
		shared = []



		for gene in glpToglc_dict:
			val = glpToglc_dict[gene][0]
			logFC = str(glpToglc_dict[gene][1])

			if val == -1:
				rightCellHigher.append(gene + "; " + logFC)
			elif val == 1:
				leftCellHigher.append(gene + "; " + logFC)
			elif val == 0:
				shared.append(gene + "; " + logFC)

		# sort the lists
		leftCellHigher.sort()
		rightCellHigher.sort()
		shared.sort()

		writeVennDiagramToFile1(title, criteria, leftCellHigher, shared, rightCellHigher)
		return

	# TO DO
	#	- make shared section

	else:
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
			valR = glpToglc_dict[gene][0]
			logFCR = str(glpToglc_dict[gene][1])

			if gene in glpTosc_dict.keys(): # the gene was statistically significant in the somatic cell too
				valL = glpTosc_dict[gene][0]
				logFCL = str(glpTosc_dict[gene][1])

				if valR == valL: # both genes increased or decreased
					# check if increase or decrease
					if valR == 1:
						if gene not in sharedIncr:
							sharedIncr.append(gene + "; " + logFCL + " - " + logFCR)
					elif valR == -1:
						if gene not in sharedDecr:
							sharedDecr.append(gene + "; " + logFCL + " - " + logFCR)
					
				elif valR > valL:
					if gene not in rightCellIncr:
						rightCellIncr.append(gene + "; " + logFCR)
					if gene not in leftCellDecr:
						leftCellDecr.append(gene + "; " + logFCL)
					
				elif valR < valL:
					if gene not in rightCellDecr:
						rightCellDecr.append(gene + "; " + logFCR)
					if gene not in leftCellIncr:
						leftCellIncr.append(gene + "; " + logFCL)

			else: # the gene wasn't statistically significant in the somatic cell
				if valR == 1:
					if gene not in rightCellIncr:
						rightCellIncr.append(gene + "; " + logFCR)
					
				elif valR == -1:
					if gene not in rightCellDecr:
						rightCellDecr.append(gene + "; " + logFCR)

		for gene in glpTosc_dict:
			valL = glpTosc_dict[gene]
			logFCL = str(glpTosc_dict[gene][1])

			if valL == 1:
				if gene not in sharedIncr and gene not in leftCellIncr:
					leftCellIncr.append(gene + "; " + logFCL)
			elif valL == -1:
				if gene not in sharedDecr and gene not in leftCellDecr:
					leftCellDecr.append(gene + "; " + logFCL)


		# sort the lists
		leftCellIncr.sort()
		sharedIncr.sort()
		rightCellIncr.sort()
		leftCellDecr.sort()
		sharedDecr.sort()
		rightCellDecr.sort()

		writeVennDiagramToFile2(title, criteria, leftCellIncr, sharedIncr, rightCellIncr, leftCellDecr, sharedDecr, rightCellDecr)
		return

def writeVennDiagramToFile1(title, criteria, leftCellHigher, shared, rightCellHigher):
	print "writing venn diagram to file"

	resultsDirName = "Results"

	if not os.path.exists(resultsDirName):
		os.makedirs(resultsDirName)


	filename = title + "_vd.csv"

	file = open(resultsDirName + "/" + filename, "w+")

	file.write(headerLine)
	file.write(criteria)

	iterator = 0

	# the increase diagram
	while iterator < len(leftCellHigher) or iterator < len(shared) or iterator < len(rightCellHigher):
		line = ""

		# left
		if iterator < len(leftCellHigher):
			line += (leftCellHigher[iterator] + comma)
		else:
			line += (" " + comma)

		# shared
		if iterator < len(shared):
			line += (shared[iterator] + comma)
		else:
			line += (" " + comma)

		# right
		if iterator < len(rightCellHigher):
			line += (rightCellHigher[iterator])
		else:
			line += " "

		file.write(newline)
		file.write(line)
		iterator += 1


	file.close()


def writeVennDiagramToFile2(title, criteria, leftCellIncr, sharedIncr, rightCellIncr, leftCellDecr, sharedDecr, rightCellDecr):
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