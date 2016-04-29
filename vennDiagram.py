# Build Venn Diagrams for Pairwise Cell Comparison
#
#	Current Interests:
#		- Increase and descrease of genes between germline cells P1-P4 and somatic siblings
#


# parameters:
#	- the 2 cell objects: germlineCell, somaticCell
#	- the genes that increased and decreased in the germ line cell
#	- the genes that increased and decreased in the somatic cell
def makeVennDiagram(germlineCell, somaticCell, increaseGLC, decreaseGLC, increaseSC, decreaseSC):
	print "\nBuilding Venn Diagrams for germ line cell {glc} and somatic cell {sc}".format(
		glc=germlineCell.getCellLineageName(), sc=somaticCell.getCellLineageName)

		# initialize lists
		leftCellIncr = []
		sharedIncr = []
		rightCellIncr = []

		leftCellDecr = []
		sharedDecr = []
		rightCellDecr = []

		# build the increase diagram
		for gene in increaseGLC:
			if gene in increaseSC: # shared
				sharedIncr.append(gene)
			else: # left only
				leftCellIncr.append(gene)

		for gene in increaseSC:
			if gene is not in sharedIncr: # right only
				rightCellIncr.append(gene)


		# build the decrease diagram
		for gene in decreaseGLC:
			if gene in decreaseSC: # shared
				sharedDecr.append(gene)
			else: # left only
				leftCellDecr.append(gene)


		for gene in decreaseSC:
			if gene is not in sharedDecr:
				rightCellDecr.append(gene)