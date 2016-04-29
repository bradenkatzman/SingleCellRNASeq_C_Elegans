import sys
import math
import mRNAEvol

class Cell(object):

	def __init__(self, lineageName, cellL, cellR):
		print "initializing {lineageName} cell".format(lineageName=lineageName)
		self.lineageName = lineageName

		self.cellL = cellL
		self.cellR = cellR

		# a list of the samples of transcriptomes for this cell
		self.samples = []

		# a list of gene RPKMs for the samples of this cell
		# self.RPKMs = []

		# a list of mRNAEvols between this cell and its entire subtree
		self.mRNAEvols = []

	def addSample(self, gene_RPKM_dict):
		self.samples.append(gene_RPKM_dict)

	def getGeneRPKMs(self, gene):
		gene_RPKMs = []

		for _dict in self.samples:
			gene_RPKMs.append(_dict[gene])

		return gene_RPKMs

	def getMedGenRPKM(self, gene):
		gene_RPKMs = self.getGeneRPKMs(gene)
		gene_RPKMs.sort()

		try:
			if len(gene_RPKMs) % 2 == 0: # even number of items in the list
				return float(gene_RPKMs[int(math.ceil(len(gene_RPKMs)/2))] + gene_RPKMs[int(math.floor(len(gene_RPKMs)/2))])
			else:
				return float(gene_RPKMs[int(math.ceil(len(gene_RPKMs)/2))])
		except:
			print "problem with {cell} on {gene}".format(cell=self.lineageName, gene=gene)
			print len(gene_RPKMs)


	def addmRNAEvol(self, childCellLineageName, numGenesIncreased, numGenesDecreased):
		evol = mRNAEvol.mRNAEvol(self.lineageName, childCellLineageName, numGenesIncreased, numGenesDecreased)
		self.mRNAEvols.append(evol)

	# search methods as necessary
	def getmRNAEvolByChildCellName(self, childCellLineageName):
		for evol in self.mRNAEvols:
			if evol.getChildCellLineageName().lower() == childCellLineageName.lower():
				return evol

		return 0

	def getCellLineageName(self):
		return self.lineageName

	def getCellR(self):
		return self.cellR

	def getCellL(self):
		return self.cellL

	def getSamples(self):
		return self.samples

	def getmRNAEvols(self):
		return self.mRNAEvols