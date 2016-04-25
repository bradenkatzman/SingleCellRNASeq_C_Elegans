import sys

class mRNAEvol(object):

	def __init__(self, parentCellLineageName, childCellLineageName, numGenesIncreased, numGenesDecreased):
		self.parentCellLineageName = parentCellLineageName
		self.childCellLineageName = childCellLineageName
		self.numGenesIncreased = numGenesIncreased
		self.numGenesDecreased = numGenesDecreased

	def getParentCellLineageName(self):
		return self.parentCellLineageName

	def getChildCellLineageName(self):
		return self.childCellLineageName

	def getNumGenesIncreased(self):
		return self.numGenesIncreased

	def getNumGenesDecreased(self):
		return self.numGenesDecreased

	def toString(self):
		return "Parent: {pC} - Child: {cC}, Genes Increased: {genesIncr}, Genes Decreased: {genesDecr}".format(pC=self.parentCellLineageName, 
			cC=self.childCellLineageName, genesIncr=self.numGenesIncreased, genesDecr=self.numGenesDecreased)