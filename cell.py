import sys
import mRNAEvol

class Cell(object):

	def __init__(self, lineageName, cellL, cellR):
		print "initializing {lineageName} cell".format(lineageName=lineageName)
		self.lineageName = lineageName

		self.cellL = cellL
		self.cellR = cellR

		# a list of mRNAEvols between this cell and its entire subtree
		self.mRNAEvols = []

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

	def getmRNAEvols(self):
		return self.mRNAEvols