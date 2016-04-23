import sys

class Cell(object):

	def __init__(self, lineageName, cellL, cellR):
		print "initializing {lineageName} cell".format(lineageName=lineageName)
		self.lineageName = lineageName

		self.cellL = cellL
		self.cellR = cellR

	def getCellLineageName(self):
		return self.lineageName

	def getCellR(self):
		return self.cellR

	def getCellL(self):
		return self.cellL