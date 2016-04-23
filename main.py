import sys
import time
import cell

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


if __name__ == '__main__':
	t0 = time.clock()
	print "start"

	P0_root = buildTree()

	

	print "\nprogram execution: {t} seconds".format(t=time.clock()-t0)
	print "exiting"