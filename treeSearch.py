import cell

def findCellInOrder(root, targetCellLineageName):
	if not root:
		return None

	result = findCellInOrder(root.getCellL(), targetCellLineageName)

	if result is not None:
		return result

	if root.getCellLineageName().lower() == targetCellLineageName.lower():
		return root

	return findCellInOrder(root.getCellR(), targetCellLineageName)