A Zhirong Bao Lab project to track the evolution of mRNAs in C. Elegans from the single cell to 16-cell stage of embryogenesis, and explore RNA seq related interestes 

In their paper, A Transcriptional Lineage of the Early C. elegans Embryo, Sophia Tintori et al profiled transcriptomes of single embryonic cells, for every cell from the 1-cell to 16-cell stage, with 5 or more replicates of each cell. 

In this project, we use the dataset to reconstruct the C. elegans lineage tree up to the 16 cell stage with mRNA evolution data, highlighting the number of genes that increase and decrease between any node and the subtree rooted at that node.

To begin, we chose a reasonable p-value cutoff and median RPKM cutoff at:
- p = 0.01
- med. rpkm = 10

Next, we traverse the lineage tree and at each step compare a node to its subtree's nodes and find the genes that meet our cutoff criteria and track the number of genes which increase and decrease in expression.

We use the data to explore the following topics:
	- Building gene modules for somatic and germline siblings, to understanding the underlying process that controls their differential gene expression
	- Expectation-Maximization of our p-value, med. RPKM, logFC parameters that best explain the biological processes behind the data 

This project is done under the guidance and supervision of Dr. Zhirong Bao and Dr. Yichi Xu and implemented by Braden Katzman (B.A. Candidate at Columbia University)


For further information and/or questions, please contact bmk2137@columbia.edu (Braden Katzman)

Resources:
Bao Lab @ MSKCC:
	- https://www.mskcc.org/research-areas/labs/zhirong-bao

Preprint: A Transcriptional Lineage of the Early C. elegans Embryo
	- http://biorxiv.org/content/early/2016/04/07/047746

Visual Interactive Tool of C. elegans mRNA exploration:
	- http://tintori.bio.unc.edu/