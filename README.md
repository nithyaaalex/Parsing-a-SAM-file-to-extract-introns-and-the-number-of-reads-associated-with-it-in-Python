Extracting introns and the number of reads associated with them from a SAM file using base Python (no existing libraries used)

Input Required:
1. Standard sam file.
2. A comma-separated file with the gene_IDs, transcript_ID and genomic location of all the genes of whatever chromosome is aligned in the SAM file.

What the code does:
1. Parses the sam file, and identifies whether each line in the file is a deletion, insertion, matches, splits or a JUNCTION(intron)
2. For each junction identified, it counts how many reads are associated with it
3. Next it checks if the junction is inside a gene from the list of genes from the second file.
4. If the junction is inside a gene, it is marked as an intron, and it is outputted to a file that tells us the gene ID, the start and end of the intron, and the number of reads in that.
