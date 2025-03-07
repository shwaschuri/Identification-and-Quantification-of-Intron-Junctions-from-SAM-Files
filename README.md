Identification and Quantification of intron junctions from SAM files.

Identification and Quantification of Intron Junctions from SAM FilesOverviewThis script processes a SAM file and a gene location file to extract and count intronic junctions. It identifies splice junctions based on the CIGAR string from the SAM file and maps them to gene locations, outputting the count of occurrences for each intron.
FeaturesParses SAM files to extract intronic junctions using CIGAR strings.
Matches intronic junctions to gene locations.
Outputs a tab-separated file with gene IDs, junction start and end positions, and counts.
Skips reads with NH:i values not equal to 1 to ensure unique mapping.
Efficiently processes large datasets with optimized regex-based parsing.
DependenciesPython 3
UsageRun the script via the command line:
python script.py <sam_file> <gene_location_file>
Input FilesSAM File: Contains alignment information, including CIGAR strings.
Gene Location File: Contains gene IDs and their genomic coordinates.
Output File3026165.txt (Generated output)
Format: GeneID    JunctionStart    JunctionEnd    Count
Returns a list of junction start and end positions.
loc_parse(loc)Extracts and formats gene locations from the input file.
seperation(start_end)Splits gene location coordinates into separate start and end positions.
Logic FlowParse the gene location file to extract gene coordinates.
Read the SAM file, extracting reads with introns (N in CIGAR string).
Match intronic junctions with gene locations.
Count occurrences of each intronic junction per gene.
Output results in a tab-separated file.
Example Runpython script.py sample.sam genes.txtNotesEnsure input files are correctly formatted.
Handles large datasets but may require optimization for very large SAM files.
LicenseThis script is open-source and free to use for research and educational purposes.