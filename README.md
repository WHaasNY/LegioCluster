# LegioCluster
NGS analysis and clustering pipeline for bacteria undergoing frequent recombination, such as Legionella


LegioCluster is a whole-genome sequencing analysis pipeline designed to determine the point source during an outbreak of Legionella pneumophila and other bacteria with high levels of genome plasticity. The pipeline optimizes genome comparisons by automatically selecting the best matching reference genome prior to mapping and variant calling. This approach reduces the number of false-positive variant calls, maximizes the fraction of all genomes that are being compared, and naturally clusters the isolates according to their reference strain. Isolates that are too distant from any genome in the database are added to the list of candidate references, thereby creating a new cluster. Short insertions and deletions are considered in addition to single nucleotide polymorphisms for increased discriminatory power. The program has been successfully deployed at the New York State Department of Health's Wadsworth Center to track outbreaks of L. pneumophila as well as other bacterial species.


### Required software:
- Linux/Unix/Mac operating system
- Docker
- Python 3.7 or higher
- the Python packages pydot, numpy and matplotlib, which can be installed with:   pip install <package-name>



