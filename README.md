# LegioCluster
NGS analysis and clustering pipeline for bacteria undergoing frequent recombination, such as Legionella

LegioCluster is a whole-genome sequencing analysis pipeline designed to determine the point source during an outbreak of Legionella pneumophila and other bacteria with high levels of genome plasticity. The pipeline optimizes genome comparisons by automatically selecting the best matching reference genome prior to mapping and variant calling. This approach reduces the number of false-positive variant calls, maximizes the fraction of all genomes that are being compared, and naturally clusters the isolates according to their reference strain. Isolates that are too distant from any genome in the database are added to the list of candidate references, thereby creating a new cluster. Short insertions and deletions are considered in addition to single nucleotide polymorphisms for increased discriminatory power. The program has been successfully deployed at the New York State Department of Health's Wadsworth Center to track outbreaks of L. pneumophila as well as other bacterial species.

### Required software:
- Linux/Unix/Mac operating system
- Docker
- Python 3.7 or higher
- the Python packages pydot, numpy and matplotlib, which can be installed with:   pip install <package-name>

### Installing LegioCluster (Ubuntu):
1) download the LegioCluster-master.zip file from GitHub: Code > Download ZIP
2) copy the LegioCluster-master.zip to a designated working directory
3) unzip the file: unzip LegioCluster-master.zip
4) run the install program, which will generate needed folders and links and change the name of the LegioCluster-master folder to LegioCluster: 
	bash LegioCluster-master/install.sh

### Tutorial:
The LegioCluster program contains a Tutorial folder. It is highly recommended that new users go through the Tutorial.pdf file to familiarize themselves with the program.

### running LegioCluster:
- LegioCluster.py  starts with an interactive menu to properly format all required input data into an input TXT file. The best way to start LegioCluster. 
- LegioCluster_main.py <input.txt>  can be used to start LegioCluster when a properly formatted input TXT file is already present in the input/ folder.

### Working with species other than Legionella:
The LegioCluster/ folder contains all files and folders needed to analyze Illumina reads from Legionella pneumophila. You can use the following python script to add other species to the pipeline (see Tutorial for details):
- add_species.py  allows users to modify the pipeline to analyze additional bacterial species. Requires at least one reference genome (FASTA file) for that species.

### Authors:
Wolfgang Haas, Pascal Lapierre, and Kimberlee A. Musser

### Author Affiliations: 
Wadsworth Center, New York State Department of Health, Albany NY, USA 

### Acknowledgements: 
The code development was partially supported by the Public Health Emergency Preparedness grant number U9OTP216988, funded by Centers for Disease Control and Prevention as well as by Wadsworth Center, New York State Department of Health. Its contents are solely the responsibility of the authors and do not necessarily represent the official views of the Wadsworth Center or the New York State Department of Health.

### License: 
GPL 3; https://www.gnu.org/licenses/gpl-3.0.en.html


