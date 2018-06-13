# Demo Workflow: from Ion torrent fastq files to SNP frequencies

2018-06-13 by staehlo  
\
For training reasons, I built a bash workflow that takes the fastq reads from a study run on an Ion Torrent *Personal Genome Machine* (PGM) and transforms the data to called variants. The fastq files were taken from study PRJEB20222 in the European Nucleotide Archive. You can find below the link to the related publication by Lugowska et al. (2018). All the software used in this workflow is open source (except the liftOver function provided by UCSC which is free for academics but billable for companies).  
\
Lugowska et al. used the *Ion AmpliSeq Cancer Hotspot Panel v2* for their study. The Browser Extensible Data (BED) file of this panel is not publicly available. So I had to use a trick: The European Nucleotide Archive of the study contains the generated BAM files, too. The genomic coordinates of the AmpliSeq Panel were approximated based on the genomic positions of one of the BAM files.


## Data source and description

Original publication:  
Lugowska I, Teterycz P, Mikula M, Kulecka M, Kluska A, Balabas A, Piatkowska M, Wagrodzki M, Pienkowski A, Rutkowski P, Ostrowski J. *IDH1/2 Mutations Predict Shorter Survival in Chondrosarcoma.* J Cancer. 2018 Feb 28;9(6):998-1005

doi: [10.7150/jca.22915](https://dx.doi.org/10.7150%2Fjca.22915)  
PubMed PMID: 29581779  
PubMed Central PMCID: PMC5868167

European Nucleotide Archive  
Accession number: PRJEB20222  
<https://www.ebi.ac.uk/ena/data/view/PRJEB20222>


## Software information

Operating system:

	$ lsb_release -a  
	No LSB modules are available.  
	Distributor ID:	Ubuntu  
	Description:	Ubuntu 16.04.4 LTS  
	Release:	16.04  
	Codename:	xenial
	
Bash:

	$ bash --version  
	GNU bash, version 4.3.48(1)-release (x86_64-pc-linux-gnu)  
	Copyright (C) 2013 Free Software Foundation, Inc.  
	License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>


## Acknowledgments

Vince Buffalo - <https://github.com/vsbuffalo> and his book *"Bioinformatics Data Skills"*
