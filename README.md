# ChIP-Seq-Peak

# Description
ChIP-Seq-Peak is a program that takes as an input ChIP-Seq data and generates 25 bp bins. Anytime a bin is not a part of an 8 consecutive bin chain (consecutive 200 bp) they are removed. Next, bins that have higher signal than both the immediately adjacent downstream and upstream bins are considered potential peaks. This creates a large list of peaks that has many false positives. To remove these false peaks, the minimal bin between two chosen bins is first found. Then, if the height of the minimum bin is less than 60% of the smaller of the two bins, the two bins are well-separated and utilized for peak calling. If not, the second peak is removed and the program searches for the next minimum between the original bin and the next bin it finds greater than its adjacent bins. This continues until peaks meet the criteria or if the distance between the two peaks is greater than 200. If so, the two bins are utilized for peak calling. Next, to best characterize the called peaks, a normal distribution is generated from the raw data utilizing peak height, position of the center of the peak, and assuming a 400 bp width. The program optimizes these parameters around the called peak and then draws a 10 bp line over the peak of this normal distribution with a height equal to the number of reads in the peak. 

# Usage
This program takes as an input a bedfile, converts to wig format, and subsequently finds peaks. 

The command line format is as follows:
  
python3 ChIP-Seq-Peak.py win genome inputfile "wig_file_name"

# Parameter description

**win**: should be set to 25. Modification of window size will require modificaiton of program.

**genome**: hg38, hg19

**inputfile**: Bedfile of fragment positions

"**wig_file_name**": name of output file




