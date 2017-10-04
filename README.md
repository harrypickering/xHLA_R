xHLA_ARRRRRR
================================================================
The Human Leukocyte Antigen (HLA) gene complex on human chromosome
6 is one of the most polymorphic regions in the human
genome, and contributes in large part to the diversity of the immune
system. Accurate typing of HLA genes with short read sequencing
data has historically been difficult due to the sequence similarity between
the polymorphic alleles.  xHLA iteratively refines the mapping results at
the amino acid level to achieve 99 to 100% 4-digit typing accuracy for both
class I and II HLA genes, taking only about 3 minutes to process a 30X
whole genome BAM file on a desktop computer.


Background to this fork
------------
We had issues getting the main repo for xHLA to work via the main docker-based repo, so we performed a quick and dirty stripdown of the code to make it work for us. The original scripts have been hacked slightly and we've added some control scripts that should be easy to run via R or R studio. 

This is intended to run as a batch mode process that will perform HLA typing on all of the fastq files in the input directory. At the end it generates a tidy table of results data. Along the way it deletes all the files except for the final json files. Obviously some of those files might be useful to some people, but for our purposes we are happy just to save disk space and keep the final typing data. To change this behaviour just edit the R scripts to take out the rm command calls to the system.

We have found that the typing algorithm works reproducibly given a sample of 15000 reads from MiSeq 2*300 data that was generated from exons amplified using the methods described by Lange, V. et al BMC Genomics. 2014; 15: 63. http://doi.org/10.1186/1471-2164-15-63

The R script calls seqtk to pull out a random sample of reads from each of the two read files. The same random seed is used to pull read 1 and read 2 so that the reads remain paired. You can change the number of reads that are pulled by changing the number in the two lines beneath "#randomly sample n reads from the fastq files" in the 0000_super_master_controller_batchmode.R file

Installation
------------

Download this repo to your computer. 

Ensure that you have installed the diamond aligner and seqtk (sequence toolkit).
Install the R packages "jsonlite" and "lpSolve"

Running 
------------
Put your paired end read fastq files in the input folder
Run the 0000_super_master_controller_batchmode.R script.


The 0000_super_master_controller_batchmode.R script will read the FASTQ files, sample 15000 reads using seqtk, then forward the data to 000_master_control_script.sh. This bash script runs bwa, samtools and the xHLA typer script. 

Output 1 is tab delimited text file *000_alldata_out.txt* of all the data you've pushed through the system. New runs are appended to this file so be aware that if you re-run the test on the same pair of fastq files, you will have two lines matching that specimen.

```
sampleID	analysis.date	report.version	report_type	hla.Ai	hla.Aii	hla.Bi	hla.Bii	hla.Ci	hla.Cii	hla.DPB1i	hla.DPB1ii	hla.DQB1i	hla.DQB1ii	hla.DRB1i	hla.DRB1ii
CD1-055_S94a	2017-10-04T09:02:08Z	1.2	hla_typing	A*01:01	A*01:01	B*58:01	B*81:01	C*03:02	C*18:01	DPB1*01:01	DPB1*02:01	DQB1*03:19	DQB1*06:02	DRB1*11:01	DRB1*15:03
CD1-056_S118	2017-10-04T09:14:50Z	1.2	hla_typing	A*24:02	A*68:02	B*18:01	B*27:03	C*02:02	C*07:01	DPB1*04:01	DPB1*04:01	DQB1*02:01	DQB1*05:01	DRB1*01:02	DRB1*03:01
CD1-057_S142	2017-10-04T10:03:12Z	1.2	hla_typing	A*02:05	A*03:01	B*35:01	B*58:01	C*07:01	C*07:05	DPB1*13:01	DPB1*17:01	DQB1*02:01	DQB1*03:01	DRB1*09:01	DRB1*11:02
CD1-058_S166	2017-10-04T10:14:17Z	1.2	hla_typing	A*24:02	A*68:02	B*18:01	B*27:03	C*02:02	C*07:01	DPB1*04:01	DPB1*04:01	DQB1*02:01	DQB1*05:01	DRB1*01:02	DRB1*03:01
```



Output is a JSON file that lists 12 HLA alleles, 2 for each of the HLA genes:

```bash
{
 "sub6444255",
 "creation_time": "2016-05-04T08:25:04Z",
 "report_version": "1.1",
 "report_type": "hla_typing",
 "sample_id": "176444255",
 "hla": {
  "alleles": [
   "A*01:01",
   "A*02:01",
   "B*13:02",
   "B*37:01",
   "C*06:02",
   "C*06:02",
   "DPB1*04:01",
   "DPB1*04:01",
   "DQB1*02:02",
   "DQB1*05:01",
   "DRB1*07:01",
   "DRB1*10:01"
  ]
 }
}
```


Citation
--------
Xie et al. (2017) Fast and accurate HLA typing from short-read next-generation
sequence data with xHLA.
***PNAS***
[doi:10.1073/pnas.1707945114](http://www.pnas.org/content/early/2017/06/27/1707945114)

License
-------
The HLA Typing Software Code (the "Code") is made available by Human
Longevity, Inc. ("HLI") on a non-exclusive, non-sublicensable,
non-transferable basis solely for non-commercial academic research use.
Commercial use of the Code is expressly prohibited.  If you would like to obtain
a license to the Code for commercial use, please contact HLI at
bizdev@humanlongevity.com.  HLI MAKES NO REPRESENTATIONS OR WARRANTIES
WHATSOEVER, EITHER EXPRESS OR IMPLIED, WITH RESPECT TO THE CODE PROVIDED
HEREUNDER. IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR
PURPOSE WITH RESPECT TO CODE ARE EXPRESSLY DISCLAIMED. THE CODE IS FURNISHED
"AS IS" AND "WITH ALL FAULTS" AND DOWNLOADING OR USING THE CODE
IS UNDERTAKEN AT YOUR OWN RISK.  TO THE FULLEST EXTENT ALLOWED BY APPLICABLE
LAW, IN NO EVENT SHALL HLI BE LIABLE, WHETHER IN CONTRACT, TORT, WARRANTY, OR
UNDER ANY STATUTE OR ON ANY OTHER BASIS FOR SPECIAL, INCIDENTAL, INDIRECT,
PUNITIVE, MULTIPLE OR CONSEQUENTIAL DAMAGES SUSTAINED BY YOU OR ANY OTHER PERSON
OR ENTITY ON ACCOUNT OF USE OR POSSESSION OF THE CODE, WHETHER OR NOT
FORESEEABLE AND WHETHER OR NOT HLI HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH
DAMAGES, INCLUDING WITHOUT LIMITATION DAMAGES ARISING FROM OR RELATED TO LOSS OF
USE, LOSS OF DATA, DOWNTIME, OR FOR LOSS OF REVENUE, PROFITS, GOODWILL, BUSINESS
OR OTHER FINANCIAL LOSS.

Any changes made on this fork inherit the terms of the above license.
