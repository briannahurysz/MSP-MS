# MSP-MS
Contained here are the scripts used by the O'Donoghue lab at UC San Diego for analysis of multiplex substrate profiling by mass spectrometry data. This data analysis is done in R studio.

Instructions on how to use these scripts are as follows:

To prepare to use the R scripts, the following packages must be downloaded:
data.table
qvalue
fitdistrplus
truncnorm
gtools

These can be installed using command: install.packages("package_name")


1) Obtain protein-peptides.csv files from the label free quantification and PEAKS 18 in the PEAKS software. Label free quantification should be edited to contain no filtering and no normalization. The PEAKS output should have the -10logp >/= adjusted such that under Table 3 the FDR (Peptide Sequences) is </= 1.0%.
2) Rename the protein-peptides.csv exported via the PEAKS 18 export "protein-peptides-FDR.csv"
3) Place both protein-peptides.csv and protein-peptides-FDR.csv into a directory of your choice
4) Open MSP_data_processing_1.R
5) Set the directory in line 3 to the directory containing your data.
6) In line 5, specify the number of mass spec timepoints and the number of replicates (the provided file is set up with 4 separate timepoints, each in quadruplicates)
7) Select the whole script and run it.
8) Access normalyzer via http://normalyzer.immunoprot.lth.se/normalize.php and input the "protein-peptides-for-normalyzer.txt" generated by the R script.
9) Download the output from normalyzer, open the file, and access the "Loess-G-normalized.txt" file. Move this to the same directory you have been working in.
10) Open MSP_data_processing_2.R and repeat steps 5-7.
11) Open MSP_data_processing_3.R and repeat steps 5-7.
12) Open MSP_data_processing_4.R and repeat steps 5-7. 
13) Your output file is protein-peptide_results.csv and can be manipulated further in excel.


Cited:
Chawade, A., Alexandersson, E., & Levander, F. (2014). Normalyzer: A Tool for Rapid Evaluation of Normalization Methods for Omics Data Sets. Journal of Proteome Research, 13(6), 3114–3120. https://doi.org/10.1021/pr401264n![image](https://user-images.githubusercontent.com/112521524/187558274-b6d7a11b-7355-40f5-8cac-c6b77084fa9c.png)

