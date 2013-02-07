These R scripts create metrics and graphs from genome assemblies.

To use, first install the dependancies with `R CMD BATCH install.R`.
After that, edit the `input.csv` file to match your situation.
Then you can run the image and stats generating script with:

Rscript AssemblyStats.R input.csv



The CSV (comma separated) columns are:
Seq_Name Seq_File

Spaces are allowed.
The order in the file will define the order in the graphics' legend.
You can comment lines with #
First line HAS to be the title:
Seq_Name,Seq_File


You can edit AssemblyStats.R to change the maximum X limit (number of contigs):
max_count <- 20000
Any assembly with more than 20000 contigs will be trimmed to this size and the total sum after this will be shown as 20001.

We have defined a new concept called "lookup value" which is defined as 75% of the harmonic mean of the sequences sizes. Once this values is set, the Nlookup (contig size) and Ilookup (contig index) is calculate for each sequence. Having the harmonic mean among all sizes as a constant lookup value seems more comparable than using 50% of the length of each assembly individualy (N50).
E.G.: assemblies size: 1Mbp, 10Mbp, 20Mbp
Hmean = 1 / (( 1/1M + 1/10M + 1/20M ) / 3) =  2.6Mbp
Mean  =      (   1M +   10M +   20M ) / 3) = 10.3Mbp

Harmonic means, usually, gives a lower value than the Mean, allowing the better analysis of discrepant values easier.
For this dataset, the lookup value would be 75% of 2.6Mbp (1.95Mbp). This fixed value will be used to calculate Nl and Il for the tree individuals.
A line is drawn showing the lookup value (Nl). In the top of the graph tick marks shows the Il. The legend has the values numerically.

The result of the script looks like this :

<img src="https://raw.github.com/jvhaarst/A50-plot/master/input.csv.png" href="https://raw.github.com/jvhaarst/A50-plot/master/input.csv.png"/>
