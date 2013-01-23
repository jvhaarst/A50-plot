These R scripts create metrics and graphs from genome assemblies.

To use, first install the dependancies with `R CMD BATCH install.R`.
After that, edit the `input.csv` file to match your situation.
Then you can run the image and stats generating script with `R CMD BATCH AssemblyStats.R`



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



The result of the script looks like this :

<img src="https://raw.github.com/jvhaarst/A50-plot/master/Rplots-1.png" href="https://raw.github.com/jvhaarst/A50-plot/master/Rplots-1.png"/>
