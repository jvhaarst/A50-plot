These R scripts create metrics and graphs from genome assemblies.

To use, first install the dependancies with `R CMD BATCH install.R`.
After that, edit the `input.csv` file to match your situation.
Then you can run the image and stats generating script with `R CMD BATCH AssemblyStats.R`

The CSV columns are:
Seq_Name Seq_File

The result of the script looks like this :

<img src="https://raw.github.com/jvhaarst/A50-plot/master/Rplots-1.png" href="https://raw.github.com/jvhaarst/A50-plot/master/Rplots-1.png"/>
