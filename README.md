These Python scripts create metrics and graphs from genome assemblies.

To use, first install the dependancies with `conda install biopython matplotlib`.
After that, edit the `input.csv` file to match your situation.
Then you can run the image and stats generating script with:

`python AssemblyStats.py input.fasta`
or
`python AssemblyStats.py input.fasta input.fasta.gz input.fasta.bz2`
or
`python AssemblyStats.py input.csv`



The CSV (comma separated) columns are:
Seq_Name Seq_File

Spaces are allowed.
The order in the file will define the order in the graphics legend.
You can comment lines with #
First line HAS to be the title:
Seq_Name,Seq_File

To change the output of the graph, one can set a couple of environment values.
These are :
```
MAX_CONTIGS
MIN_LENGTH
TITLE
TYPE
EXPECTED_GENOME_SIZE

```
The result of the script looks like this :

<img src="https://raw.github.com/jvhaarst/A50-plot/master/input.csv.png" href="https://raw.github.com/jvhaarst/A50-plot/master/input.csv.png"/>
