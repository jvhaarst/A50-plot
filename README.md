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
`Seq_Name,Seq_File`

To change the output of the graph, one can set a couple of environment values.
These are :
```
MAX_CONTIGS
MIN_LENGTH
TITLE
TYPE
EXPECTED_GENOME_SIZE
```

So an example looks like this :
```
export TITLE="103"
export MAX_CONTIGS=4000
export MIN_LENGTH=100
python ~/A50-plot/AssemblyStats.py 103_assemblies.csv > 103.stats
```
To parse the output to human readable, one can use this :
`
cat 103.stats | tr ',' '\t' | numfmt --header --field 2-5,7-18,20 --grouping | numfmt --header --field 6,19,21 --format '%.1f' | column -t`

This then results in:

```
Name                      Count  Sum            Max         Min  Average   Median  N50         L50  NG50  LG50  N90        L90  N95        L95  Count>1000  Count>10000  #GC            GC    #N     N
103/assembly.fasta        3,467  3,002,559,922  53,397,245  104  866039.8  6,152   10,128,795  77   0     0     2,373,600  302  1,366,482  382  2,897       1,327        1,045,685,741  34.9  5,100  0.1
103_final/assembly.fasta  3,636  3,006,248,922  40,557,336  104  826801.2  5,803   11,571,573  76   0     0     2,427,154  296  1,406,690  374  3,019       1,323        1,047,101,522  34.9  4,800  0.1
103_pbont/assembly.fasta  3,547  3,007,150,960  38,923,932  100  847801.3  5,785   11,975,311  76   0     0     2,910,579  273  1,394,760  344  2,870       1,297        1,047,384,588  34.9  5,500  0.1

```
As I had not entered the expected genome size, the NG50 and LG50 are zero.

The image  output of the script looks like this :

<img src="https://raw.github.com/jvhaarst/A50-plot/master/input.csv.png" href="https://raw.github.com/jvhaarst/A50-plot/master/input.csv.png"/>
