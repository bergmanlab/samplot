<center><img src="/doc/imgs/samplot_logo_v5.png" width="300"/></center>

<center><img src="/doc/imgs/montage.jpg" width="100%"/></center>

`samplot` is a command line tool for rapid, multi-sample structural variant
visualization. `samplot` takes SV coordinates and bam files and produces
high-quality images that highlight any alignment and depth signals that
substantiate the SV.

## Usage
```
Usage: samplot.py [options]

the following arguments are required: -b/--bams, -o/--output_file, -s/--start, -e/--end, -c/--chrom

Options:
  -h, --help            show this help message and exit
  --marker_size MARKER_SIZE
                        Size of marks on pairs and splits (default 3)
  -n TITLES [TITLES ...], --titles TITLES [TITLES ...]
                        Space-delimited list of plot titles. Use quote marks
                        to include spaces (i.e. "plot 1" "plot 2")
  -r REFERENCE, --reference REFERENCE
                        Reference file for CRAM, required if CRAM files used
  -z Z, --z Z           Number of stdevs from the mean (default 4)
  -b BAMS [BAMS ...], --bams BAMS [BAMS ...]
                        Space-delimited list of BAM/CRAM file names
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        Output file name
  -s START, --start START
                        Start position of region/variant
  -e END, --end END     End position of region/variant
  -c CHROM, --chrom CHROM
                        Chromosome
  -w WINDOW, --window WINDOW
                        Window size (count of bases to include in view),
                        default(0.5 * len)
  -d MAX_DEPTH, --max_depth MAX_DEPTH
                        Max number of normal pairs to plot
  --minq MINQ           coverage from reads with MAPQ <= minq plotted in
                        lighter grey. To disable, pass in negative value
  -t SV_TYPE, --sv_type SV_TYPE
                        SV type. If omitted, plot is created without variant
                        bar
  -T TRANSCRIPT_FILE, --transcript_file TRANSCRIPT_FILE
                        GFF of transcripts
  -A ANNOTATION_FILES [ANNOTATION_FILES ...], --annotation_files ANNOTATION_FILES [ANNOTATION_FILES ...]
                        Space-delimited list of bed.gz tabixed files of
                        annotations (such as repeats, mappability, etc.)
  --coverage_tracktype {stack,superimpose}
                        type of track to use for low MAPQ coverage plot.
  -a, --print_args      Print commandline arguments
  -H PLOT_HEIGHT, --plot_height PLOT_HEIGHT
                        Plot height
  -W PLOT_WIDTH, --plot_width PLOT_WIDTH
                        Plot width
  -q MIN_MQUAL, --min_mqual MIN_MQUAL
                        Min mapping quality of reads to be included in plot
  -j, --json_only       Create only the json file, not the image plot
  --start_ci START_CI   confidence intervals of SV first breakpoint (distance
                        from the breakpoint). Must be a comma-separated pair
                        of ints (i.e. 20,40)
  --end_ci END_CI       confidence intervals of SV end breakpoint (distance
                        from the breakpoint). Must be a comma-separated pair
                        of ints (i.e. 20,40)
  --long_read LONG_READ
                        Min length of a read to be treated as a long-read
                        (default 1000)
  --min_event_size MIN_EVENT_SIZE
                        Min size of an event in long-read CIGAR to include
                        (default 100)
  --xaxis_label_fontsize XAXIS_LABEL_FONTSIZE
                        Font size for X-axis labels (default 6)
  --yaxis_label_fontsize YAXIS_LABEL_FONTSIZE
                        Font size for Y-axis labels (default 6)
  --legend_fontsize LEGEND_FONTSIZE
                        Font size for legend labels (default 6)
  --annotation_fontsize ANNOTATION_FONTSIZE
                        Font size for annotation labels (default 6)
  --common_insert_size  Set common insert size for all plots
  --hide_annotation_labels
                        Hide the label (fourth column text) from annotation
                        files, useful for region with many annotations
  --coverage_only       Hide all reads and show only coverage
  --same_yaxis_scales   Set the scales of the Y axes to the max of all
  --multifig_vertical   when making multiple plots, combine plots vertically
                        rather than horizontally

```

## Dependencies

* numpy
* matplotlib
* pysam
* statistics
* Pillow

All of these are available from [pip](https://pypi.python.org/pypi/pip).

## Examples: 

Samplot requires either BAM files or CRAM files as primary input. If you use
CRAM, you'll also need a reference genome like one used the the 1000 Genomes
Project
(ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz).

### Basic use case
Using data from NA12878, NA12889, and NA12890 in the 
[1000 Genomes Project](http://www.internationalgenome.org/about), we will
inspect a possible deletion in NA12878 at 4:115928726-115931880 with respect
to that same region in two unrelated samples NA12889 and NA12890.

The following command will create an image of that region:
```
time samplot/src/samplot.py \
    -n NA12878 NA12889 NA12890 \
    -b samplot/test/data/NA12878_restricted.bam \
      samplot/test/data/NA12889_restricted.bam \
      samplot/test/data/NA12890_restricted.bam \
    -o 4_115928726_115931880.png \
    -c chr4 \
    -s 115928726 \
    -e 115931880 \
    -t DEL

real    0m9.450s
user    0m9.199s
sys     0m0.217s
```

The arguments used above are:

`-n` The names to be shown for each sample in the plot

`-b` The BAM/CRAM files of the samples (space-delimited)

`-o` The name of the output file containing the plot

`-c` The chromosome of the region of interest

`-s` The start location of the region of interest

`-e` The end location of the region of interest

`-t` The type of the variant of interest

This will create an image file named `4_115928726_115931880.png`, shown below:

<img src="/doc/imgs/4_115928726_115931880.png">

### Downsampling "normal" pairs

The runtime of `samplot` can be reduced by only plotting a portion of the concordant 
pair-end reads (+/- strand orientation, within z s.d. of the mean insert size where z 
is a command line option the defaults to 4). If we rerun the prior example, but only plot
a random sampling of 100 normal pairs we get a similar result 3.6X faster.

```
time samplot/src/samplot.py \
    -n NA12878 NA12889 NA12890 \
    -b samplot/test/data/NA12878_restricted.bam \
      samplot/test/data/NA12889_restricted.bam \
      samplot/test/data/NA12890_restricted.bam \
    -o 4_115928726_115931880.d100.png \
    -c chr4 \
    -s 115928726 \
    -e 115931880 \
    -t DEL \
    -d 100

real    0m2.621s
user    0m2.466s
sys     0m0.124s
```

<img src="/doc/imgs/4_115928726_115931880.d100.png">


### Gene and other genomic feature annotations

Gene annotations (tabixed, gff3 file) and genome features (tabixed, bgzipped, bed file) can be 
included in the plots.

Get the gene annotations:
```
wget ftp://ftp.ensembl.org/pub/grch37/release-84/gff3/homo_sapiens/Homo_sapiens.GRCh37.82.gff3.gz
bedtools sort -i Homo_sapiens.GRCh37.82.gff3.gz \
| bgzip -c > Homo_sapiens.GRCh37.82.sort.gff3.gz
tabix Homo_sapiens.GRCh37.82.sort.gff3.gz
```

Get genome annotations, in this case Repeat Masker tracks and a mappability track:
```
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDukeMapabilityUniqueness35bp.bigWig
bigWigToBedGraph wgEncodeDukeMapabilityUniqueness35bp.bigWig wgEncodeDukeMapabilityUniqueness35bp.bed
bgzip wgEncodeDukeMapabilityUniqueness35bp.bed
tabix wgEncodeDukeMapabilityUniqueness35bp.bed.gz

curl http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz \
| bgzip -d -c \
| cut -f 6,7,8,13 \
| bedtools sort -i stdin \
| bgzip -c > rmsk.bed.gz
tabix rmsk.bed.gz
```

Plot:
```
samplot/src/samplot.py \
    -n NA12878 NA12889 NA12890 \
    -b samplot/test/data/NA12878_restricted.bam \
      samplot/test/data/NA12889_restricted.bam \
      samplot/test/data/NA12890_restricted.bam \
    -o 4_115928726_115931880.d100.genes_reps_map.png \
    -c chr4 \
    -s 115928726 \
    -e 115931880 \
    -t DEL \
    -d 100 \
    -T Homo_sapiens.GRCh37.82.sort.gff3.gz \
    -A rmsk.bed.gz wgEncodeDukeMapabilityUniqueness35bp.bed.gz

real    0m2.784s
user    0m2.633s
sys 0m0.129s
```

<img src="/doc/imgs/4_115928726_115931880.d100.genes_reps_map.png">

### Generating images from a VCF file
To plot images from all structural variants in a VCF file, use samplot's
`samplot_vcf.sh` script. This accepts a VCF file and the BAM files of samples
you wish to plot, outputting images and related metadata to a directory of your
choosing.

This script is especially useful as part of the 
[SV-plaudit pipeline](https://github.com/jbelyeu/SV-plaudit) and creates
metadata files for
all images which SV-plaudit requires.

```
samplot/src/samplot_vcf.sh \
    -o output_dir \
    -B $HOME/bin/bcftools \
    -S samplot/src/samplot.py \
    -v samplot/test/data/NA12878.trio.svt.subset.vcf \
    samplot/test/data/NA12878_restricted.bam \
    samplot/test/data/NA12889_restricted.bam \
    samplot/test/data/NA12890_restricted.bam
```
The arguments used above are:

`-o` output directory (make this directory before executing)

`-B` Executable file of [bcftools](https://samtools.github.io/bcftools/)

`-S` samplot.py script

`-v` VCF file with variants to plot


#### CRAM inputs
Samplot also support CRAM input, which requires a reference fasta file for
reading as noted above. Notice that the reference file is not included in this
repository due to size. This time we'll plot an interesting duplication at
X:101055330-101067156.

```
samplot/src/samplot.py \
    -n NA12878 NA12889 NA12890 \
    -b samplot/test/data/NA12878_restricted.cram \
      samplot/test/data/NA12889_restricted.cram \
      samplot/test/data/NA12890_restricted.cram \
    -o cramX_101055330_101067156.png 
    -c chrX \
    -s 101055330 \
    -e 101067156 \
    -t DUP \
    -r hg19.fa
```

The arguments used above are the same as those used for the basic use case, with the addition of the following:

`-r` The reference file used for reading CRAM files

#### Plotting without the SV 
Samplot can also plot genomic regions that are unrelated to an SV. If you do
not pass the SV type option (`-t`) then the top SV bar will go away and only
the region that is given by `-c` `-s` and `-e` will be displayed.

#### Long read (Oxford nanopore and PacBio) and linked read support
Any alignment that is longer than 1000 bp are treated as a longread, and
the plot design will focus on aligned regions and gaps. Aligned regions are in orange, and gaps follow the same DEL/DUP/INV color code used for short reads. The height of the alignment is based on the size of its largest gap.

<img src="/doc/imgs/longread_del.png">

If the bam file has an MI tag, then the reads will be treated as linked reads.
The plots will be similar to short read plots, but all alignments with the same MI is plotted at the same height according to alignment with the largest gap in the group. A green line connects all alignments in a group.

<img src="/doc/imgs/linkedread_del.png">
