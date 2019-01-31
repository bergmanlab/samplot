#!/usr/bin/env python
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.gridspec as gridspec
import pylab
import pysam
import os
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import argparse
from matplotlib.offsetbox import AnchoredText
import matplotlib.ticker as ticker

#pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples
CIGAR_MAP = { 
    'M' : 0,
    'I' : 1,
    'D' : 2,
    'N' : 3,
    'S' : 4,
    'H' : 5,
    'P' : 6,    
    '=' : 7,
    'X' : 8,
    'B' : 9
}


def main():
    options = setup_arguments()
    
    # set up plot 
    plot_height,plot_width,window = set_plot_dimensions(options.start, 
                                                        options.end, 
                                                        options.plot_height, 
                                                        options.plot_width, 
                                                        options.bams, 
                                                        options.window)
    
    range_min = max(0,int(options.start) - window)
    range_max = int(options.end) + window

    # set up sub plots
    matplotlib.rcParams.update({'font.size': 12})
    fig = matplotlib.pyplot.figure(figsize=(plot_width, plot_height), dpi=300)

    # read alignment data
    read_data = get_read_data(  options.chrom, 
                                options.start, 
                                options.end, 
                                options.bams, 
                                options.reference, 
                                options.min_mqual, 
                                window,
                                options.minq)

    # set up grid organizer
    grid,num_ax = create_gridspec(options.bams, read_data)
    current_axis_idx = 0

    # Plot each sample
    current_axis_idx = plot_samples(read_data, 
                                    grid, 
                                    current_axis_idx, 
                                    num_ax, 
                                    options.bams, 
                                    options.chrom,  
                                    options.coverage_tracktype, 
                                    options.titles,  
                                    options.xaxis_label_fontsize, 
                                    options.yaxis_label_fontsize, 
                                    range_min,
                                    range_max,
                                    options.normalize_cov,
                                    options.add_hline)

    plot_legend(fig, options.legend_fontsize, options.minq, options.add_hline)

    # save
    matplotlib.rcParams['agg.path.chunksize'] = 100000
    matplotlib.pyplot.tight_layout(pad=0.8,h_pad=.00001, w_pad=.00001)
    matplotlib.pyplot.savefig(options.output_file)


    


def setup_arguments():
    """Defines the allowed arguments for samplot
    """
    parser = argparse.ArgumentParser(
            description="SAMPLOT creates images of genome regions from CRAM/SAM alignments, "+\
                    "optimized for structural variant call review")

    parser.add_argument("-n",
                      "--titles",
                      help="Space-delimited list of plot titles. "+\
                              "Use quote marks to include spaces (i.e. \"plot 1\" \"plot 2\")",
                      type=str,
                      nargs="+",
                      required=False);

    parser.add_argument("-r",
                      "--reference",
                      help="Reference file for CRAM, required if CRAM files used",
                      type=str,
                      required=False);

    parser.add_argument("-z",
                      "--z",
                      type=int,
                      default=4,
                      help="Number of stdevs from the mean (default 4)",
                      required=False)

    parser.add_argument("-b",
                      "--bams",
                      type=str,
                      nargs="+",
                      help="Space-delimited list of BAM/CRAM file names",
                      required=True)

    parser.add_argument("-o",
                      "--output_file",
                      type=str,
                      help="Output file name",
                      required=True)

    parser.add_argument("-s",
                      "--start",
                      type=int,
                      help="Start position of region/variant",
                      required=True)

    parser.add_argument("-e",
                      "--end",
                      type=int,
                      help="End position of region/variant",
                      required=True)

    parser.add_argument("-c",
                      "--chrom",
                      type=str,
                      help="Chromosome",
                      required=True)

    parser.add_argument("-w",
                      "--window",
                      type=int,
                      help="Window size (count of bases to include in view), default(0.5 * len)",
                      required=False)

    parser.add_argument("--minq",
                      type=int,
                      help="coverage from reads with MAPQ <= minq plotted in lighter grey."+\
                              " To disable, pass in negative value",
                      default=0,
                      required=False)



    parser.add_argument("--coverage_tracktype",
                      type=str,
                      help="type of track to use for low MAPQ coverage plot.",
                      choices=['stack','superimpose'],
                      default="stack",
                      required=False)


    parser.add_argument("-H",
                      "--plot_height",
                      type=int,
                      help="Plot height",
                      required=False)

    parser.add_argument("-W",
                      "--plot_width",
                      type=int,
                      help="Plot width",
                      required=False)

    parser.add_argument("-q",
                      "--min_mqual",
                      type=int,
                      help="Min mapping quality of reads to be included in plot",
                      required=False)



    parser.add_argument("--xaxis_label_fontsize",
                      type=int,
                      default=6,
                      help="Font size for X-axis labels (default 6)",
                      required=False)

    parser.add_argument("--yaxis_label_fontsize",
                      type=int,
                      default=6,
                      help="Font size for Y-axis labels (default 6)",
                      required=False)

    parser.add_argument("--legend_fontsize",
                      type=int,
                      default=6,
                      help="Font size for legend labels (default 6)",
                      required=False)

    parser.add_argument("--normalize_cov",
                type=float,
                help="Divides raw coverage by given value to plot relative coverage",
                required=False)

    parser.add_argument("--add_hline",
                type=float,
                help="Adds a horizontal line to plot at given y-axis value",
                required=False)

    options = parser.parse_args()

    for bam in options.bams:
        if ".cram" in bam:
            if not options.reference:
                parser.print_help(sys.stderr)
                sys.exit("Error: Missing reference for CRAM")
    return options

def set_plot_dimensions(start, end, arg_plot_height, arg_plot_width, bams, arg_window):
    """Chooses appropriate dimensions for the plot

    Includes the number of samples, whether a variant type is included, and any annotations in height
    Includes the start, end, and window argument in width
    If height and width are chosen by used, these are used instead

    Return plot height, width, and window as integers
    """

    plot_height = 5
    plot_width = 8
    if arg_plot_height:
        plot_height = arg_plot_height
    else:
        num_subplots = len(bams)
        plot_height = 2 + num_subplots

    if arg_plot_width:
        plot_width = arg_plot_width

    # if an SV type is given, then expand the window around its bounds
    window = 0
    if arg_window:
        window = arg_window

    return plot_height,plot_width,window

def add_coverage(read, coverage, minq):
    """Adds a read to the known coverage 
    
    Coverage from Pysam read is added to the coverage list
    Coverage list is a pair of high- and low-quality lists
    Quality is determined by minq, which is min quality
    """

    hp = 0

    if read.has_tag('HP'):
        hp = int(read.get_tag('HP'))

    if hp not in coverage:
        coverage[hp] = {}

    curr_pos = read.reference_start
    if not read.cigartuples: return

    for op,length in read.cigartuples:
        if op in [CIGAR_MAP['M'], CIGAR_MAP['='], CIGAR_MAP['X']]:
            for pos in range(curr_pos, curr_pos+length+1):
                if pos not in coverage[hp]:
                    coverage[hp][pos] = [0,0]

                #the two coverage tracks are [0] high-quality and [1] low-quality
                if read.mapping_quality > minq:
                    coverage[hp][pos][0] += 1
                else:
                    coverage[hp][pos][1] += 1
            curr_pos += length
        elif op == CIGAR_MAP['I']:
            curr_pos = curr_pos
        elif op == CIGAR_MAP['D']:
            curr_pos += length
        elif op == CIGAR_MAP['N']:
            curr_pos = length
        elif op == CIGAR_MAP['S']:
            curr_pos = curr_pos
        elif op == CIGAR_MAP['H']:
            curr_pos = curr_pos
        else:
            curr_pos += length

def get_read_data(chrom, start, end, bams, reference, min_mqual, window, minq):
    """Reads alignment files to extract reads for the region

    Region and alignment files given with chrom, start, end, bams
    If CRAM files are used, reference must be provided
    Reads with mapping quality below min_mqual will not be retrieved
    If coverage_only, reads are not kept and used only for checking coverage
    Reads longer than long_read_length will be treated as long reads
    Max coverages values will be set to same value for all samples if same_yaxis_scales
    If max_depth, only max_depth reads will be retrieved, although all will be included in coverage
    If PairedEnd read insert size is greater than z_score standard deviations from mean, read will be treated as discordant
    """

    all_coverages = []

    range_min = max(0,int(start) - window)
    range_max = int(end) + window

    for bam_file_name in bams:
        bam_file = None
        if not reference:
            bam_file = pysam.AlignmentFile(bam_file_name, "rb")
        else:
            bam_file = pysam.AlignmentFile(bam_file_name, \
                                           "rc", \
                                           reference_filename=reference)

        coverage = {}
        for read in bam_file.fetch(chrom,
                                   max(0,range_min-1000), 
                                   range_max+1000):
            if min_mqual and int(read.mapping_quality) < min_mqual:
                continue
            
            add_coverage(read, coverage, minq)

        all_coverages.append(coverage)
    
    read_data = {
            "all_coverages":    all_coverages,
        }

    
    return read_data

def create_gridspec(bams, read_data):
    """Helper function for creation of a correctly-sized GridSpec instance
    """
    # give one axis to display each sample
    num_ax = len(bams)

    # set the relative sizes for each
    ratios = []
    
    for i in range(len(bams)):
        ratios.append( len(read_data['all_coverages'][i]) * 3 )
    
    return gridspec.GridSpec(num_ax, 1, height_ratios = ratios), num_ax

def set_haplotypes(curr_coverage):
    """Creates a list to manage counting haplotypes for subplots
    """
    hps = sorted(curr_coverage.keys(), reverse=True)
    #if there are multiple haplotypes, must have 0,1,2
    if len(hps) > 1 or (len(hps) == 1 and hps[0] != 0):
        if 0 not in hps:
            hps.append(0)
        if 1 not in hps:
            hps.append(1)
        if 2 not in hps:
            hps.append(2)
    elif 0 not in hps:
        hps.append(0)
    hps.sort(reverse=True)
    return hps

def plot_coverage(coverage, 
                ax, 
                range_min, 
                range_max, 
                hp_count,
                max_coverage, 
                tracktype,
                yaxis_label_fontsize, normalize, hline):
    """Plots high and low quality coverage for the region

    User may specify a preference between stacked and superimposed 
        superimposed may cause unexpected behavior if low-quality depth is greater than high 
    """

    cover_x = []
    cover_y_lowqual = []
    cover_y_highqual = []
    cover_y_all = []

    for pos in range(range_min,range_max+1):
        if pos in coverage:
            cover_x.append(\
                    float(pos-range_min)/float(range_max - range_min))
            cover_y_all.append(coverage[pos][0] + coverage[pos][1])
            cover_y_highqual.append(coverage[pos][0])
            cover_y_lowqual.append(coverage[pos][1])
        else:
            cover_x.append(\
                    float(pos-range_min)/float(range_max - range_min))
            cover_y_lowqual.append(0)
            cover_y_highqual.append(0)
            cover_y_all.append(0)

        print(pos,cover_y_all[len(cover_y_all)-1], sep="\t")

    cover_y_lowqual = np.array(cover_y_lowqual)
    cover_y_highqual = np.array(cover_y_highqual)
    cover_y_all = np.array(cover_y_all)

    if normalize is not None:
        cover_y_lowqual = cover_y_lowqual/normalize
        cover_y_highqual = cover_y_highqual/normalize
        cover_y_all = cover_y_all/normalize

    if max_coverage > 0:
        max_plot_depth = max_coverage
    elif cover_y_all.max() > 3 * cover_y_all.mean():
        max_plot_depth = max(np.percentile(cover_y_all, 99.5), np.percentile(cover_y_all, 99.5))
    else:
        max_plot_depth = max(cover_y_all.max(), cover_y_all.max())
    ax2 = ax.twinx()
    ax2.set_xlim([0,1])
    
    if 0 == max_plot_depth:
        max_plot_depth = 0.01

    ax2.set_ylim([0,max_plot_depth])
    bottom_fill = np.zeros(len(cover_y_all))
    if tracktype == "stack":
        ax2.fill_between(cover_x, \
                        cover_y_highqual, \
                        bottom_fill,\
                        color='darkgrey',
                        step="pre",
                        alpha=.4)

        ax2.fill_between(cover_x, \
                        cover_y_all, \
                        cover_y_highqual,
                        color='grey',
                        step="pre",
                        alpha=0.15)

        
    elif tracktype == "superimpose": 
        ax2.fill_between(cover_x, \
                        cover_y_lowqual, \
                        bottom_fill,\
                        color='grey',
                        step="pre",
                        alpha=.15)


        ax2.fill_between(cover_x, \
                        cover_y_highqual, \
                        cover_y_lowqual,\
                        color='darkgrey',
                        step="pre",
                        alpha=.4)

        ax2.fill_between(cover_x, \
                        cover_y_lowqual, \
                        bottom_fill,
                        color='grey',
                        step="pre",
                        alpha=0.15)
    
    #number of ticks should be 6 if there's one hp, 3 otherwise
    tick_count = 5 if hp_count==1 else 2
    tick_count = max(int(max_plot_depth/tick_count), 1)

    # set axis parameters
    ax2.yaxis.set_major_locator(ticker.MultipleLocator(tick_count))
    ax2.tick_params(axis='y', colors='grey', labelsize=yaxis_label_fontsize)
    ax2.spines['top'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.tick_params(axis='x',length=0)
    ax2.tick_params(axis='y',length=0)

    if hline is not None:
        ax2.axhline(y=hline, color='black',alpha=0.5)

    return ax2

def plot_samples(   read_data, 
                    grid, 
                    ax_i, 
                    number_of_axes, 
                    bams, 
                    chrom, 
                    coverage_tracktype, 
                    titles, 
                    xaxis_label_fontsize, 
                    yaxis_label_fontsize,
                    range_min,
                    range_max, 
                    normalize,
                    hline):



    ### SUB ROUTINES ####

    """Plots all samples
    """
    # max_insert_size = 0
    for i in range(len(bams)):
        ax =  matplotlib.pyplot.subplot(grid[ax_i])
        hps = set_haplotypes(read_data['all_coverages'][i])
        inner_axs = gridspec.GridSpecFromSubplotSpec(len(hps), 
                                                     1,
                                                     subplot_spec=grid[ax_i],
                                                     wspace=0.0,
                                                     hspace=0.5)
        max_coverage = 0
        axs = {}
        for j in range(len(hps)):
            axs[j] = matplotlib.pyplot.subplot(inner_axs[hps[j]])

        cover_axs = {}
        # curr_axs = ''
        # overall_coverage_range = []
        for hp in hps:
            curr_ax = axs[hp]

            curr_coverage = []
            if hp in read_data['all_coverages'][i]:
                curr_coverage = read_data['all_coverages'][i][hp]

            cover_ax = plot_coverage(curr_coverage, 
                    curr_ax, 
                    range_min, 
                    range_max,
                    len(hps), 
                    max_coverage, 
                    coverage_tracktype, 
                    yaxis_label_fontsize, normalize, hline)

            cover_axs[hp] = cover_ax

        #{{{ set axis parameters
        #set the axis title to be either one passed in or filename
        curr_ax = axs[hps[0]]

        if titles and len(titles) == len(bams):
            curr_ax.set_title(titles[ax_i], fontsize=8, loc='left')

        else:
            curr_ax.set_title(os.path.basename(bams[ax_i]), \
                        fontsize=8, loc='left')


        if len(axs) > 1:
            for j in axs:
                curr_ax = axs[j]
                fp = dict(size=8, backgroundcolor='white')
                text = 'HP: '
                if j == 0:
                    text += 'Undef'
                else:
                    text += str(j)
                at = AnchoredText(text,
                                  loc=2,
                                  prop=fp,
                                  borderpad=0,
                                  pad=0,
                                  frameon=False)
                curr_ax.add_artist(at)
        
        for j in hps:
            curr_ax = axs[j]
            curr_ax.set_xlim([0,1])
            curr_ax.spines['top'].set_visible(False)
            curr_ax.spines['bottom'].set_visible(False)
            curr_ax.spines['left'].set_visible(False)
            curr_ax.spines['right'].set_visible(False)
            curr_ax.tick_params(axis='y', labelsize=yaxis_label_fontsize)
            #if there's one hp, 6 ticks fit. Otherwise, do 3
            tick_count = 6 if len(hps)==1 else 3
            curr_ax.yaxis.set_major_locator(ticker.LinearLocator(tick_count))
            curr_ax.tick_params(axis='both', length=0)
            curr_ax.set_xticklabels([])
            curr_ax.set_yticklabels([])

        last_sample_num = number_of_axes - 1
        
        if (ax_i == last_sample_num):
            curr_ax = axs[ hps[-1] ]
            labels = [int(range_min + l*(range_max-range_min)) \
                    for l in curr_ax.xaxis.get_majorticklocs()]
            curr_ax.set_xticklabels(labels, fontsize=xaxis_label_fontsize)
            curr_ax.set_xlabel('Chromosomal position on ' + chrom, fontsize=8)

        curr_ax = axs[hps[ int(len(hps)/2)    ]]

        curr_ax.set_yticklabels([])
        
        cover_ax = cover_axs[hps[ int(len(hps)/2)    ]]
        if normalize is None:
            cover_ax.set_ylabel('Coverage', fontsize=8)
        else:
            cover_ax.set_ylabel('Normalized Coverage', fontsize=8)
        cover_ax.yaxis.tick_left()
        cover_ax.yaxis.set_label_position("left")
        #}}}

        ax_i += 1
    return ax_i

def plot_legend(fig, legend_fontsize, minq, hline):
    """Plots the figure legend
    """

    marker_labels = []
    legend_elements = []

    if minq > 0:
        legend_elements += [mpatches.Patch(color='darkgrey',alpha=.4)]
        marker_labels.append("Coverage (MAPQ > "+str(minq)+")")
        legend_elements += [mpatches.Patch(color='grey',alpha=.15)]
        marker_labels.append("Coverage (MAPQ <= "+str(minq)+")")

    elif minq == 0:
        legend_elements += [mpatches.Patch(color='darkgrey',alpha=.4)]
        marker_labels.append("Coverage (MAPQ > "+str(minq)+")")
        legend_elements += [mpatches.Patch(color='grey',alpha=.15)]
        marker_labels.append("Coverage (MAPQ = "+str(minq)+")")

    if hline is not None:
        legend_elements += [matplotlib.pyplot.Line2D([0,0],[0,1],color='black', alpha=0.5)]
        marker_labels.append("Avg. Norm. Cov. ("+format(hline, '.2f')+")")
    
    else:
        legend_elements += [mpatches.Patch(color='darkgrey',alpha=.4)]
        marker_labels.append("Coverage")





    fig.legend( legend_elements ,
                marker_labels, 
                loc=1,
                fontsize = legend_fontsize,
                frameon=False)


if __name__ == '__main__':
    main()