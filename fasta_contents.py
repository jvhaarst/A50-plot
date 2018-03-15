#!/usr/bin/env python
#import sys
#sys.path.append('/home/assembly/lib/Cogent/build/lib.linux-x86_64-2.7/cogent/parse')
#sys.path.append('/home/assembly/lib/Cogent/build/lib.linux-x86_64-2.7/')

from __future__ import division
from sys import argv, exit
import sys
sys.path.append('/home/assembly/lib/Cogent/build/lib.linux-x86_64-2.7/cogent/parse')
sys.path.append('/home/assembly/lib/Cogent/build/lib.linux-x86_64-2.7/')

from cogent.parse.fasta import *
from numpy import sum, array, cumsum, concatenate, median
#from pylab import *


def calc_N_values(lengths):

    num_bases = sum(lengths)
    data = array(lengths)
    data.sort()
    data = data[::-1]
    cum = cumsum(data)

    N25_lookup = num_bases*0.25
    N50_lookup = num_bases*0.5
    N75_lookup = num_bases*0.75
    N95_lookup = num_bases*0.95

    lookups = [('N25', num_bases*0.25),\
        ('N50', num_bases*0.5),\
        ('N75', num_bases*0.75),\
        ('N95', num_bases*0.95)]

    res = []
    for lookup_label, lookup_val in lookups:
        sub_total = 0
        for idx, sub_total in enumerate(cum):
            if sub_total >= lookup_val:
                res.append((lookup_label, lookup_val, idx, data[idx]))
                #res.append((idx, data[idx], sub_total, lookup_val))
                #print 'IDX', idx
                #print 'LEN', data[idx]
                break
    return res

def get_len_stats(lengths):
    res = {}
    res['numseq'] = len(lengths)
    res['numbps'] = sum(lengths)
    res['avglen'] = lengths.mean()
    res['stdlen'] = lengths.std()
    res['maxlen'] = max(lengths)
    res['minlen'] = min(lengths)
    res['medlen'] = float(median(lengths))
    return res

"""
def plot_lengths(data, one_bar_per=1000, n_vals=None):
    if n_vals is None:
        n_vals = calc_N_values(data)

    clf()
    fig = figure()
    data = array(data)/1000000 # Turn into Mb
    data.sort()
    data = data[::-1]
    cum = cumsum(data)
    xs = concatenate([arange(1,50),arange(50,len(cum),one_bar_per)])
    #xs_plot = []
    #for i in xs:
    #    if i == 0:
    #        xs_plot.append(i)
    #    else:
    #        xs_plot.append(log10(i))
    #xs_plot = array(xs_plot)
    ys_cum = []
    for idx in xs:
        ys_cum.append(cum[idx])
    plot(xs, ys_cum, color='k', linewidth=4)
    #bar(array(xs)-(.3*one_bar_per), ys_cum, width=0.6*one_bar_per,\
    #    color='0.75')

    text_items = []
    text_labels = ['N25','N50','N75','N95']
    colors=['c','b','g','r']
    x=0
    for n_idx, n_len, n_cum, n_lookup in n_vals:
        n_idx = n_idx
        y_val = n_lookup/1000000
        if x == 3:
            plot([n_idx,n_idx],[0,y_val], color=colors[x], linewidth=3)
            plot([0,n_idx],[y_val,y_val], color=colors[x], linewidth=3)
        text_str = "%s (%.1f Mb): %s, length: %.1f Kb"\
            %(text_labels[x], y_val, n_idx, n_len/1000)
        text_items.append(text_str)
        #axvline(n_idx, color='r')
        #axhline(n_cum/1000000, color='r')
        x += 1
    
    ylim(0,950)
    yticks(arange(50,1000,100))
    #xlim(xmin=0)
    xlim(0,2000)
    ylabel("Cumulative number of bases (Mb)")
    xlabel("Contig/Scaffold")

    for x in range(4):
        text(len(data)/3, 200-(x/20*900), "%s"%\
            (text_items[x]), color=colors[x])
    title("Length distribution of contigs/scaffolds")

    return fig
"""

if __name__ == "__main__":

    KEYS = ['numseq', 'numbps','avglen','stdlen','minlen','maxlen','medlen','n25idx','n25len','n50idx','n50len','n75idx','n75len','n95idx','n95len','counts']

    result = {}
    file_storage = []
    len_storage = []
    counts_storage = []

    for fna_file in argv[1:]:
        file_storage.append(fna_file)

        lengths = []
        raw_counts = {}
	raw_counts['N']=0
        for label, read in MinimalFastaParser(open(fna_file)):
            lengths.append(len(read))
            for s in read:
                try:
                    raw_counts[s] +=1
                except KeyError:
                    raw_counts[s] = 1
                    
                #s = s.upper()
                #if s not in counts:
                #    counts[s] = 0
                #counts[s] += 1
            
            #if len(lengths) == 100:
            #    break
        counts = {}
        for k in raw_counts:
            if k.upper() not in counts:
                counts[k.upper()] = 0
            counts[k.upper()] += raw_counts[k]
        lengths = array(lengths)
        
        len_storage.append(lengths)
        counts_storage.append(counts)

    for fn, lengths, counts in zip(file_storage, len_storage, counts_storage):
        len_stats = get_len_stats(lengths)
        N_vals = calc_N_values(lengths)
        result[fn] = {}
        result[fn].update(len_stats)
        for lookup_label, lookup_val, item_idx, item_len in N_vals:
            result[fn][lookup_label+'idx'] = item_idx
            result[fn][lookup_label+'len'] = item_len
        result[fn]['counts'] = counts


    
    file_items = [""]
    for idx, fn in enumerate(file_storage):
        file_items.append(fn)
    print '\t'.join(file_items)
    #for idx, fn in enumerate(file_storage):
    #    print "Column %s: %s"%(idx+1,fn)
    
    #print ""
    print 'Number of sequences:\t%s'%('\t'.join(["%s"%(result[fn]['numseq']) for fn in file_storage]))
    print 'Total sequence length:\t%s'%('\t'.join(["%s"%(result[fn]['numbps']) for fn in file_storage]))
    print 'Average sequence length:\t%s'\
        %('\t'.join(["%.2f"%(result[fn]['avglen']) for fn in file_storage]))
    print 'Std. dev. sequence length:\t%s'\
        %('\t'.join(["%.2f"%(result[fn]['stdlen']) for fn in file_storage]))
    print 'Minimum sequence length:\t%s'\
        %('\t'.join(["%s"%(result[fn]['minlen']) for fn in file_storage]))
    print 'Maximum sequence length:\t%s'\
        %('\t'.join(["%s"%(result[fn]['maxlen']) for fn in file_storage]))
    print 'Median sequence length:\t%s'\
        %('\t'.join(["%.2f"%(result[fn]['medlen']) for fn in file_storage]))
    print 'N25 sequence index:\t%s'\
        %('\t'.join(["%s"%(result[fn]['N25idx']) for fn in file_storage]))
    print 'N25 sequence length:\t%s'\
        %('\t'.join(["%s"%(result[fn]['N25len']) for fn in file_storage]))
    print 'N50 sequence index:\t%s'\
        %('\t'.join(["%s"%(result[fn]['N50idx']) for fn in file_storage]))
    print 'N50 sequence length:\t%s'\
        %('\t'.join(["%s"%(result[fn]['N50len']) for fn in file_storage]))
    print 'N75 sequence index:\t%s'\
        %('\t'.join(["%s"%(result[fn]['N75idx']) for fn in file_storage]))
    print 'N75 sequence length:\t%s'\
        %('\t'.join(["%s"%(result[fn]['N75len']) for fn in file_storage]))
    print 'N95 sequence index:\t%s'\
        %('\t'.join(["%s"%(result[fn]['N95idx']) for fn in file_storage]))
    print 'N95 sequence length:\t%s'\
        %('\t'.join(["%s"%(result[fn]['N95len']) for fn in file_storage]))

    all_count_keys = []
    for fn in file_storage:
        all_count_keys.extend(result[fn]['counts'].keys())
    all_count_keys = list(set(all_count_keys))
    for k in all_count_keys:
        items = []
        for fn in file_storage:
            try:
                items.append("%s (%.2f%%)"%(result[fn]['counts'][k],\
                    result[fn]['counts'][k]/result[fn]['numbps']*100))
            except KeyError:
                items.append("--")
        print '%s content:\t%s'%(k, '\t'.join(items))
    gc_items = []
    for fn in file_storage:
        count_gc = result[fn]['counts']['G'] + result[fn]['counts']['C']
        count_tgacN = result[fn]['counts']['G'] + result[fn]['counts']['C'] +\
            result[fn]['counts']['A'] + result[fn]['counts']['T'] + result[fn]['counts']['N']
        perc_gc = (count_gc/count_tgacN)*100
        gc_items.append("%.2f%%"%(perc_gc))
    print 'GC content:\t%s'%('\t'.join(gc_items))


