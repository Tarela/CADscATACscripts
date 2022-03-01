#!/usr/bin/env python2.7
import regions as R
from multiprocessing import Pool
import subprocess
import sys,os
import numpy as np
import math
import argparse
import time

def extractBigwig(cmd,bins):
    result = os.popen(cmd)
    content = result.readlines()

    if content:
        temp = content[0].strip().replace('n/a','0')
        temp = temp.split('\t')
        return temp
    else:
        temp=[0]*int(bins)
        return temp

def getInfo(genome, name, alpha, bwfilename, tss, chromesize, bwsum):
    
    Infos = []
    chromeinfo = {}
    inf = open(chromesize)
    for line in inf:
        line = line.strip().split('\t')
        chrome = line[0]
        size = line[1]
        chromeinfo[chrome] = size

    G = R.interval(genome=genome)
    G.read_bed(tss)
    
    padding = int(1e5)
  
    weight  = np.array( [ 2.0*math.exp(-alpha*math.fabs(z)/1e5)/(1.0+math.exp(-alpha*math.fabs(z)/1e5))  for z in range( -padding,padding+1) ] )
    #print weight
    #print np.sum(weight)    
    #chrom_p = None
    for i,chrom in enumerate(G.chrom):
        center = G.start[i]
        chrom = G.chrom[i]
        bwstart = center - padding
        bwend = center + padding
        if bwstart < 0:
            bwstart = 0
        if bwend > int(chromeinfo[G.chrom[i]]):
            bwend = int(chromeinfo[G.chrom[i]])
        bins = bwend - bwstart + 1
        #print np.shape(weight)
        Infos.append([bwfilename,chrom,bwstart,bwend,bins,weight,name,G.start[i],G.end[i], G.name[i].split(':')[0], G.name[i].split(':')[1],G.strand[i], bwsum])
    #print Infos
    return Infos

def getRP(Info):
    bwfilename = Info[0]
    chrom = Info[1]
    bwstart = Info[2]
    bwend = Info[3]
    bins = Info[4]
    weight = Info[5]
    name = Info[6]
    Gstart = Info[7]
    Gend = Info[8]
    Gname1 = Info[9]
    Gname2 = Info[10]
    Gstrand = Info[11]

    bwsum = Info[12]

    padding = int(1e5)
    center = Gstart
    cmd = '%s {0} {1} {2} {3} {4}'%bwsum
    try:
        cmd = cmd.format(bwfilename,chrom,bwstart,bwend + 1, bins)
        #print cmd
        values = extractBigwig(cmd,bins)
        #print len(values)
        values = np.array(values,dtype=np.float64)
   
        values2 = np.hstack((np.zeros(bwstart - center + padding),values,np.zeros(center + padding - bwend )))
        invalid = np.isnan( values2 )
        values2[ invalid ]   = 0
        #print np.sum(values2)

        return chrom + '\t' + str(Gstart) + '\t' + str(Gend) + '\t' + Gname1 + '\t' + str(np.dot( values2, weight ))  + '\t' + Gname2 + '\t' + Gstrand + '\n'
    except:
        return chrom + '\t' + str(Gstart) + '\t' + str(Gend) + '\t' + Gname1 + '\t' + str(0)  + '\t' + Gname2 + '\t' + Gstrand + '\n'


if __name__ == "__main__":
    start = time.time()
    results = []
    try:
        parser = argparse.ArgumentParser(description="""Map TSS to the nearest regions.""")
        parser.add_argument( '-g', dest='genome', default='hg38', choices=['mm9','hg19','mm10','hg38'], required=False, help='genome' )
        parser.add_argument( '-n', dest='name',   type=str, required=True, help='name to associate with input bed file' )
        parser.add_argument( '-a', dest='alpha',  type=float, default=1e4, required=False, help='effect of distance on regulatory potential. Distance at which regulatory potential is 1/2, (default=10kb)' )
        parser.add_argument( '-b', dest='bw',     type=str, required=True, help='Bigwig file name' )
        parser.add_argument( '--cs', dest='chromesize',  type=str, required=False, default="/nv/vol190/zanglab/sh8tv/Data/Genome/hg38/chromInfo_hg38.txt", help='chrom.sizes is two columns: <chromosome name> <size in bases>' )
#        parser.add_argument( '--tss', dest='refseqTSS',  type=str, required=False, default="/scratch/sh8tv/Project/regulation_network/Data/annotation/hg38_uniq_symbol_withExonArrayExp_TSS.bed", help='refseqTSS is six columns: <chromosome name> <TSS> <TSS+1> <refseq:genesymbok> <score> <strand>' )
#        parser.add_argument( '--tss', dest='refseqTSS',  type=str, required=False, default="/nv/vol190/zanglab/sh8tv/Data/refgenes/hg38/hg38_refseq_clean_tss_forRP.bed", help='refseqTSS is six columns: <chromosome name> <TSS> <TSS+1> <refseq:genesymbok> <score> <strand>' )
        parser.add_argument( '--tss', dest='refseqTSS',  type=str, required=False, default="/nv/vol190/zanglab/sh8tv/Data/refgenes/geneID_annotation/hg38_gene_annotation_geneID_LenOrder_TSS_forRP.bed", help='refseqTSS is six columns: <chromosome name> <TSS> <TSS+1> <refseq:genesymbok> <score> <strand>' )
        parser.add_argument( '--thread', dest='threads',  type=int, required=False, default=8, help='Number of threads to calcuate the Regulatory Potential, DEFAULT=8' )
        parser.add_argument( '--bwsum', dest='bwsum',  type=str, required=False, help='Path fot the bigWigSummary script', default='/nv/vol190/zanglab/sh8tv/bin/UCSCTools/bigWigSummary' )
        args = parser.parse_args()
        #chromeinfo = chromesizeinfo(args.chromesize)
        output = args.name
        outf = open(output,'w+')
        alpha = -math.log(1.0/3.0)*1e5/args.alpha
        Infos = getInfo(args.genome, args.name, alpha, args.bw, args.refseqTSS, args.chromesize, args.bwsum)
        #print Infos[0]
        p = Pool(args.threads)
        result = p.map_async(getRP, Infos, callback=results.append)
        p.close()
        p.join()
       # print results
        for line in results:
            for element in line:
                outf.write(element)

    except KeyboardInterrupt:
        sys.stderr.write("User interrunpt me! ;-) Bye!\n")
        sys.exit(0)
    outf.close()       
    end = time.time()
    total = end - start
    hour = int(total/3600)
    minite = int(total - hour*3600)/60
    second = int(total - hour*3600 - minite*60)
    print 'total time: %s:%s:%s'%(hour, minite, second)


