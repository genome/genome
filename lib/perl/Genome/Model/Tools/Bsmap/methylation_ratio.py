'''
Created on Jan 26, 2015

@author: ygindin
'''
import optparse
import logging
import time
import tempfile
import pysam
import os
import pdb
from multiprocessing import Process, Queue

def write_data(result_queue, out_file, work_units):
    start = time.clock()
    fout = open(out_file, 'w')
    fout.write('chr\tpos\tratio\ttotal_C\tmethy_C\tunmeth_C\n')
    tab = '\t'
    new_line = '\n'
    print 'writer opens shop. waiting for ', work_units, ' work units'
    work_written = 0
    logging.basicConfig(level=logging.DEBUG, format='(%(threadName)-10s) %(message)s')

    while work_written < work_units:
        counter = result_queue.get()
        work_written += 1
        for coordinate in counter.keys():
            depth = float(counter[coordinate].getCoverage())
            meth_depth = float(counter[coordinate].getMethylationCoverage())
            un_methylated = float(counter[coordinate].get_un_methylation_coverage())
            if meth_depth + un_methylated == 0:
                # it's a mismatch, not a methylation eligible event
                continue
            ratio = meth_depth / (meth_depth + un_methylated)
            chromosome = counter[coordinate].get_chromosome()
            
            fout.write(tab.join([chromosome, str(coordinate+1), str(ratio), str(depth), str(meth_depth), 
                                 str(un_methylated) + new_line]))            
        
        logging.debug('wrote %s', ''.join(['work unit ', str(work_written), ' out of ', str(work_units)]))
    print 'writer calling quits'
    fout.close()
    end = time.clock()
    run_time = end - start
    print 'Running time: ', run_time, ' seconds'

def isOnWatsonStrand(tags):
    """ Parses the ZS field to determine if the read is on the Crick strand """
    for tag in reversed(tags):
        if 'ZS' in tag:
            if tag[1].startswith('+'):
                return True
            else:
                return False
    return None


def process_region_pileup(i, bam_file_path, fasta_path, work_queue, out_queue):
    bam_file = pysam.AlignmentFile(bam_file_path, 'rb')
    logging.basicConfig(level=logging.DEBUG,
                        format='(%(threadName)-10s) %(message)s')

    genome = pysam.FastaFile(fasta_path)
    nuc_A = 'A'
    nuc_T = 'T'
    nuc_C = 'C'
    nuc_G = 'G'
    cg = ['C', 'G']
    while True:
#     while work_queue.qsize() > 0:
        chromosome, start, end = work_queue.get()      
        counter = {}
        ref_seq = None
        # stepper=all does basic read filtering
        for pileupcolumn in bam_file.pileup(chromosome, start, end, stepper="all"):
            if ref_seq is None: 
                ref_seq = genome.fetch(chromosome, start, end+3)
            position = pileupcolumn.pos
            if position < start or position >= end or pileupcolumn.n == 0:
                continue
            ref_seq_offset = position-start
            if ref_seq[ref_seq_offset] in cg:
                counter[position] = Counter()
                counter[position].set_coverage(pileupcolumn.n)
                counter[position].set_chromosome(chromosome)
                meth_cov = 0
                un_meth_cov = 0
                for pileupread in pileupcolumn.pileups:          
                    read = pileupread.alignment
                    watson_strand = isOnWatsonStrand(read.tags)
                    if read.query_sequence[pileupread.query_position] == ref_seq[ref_seq_offset]:
                        if ref_seq[ref_seq_offset] == nuc_C and watson_strand:
                            meth_cov += 1
                        elif ref_seq[ref_seq_offset] == nuc_G and watson_strand is False:
                            meth_cov += 1
                    elif ref_seq[ref_seq_offset] == nuc_C and watson_strand and read.query_sequence[pileupread.query_position] == nuc_T:
                        un_meth_cov += 1
                    elif ref_seq[ref_seq_offset] == nuc_G and watson_strand is False and read.query_sequence[pileupread.query_position] == nuc_A:
                        un_meth_cov += 1           
                counter[position].set_methylation_coverage(meth_cov)
                counter[position].set_un_methylation_coverage(un_meth_cov)
        out_queue.put(counter)
    print 'Quitting work'

class Counter:
    def __init__(self):
        self.coverage = 0
        self.methylation_coverage = 0
        self.un_methylation_coverage = 0
    def set_chromosome(self, chr):
        self.chromosome = chr
    def set_coverage(self, cov):
        self.coverage = cov
    def set_context(self, context):
        self.context = context
    def set_methylation_coverage(self, cov):
        self.methylation_coverage = cov
    def set_un_methylation_coverage(self, cov):
        self.un_methylation_coverage = cov
    def getCoverage(self):
        return self.coverage
    def getMethylationCoverage(self):
        return self.methylation_coverage
    def get_un_methylation_coverage(self):
        return self.un_methylation_coverage
    def get_chromosome(self):
        return self.chromosome
    def get_context(self):
        return self.context

if __name__ == '__main__':
     
    do_debug = False
    usage = "usage: %prog [options] BSMAP_MAPPING_FILES (sorted BAM files)"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-o", "--out", dest="outfile", metavar="FILE", help="output file name. (required)", default="")
    parser.add_option("-d", "--ref", dest="reffile", metavar="FILE", help="reference genome fasta file. (required)", default="")
    parser.add_option("-t", "--threads", dest="threads", metavar="THREADS", help="number of threads", default="1")
    options, infile = parser.parse_args()
    
    if len(options.reffile) == 0: parser.error("Missing reference file, use -d or --ref option.")
    if len(options.outfile) == 0: parser.error("Missing output file name, use -o or --out option.")
    if len(infile) == 0: parser.error("Require at least one BSMAP_MAPPING_FILE.")
    
    
    step = 100000
    genome = pysam.FastaFile(options.reffile)
    work_queue = Queue()
    work_units = 0
    for chromosome in genome.references:
        chr_length = genome.get_reference_length(chromosome)
        start = 0
        end = min(start + step, chr_length)
        while True:
            work_queue.put([chromosome, start, end])
            work_units += 1
            if end == chr_length: break
            start += step
            end = min(start + step, chr_length)
    
    num_threads = int(options.threads)
    write_queue = Queue(100)
    for i in range(num_threads):
        worker = Process(target=process_region_pileup, args=(i, infile[0], options.reffile, work_queue, write_queue,))
        worker.daemon = True
        worker.start()
    
    raw_file = tempfile.NamedTemporaryFile()
    writer = Process(target=write_data, args=(write_queue, raw_file.name, work_units,))
    writer.start()
    writer.join()
    print "Sorting and compressing reads..."
    command = ''.join(['tail -n +2 ', raw_file.name, ' | sort -k1V -k2n | bgzip > ', options.outfile, '.bgz'])
    print(command)
    os.system(command)
    
    print "Indexing..."
    command = ''.join(['tabix -f -s 1 -b 2 -e 2 ', options.outfile, '.bgz'])
    print(command)
    os.system(command)