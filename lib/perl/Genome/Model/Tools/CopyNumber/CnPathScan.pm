package Genome::Model::Tools::CopyNumber::CnPathScan;

##############################################################################
#
#
#       AUTHOR:         Chris Miller (cmiller@genome.wustl.edu)
#
#       CREATED:        07/25/11
#
#
##############################################################################

use strict;
use Genome;
use IO::File;
use File::Basename;
use warnings;
require Genome::Sys;
use FileHandle;
use Cwd;


class Genome::Model::Tools::CopyNumber::CnPathScan {
    is => 'Command',
    has => [
        cbs_files => {
            is => 'String',
            is_optional => 1,
            doc => 'comma separated list of files containing CBS output (6 col - sample, chr, st, sp, nbins, log2). Used if you do not already have a copy number matrix',
        },

        zero_based_cbs => {
            is => 'Boolean',
            is_optional => 1,
            doc => 'Specifies whether the input cbs_files are zero-based (alternative is one-based)',
            default => 0,
        },

        qval_threshold => {
            is => 'Float',
            is_optional => 1,
            doc => 'CN Regions with FDR q-value bellow the threshold are considered significant by Jistic',
            default => 0.20,
        },

        output_directory => {
            is => 'String',
            is_optional => 0,
            doc => 'output directory that will contain results',
        },
        peak_gene_matrix => {
            is => 'String',
            is_optional => 1,
            doc => 'path to a matrix of genes in peaks and which samples they are altered in. (causes script to skip Jistic steps)',
        },

        ]
};






sub help_brief {
    "create a copy-number matrix (if necessary), apply Jistic to find recurrent peaks, then determine the significance of those peaks in the given pathway",
}

sub help_detail {
    "create a copy-number matrix (if necessary), apply Jistic to find recurrent peaks, then determine the significance of those peaks in the given pathway",
}


#########################################################################


sub execute {
    my $self = shift;
    my $qval_threshold = $self->qval_threshold;
    my $output_directory = $self->output_directory;
    my $peak_gene_matrix = $self->peak_gene_matrix;
    my $cbs_files = $self->cbs_files;
    my $zero_based_cbs = $self->zero_based_cbs;

    #sanity checks
    unless(defined($peak_gene_matrix) xor defined($cbs_files)){
        die("need either a peak-gene-matrix or a list of cbs-files");
    }

    #get to the output directory
    my $dir = cwd();
    chdir($output_directory);

    my %samples;
    my $ret;

    if(defined($cbs_files)){
        print STDERR "concatenating segments \n";

        my @files = split(/,/,$cbs_files);

        my $outFh = open (OUTFILE, ">all.cbs.segs") || die "Can't open output file.\n";
        my $ampbed = open (OUTFILE2, ">all.cbs.amp.bed") || die "Can't open output file.\n";
        my $delbed = open (OUTFILE3, ">all.cbs.del.bed") || die "Can't open output file.\n";
        foreach my $file (@files){

            next if ($file eq "");

            my $inFh = IO::File->new( $file ) || die "can't open cbs file: $file\n";

            while( my $line = $inFh->getline )
            {
                next if $line =~ /^#/;

                chomp($line);
                my @F = split("\t",$line);

                #do the output for the matrix
                if($zero_based_cbs){
                    $F[2] = $F[2] + 1;
                    $F[3] = $F[3] + 1;
                    print OUTFILE join("\t",@F) . "\n";
                } else {
                    print OUTFILE $line . "\n";
                }

                #do conversion to amp/del bed files
                #use 1.5 and 2.5 as thresholds for loss/gain
                if($F[5] >= 0.3219281){
                    print OUTFILE2 join("\t",(@F[1..3],$F[0],$F[5])) . "\n";
                } elsif ($F[5] <= -0.4150375){
                    print OUTFILE3 join("\t",(@F[1..3],$F[0],$F[5])) . "\n";
                }

                #store sample names for use later
                $samples{$F[0]} = 1;
            }
            close($inFh);
        }
        close(OUTFILE);
        close(OUTFILE2);
        close(OUTFILE3);

        print STDERR "running JISTIC\n";
        print STDERR "you can typically ignore these duplicate gene and cytoband errors\n";
        
        #now run Jistic on the matrix to get the recurrent peaks
        my $cmd = "perl -I ~/gscCode/genome/lib/perl `which gmt` copy-number jistic --segment-file all.cbs.segs --gene-annotations-file /gscuser/cmiller/sata921/annotations/jistic/hg18_Gene_Info.txt --mirna-annotations-file /gscuser/cmiller/sata921/annotations/jistic/hg18_miRNA_Info.txt --histogram-bin-size 0.01 --by-chromosome --qval-threshold $qval_threshold --output-dir $output_directory --cytoband-file /gscuser/cmiller/sata921/annotations/jistic/Human_cytoBand.txt";

        $ret = Genome::Sys->shellcmd(
            cmd => "$cmd",
            );
        unless($ret) {
            $self->error_message("Failed to execute: Returned $ret");
            die $self->error_message;
        }



        # } else { #already have the input matrix

        #     #run Jistic on the matrix to get the recurrent peaks
        #     my $cmd = "perl -I ~/gscCode/genome/lib/perl `which gmt` copy-number jistic --copy-number-matrix $copy_number_matrix --gene-annotations-file /gscuser/cmiller/sata921/annotations/jistic/hg18_Gene_Info.txt --mirna-annotations-file /gscuser/cmiller/sata921/annotations/jistic/hg18_miRNA_Info.txt --histogram-bin-size 0.01 --by-chromosome --qval-threshold $qval_threshold --output-dir $output_directory --cytoband-file /gscuser/cmiller/sata921/annotations/jistic/Human_cytoBand.txt";
        #     my $ret = Genome::Sys->shellcmd(
        #         cmd => "$cmd",
        #         );
        #     unless($ret) {
        #         $self->error_message("Failed to execute: Returned $ret");
        #         die $self->error_message;
        #     }
        # }


        #temporary solution for munging some files - these will be cleaned up eventually
        `bash /gscuser/cmiller/oneoffs/getPeakGenes.sh`;


        
        #now construct a matrix
        
        my %genes;
        #store all the peak genes and samples
        my %sampHash;
        my $inFh = IO::File->new( "peak.amp.int" ) || die "can't open peak amp file\n";
        while( my $line = $inFh->getline )
        {
            chomp($line);
            my @F = split("\t",$line);
            $sampHash{$F[3]}{$F[8]} = "amp";
            $genes{$F[8]} = 1;
            
        }
        close($inFh);
        
        $inFh = IO::File->new( "peak.del.int" ) || die "can't open peak del file\n";
        while( my $line = $inFh->getline )
        {
            chomp($line);
            my @F = split("\t",$line);
            $sampHash{$F[3]}{$F[8]} = "del";
            $genes{$F[8]} = 1;
        }
        close($inFh);
        
        $peak_gene_matrix = "peak.matrix.dat";
        
        $outFh = open (OUTFILE, ">$peak_gene_matrix") || die "Can't open output file.\n";
        #print header
        foreach my $sample (sort(keys(%samples))){
            print OUTFILE "\t" . $sample;
        }
        print OUTFILE "\n";


        foreach my $gene (sort(keys(%genes))){
            print OUTFILE $gene;
            foreach my $sample (sort(keys(%samples))){
                print OUTFILE "\t";
                if(!(exists($sampHash{$sample}{$gene}))){
                    print OUTFILE "0";
                } else {
                    if(($sampHash{$sample}{$gene} eq "amp") && ($sampHash{$sample}{$gene} eq "del")){
                        print OUTFILE "B";
                    } elsif ($sampHash{$sample}{$gene} eq "amp"){
                        print OUTFILE "A";
                    } elsif ($sampHash{$sample}{$gene} eq "del"){
                        print OUTFILE "D";
                    } else {
                        print STDERR "unknown value: $sampHash{$sample}{$gene}\n";
                    }
                }
            }
            print OUTFILE "\n";
        }

        #we also need to add in all of the genes that aren't in a peak
        $inFh = IO::File->new( "/gscuser/cmiller/sata921/annotations/genes.36.names" ) || die "can't open peak del file\n";
        while( my $line = $inFh->getline )
        {
            chomp($line);
            unless(exists($genes{$line})){
                print OUTFILE $line;

                foreach my $sample (sort(keys(%samples))){
                    print OUTFILE "\t0";
                }
                print OUTFILE "\n"
            }            
        }
        close($inFh);

        close(OUTFILE);
    }

    #create an R file
    my $outFh = open (OUTFILE, ">pathTest.R") || die "Can't open output file.\n";

    print OUTFILE 'source("/gscuser/qzhang/gstat/stat.lib")' . "\n";
    print OUTFILE 'cn.status.file="' . $peak_gene_matrix . '"' . "\n";
    print OUTFILE 'pathway.file="~/brc/pathtest/keggpath.byid"' . "\n";
    print OUTFILE 'gene.col="GENE"   # gene column in pathway.file' . "\n";
    print OUTFILE 'path.col="ID" # pathway column in pathway file' . "\n";
    print OUTFILE 'out.file="pathwayTestOutput.dat"' . "\n";
    print OUTFILE '' . "\n";
    print OUTFILE 'pa <- read.table(pathway.file,sep="\t",header=T)' . "\n";
    print OUTFILE 'pa=unique(pa[,c(gene.col,path.col)])' . "\n";
    print OUTFILE 'cn <- read.table(cn.status.file,stringsAsFactors=F)' . "\n";
    print OUTFILE 'cn=cn[rownames(cn) %in% pa[,gene.col],]' . "\n";
    print OUTFILE 'cn[cn!="0"]="1"' . "\n";
    print OUTFILE 'cn <- as.numeric.frame(cn)' . "\n";
    print OUTFILE '' . "\n";
    print OUTFILE 'cn.pathway.test(cn,pa,out.file)' . "\n";
    close(OUTFILE);

    `Rscript pathTest.R`;


    #return to the directory that launched the script
    chdir($dir);

}

1;
