package Genome::Model::Tools::CopyNumber::Jistic;

##############################################################################
#
#
#	AUTHOR:		Chris Miller (cmiller@genome.wustl.edu)
#
#	CREATED:	07/25/11
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


class Genome::Model::Tools::CopyNumber::Jistic {
    is => 'Command',
    has => [
	histogram_bin_size => {
	    is => 'Float',
	    is_optional => 1,
	    doc => 'Defines the bin size used to quantize the copy number histogram',
            default => 0.001,
	},

	qval_threshold => {
	    is => 'Float',
	    is_optional => 1,
	    doc => 'Regions with FDR q-value bellow the threshold are considered significant. Default:
0.25',
            default => 0.25,
	},

        gene_distance_max => {
            is => 'Integer',
	    is_optional => 1,
            doc => 'Each aberrant region is associated with all genes overlapping it, and also with the closest gene in each direction, as long as the distance is less than $distance base pairs. This feature tries to model cases in which the aberration is affecting the gene‘s control regions (promoter and regulatory sequences).',
            default => 10000,
        },
        
        by_chromosome => {
            is => 'Boolean',
	    is_optional => 1,
            doc => 'When this flag is set, each chromosome will be processed separately after calculatingthe null distribution.',
            default => 1,
        },
        
        fix_for_single_aberrantSample => {
            is => 'Boolean',
            is_optional => 1,
            doc => 'Used for data that has peaks within a single aberrant sample. Normally, if the data presents this type of behavior, it is probably unsuitable for the GISTIC algorithm either because there are too few samples or because the samples have too few aberrations. If the user still wants to run the algorithm, this parameter will allow cause these peaks to be ignored (they will simply be peeled-off and will not be presented as significant peaks).',
            default => 0,
        },

        skip_all_normal => {
            is => 'Boolean',
            is_optional => 1,
            doc => 'Normally, all genes affected by the aberration (as defined by the genedistancemax) will be outputted. This flag is used when the user does not want to output genes for which the markers in their actual location do not pass the aberration significant threshold for any sample.',
            default => 0,
        },

        mixed_aberrations => {
            is => 'Boolean',
            is_optional => 1,
            doc => 'Flag used when the user wants to consider aberrations in absolute value instead of considering amplifications and deletions separately (this will rarely happen)',
            default => 0,
        },

        debug => {
            is => 'Boolean',
            is_optional => 1,
            doc => 'Use this flag to output debugging information to the standard error.',
            default => 0,
        },

        limted_peeloff => {
            is => 'Boolean',
            is_optional => 1,
            doc => 'Selects the limited peel-off variant of JISTIC (include reference when available). If a distribution file is not specified, the qvalue distribution created by this run is used to determine the peel-off limits.',
            default => 0,
        },

        aberration_gain_threshold => {
            is => 'Float',
            is_optional => 1,
            doc => 'Threshold copy number values to be considered aberrant for gains',
            default => 0.32,
        },

        aberration_loss_threshold => {
            is => 'Float',
            is_optional => 1,
            doc => 'Threshold copy number values to be considered aberrant for loss',
            default => -0.41,
        },

        focal => {
            is => 'Boolean',
            is_optional => 1,
            doc => 'use the focal variant of JISTIC',
            default => 0,
        },

        arm_peeloff => {
            is => 'Boolean',
            is_optional => 1,
            doc => 'Selects the arm peel-off variant of GISTIC. Note that arm peel-off obtains better results when using focal GISTIC',
            default => 0,
        },
        
        high_level_loss_threshold => {
            is => 'Float',
            is_optional => 1,
            doc => 'Threshold for copy number deletions',
            default => -2.0,
        },

        high_level_gain_threshold => {
            is => 'Float',
            is_optional => 1,
            doc => 'Threshold for copy number amplifications',
            default => 0.807,
        },

        cap_gain => {
            is => 'Float',
            is_optional => 1,
            doc => 'Values used to cap copy number values before applying GISTIC. Copy numbers above cap-gain will be treated as cap-gain',
            default => 2.32,
        },

        cap_loss => {
            is => 'Float',
            is_optional => 1,
            doc => 'Values used to cap copy number values before applying GISTIC. Copy numbers above cap-loss will be treated as cap-loss',
            default => -2.0,
        },

        copy_number_matrix => {
            is => 'String',
            is_optional => 1,
            doc => 'path to copy-number matrix file - format: marker, chromosome, position, Sample1, Sample2, Sample3, etc. Expects log2 tumor/normal ratios as input',
        }, 

        segment_file => {
            is => 'String',
            is_optional => 1,
            doc => 'path to segment file. Format: Sample, chr, st, sp, nmarkers, log2ratio (native CBS output)',
        },

        gene_annotations_file => {
            is => 'String',
            is_optional => 0,
            doc => 'path to gene-location file - format: uniqueId, geneName, Alias|Alias, Chr, strand(W/C), st, sp, cytoband, description',
        },

        mirna_annotations_file => {
            is => 'String',
            is_optional => 0,
            doc => 'path to mirna-location file (same format as gene-annotation file)',
        },

        cytoband_file => {
            is => 'String',
            is_optional => 1,
            doc => 'path to chromosomal bands file. Required if using focal or arm-peeloff methods. Format: chr, start, stop, name, label',
        },

        output_dir => {
            is => 'String',
            is_optional => 0,
            doc => 'output directory',
        },

        probe_sites => {
            is => 'String',
            is_optional => 0,
            doc => '3-column file containing location of probes (microarray), targeted regions (exome), or tiled windows (WGS). For WGS hg18 female, you can use /gscmnt/sata921/info/medseq/cmiller/annotations/jistic/probes.10k.hg18.female.dat. File format: ProbeName\tChr\tPosition (1-based)',
        },

        filter_dgv_sites => {
            is => 'Boolean',
            is_optional => 1,
            doc => 'Filter input to remove sites that intersect with a high-resolution set of CNVs from the Database of Genomic Variants',
            default => 0,
        },

        numchrom => {
            is => 'Integer',
            is_optional => 1,
            doc => 'The number of chromosomes in your file. For female-only sample sets, this will be 23, for sets including males, the default of 24 will work',
            default => 24,
        }
    ]
};






sub help_brief {
    "get recurrent peaks of amplification and deletion from copy number data using the JISTIC algorithm"
}

sub help_detail {
    "note that a few little-used params aren't properly implemented."
}


#########################################################################


sub execute {
    my $self = shift;
    my $histogram_bin_size = $self->histogram_bin_size;
    my $qval_threshold = $self->qval_threshold;
    my $gene_distance_max = $self->gene_distance_max;
    my $by_chromosome = $self->by_chromosome;
    my $fix_for_single_aberrantSample = $self->fix_for_single_aberrantSample;
    my $skip_all_normal = $self->skip_all_normal;
    my $mixed_aberrations = $self->mixed_aberrations;
    my $debug = $self->debug;
    my $limted_peeloff = $self->limted_peeloff;
    my $aberration_gain_threshold = $self->aberration_gain_threshold;
    my $aberration_loss_threshold = $self->aberration_loss_threshold;
    my $focal = $self->focal;
    my $arm_peeloff = $self->arm_peeloff;
    my $high_level_loss_threshold = $self->high_level_loss_threshold;
    my $high_level_gain_threshold = $self->high_level_gain_threshold;
    my $cap_gain = $self->cap_gain;
    my $cap_loss = $self->cap_loss;
    my $copy_number_matrix = $self->copy_number_matrix;
    my $segment_file = $self->segment_file;
    my $gene_annotations_file = $self->gene_annotations_file;
    my $mirna_annotations_file = $self->mirna_annotations_file;
    my $cytoband_file = $self->cytoband_file;
    my $output_dir = $self->output_dir;
    my $filter_dgv_sites  = $self->filter_dgv_sites;
    my $probe_sites = $self->probe_sites;
    my $numchrom = $self->numchrom;


    #sanity checks
    if($focal || $arm_peeloff){
        unless(defined($cytoband_file)){
            die("Need a cytoband-file")
        }
    }
    
    unless(defined($copy_number_matrix) xor defined($segment_file)){
        die("Either --segment-file or --copy-number-matrix must be defined, but not both");
    }


    #get to the output directory
    my $dir = cwd();
    print STDERR "current dir: " . $dir . "\n";
    chdir($output_dir);

    print STDERR "now in " . cwd() . "\n";

    if(defined($segment_file)){
        my $cmd = "java -Xmx1500m -classpath /gscuser/cmiller/usr/src/JISTIC/JISTIC.jar JISTIC.convertSEG $segment_file $probe_sites";
        if ($filter_dgv_sites){
            $cmd = $cmd . " excludedregions=/gscmnt/sata809/info/medseq/dkoboldt/SNPseek2/ucsc/hg19/dgv.nochr.pop-cnv-indel.excludedregions"; #/gscuser/cmiller/annotations/jistic/dgvExclude.reg
        }
        $cmd = $cmd . " numchrom=$numchrom";
        $cmd = $cmd . " >cnaMatrix.dat";

        my $return = Genome::Sys->shellcmd(
            cmd => "$cmd",
            );
        unless($return) {
            $self->error_message("Failed to execute: Returned $return");
            die $self->error_message;
        }
        $copy_number_matrix = "cnaMatrix.dat";
    }


    #print out the specs file
    open(OUTFILE,">specs") || die "can't open specs file for writing\n";
    print OUTFILE "histogramBinsize	$histogram_bin_size\n";
    print OUTFILE "qvalthreshold	$qval_threshold\n";
    if($focal){
        print OUTFILE "aberrationThresholds	focal $aberration_gain_threshold $aberration_loss_threshold\n";
    } else {
        print OUTFILE "aberrationThresholds	$aberration_gain_threshold $aberration_loss_threshold\n";
    }
    print OUTFILE "highLevelThresholds	 $high_level_gain_threshold $high_level_loss_threshold \n";

    print OUTFILE "caps $cap_gain $cap_loss\n";
    print OUTFILE "genedistancemax	$gene_distance_max\n";

    if($by_chromosome){
        print OUTFILE "ByChromosome\n";
    }
    if($arm_peeloff){
        print OUTFILE "ArmPeeloff\n"
    }
    if($skip_all_normal){
        print OUTFILE "SkipAllNormal\n";
    }
    if($mixed_aberrations){
        print OUTFILE "MixedAberrations\n";
    }
    if($debug){
        print OUTFILE "debug\n";
    }
    if($limted_peeloff){
        print OUTFILE "LimitedPeeloff\n";
    }
    if ($fix_for_single_aberrantSample){
        print OUTFILE "FixForSingleAberrantSample\n";
    }

    close(OUTFILE);
    
    
    #construct the command
    my $cmd = "java -Xmx5000m -Xms5000m -classpath /gscuser/cmiller/usr/src/JISTIC/JISTIC.jar JISTIC.Distribution spec=specs copynumber=" . $copy_number_matrix;
    if(defined($gene_annotations_file)){
        $cmd = $cmd . " locations=" . $gene_annotations_file;
    }
    if(defined($mirna_annotations_file)){
        $cmd = $cmd . " mirna=" . $mirna_annotations_file;
    }
    if(defined($cytoband_file)){
        $cmd = $cmd . " bands=" . $cytoband_file;
    }

    print STDERR $cmd . "\n";

    my $return = Genome::Sys->shellcmd(
	cmd => "$cmd",
        );
    unless($return) {
	$self->error_message("Failed to execute: Returned $return");
	die $self->error_message;
    }

    `mkdir matrices`;
    `mv AMP*.matrix matrices`;
    `mv DEL*.matrix matrices`;
    `mv miRNA*.matrix matrices`;
    `mv gene*.matrix matrices`;
    `mkdir logs`;
    `mv *.log logs`;


    # #do some conversion of the output    
    # my $inFh = IO::File->new( "peaks.txt" ) || die "can't open peaks file\n";
    # my $outFh = open (OUTFILE, ">peaks.bed") || die "Can't open output file.\n";

    # while( my $line = $inFh->getline )
    # {
    #     chomp($line);
    #     my @F = split("\t",$line);

    #     #only keep focal peaks
    #     if($F[3] eq "Y"){
    #         if ($line =~ /PEAK:\w+:chr(\d+):(\d+)-(\d+)\(/){
    #             print OUTFILE join("\t",($1,$2,$3));
    #         }
    #     }

    # }
    # close($inFh);
    # close($outFh);

    # #quick intersection
    # `intersectBed -a peaks.bed -b ~/annotations/genes.36.bed -wao >zz`;
    # $inFh = IO::File->new( "peaks.bed" ) || die "can't open peaks file\n";
    # while( my $line = $inFh->getline )
    # {
    #     chomp($line);
    #     my @F = split("\t",$line);
    #     next if $F[4] eq "-1";
    #     @a=split(/:/,$F[6]);
    #     print $a[0] . "\n"
    # }

    
    #`bash /gscuser/cmiller/oneoffs/getPeakGenes.sh $copy_number_matrix`;


    chdir($dir);
    

    return $return;
}

1;
