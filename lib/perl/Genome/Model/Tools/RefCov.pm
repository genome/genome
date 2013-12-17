package Genome::Model::Tools::RefCov;

use strict;
use warnings;

use Genome;
use Devel::Size qw/size total_size/;

our $VERSION = '0.003';

# Input Parameters
my $DEFAULT_ALIGNMENT_FILE_FORMAT = 'bam';
my $DEFAULT_ROI_FILE_FORMAT = 'bed';
my $DEFAULT_VALIDATE_CHROMOSOMES = 1;

# Output Parameters
my $DEFAULT_EVALUATE_GC_CONTENT = 0;
my $DEFAULT_GENOME_NORMALIZED_COVERAGE = 0;
my $DEFAULT_ROI_NORMALIZED_COVERAGE = 0;
my $DEFAULT_ALIGNMENT_COUNT = 0;
my $DEFAULT_MAXIMUM_ALIGNMENT_COUNT = 0;
my $DEFAULT_PRINT_MIN_MAX = 0;
my $DEFAULT_PRINT_HEADERS = 0;
my $DEFAULT_RELATIVE_COVERAGE = 0;
my $DEFAULT_EMBED_BED = 0;
my $DEFAULT_NORMALIZE_WITH_FORMULA = 0;

# Coverage Parameters
my $DEFAULT_MINIMUM_MAPPING_QUALITY = 0;
my $DEFAULT_MINIMUM_BASE_QUALITY = 0;
my $DEFAULT_MINIMUM_DEPTH = 1;
#my $DEFAULT_WINGSPAN = 0;
my $DEFAULT_RELATIVE_COVERAGE_BINS = [500,500,1000,1000,1000,1000,5000,5000,5000];

# Memory/CPU Parameters
my $DEFAULT_MAXIMUM_DEPTH = 100_000_000;
my $DEFAULT_MAXIMUM_ROI_LENGTH = 10_000_000;
my $DEFAULT_WINDOW_SIZE = 10_000_000;

my @GC_HEADERS = qw/
                       gc_reflen_bp
                       gc_reflen_percent
                       gc_covlen_bp
                       gc_covlen_percent
                       gc_uncovlen_bp
                       gc_uncovlen_percent
                   /;

#Possibly replace with subroutine/CODEREF?
# Note: ROI chrom|start|stop not supported in merging. (T. Wylie)
# Wed Mar  7 15:49:02 CST 2012
my %MERGE_STATS_OPERATION = (
    name => undef,
    percent_ref_bases_covered => undef,
    total_ref_bases => '+',
    total_covered_bases => '+',
    missing_bases => '+',
    ave_cov_depth => '* total_covered_bases',
    sdev_ave_cov_depth => 'weighted_mean',
    med_cov_depth => 'weighted_mean',
    gap_number => '+',
    ave_gap_length => '* gap_number',
    sdev_ave_gap_length => 'weighted_mean',
    med_gap_length => 'weighted_mean',
    min_depth_filter => 'min_depth_filter',
    min_depth_discarded_bases => '+',
    percent_min_depth_discarded => undef,
    gc_reflen_bp => '+',
    gc_reflen_percent => undef,
    gc_covlen_bp => '+',
    gc_covlen_percent => undef,
    gc_uncovlen_bp => '+',
    gc_uncovlen_percent => undef,
    roi_normalized_depth => 'weighted_mean',
    genome_normalized_depth => 'weighted_mean',
    intervals => 'intervals',
    alignment_count => '+',
    fwd_strand => '+',
    rev_strand => '+',
    minimum_coverage_depth => '<',
    maximum_coverage_depth => '>',
);


class Genome::Model::Tools::RefCov {
    is => ['Command'],
    has_input => [
        alignment_file_path => {
            doc => 'The path to the alignment file path.',
        },
        roi_file_path => {
            doc => 'The format of the region-of-interest file.',
        },
        alignment_file_format => {
            is => 'String',
            doc => 'The format of the alignment file.',
            default_value => $DEFAULT_ALIGNMENT_FILE_FORMAT,
            valid_values => ['bam'],
            is_optional => 1,
        },
        roi_file_format => {
            is => 'String',
            doc => 'The format of the region-of-interest file.',
            default_value => $DEFAULT_ROI_FILE_FORMAT,
            valid_values => ['bed','bam'],
            is_optional => 1,
        },
        min_depth_filter => {
            is => 'String',
            doc => 'The minimum depth at each position to consider coverage.  For more than one, supply a comma delimited list(ie. 1,5,10,15,20)',
            default_value => $DEFAULT_MINIMUM_DEPTH,
            is_optional => 1,
        },
        maximum_depth => {
            is => 'Integer',
            doc => 'The maximum depth to evaluate by samtools pileup.',
            default_value => $DEFAULT_MAXIMUM_DEPTH,
            is_optional => 1,
        },
        wingspan => {
            is => 'Integer',
            doc => 'A base pair wingspan value to add +/- of the input regions',
            is_optional => 1,
        },
        min_base_quality => {
            is => 'Integer',
            doc => 'only consider bases with a minimum phred quality',
            default_value => $DEFAULT_MINIMUM_BASE_QUALITY,
            is_optional => 1,
        },
        min_mapping_quality => {
            is => 'Integer',
            doc => 'only consider alignments with minimum mapping quality',
            default_value => $DEFAULT_MINIMUM_MAPPING_QUALITY,
            is_optional => 1,
        },
        genome_normalized_coverage => {
            is => 'Boolean',
            doc => 'Normalized coverage based on average depth across entire reference genome.',
            default_value => $DEFAULT_GENOME_NORMALIZED_COVERAGE,
            is_optional => 1,
        },
        roi_normalized_coverage => {
            is => 'Boolean',
            doc => 'Normalized coverage based on average depth across supplied regions-of-interest.',
            default_value => $DEFAULT_ROI_NORMALIZED_COVERAGE,
            is_optional => 1,
        },
        evaluate_gc_content => {
            is => 'Boolean',
            doc => 'Evaluate G+C content of the defined regions-of-interest.',
            default_value => $DEFAULT_EVALUATE_GC_CONTENT,
            is_optional => 1,
        },
        alignment_count => {
            is => 'Boolean',
            doc => 'Calculate the number of alignments that overlap region.',
            default_value => $DEFAULT_ALIGNMENT_COUNT,
            is_optional => 1,
        },
        maximum_alignment_count => {
            is => 'Integer',
            doc => 'The maximum number of alignments to count for a region.',
            default_value => $DEFAULT_MAXIMUM_ALIGNMENT_COUNT,
            is_optional => 1,
        },
        relative_coverage => {
            is => 'Boolean',
            doc => 'Determine relative coverage by position along the length of the ROI',
            default_value => $DEFAULT_RELATIVE_COVERAGE,
            is_optional => 1,
        },
        validate_chromosomes => {
            is => 'Boolean',
            doc => 'Validate the chromosome names from the ROI and that they exist in the BAM.',
            default_value => $DEFAULT_VALIDATE_CHROMOSOMES,
            is_optional => 1,
        },
        relative_coverage_bins => {
            doc => 'The size in bp of each bin increment.  The default settings result in bins <= 500, 1000, 2000, 3000, 4000, 5000, 10000, 15000, 20000,',
            default_value => $DEFAULT_RELATIVE_COVERAGE_BINS,
            is_optional => 1,
        },
        relative_coverage_file_basename => {
            is => 'Text',
            doc => 'An output file to print a summary of relative coverage by size bin',
            is_optional => 1,
        },
        print_min_max => {
            is => 'Boolean',
            doc => 'Print the minimum and maximum depth of coverage.',
            default_value => $DEFAULT_PRINT_MIN_MAX,
            is_optional => 1,
        },
        reference_fasta => {
            is => 'String',
            doc => 'The path to the reference genome fasta file',
            is_optional => 1,
        },
        output_directory => {
            doc => 'When run in parallel, this directory will contain all output and intermediate STATS files. Sub-directories will be made for wingspan and min_depth_filter params. Do not define if stats_file is defined.',
            is_optional => 1,
        },
        print_headers => {
            is => 'Boolean',
            doc => 'Print a header describing the output including column headers.',
            is_optional => 1,
            default_value => $DEFAULT_PRINT_HEADERS,
        },
        embed_bed => {
            is		  => 'Boolean',
            doc		  => 'Include associated BED information (REF:START-STOP) per target ID when reporting. Returns 0-based BED coords.',
            is_optional	  => 1,
            default_value => $DEFAULT_EMBED_BED,
        },
        normalize_with_formula => {
            is		  => 'String',
            doc		  => 'Runs normalization on AVE DEPTH and MIN & MAX DEPTHs (opt) given a Perl compatible formula, replacing $X with depth value per ROI. EXAMPLE: --normalize-with-formula=\'log( ( $X / 250_000_000 ) * 1_000_000 ) / log( 2 )\'',
            is_optional	  => 1,
            default_value => $DEFAULT_NORMALIZE_WITH_FORMULA,
        },
        merged_stats_file => {
            is => 'Text',
            doc => 'The final merged stats file only created if merge_by parameter defined',
            is_optional => 1,
        },
        merge_by => {
            is => 'Text',
            doc => 'The level of granularity to merge coverage statistics.  Requires ROI file uses interval names like $GENE:$TRANSCRIPT:$TYPE:$ORDINAL:$DIRECTION',
            is_optional => 1,
            valid_values => ['exome','gene','transcript'],
        },
        window_size => {
            is => 'Number',
            doc => 'The size of a reference window when calculating coverage of ROI greater than maximum_roi_length.',
            default_value => $DEFAULT_WINDOW_SIZE,
            is_optional => 1,
        },
        maximum_roi_length => {
            is => 'Number',
            doc => 'The maximum length of an ROI(region-of-interest) before using a window approach',
            default_value => $DEFAULT_MAXIMUM_ROI_LENGTH,
            is_optional => 1,
        },
    ],
    has_output => [
        stats_file => {
            doc => 'When run in parallel, do not define.  From the command line this file will contain the output metrics for each region.',
            is_input => 1,
            is_optional => 1,
        },
        final_directory => {
            doc => 'The directory where parallel output is written to when wingspan is defined in parallel fashion',
            is_optional => 1,
        },
    ],
    has_param => [
        lsf_queue => {
            doc => 'When run in parallel, the LSF queue to submit jobs to.',
            is_optional => 1,
            default_value => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER_ALT},
        },
        lsf_resource => {
            doc => 'When run in parallel, the resource request necessary to run jobs on LSF.',
            is_optional => 1,
            default_value => "-R 'select[type==LINUX64]'",
        },
    ],
    has_optional => [
        _alignments => {},
        _roi => {},
        _fai => {},
        _roi_stats => {},
        _genome_stats => {},
        _nuc_cov => {},
        _roi_cov => {},
        _brief_roi_cov => {},
        _relative_coverage_hash_ref => {},
        _relative_coverage_object => {},
    ],
};

sub help_brief {
    "Tools to run the RefCov tookit.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt ref-cov ...
EOS
}

sub help_detail {
'
 Output file format(STANDARD):
 [1] Region Name (column 4 of BED file)
 [2] Percent of Reference Bases Covered
 [3] Total Number of Reference Bases
 [4] Total Number of Covered Bases
 [5] Number of Missing Bases
 [6] Average Coverage Depth
 [7] Standard Deviation Average Coverage Depth
 [8] Median Coverage Depth
 [9] Number of Gaps
 [10] Average Gap Length
 [11] Standard Deviation Average Gap Length
 [12] Median Gap Length
 [13] Min. Depth Filter
 [14] Discarded Bases (Min. Depth Filter)
 [15] Percent Discarded Bases (Min. Depth Filter)

 OPTIONAL GC FIELDS:
 [1] G+C Reference Base Pair
 [2] G+C Percent of Reference
 [3] G+C Covered Base Pair
 [4] G+C Percent of Reference Covered
 [5] G+C Uncovered Base Pair
 [6] G+C Percent of Reference Uncovered

 OPTIONAL ROI NORMALIZED COVERAGE FIELD:
 [1] ROI Normalized Depth

 OPTIONAL GENOME NORMALIZED COVERAGE FIELD:
 [1] Genome Normalized Depth

 OPTIONAL ALIGNMENT COUNT FIELD:
 [1] Count of Overlapping Alignments

 OPTIONAL MIN/MAX FIELD:
 [1] Minimum Coverage Depth
 [2] Maximum Coverage Depth

 OPTIONAL BED EMBED FIELDS:
 [1] ROI reference (chromosome/scaffold)
 [2] ROI start
 [3] ROI stop

 ------------------
 Useful Definitions
 ------------------

 Target ID or Name:
 ID for the target reference space (region of interest) being
 evaluated, as determined by a corresponding BED file of target
 regions.
 Example Value:  NOTCH1.23

 Percent of Reference Bases Covered:
 Breadth-of-coverage as a percent of total target length (total covered
 bases divided by target length).
 Example Value:  78.57%

 Total Number of Reference Bases:
 Total nt bases in the target reference space.
 Example Value:  112

 Total Number of Covered Bases:
 Total covered bases in the target reference space.
 Example Value:  88

 Number of Missing Bases:
 Total uncovered bases in the target reference space (non-contiguous).
 Example Value:  24

 Average Coverage Depth:
 Mean depth-of-coverage for the target reference space; calculated
 across all positions in the target reference space whether covered or
 uncovered.
 Example Value:  2.36x

 Standard Deviation Average Coverage Depth:
 Standard deviation for mean depth-of-coverage; calculated across all
 positions in the target reference space whether covered or uncovered.
 Example Value:  1.62x

 Median Coverage Depth:
 Median depth-of-coverage for the target reference space, non-contiguous in nature.
 Example Value:  3x

 Number of Gaps:
 Number of gapped (uncovered) areas in the target reference space.
 Example Value:  1

 Average Gap Length:
 Mean gap length for areas uncovered in the target reference space.
 Example Value:  24

 Standard Deviation Average Gap Length:
 Standard deviation for mean gap length in the target reference space.
 Example Value:  0

 Median Gap Length:
 Median gap length in context of all gap lengths in the target reference space.
 Example Value:  24

 Min. Depth Filter:
 Minimum depth filter is a value passed to the coverage evaluation
 software which will limit the target reference postions contributing
 to coverage calculations--default is 0 (no filtering) or >= 1x
 requirement. A min. depth filter of 10 would require any position
 included in coverage calculations to be of depth >= 10x at all
 reference positions; any position < 10x would be considered 0x
 coverage.
 Example Value:  0

 Discarded Bases (Min. Depth Filter):
 Number of target reference space positions affected by the minimum
 depth filter--i.e., the positions set to 0x coverage in subsequent
 calculations.
 Example Value:  0

 Percent Discarded Bases (Min. Depth Filter):
 Percent of target reference space positions affected by the minimum
 depth filter.
 Example Value:  0

 GC bp (reference length):
 A hard, non-contiguous value associated with the number of G and C
 nucleotides in the target reference space.
 Example Value:  69

 GC percent (reference length):
 Percent of the target reference space comprised of G and C
 nucleotides.
 Example Value:  61.61%

 GC bp (coverage length):
 Number of covered G and C nucleotides in the target reference space.
 Example Value:  55

 GC percent (coverage length):
 Percent of all covered positions in the reference target space that
 are G and C nucleotides. Describes how much of covered sequence is G
 or C in nature.
 Example Value:  62.50%

 GC bp (uncovered length):
 Number of uncovered G and C nucleotides in the target reference space.
 Example Value:  14

 GC percent (uncovered length):
 Percent of all uncovered positions in the reference target space that
 are G and C nucleotides. Describes how much of uncovered sequence is G
 or C in nature.
 Example Value:  58.33%
';
}

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    unless ($self) { return; }

    require Bio::DB::Sam;

    if ($self->evaluate_gc_content) {
        unless ($self->reference_fasta) {
            die('In order to evaluate_gc_content a FASTA file of the reference genome must be provided');
        }
    }
    if ($self->merge_by) {
        unless ($self->merged_stats_file) {
            die('Please define a merged_stats_file in order to merge by '. $self->merge_by);
        }
    }
    # This is only necessary when run in parallel
    $self->resolve_final_directory;
    $self->resolve_stats_file;
    return $self;
}


# NOTE: I doubt we ever support anything but BAM.  If we do, some common adaptor/iterator will be necessary to do something like $alignments->next_alignment
sub _load_alignments {
    my $self = shift;
    my $alignments;
    if ($self->alignment_file_format eq 'bam') {
        $alignments  = Genome::Model::Tools::RefCov::Bam->create(bam_file => $self->alignment_file_path );
        unless ($alignments) {
            die('Failed to load alignment file '. $self->alignment_file_path);
        }
    } else {
        die('Failed to load '. $self->alignment_file_format .' file '. $self->alignment_file_path);
    }
    $self->_alignments($alignments);
    return $alignments;
}

sub alignments {
    my $self = shift;
    unless ($self->_alignments) {
        $self->_load_alignments;
    }
    return $self->_alignments;
}

sub _load_roi {
    my $self = shift;
    # TODO: Can the class Genome::Model::Tools::RefCov::ROI or a new class Genome::Model::Tools::RefCov::ROI::File resolve the appropriate adaptor based on the file type?
    my $format = $self->roi_file_format;
    my $subclass = ucfirst($format);
    my $class = 'Genome::Model::Tools::RefCov::ROI::'. $subclass;
    my %roi_params = (file => $self->roi_file_path);
    if (defined($self->wingspan)) {
        $roi_params{wingspan} = $self->wingspan;
    }
    my $regions = $class->create(%roi_params);
    unless ($regions) {
        die('Failed to load '. $self->roi_file_format .' regions-of-interest file '. $self->roi_file_path );
    }
    $self->_roi($regions);
    return $regions;
}

sub roi {
    my $self = shift;
    unless ($self->_roi) {
        $self->_load_roi;
    }
    return $self->_roi;
}

sub nucleotide_coverage {
    my $self = shift;
    unless ($self->_nuc_cov) {
        my $gc = Genome::Model::Tools::RefCov::Reference::GC->create();
        $self->_nuc_cov($gc);
    }
    return $self->_nuc_cov;
}

sub region_coverage_stat {
    my $self = shift;
    unless ($self->_roi_cov) {
        my $stat = Genome::Model::Tools::RefCov::Stats->create();
        $self->_roi_cov($stat);
    }
    return $self->_roi_cov;
}

sub brief_region_coverage_stat {
    my $self = shift;
    unless ($self->_brief_roi_cov) {
        my $stat = Genome::Model::Tools::RefCov::Stats->create( stats_mode => 'brief' );
        $self->_brief_roi_cov($stat);
    }
    return $self->_brief_roi_cov;
}

sub _load_fai {
    my $self = shift;
    my $fasta_file = $self->reference_fasta;
    unless ($fasta_file) { return; }
    unless (-f $fasta_file .'.fai') {
        # TODO: We chould try to create the fasta index
        die('Failed to find fai index for fasta file '. $fasta_file);
    }
    my $fai = Bio::DB::Sam::Fai->load($fasta_file);
    unless ($fai) {
        die('Failed to load fai index for fasta file '. $fasta_file);
    }
    $self->_fai($fai);
    return $fai;
}

sub fai {
    my $self = shift;
    unless ($self->_fai) {
        $self->_load_fai;
    }
    return $self->_fai;
}

sub _load_roi_stats {
    my $self = shift;
    my $alignments = $self->alignments;
    $self->status_message('Loading ROI Reference Stats...');
    my $roi_stats = Genome::Model::Tools::RefCov::Reference::Stats->create(
        bam => $alignments->bio_db_bam,
        bam_index => $alignments->bio_db_index,
        bed_file => $self->roi_file_path,
    );
    $self->_roi_stats($roi_stats);
    $self->status_message('Finished loading ROI Reference Stats!');
    return $roi_stats;
}

sub roi_stats {
    my $self = shift;
    unless ($self->_roi_stats) {
        $self->_load_roi_stats;
    }
    return $self->_roi_stats;
}

sub _load_genome_stats {
    my $self = shift;
    my $alignments = $self->alignments;
    $self->status_message('Loading Genome Reference Stats...');
    my $genome_stats = Genome::Model::Tools::RefCov::Reference::Stats->create(
        bam => $alignments->bio_db_bam,
        bam_index => $alignments->bio_db_index,
    );
    $self->_genome_stats($genome_stats);
    $self->status_message('Finished loading Genome Reference Stats!');
    return $genome_stats;
}

sub genome_stats {
    my $self = shift;
    unless ($self->_genome_stats) {
        $self->_load_genome_stats;
    }
    return $self->_genome_stats;
}


sub _load_relative_coverage_object {
    my $self = shift;
    my $relative_coverage = Genome::Model::Tools::RefCov::RelativeCoverage->create();
    $self->_relative_coverage_object($relative_coverage);
    return $relative_coverage;
}

sub relative_coverage_object {
    my $self = shift;
    unless ($self->_relative_coverage_object) {
        $self->_load_relative_coverage_object;
    }
    return $self->_relative_coverage_object;
}

# This is only necessary when running in parallel as a part of a workflow
# There is probably a better way of doing this
sub resolve_final_directory {
    my $self = shift;

    my $output_directory = $self->output_directory;
    if ($output_directory) {
        my $wingspan = $self->wingspan;
        if (defined($wingspan)) {
            $output_directory .= '/wingspan_'. $wingspan;
        }
        unless (-d $output_directory){
            unless (Genome::Sys->create_directory($output_directory)) {
                die('Failed to create output directory '. $output_directory);
            }
        }
        $self->final_directory($output_directory);
    }
    return 1;
}

# This is only necessary when running in parallel as a part of a workflow
# There is probably a better way of doing this
sub resolve_stats_file {
    my $self = shift;
    unless (defined($self->stats_file)) {
        unless (defined($self->output_directory)) {
            die('Failed to define output_directory or stats_file!');
        }
        my $alignment_file_expected_suffix = '.'. $self->alignment_file_format;
        my ($alignment_basename,$alignment_dirname,$alignment_suffix) = File::Basename::fileparse($self->alignment_file_path,($alignment_file_expected_suffix));
        unless ($alignment_suffix) {
            die('Failed to recognize file '. $self->alignment_file_path .' without expected suffix '. $alignment_file_expected_suffix);
        }
        my $roi_file_expected_suffix = '.'. $self->roi_file_format;
        my ($regions_basename,$roi_dirname,$roi_suffix) = File::Basename::fileparse($self->roi_file_path,($roi_file_expected_suffix));
        unless ($roi_suffix) {
            die('Failed to recognize file '. $self->roi_file_path .' without bed suffix');
        }
        $self->stats_file($self->final_directory .'/'. $alignment_basename .'_'. $regions_basename .'_STATS.tsv');
    }
    return 1;
}

sub region_sequence_array_ref {
    my $self = shift;
    my $region = shift;

    my $fai = $self->fai;
    my $id = $region->{id};
    #$self->status_message('Fetching sequence for region '. $id);
    my $dna_string = $fai->fetch($id);
    if ( !defined($dna_string) || ( length($dna_string) != $region->{length} )  ) {
        die('Failed to fetch the proper length ('. $region->{length} .') dna.  Fetched '. length($dna_string) .' bases for region: '. $region->{id});
    }
    my @dna = split('',$dna_string);
    return \@dna;
}

sub region_alignment_count {
    my $self = shift;
    my $region = shift;

    my $alignments = $self->alignments;
    my $bam = $alignments->bio_db_bam;
    my $index = $alignments->bio_db_index;
    my $tid = $alignments->tid_for_chr($region->{chrom});
    my $alignment_count = 0;
    my $forward_strand = 0;
    my $reverse_strand = 0;
    my $max_alignment_count = $self->maximum_alignment_count;
    my $alignment_count_callback = sub {
        my ($alignment,$data) = @_;
        my $flag = $alignment->flag;
        # Only count mapped reads
        unless ($flag & 4) {
            if ($flag & 16) {
                # Reverse strand
                $reverse_strand++;
            } else {
                # Positive Strand
                $forward_strand++;
            }
            $alignment_count++;
            if ($max_alignment_count && ($alignment_count >= $max_alignment_count) ) {
                die('Exceeded maximum alignment count '. $max_alignment_count .' for region '. $region->{id});
            }
        }
    };
    # The callback may die if the alignment count exceeds the maximum_alignment_count
    # If it does die, just ignore it.  die was used since I'm unsure how to return from the callback gracefully
    eval {
        $index->fetch( $bam, $tid, ($region->{start} - 1), $region->{end},$alignment_count_callback);
    };
    if ($@) {
        $self->status_message($@);
    }
    return ($alignment_count,$forward_strand,$reverse_strand);
}

sub region_coverage_array_ref {
    my $self = shift;
    my $region = shift;

    my $alignments = $self->alignments;
    my $bam = $alignments->bio_db_bam;
    my $index = $alignments->bio_db_index;
    my $tid = $alignments->tid_for_chr($region->{chrom});
    #$self->status_message('Fetching coverage for region '. $region->{id});
    my $coverage;
    if ($self->min_base_quality || $self->min_mapping_quality) {
        $coverage = $self->region_coverage_with_quality_filter($region);
    } else {
        # low-level API uses zero based coordinates
        # all regions should be zero based, but for some reason the correct length is never returned
        # the API must be expecting BED like inputs where the start is zero based and the end is 1-based
        # you can see in the docs for the low-level Bio::DB::BAM::Alignment class that start 'pos' is 0-based,but calend really returns 1-based
        # zero bins returns each position, and the max depth is required to avoid samtools 8000 read depth max
        $coverage = $index->coverage( $bam, $tid, ($region->{start} - 1), $region->{end},0,$self->maximum_depth);
    }
    unless (scalar( @{ $coverage } ) == $region->{length}) {
        die('The length of region '. $region->{name} .' '. $region->{id}
                .'('. $region->{length} .') does not match the coverage array length '. scalar( @{ $coverage }));
    }
    return $coverage;
}

sub region_coverage_with_quality_filter {
    my $self = shift;
    my $region = shift;

    my $alignments = $self->alignments;
    my $bam = $alignments->bio_db_bam;
    my $index = $alignments->bio_db_index;
    my $tid = $alignments->tid_for_chr($region->{chrom});

    my $min_mapping_quality = $self->min_mapping_quality;
    my $min_base_quality = $self->min_base_quality;
    my $quality_coverage_callback = sub {
        my ($tid,$pos,$pileups,$data) = @_;
        my ($start,$end,$coverage) = @$data;
        #Here the position $pos is always zero-based, but the end position has to be 1-based in the coverage function
        if ($pos < $start || $pos >= $end) { return; }
        my $index = $pos - $start;
        for my $pileup (@$pileups) {
            my $base_position = $pileup->qpos;
            my $alignment = $pileup->alignment;
            if ($min_mapping_quality) {
                unless ($alignment->qual >= $min_mapping_quality) {
                    next;
                }
            }
            my @base_qualities = $alignment->qscore;
            my $quality = $base_qualities[$base_position];
            if ($quality >= $min_base_quality) {
                @$coverage[$index]++;
            }
        }
    };
    #Start with an empty array of zeros
    my @coverage = map { 0 } (1 .. $region->{length});
    my $coverage = \@coverage;
    # the pileup callback will add each base gt or eq to the quality_filter to the index position in the array ref
    $index->pileup($bam,$tid,($region->{start} - 1),$region->{end},$quality_coverage_callback,[($region->{start} - 1),$region->{end},$coverage]);
    return $coverage;
}

#sub resolve_stats_class {
#    my $self = shift;
#    #TODO: this is suboptimal there is a better way to do this through delegation, inheritance, factory... something
#    if ($self->evaluate_gc_content) {
#        if ($self->roi_normalized_coverage) {
#            return Genome::Model::Tools::RefCov::ExomeCaptureStats;
#        }
#        if ($self->genome_normalized_coverage) {
#            return Genome::Model::Tools::RefCov::WholeGenomeStats;
#        }
#    } else {
#        if (!$self->roi_normalized_coverage && !$self->genome_normalized_coverage) {
#            return Genome::Model::Tools::RefCov::Stats;
#        }
#    }
#    die('No class implemented for combination of parameters!');
#}

sub merge_stats_by {
    my $self = shift;
    my $merge_by = shift;
    my $file = shift;

    unless (defined($merge_by)) {
        die('Must provide a merge_by option!');
    }

    # **NOTE**
    # Operation should be performed post print_standard_roi_coverage()
    # execution.
    my %params = (
        input   => $self->stats_file,
        separator => "\t",
    );
    my @headers;
    unless ($self->print_headers){
        @headers = $self->resolve_stats_file_headers;
        $params{headers} = \@headers;
    }
    my $stats_reader = Genome::Utility::IO::SeparatedValueReader->new(%params);
    unless (@headers) {
        @headers = @{$stats_reader->headers};
    }
    my %merge_by_stats;
    while (my $data = $stats_reader->next) {
        my $name = $data->{name};
        unless ($name) {
            die('Failed to find name for stats region: '. Data::Dumper::Dumper($data));
        }
        my ($gene,$transcript,$type,$ordinal,$strand) = split(':',$name);
        if (!defined($gene) || $gene eq '') {
            die('Failed to parse gene from ROI name:  '. $name);
        }
        unless ($merge_by eq 'gene') {
            if (!defined($transcript) || $transcript eq '') {
                die('Failed to parse transcript from ROI name:  '. $name);
            }
            unless ($merge_by eq 'transcript') {
                if (!defined($type) || $type eq '') {
                    die('Failed to parse type from ROI name:  '. $name);
                }
            }
        }
        my $merge_key = $merge_by;
        if ($merge_by eq 'gene') {
            $merge_key = $gene;
        } elsif ($merge_by eq 'transcript') {
            $merge_key = $transcript;
        } elsif ($merge_by eq 'exon') {
            die('The ROI should be at the exon level.  Why would it not?');
            $merge_key = $gene .':'. $transcript .':'. $type;
            if (defined $ordinal) {
                $merge_key .= ':' . $ordinal;
            }
        }
        $merge_by_stats{$merge_key}{intervals}++;
        for my $data_key (keys %{$data}) {
            if ($data_key eq 'name') {
                next;
            }
            my $data_value = $data->{$data_key};
            my $operation = $MERGE_STATS_OPERATION{$data_key};
            if (defined($operation)) {
                if ($operation eq '+') {
                    $merge_by_stats{$merge_key}{$data_key} += $data_value;
                } elsif ($operation =~ /^\*\s+(\S+)/) {
                    my $multiplier_key = $1;
                    my $multiplier_value = $data->{$multiplier_key};
                    $merge_by_stats{$merge_key}{$data_key} += ($data_value * $multiplier_value);
                } elsif ($operation eq 'weighted_mean') {
                    $merge_by_stats{$merge_key}{$data_key} += ($data_value * $data->{total_ref_bases});
                } elsif ($operation) {
                    $merge_by_stats{$merge_key}{$data_key} = $data_value;
                } elsif ($operation eq '<') {
                    unless (defined($merge_by_stats{$merge_key}{$data_key})) {
                        $merge_by_stats{$merge_key}{$data_key} = $data_value;
                    } elsif ($data_value < $merge_by_stats{$merge_key}{$data_key}) {
                        $merge_by_stats{$merge_key}{$data_key} = $data_value;
                    }
                } elsif ($operation eq '>') {
                    unless (defined($merge_by_stats{$merge_key}{$data_key})) {
                        $merge_by_stats{$merge_key}{$data_key} = $data_value;
                    } elsif ($data_value > $merge_by_stats{$merge_key}{$data_key}) {
                        $merge_by_stats{$merge_key}{$data_key} = $data_value;
                    }
                }
            }
        }
    }
    my $writer = Genome::Utility::IO::SeparatedValueWriter->new(
        output => $file,
        separator => "\t",
        headers => \@headers,
    );
    for my $merge_key (keys %merge_by_stats) {
        my %data;
        my $length = $merge_by_stats{$merge_key}{total_ref_bases};
        my $covered = $merge_by_stats{$merge_key}{total_covered_bases};
        my $uncovered = $merge_by_stats{$merge_key}{missing_bases};
        for my $header (@headers) {
            if (defined $merge_by_stats{$merge_key}{$header}) {
                my $data_value = $merge_by_stats{$merge_key}{$header};
                my $operation = $MERGE_STATS_OPERATION{$header};
                if ($operation =~ /^\+$/) {
                    $data{$header} = $data_value;
                } elsif ($operation =~ /^\*\s+(\S+)/) {
                    my $multiplier_key = $1;
                    my $multiplier_value = $merge_by_stats{$merge_key}{$multiplier_key};
                    if ($multiplier_value) {
                        $data{$header} = $self->_round(($data_value / $multiplier_value));
                    } elsif ($data_value) {
                        $self->error_message('For header '. $header .' found value '. $data_value .' but no denominator '. $multiplier_value);
                        die($self->error_message);
                    } else {
                        $data{$header} = 0;
                    }
                } elsif ($operation eq 'weighted_mean') {
                    $data{$header} = $self->_round(($data_value / $length));
                } elsif ($operation) {
                    $data{$header} = $data_value;
                } else {
                    die('Not sure what to do with '. $header);
                }
            } elsif ($header =~ /^gc_(\S+)_percent$/) {
                my $type = $1;
                my $method = 'gc_'. $type .'_bp';
                if ($type eq 'reflen') {
                    $data{$header} = $self->_round((($merge_by_stats{$merge_key}{$method} / $length )* 100));
                } elsif ($type eq 'covlen') {
                    if ($covered) {
                        $data{$header} = $self->_round((($merge_by_stats{$merge_key}{$method} / $covered )* 100));
                    } else {
                        $data{$header} = 0;
                    }
                } elsif ($type eq 'uncovlen') {
                    if ($uncovered) {
                        $data{$header} = $self->_round((($merge_by_stats{$merge_key}{$method} / $uncovered )* 100));
                    } else {
                        $data{$header} = 0;
                    }
                }
            } elsif ($header eq 'percent_min_depth_discarded') {
                $data{$header} = $self->_round((($merge_by_stats{$merge_key}{'min_depth_discarded_bases'} / $merge_by_stats{$merge_key}{total_ref_bases}) * 100));
            } elsif ($header eq 'percent_ref_bases_covered') {
                $data{$header} = $self->_round( ( ($covered / $length) * 100 ) );
            } elsif ($header eq 'name') {
                $data{$header} = $merge_key;
            } else {
                die('Please implement condition to deal with header: '. $header);
            }
        }
        unless ($writer->write_one(\%data)) {
            die($writer->error_message);
        }
    }
    $writer->output->close;

    return 1;
}

sub resolve_stats_file_headers {
    my $self = shift;
    my @headers = Genome::Model::Tools::RefCov::Stats->headers;
    if ($self->evaluate_gc_content) {
        push @headers, $self->gc_headers;
    }
    if ($self->roi_normalized_coverage) {
        push @headers, 'roi_normalized_depth';
    }
    if ($self->genome_normalized_coverage) {
        push @headers, 'genome_normalized_depth';
    }
    if ($self->alignment_count) {
        push @headers, 'alignment_count';
        push @headers, 'fwd_strand';
        push @headers, 'rev_strand';
    }
    if ($self->print_min_max) {
        push @headers, 'minimum_coverage_depth';
        push @headers, 'maximum_coverage_depth';
    }
    if ($self->embed_bed()) {
        push( @headers, 'ROI_ref', 'ROI_start', 'ROI_stop' );
    }
    if ($self->normalize_with_formula()) {
        push( @headers, 'norm_ave_cov_depth' );
	if ($self->normalize_with_formula() && $self->print_min_max()) {
	    push( @headers, 'norm_minimum_coverage_depth', 'norm_maximum_coverage_depth' );
	}
    }
    return @headers;
}

sub validate_chromosome_names {
    my $self = shift;
    $self->status_message('Validate chromosomes...');
    my $roi = $self->roi;
    my $refcov_bam = $self->alignments;
    for my $chr ($roi->chromosomes) {
        eval {
            my $tid = $refcov_bam->tid_for_chr($chr);
        };
        if ($@) {
            my $err = $@;
            die('Failed to validate chromosomes in ROI '. $self->roi_file_format .' file '. $self->roi_file_path .' with alignment '. $self->alignment_file_format .' file '. $self->alignment_file_path .' with error:' ."\n". $err);
        }
    }
    $self->status_message('Validate chromosomes...OK');
    return 1;
}

sub generate_coverage_stats {
    my $self = shift;
    my $region = shift;
    my $min_depth = shift;

    my $stat = $self->region_coverage_stat;
    my $length = $region->{length};
    if ($length < $self->maximum_roi_length) {
        my $coverage_array_ref = $self->region_coverage_array_ref($region);
        $stat->calculate_coverage_stats(
            coverage => $coverage_array_ref,
            min_depth => $min_depth,
            name => $region->{name},
        );
        if ($self->relative_coverage) {
            my $relative_coverage = $self->relative_coverage_object;
            $relative_coverage->min_depth($min_depth);
            $relative_coverage->coverage($coverage_array_ref);
            $relative_coverage->_revise_relative_coverage;

            my $size_bin;
            my $total_size = 0;
            for my $bin (@{$self->relative_coverage_bins}) {
                $total_size += $bin;
                if ($length <= $total_size) {
                    $size_bin = $total_size;
                    last;
                }
            }
            unless ($size_bin) {
                $size_bin = 'gt_'. $total_size;
            }

            my $relative_coverage_hash_ref = $self->_relative_coverage_hash_ref;
            unless (defined($relative_coverage_hash_ref)) {
                $relative_coverage_hash_ref = {};
            }
            my $hash_ref = $relative_coverage->relative_coverage;
            for my $relative_position (keys %{$hash_ref}) {
                $relative_coverage_hash_ref->{$size_bin}->{$relative_position} += $hash_ref->{$relative_position};
            }
            $self->_relative_coverage_hash_ref($relative_coverage_hash_ref);
        }
    } else {
        $self->status_message('Region '. $region->{name} .' for chr/reference '. $region->{chrom} .' is longer than the defined maximum_roi_length '. $self->maximum_roi_length. '!  A window approach will be used with a '. $self->window_size .' offset.');
        my @region_gaps;
        my $offset = $self->window_size;

        # 1-based coordinates
        my $start = $region->{start};
        my $end = $region->{end};

        my $seq_id = $region->{chrom};
        my $bases_covered = 0;
        my $discarded_bases = 0;
        my $coverage_stats = Statistics::Descriptive::Sparse->new();
        # This is only used when a region is larger than the maximum_roi_length
        my $brief_stat = $self->brief_region_coverage_stat;
        for (my $new_start = $start; $new_start <= $end; $new_start += $offset) {
            my $new_end = $new_start + $offset - 1;
            if ($new_end > $end) {
                $new_end = $end;
            }
            my $id = $seq_id .':'. $new_start .'-'. $new_end;
            $self->status_message('Evaluating new sub-region: '. $id);
            #$self->status_message(arena_table());
            #$self->status_message('Total Coverage Stats Size: '. total_size($coverage_stats));
            #$self->status_message('Total Region Gaps Size: '. total_size(\@region_gaps));
            my $new_region = {
                name => $id,
                start => $new_start,
                end => $new_end,
                chrom => $seq_id,
                id => $id,
                length => (($new_end - $new_start) + 1),
            };
            my $coverage_array_ref = $self->region_coverage_array_ref($new_region);
            #$self->status_message('Total Coverage Size: '. total_size($coverage_array_ref));
            $brief_stat->calculate_coverage_stats(
                min_depth => $min_depth,
                coverage => $coverage_array_ref,
                name => $new_region->{name},
            );
            #$self->status_message('Total Window Stat Size: '. total_size($brief_stat));
            $bases_covered += $brief_stat->total_covered_bases;
            $discarded_bases += $brief_stat->min_depth_discarded_bases;
            my $filtered_coverage = $brief_stat->min_depth_filtered_coverage;
            #$self->status_message('Total Window Stat with Min Depth Filter Size: '. total_size($brief_stat));
            $coverage_stats->add_data($filtered_coverage);

            #The offset should be zero-based start coordinate
            my @window_gaps = $brief_stat->generate_gaps( ($new_start - 1) );
            #$self->status_message('Total Window Stat with Gaps Size: '. total_size($brief_stat));
            if (@window_gaps) {
                if (@region_gaps) {
                    my $last_gap = $region_gaps[-1];
                    my $first_gap = $window_gaps[0];
                    if (($first_gap->[0] == $last_gap->[1])) {
                        #Blunt end clusters: merge
                        $last_gap->[1] = $first_gap->[1];
                        shift(@window_gaps);
                    }
                }
                push @region_gaps, @window_gaps;
            }
        }
        my $gap_stats = Statistics::Descriptive::Full->new();
        for my $gap (@region_gaps) {
            my $gap_length = $gap->[1] - $gap->[0];
            $gap_stats->add_data($gap_length);
        }
        #$self->status_message('Total Gap Stats Size: '. total_size($gap_stats));
        $stat->{_coverage_pdl} = undef;
        $stat->name($region->{name});
        # [0] Percent of Reference Bases Covered
        $stat->percent_ref_bases_covered( $self->_round( ( $bases_covered / $coverage_stats->count ) * 100 ) );

        # [1] Total Number of Reference Bases
        $stat->total_ref_bases( $coverage_stats->count );

        # [2] Total Number of Covered Bases
        $stat->total_covered_bases( $bases_covered ) ;

        # [3] Number of Missing Bases
        $stat->missing_bases( $stat->total_ref_bases - $stat->total_covered_bases() );

        # [4] Average Coverage Depth
        $stat->ave_cov_depth( $self->_round( $coverage_stats->mean ) );

        # [5] Standard Deviation Average Coverage Depth
        $stat->sdev_ave_cov_depth( $self->_round( $coverage_stats->standard_deviation ) );

        # [6] Median Coverage Depth
        $stat->med_cov_depth( 'n/a' );

        # [7] Number of Gaps
        if ($gap_stats->count) {
            $stat->gap_number( $gap_stats->count() );
        } else {
            $stat->gap_number( '0' );
        }

        # [8] Average Gap Length
        $stat->ave_gap_length( $self->_round( $gap_stats->mean ) );

        # [9] Standard Deviation Average Gap Length
        $stat->sdev_ave_gap_length( $self->_round( $gap_stats->standard_deviation ) );

        # [10] Median Gap Length
        $stat->med_gap_length( $self->_round( $gap_stats->median ) );

        # [11] Min. Depth Filter
        $stat->min_depth_filter( $min_depth );

        # [12] Discarded Bases (Min. Depth Filter)
        $stat->min_depth_discarded_bases( $discarded_bases );

        # [13] Percent Discarded Bases (Min. Depth Filter)
        $stat->percent_min_depth_discarded( $self->_round( ($stat->min_depth_discarded_bases() / $stat->total_ref_bases()) * 100 ) );

        # OPTIONAL COLUMNS
        $stat->minimum_coverage_depth($coverage_stats->min);
        $stat->maximum_coverage_depth($coverage_stats->max);
    }
    return $stat;
}

sub print_roi_coverage {
    my $self = shift;
    $self->status_message('Printing ROI Coverage...');

    if ($self->validate_chromosomes) {
        $self->validate_chromosome_names;
    }

    my $temp_stats_file = Genome::Sys->create_temp_file_path;
    my @headers = $self->resolve_stats_file_headers;
    my $writer = Genome::Utility::IO::SeparatedValueWriter->create(
        separator => "\t",
        headers => \@headers,
        output => $temp_stats_file,
        print_headers => $self->print_headers,
    );
    unless ($writer) {
        die 'Failed to open stats file for writing '. $temp_stats_file;
    }
    my $roi = $self->roi;
    my @min_depths = split(',',$self->min_depth_filter);
    my $region_count;
    while (my $region = $roi->next_region) {
        unless ($region_count++ % 10000) {
            $self->status_message('Processed '. $region_count .' ROI...');
        }
        my $sequence_array_ref;
        # TODO: Add a windowing method for gc content, G+C evaluation will fail at ~50Mb regions..
        if ($self->evaluate_gc_content) {
            if ($region->{length} >= $self->maximum_roi_length) {
                die('Evaluate G+C content is not supported for regions '. $region->{length} .' and a maximum_roi_length of '. $self->maximum_roi_length);
            }
            $sequence_array_ref = $self->region_sequence_array_ref($region);
        }
        for my $min_depth (@min_depths) {
            my $stat = $self->generate_coverage_stats($region,$min_depth);
            my $data = $stat->stats_hash_ref;
            if ($self->evaluate_gc_content) {
                #$self->status_message('Evaluating GC content of '. $data->{name} .' '. $region->{id});
                # TODO: getting the coverage array ref back from the stat object will not work with windows
                my $coverage_array_ref = $stat->min_depth_filtered_coverage;
                my $gc_data = $self->evaluate_region_gc_content($sequence_array_ref,$coverage_array_ref);
                for my $key (keys %{$gc_data}) {
                    $data->{$key} = $gc_data->{$key};
                }
            }
            if ($self->roi_normalized_coverage) {
                my $roi_stats = $self->roi_stats;
                my $roi_normalized_depth;
                if ($roi_stats->mean_coverage) {
                    $roi_normalized_depth = $self->_round( ($stat->ave_cov_depth / $roi_stats->mean_coverage) );
                } else {
                    # The ROI coverage was zero for all probes
                    $roi_normalized_depth = 0;
                }
                $data->{'roi_normalized_depth'} = $roi_normalized_depth;
            }
            if ($self->genome_normalized_coverage) {
                my $genome_stats = $self->genome_stats;
                my $genome_normalized_depth;
                if ($genome_stats->mean_coverage) {
                    $genome_normalized_depth = $self->_round( ($stat->ave_cov_depth / $genome_stats->mean_coverage) );
                } else {
                    # The genome coverage was zero... I hope this doesn't happen very often.
                    $genome_normalized_depth = 0;
                }
                $data->{'genome_normalized_depth'} = $genome_normalized_depth;
            }
            if ($self->alignment_count) {
                my @values = $self->region_alignment_count($region);
                 $data->{'alignment_count'} = $values[0];
                 $data->{'fwd_strand'} = $values[1];
                 $data->{'rev_strand'} = $values[2];
            }
            if ($self->print_min_max) {
                $data->{'minimum_coverage_depth'} = $stat->minimum_coverage_depth;
                $data->{'maximum_coverage_depth'} = $stat->maximum_coverage_depth;
            }
	    if ($self->embed_bed()) {
		# REGION getter always returns 1-based BED. We will
		# return 0-based BED by default.
		$data->{'ROI_ref'}   = $region->{'chrom'};
		$data->{'ROI_start'} = $region->{'start'} - 1; # TRUE BED
		$data->{'ROI_stop'}  = $region->{'end'};       # TRUE BED
	    }
	    if ($self->normalize_with_formula()) {
		# Ave Depth (normalization)
		my $ave_formula = $self->normalize_with_formula();
		my $ave_depth_x = $stat->ave_cov_depth();
		$ave_formula =~ s/\$X/$ave_depth_x/g;
		$data->{'norm_ave_cov_depth'} = eval( $ave_formula );
		if ($self->normalize_with_formula() && $self->print_min_max()) {
		    # Min Depth (normalization)
		    my $min_formula = $self->normalize_with_formula();
		    my $min_depth_x = $stat->minimum_coverage_depth();
		    $min_formula =~ s/\$X/$min_depth_x/g;
		    $data->{'norm_minimum_coverage_depth'} = eval( $min_formula );
		    # Max Depth (normalization)
		    my $max_formula = $self->normalize_with_formula();
		    my $max_depth_x = $stat->maximum_coverage_depth();
		    $max_formula =~ s/\$X/$max_depth_x/g;
		    $data->{'norm_maximum_coverage_depth'} = eval( $max_formula );
		}
	    }
            unless ($writer->write_one($data)) {
                die($writer->error_message);
            }
        }
    }
    $writer->output->close;
    $self->status_message('Copying Temporary Coverage File...');
    Genome::Sys->copy_file($temp_stats_file, $self->stats_file);

    if ($self->merge_by && $self->merged_stats_file) {
        $self->merge_stats_by($self->merge_by,$self->merged_stats_file);
    }
    if ($self->relative_coverage) {
        my $relative_coverage_hash_ref = $self->_relative_coverage_hash_ref();
        my %total;
        for my $bin (keys %{$relative_coverage_hash_ref}) {
            for my $position (sort {$a <=> $b} keys %{$relative_coverage_hash_ref->{$bin}}) {
                $total{$position} += $relative_coverage_hash_ref->{$bin}->{$position};
            }
        }
        for my $bin (keys %{$relative_coverage_hash_ref}) {
            my $relative_coverage_file = $self->relative_coverage_file_basename .'_'. $bin .'.tsv';
            my $relative_coverage_fh = Genome::Sys->open_file_for_writing($relative_coverage_file);
            unless ($relative_coverage_fh) {
                die('Failed to open output file handle for writing: '. $relative_coverage_file);
            }
            print $relative_coverage_fh "relative_position\trelative_depth\ttotal_depth\n";
            for my $position (sort {$a <=> $b} keys %{$relative_coverage_hash_ref->{$bin}}) {
                my $relative_depth;
                if ($relative_coverage_hash_ref->{$bin}->{$position}) {
                    # TODO: This may or may not be useful....
                    $relative_depth = $relative_coverage_hash_ref->{$bin}->{$position} / $total{$position};
                }
                print $relative_coverage_fh $position ."\t". $relative_depth ."\t". $relative_coverage_hash_ref->{$bin}->{$position} ."\n";
            }
            $relative_coverage_fh->close;
        }
        # TODO: print a total depth per relative_position file
    }
    $self->status_message('Print ROI Coverage...OK');
    return 1;
}


sub evaluate_region_gc_content {
    my $self = shift;
    my $sequence = shift;
    my $coverage = shift;

    #$self->status_message('Loading GC Reference Coverage...');
    my $nucleotide_coverage = $self->nucleotide_coverage;
    $nucleotide_coverage->calculate_nucleotide_coverage(
        sequence => $sequence,
        coverage => $coverage,
    );
    unless ($nucleotide_coverage) {
        die('Failed to create GC coverage!');
    }
    #$self->status_message('Finished loading GC Reference Coverage...');
    my $gc_hash_ref = $nucleotide_coverage->gc_hash_ref;
    return $gc_hash_ref;
}

sub gc_headers {
    return @GC_HEADERS;
}

sub _round {
    my $self = shift;
    my $value = shift;
    return sprintf( "%.2f", $value );
}

1;
