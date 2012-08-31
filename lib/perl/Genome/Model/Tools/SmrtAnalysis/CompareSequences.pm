package Genome::Model::Tools::SmrtAnalysis::CompareSequences;

use strict;
use warnings;

use Genome;

use File::Temp qw/tempfile/;

my $DEFAULT_LSF_RESOURCE = "-g /pacbio/smrtanalysis -M 4000000 -R 'select[type==LINUX64 && mem>=4000 && tmp>=40000] rusage[mem=4000,tmp=20000] span[hosts=1]' -n 4";

class Genome::Model::Tools::SmrtAnalysis::CompareSequences {
    is  => 'Genome::Model::Tools::SmrtAnalysis::Base',
    has_input => [
        query => {
            doc => 'A file, fofn, or directory.',
        },
        target => {
            doc => 'A file, fofn, or directory.',
        },
        algorithm => {
            is => 'Text',
            valid_values => ['exonerate','yass','blasr','lastz'],
        },
    ],
    has_optional_input => [
        base_output_directory => {
            is => 'Text',
            doc => 'A base output directory used for parallel processing',
        },
        aligner_output_file => {
            is => 'Text',
            doc => 'The STDOUT from the aligner.',
        },
        debug => {
            is => 'Boolean',
            doc => 'Turn on debugging output',
            default_value => 0,
        },
        local => {
            is => 'Boolean',
            doc => 'Use local alignment instead of global alignment',
            default_value => 0,
        },
        dirty => {
            is => 'Boolean',
            doc => 'Use settings that are good for aligning dirty reads to a reference',
            default_value => 0,
        },
        algorithm_params => {
            is => 'Text',
            doc => 'Pass through options for algorithm (e.g. exonerate)',
        },
        min_accuracy => {
            is => 'Number',
            doc => 'Min accuracy to output a hit',
        },
        min_length => {
            is => 'Number',
            doc => 'Min length to output a hit',
        },
        min_z => {
            is => 'Number',
            doc => 'Min Z score to output a hit',
        },
        unique_ids => {
            is => 'Boolean',
            doc => 'Modify query ids so that each hit has a unique id',
            default_value => 0,
        },
        nproc => {
            is => 'Number',
            doc => 'Number of processors to use for alignment computation',
        },
        delta => {
            is => 'Boolean',
            doc => 'Output hits in the delta format from MUMmer',
            default_value => 0,
        },
        show_alignment => {
            is => 'Boolean',
            doc => 'Whether or not to show the alignment in the XML output.  This option must be specified for some versions of the Z calculation.',
            default_value => 0,
        },
        noise_data => {
            is => 'Text',
            doc => 'noise triplet, .xy or compare XML file containing noise data for estimating a z-score',
        },
        trim_window => {
            is => 'Number',
            doc => 'Window size for trimming ends [not implemented]',
        },
        trim_errors => {
            is => 'Number',
            doc => 'Trim ends until there are less than trimErrors in the window [not implemented]',
        },
        multiple => {
            is => 'Text',
            doc => 'Specify a policy for how to treat multiple hits. bestscore returns the best alignment score hit.  All returns all hits, limited by the alignment routine, and deltaz limits by z value.',
            valid_values => ['random','all','deltaz','bestscore','leftpositive'],
            default_value => 'bestscore',
        },
        hdf5_file => {
            is => 'Text',
            doc => 'Write to specified file in HDF5 format, set showAlignmet to be true.',
            is_output => 1,
        },
        hdf5_mode => {
            is => 'Text',
            doc => "'w' for creating a new hdf5 file, 'a' for appending",
            valid_values => ['w','a'],
        },
        hdf5_pbi => {
            is => 'Boolean',
            doc => 'Create a PBCmpH5 file. Only works with .fofn input and the --hd5-file option.',
            default_value => 0,
        },
        ref_seq_name => {
            is => 'Text',
            doc => 'string for the name of the target sequence',
        },
        keep_temp => {
            is => 'Boolean',
            doc => 'Do not delete the temporary output file.',
            default_value => 0,
        },
        use_temp => {
            is => 'Text',
            doc => 'Specify a temporary output name, rather than using an auto-generated one in /tmp.',
        },
        pls2fasta => {
            is => 'Boolean',
            doc => 'Convert pls files into fasta before alignment.',
            default_value => 0,
        },
        split_subreads => {
            is => 'Boolean',
            doc => 'Split reads into subreads if subread regions are available.',
            default_value => 1,
        },
        advance_half => {
            is => 'Boolean',
            doc => 'Hack for speeding up blasr.',
            default_value => 0,
        },
        region_table => {
            is => 'Text',
            doc => 'Specify a regions table for filtering reads.',
            is_output => 1,
        },
        lookup_region_table => {
            is => 'Text',
            doc => 'Lookup the region_table based on the input.fofn file and the extenstion defined by this variable.  Used only for parallel processing',
        },
        ignore_quality => {
            is => 'Boolean',
            doc => 'Force blasr to ignore quality values when computing a local alignment.',
            default_value => 0,
        },
        ignore_regions => {
            is => 'Boolean',
            doc => 'Ignore all information in a regions table, even if one exists.',
            default_value => 0,
        },
        ignore_hq_regions => {
            is => 'Boolean',
            doc => 'Do not use high-quality region information, even if it exists (no masking is done).',
            default_value => 0,
        },
        xml => {
            is => 'Boolean',
            doc => 'Generate XML output.',
            default_value => 1,
        },
        single_read_group => {
            is => 'Boolean',
            doc => 'When writing cmpH5 store reads in a single default read group.',
            default_value => 0,
        },
        subsample => {
            is => 'Number',
            doc => 'Makes blasr subsample and align reads at a desired fraction of input.',
        },
        randomize_deletions => {
            is => 'Boolean',
            doc => 'Post-process alignment using a gap randomizer for deletions in homopolymers',
            default_value => 0,
        },
        filter_adapter_only => {
            is => 'Boolean',
            doc => 'If specified, do not report adapter-only hits using annotations associated with the reference entry.',
            default_value => 0,
        },
        respect_fasta_given_subread_location => {
            is => 'Boolean',
            doc => 'If specified, add the given subread start to the hit coordinates of the query_id in the cmp.h5 file giving global read location',
            default_value => 0,
        },
        use_guided_align => {
            is => 'Boolean',
            doc => 'Pass the useGuidedAlign option to BLASR',
            default_value => 0,
        },
        use_ccs => {
            is => 'Text',
            doc => 'Map the ccsSequence to the genome first, then align subreads to the interval that the CCS read mapped to.  fullpass only aligns subreads that span the length of the template.  Specifying allpass maps all subreads.',
            valid_values => ['fullpass','allpass','denovo'],
        },
        profile => {
            is => 'Boolean',
            doc => 'Use the python cProfile module to profile the running of this script.',
            default_value => 0,
        },
    ],
    has_optional_param => [
        lsf_resource => { default_value => $DEFAULT_LSF_RESOURCE },
    ],
};

sub help_brief {
    'Run Pacific Biosciences aligners.'
}


sub help_detail {
    return <<EOS 
Compares sequences to each other using an algorithm selected from a
selection of supported command-line alignment algorithms.  Supports
output in an XML format, MUMmer delta format, or cmp.h5.
Miscellaneous options for parallel processing, calculating read
significance and filtering alternative hits.
EOS
}

sub execute {
    my $self = shift;

    my @input_files;
    my @output_files;
    if (-f $self->query) {
        push @input_files, $self->query;
    }
    if (-f $self->target) {
        push @input_files, $self->target;
    }
    my $cmd = $self->analysis_bin .'/compareSequences.py --algorithm='. $self->algorithm;
    if ($self->debug) {
        $cmd .= ' --debug';
    }
    if (defined($self->algorithm_params)) {
        # TODO: The split on - is dangerous(the dash character may be in the file path, etc.) and presumes all command line args start with -(blasr)
        my @params = split('-',$self->algorithm_params);
        for my $param (@params) {
            if ($param =~ /^\s*$/) { next; }
            #This will not work for boolean or flags
            #unless ($param =~ /^\S+\s+\S+$/) {next;}
            $cmd .= ' -x -'. $param;
        }
    }
    if ($self->local) {
        $cmd .= ' --local';
    }
    if ($self->dirty) {
        $cmd .= ' --dirty';
    }
    if (defined($self->min_accuracy)) {
        $cmd .= ' --minAccuracy='. $self->min_accuracy;
    }
    if (defined($self->min_length)) {
        $cmd .= ' --minLength='. $self->min_length;
    }
    if (defined($self->min_z)) {
        $cmd .= ' --minZ='. $self->min_z;
    }
    if ($self->unique_ids) {
        $cmd .= ' --uniqueIds';
    }
    if (defined($self->nproc)) {
        $cmd .= ' --nproc='. $self->nproc;
    }
    if ($self->delta) {
        $cmd .= ' --delta';
    }
    if ($self->show_alignment) {
        $cmd .= ' --showAlignment';
    }
    if (defined($self->noise_data)) {
        $cmd .= ' --noiseData='. $self->noise_data;
    }
    if (defined($self->trim_window)) {
        $cmd .= ' --trimWindow='. $self->trim_window;
    }
    if (defined($self->trim_errors)) {
        $cmd .= ' --trimErrors='. $self->trim_errors;
    }
    if (defined($self->multiple)) {
        $cmd .= ' --multiple='. $self->multiple;
    }
    if (defined($self->lookup_region_table)) {
        if ($self->region_table) {
            die('Do not define a region_table when attemtping to lookup one!');
        }
        my $query = $self->query;
        my ($basename,$dirname,$suffix) = File::Basename::fileparse($query,qw/.fofn/);
        unless ($basename =~ /(\S+)_input/) {
            die('Failed to lookup the region table from query '. $query);
        } else {
            my $region_table = $dirname .'/'. $1 . $self->lookup_region_table;
            unless (-e $region_table) {
                die('Expected region_table '. $region_table .' not found!');
            }
            $self->region_table($region_table);
        }
    }
    if (defined($self->base_output_directory)) {
        my (undef, $hdf5_file) = tempfile('XXXXXXX', DIR => $self->base_output_directory, SUFFIX => '.cmp.h5');
        $self->hdf5_file($hdf5_file);
    }
    if (defined($self->hdf5_file)) {
        unless (defined($self->hdf5_mode)) {
            die('Failed to define HDF5 mode while defining file.');
        }
        push @output_files, $self->hdf5_file;
        $cmd .= ' --h5fn='. $self->hdf5_file;
        unless (defined($self->hdf5_mode)) {
            die('Failed to define an hdf5_mode!');
        }
        $cmd .= ' --h5mode='. $self->hdf5_mode;
        if ($self->hdf5_pbi) {
            $cmd .= ' --h5pbi';
        }
    } else {
        if ($self->hdf5_pbi || defined($self->hdf5_mode)) {
            die('Defined HDF5 params, but no HDF5 file defined.');
        }
    }
    if ($self->ref_seq_name) {
        $cmd .= ' --refSeqName='. $self->ref_seq_name;
    }
    if ($self->keep_temp) {
        $cmd .= ' --keepTemp';
    }
    if (defined($self->use_temp)) {
        $cmd .= ' --useTemp='. $self->use_temp;
    }
    if ($self->pls2fasta) {
        $cmd .= ' --pls2fasta';
    }
    unless ($self->split_subreads) {
        $cmd .= ' --noSplitSubreads';
    }
    if ($self->advance_half) {
        $cmd .= ' --advanceHalf';
    }
    if ($self->region_table) {
        push @input_files, $self->region_table;
        $cmd .= ' --regionTable='. $self->region_table;
    }
    if ($self->ignore_quality) {
        $cmd .= ' --ignoreQuality';
    }
    if ($self->ignore_regions) {
        $cmd .= ' --ignoreRegions';
    }
    if ($self->ignore_hq_regions) {
        $cmd .= ' --ignoreHQRegions';
    }
    unless ($self->xml) {
        $cmd .= ' --noXML',
    }
    if ($self->single_read_group) {
        $cmd .= ' --singleReadGroup';
    }
    if (defined($self->subsample)) {
        $cmd .= ' --subsample='. $self->subsample;
    }
    if ($self->randomize_deletions) {
        $cmd .= ' --randomizeDeletions';
    }
    if ($self->filter_adapter_only) {
        $cmd .= ' --filterAdapterOnly';
    }
    if ($self->respect_fasta_given_subread_location) {
        $cmd .= ' --respectFastaGivenSubreadLocation';
    }
    if ($self->use_guided_align) {
        $cmd .= ' --useGuidedAlign';
    }
    if (defined($self->use_ccs)) {
        $cmd .= ' --useCcs='. $self->use_ccs;
    }
    if ($self->profile) {
        $cmd .= ' --profile';
    }
    $cmd .= ' '. $self->query .' '. $self->target;
    if ($self->aligner_output_file) {
        push @output_files, $self->aligner_output_file;
        $cmd .= ' > ', $self->aligner_output_file;
    }
    $self->shellcmd(
        cmd => $cmd,
        input_files => \@input_files,
        output_files => \@output_files,
        skip_if_output_is_present => 0,
    );
    return 1;
}


1;
