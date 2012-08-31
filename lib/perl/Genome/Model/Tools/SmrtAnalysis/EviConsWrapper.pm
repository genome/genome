package Genome::Model::Tools::SmrtAnalysis::EviConsWrapper;

use strict;
use warnings;

use Genome;

use File::Temp qw/tempfile/;

my $DEFAULT_LSF_QUEUE = 'pacbio';
my $DEFAULT_LSF_RESOURCE = "-g /pacbio/smrtanalysis -M 16000000 -R 'select[type==LINUX64 && mem>=16000 && tmp>=40000] rusage[mem=16000,tmp=20000]'";

class Genome::Model::Tools::SmrtAnalysis::EviConsWrapper {
    is => ['Genome::Model::Tools::SmrtAnalysis::Base'],
    has_input => [
        cmp_hdf5_file => {
            is => 'Text',
            doc => 'The aligned reads cmp.h5 file to call consensus with',
        },
    ],
    has_optional_input => [
        base_output_directory => {
            doc => ''
         },
        base_temp_output_directory => {
            is_output => 1,
        },
        evi_cons_params => {},
        consensus_file => {
            is_output => 1,
        },
        # TODO: Add ability to run from command line with each option below or with above "params" string
        result_directory => {
            is => 'Text',
            doc => 'Absolute path of the specified directory under which to write consensus results.',
            is_output => 1,
        },
        # threshold => {
        #     doc => 'Threshold for defining a post',
        # },
        # debug => {
        #     is => 'Boolean',
        #     doc => 'Enable debugging output',
        #     default_value => 0,
        # },
        # keep_annotation => {
        #     is => 'Boolean',
        #     doc => 'Whether or not to keep the HMM (lower-case) basecalls',
        #     default_value => 0,
        # },
        # hmm_param_file => {
        #     is => 'Text',
        #     doc => 'Pre-trained HMM parameters file for initializing HMM model',
        # },
        # base_map => {
        #     is => 'Text',
        #     doc => 'Comma-separted basemap string: dye numbers for (A,C,G,T)',
        #     default_value => '3,4,2,1',
        # },
        # no_gap_posts => {
        #     is => 'Boolean',
        #     doc => 'Whether to turn off calling gap columns as potential posts',
        #     default_value => 0,
        # },
        # insertion_prior_probability => {
        #     is => 'Number',
        #     doc => 'Prior probability for insertion',
        # },
        # deletion_prior_probability => {
        #    is => 'Number',
        #    doc => 'Prior probability for deletion',
        # },
        # call_prior_probability => {
        #    is => 'Number',
        #    doc => 'Prior probability for aligning a base to reference',
        # },
        # gapped_output => {
        #     is => 'Boolean',
        #     doc => 'Whether or not to retain gaps ('-') in the consensus output',
        #     default_value => 0,
        # },
        # nproc => {
        #     is => 'Number',
        #     doc => 'Number of processors to use for parallel execution',
        #     default_value => 1,
        # },
        # run_decode => {
        #     is => 'Boolean',
        #     doc => 'Whether or not to use Decode for post identification',
        #     default_value => 1,
        # },
        # decode_file => {
        #     is => 'Text',
        #     doc => 'File containing Decode matrix',
        #     default_value => Genome::Model::Tools::SmrtAnalysis::Base->seymour_home .'/analysis/etc/defaultDecode.params',
        # },
        # fast_mode => {
        #     is => 'Boolean',
        #     doc => 'Whether or not to run in Fast mode',
        #     default_value => 1,
        # },
        variants_file => {
            is => 'Text',
            doc => 'Specify name of variants GFF file outputted after consensus calling.',
            is_output => 1,
        },
        # preserveAllGaps => {
        #     is => 'Boolean',
        #     doc => 'Whether or not to keep all gaps present in aligned sequences',
        #     default_value => 0,
        # },
        # subAlignment => {
        #     is => 'Boolean',
        #     doc => 'Whether to compute consensus for partial alignment',
        #     default_value => 1,
        # },
        # ref_start => {
        #     is => 'Number',
        #     doc => 'Start position in reference (unaligned coordinates, 0-indexing) if doing partial alignment consensus',
        # },
        # ref_end => {
        #     is => 'Number',
        #     doc => 'Inclusive end position in reference (unaligned coordinates, 0-indexing) if doing partial alignment consensus',
        # },
        # refine_alignment => {
        #     is => 'Boolean',
        #     doc => 'Whether to perform MSA refinement prior to consensus calling',
        # },
        confidence_file => {
            is => 'Text',
            doc => 'Specify name of consensus confidence file.',
            is_output => 1,
        },
        # hdf5_reference => {
        #     doc => 'Which reference group to call consensus on, if using cmpH5 as input.',
        # },
        hdf5_output => {
            is => 'Text',
            doc => 'Actively specify a HDF5 file path for writing consensus results to, as opposed to appending to the CMP H5 input.',
            is_output => 1,
        },
        # min_confidence => {
        #     doc => 'Only output variants above this confidence threshold.',
        # },
    ],
    has_optional_param => [
        lsf_queue => { default_value => $DEFAULT_LSF_QUEUE },
        lsf_resource => { default_value => $DEFAULT_LSF_RESOURCE },
    ],
};


sub execute {
    my $self = shift;
    if ($self->base_output_directory) {
        my $tmp_dir = File::Temp::tempdir('EviConsWrapper-XXXXXX',DIR => $self->base_output_directory);
        $self->base_temp_output_directory($tmp_dir);
        unless ($self->result_directory) {
            my $results_dir = File::Temp::tempdir('Results-XXXXXX',DIR => $tmp_dir);
            $self->result_directory($results_dir);
        }
        unless ($self->consensus_file) {
            my (undef, $consensus_file) = tempfile('XXXXXXX', DIR => $tmp_dir, SUFFIX => '_consensus.txt');
            $self->consensus_file($consensus_file);
        }
        unless ($self->hdf5_output) {
            my (undef, $hdf5_output) = tempfile('XXXXXXX', DIR => $tmp_dir, SUFFIX => '_cons.h5');
            $self->hdf5_output($hdf5_output);
        }
        unless ($self->confidence_file) {
            my (undef, $confidence_file) = tempfile('XXXXXXX', DIR => $tmp_dir, SUFFIX => '_consensus.conf');
            $self->confidence_file($confidence_file);
        }
        unless ($self->variants_file) {
            my (undef, $variants_file) = tempfile('XXXXXXX', DIR => $tmp_dir, SUFFIX => '_variants.gff');
            $self->variants_file($variants_file);
        }
    }
    my $cmd = $self->analysis_bin .'/eviConsWrapper.py '. $self->evi_cons_params;
    if ($self->result_directory) {   
        $cmd .= ' --resultDir='. $self->result_directory;
    }
    my @output_files;
    if ($self->variants_file) {
        $cmd .= ' --variantsFile='. $self->variants_file;
        push @output_files, $self->variants_file;
    }
    if ($self->confidence_file) { 
        $cmd .= ' --confFile='. $self->confidence_file;
        push @output_files, $self->confidence_file;
    }
    if ($self->hdf5_output) {
        $cmd .= ' --hdf5Output='. $self->hdf5_output;
        push @output_files, $self->hdf5_output;
    }
    $cmd .= ' '. $self->cmp_hdf5_file;
 
    if ($self->consensus_file) {
        $cmd .= ' > '. $self->consensus_file;
        push @output_files, $self->consensus_file;
    }
    $self->shellcmd(
        cmd => $cmd,
        input_files => [$self->cmp_hdf5_file],
        output_files => \@output_files,
        skip_if_output_is_present => 0,
    );
    return 1;
}


1;
