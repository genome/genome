package Genome::ProcessingProfile::SimpleAlignment;

use strict;
use warnings;

use Genome;

class Genome::ProcessingProfile::SimpleAlignment {
    is => 'Genome::ProcessingProfile::Staged',
    has_param => [
        reference_sequence_name => {
                                    doc => 'identifies the reference sequence used in the model(required if no prior_ref_seq)',
                                    is_optional => 1,
        },
        #picard_version => {
        #                      doc => 'picard version for MarkDuplicates, MergeSamfiles, CreateSequenceDictionary...',
        #                      is_optional => 1,
        #},
        #samtools_version => {
        #                      doc => 'samtools version for SamToBam, samtools merge, etc...',
        #                      is_optional => 1,
        #},
        #read_aligner_name => {
        #                      doc => 'alignment algorithm/software used for this model',
        #},
        #read_aligner_version => {
        #                         doc => 'the aligner version used for this model',
        #                         is_optional => 1,
        #},
        #read_aligner_params => {
        #                        doc => 'command line args for the aligner',
        #                        is_optional => 1,
        #},
    ],
};

sub stages {
    return (qw/
            shotgun
            /);
}

sub shotgun_job_classes {
    return (qw/
            Genome::Model::Event::Build::SimpleAlignmentWorkflow
        /);
}

sub shotgun_objects {
    return 1;
}


1;
