package Genome::Model::Tools::TechD::AlignmentSummary;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::TechD::AlignmentSummary {
    is => ['Command'],
    has => {
        build_id => {
            is => 'Number',
            doc => 'The build id used to resolve the path to the alignment metrics file of interest.',
        },
    },
};

sub help_brief {
    'A temporary tool to print alignment metrics for a build'
}
sub help_synopsis {
    'Print build alignment metrics to terminal'
}
sub help_detail{
    return <<"EOS"
Each capture build generates an alignment summary file per wingspan value.
The alignment summary file contains basic alignment metrics as well as on/off target alignments.
Duplication rates and paired-end metrics are also included.  This tool simply parses the tsv file
and dumps the metrics to the a terminal.
EOS
}

sub execute {
    my $self = shift;
    my $build = Genome::Model::Build->get($self->build_id);
    unless ($build) {
        die('Failed to find build by id '. $self->build_id);
    }
    my $alignment_summary = $build->alignment_summary_hash_ref;
    print Data::Dumper::Dumper($alignment_summary);
    return 1;
}
