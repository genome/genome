package Genome::Model::Tools::TechD::PicardDuplicationRatios;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::TechD::PicardDuplicationRatios{
    is => ['Command'],
    has => {
        build_id => {
            is => 'Number',
            doc => 'The build id used to resolve the path to the Picard MarkDuplicates metrics file of interest.',
        },
    },
};

sub help_brief {
    'A temporary tool to print Picard MarkDuplicates duplication metrics for a build'
}
sub help_synopsis {
    'Print build duplication metrics to terminal'
}
sub help_detail{
    return <<"EOS"
All builds currently use Picard MarkDuplicates to flag PCR duplicates after merging the individual instrument data alignments.
One output from Picard MarkDuplicates is a metrics file that provides a per library breakdown of PCR and optical duplicate reads.  This tool simply parses the metrics file and dumps a hash reference of the per library breakdown.
EOS
}

sub execute {
    my $self = shift;
    my $build = Genome::Model::Build->get($self->build_id);
    unless ($build) {
        die('Failed to find build by id '. $self->build_id);
    }
    my $subject = $build->mark_duplicates_library_metrics_hash_ref;
    print Data::Dumper::Dumper($subject);
    return 1;
}

1;
