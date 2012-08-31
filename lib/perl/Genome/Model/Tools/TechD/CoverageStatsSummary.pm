package Genome::Model::Tools::TechD::CoverageStatsSummary;

use strict;
use warnings;

use Genome;


class Genome::Model::Tools::TechD::CoverageStatsSummary{
    is => ['Command'],
    has => {
        build_id => { 
            is => 'Number',
            doc => 'The build id used to resolve the path to the RefCov stats summary file of interest.',},
    },
};

sub help_brief {
    'A temporary tool to print a RefCov stats summary for a build'
}
sub help_synopsis {
    'Print RefCov stats summary to terminal'
}
sub help_detail{
    return <<"EOS"
All capture builds currently use RefCov to evaluate coverage metrics including depth, breadth, contiguity, etc. of eacth region of interest.
RefCov will output a STATS file that contains one line per region of interest and minimum depth filter.
Typically coverage of a region is assessed at 1X, 5X, 10X, etc. minimum depth filters.  A summary file with means, medians and variance are produced as well.  This tool simply parses the summary file and dumps a hash reference of metrics to the terminal.
EOS
}

sub execute {
    my $self = shift;
    my $build = Genome::Model::Build->get($self->build_id);
    unless ($build) {
        die('Failed to find build by id '. $self->build_id);
    }
    my $stats_summary = $build->coverage_stats_summary_hash_ref;
    print Data::Dumper::Dumper($stats_summary);
    return 1;
}
