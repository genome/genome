package Genome::VariantReporting::Process::Trio;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Process::Trio {
    is => 'Genome::Process',
    has_input => [
        builds => {
            is => 'Genome::Model::Build::SomaticValidation',
            is_many => 1,
        },
        coverage_builds => {
            is => 'Genome::Model::Build::SomaticValidation',
            is_many => 1,
            is_optional => 1,
        },
        tumor_sample => {
            is => 'Genome::Sample',
        },
        followup_sample => {
            is => 'Genome::Sample',
        },
        normal_sample => {
            is => 'Genome::Sample',
        },
    ],
};

sub get_reports_structure {
    my $self = shift;

    my @report_users = $self->result_users('label like' => '%.%.%');
    use Data::Dump qw(pp);
    print "#\e[34mDEBUG:\e[0m \e[31m" . '\@report_users' . "\e[0m " . pp(\@report_users) . "\n";
    return 1;
}

1;
