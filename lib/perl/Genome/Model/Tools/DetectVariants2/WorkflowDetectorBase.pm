package Genome::Model::Tools::DetectVariants2::WorkflowDetectorBase;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::DetectVariants2::WorkflowDetectorBase {
    doc => "A base class for detectors that run workflows",
    is_abstract => 1,
    is => ['Genome::Model::Tools::DetectVariants2::Detector'],
    has => [
        chromosome_list => {
            is => 'ARRAY',
            is_optional => 1,
            doc => 'list of chromosomes to run on. This will run on a default set of chromosomes from the reference sequence if not set.',
        },
    ],
    has_param => [
        lsf_queue => {
            default_value => 'workflow'
        },
    ],
    has_transient_optional => [
        _workflow_result => {
            doc => 'Result of the workflow',
        },
    ],
};

sub chromosome_sort {
    # numeric sort if chrom starts with number
    # alphabetic sort if chrom starts with non-number
    # numeric before alphabetic
    # Not perfect, NT_1 will sort the same as NT_100
    my ($a_chrom) = $a =~ /^(\S+)/;
    my ($b_chrom) = $b =~ /^(\S+)/;
    ($a_chrom, $b_chrom) = (uc($a_chrom), uc($b_chrom));
    my $a_is_numeric = ($a_chrom =~ /^\d+$/ ? 1 : 0);
    my $b_is_numeric = ($b_chrom =~ /^\d+$/ ? 1 : 0);

    my $rv;
    if ($a_chrom eq $b_chrom) {
        $rv = 0;
    }
    elsif ($a_is_numeric && $b_is_numeric) {
        $rv = ($a_chrom <=> $b_chrom);
    }
    elsif ($a_is_numeric && !$b_is_numeric) {
        $rv = -1;
    }
    elsif ($b_is_numeric && !$a_is_numeric) {
        $rv = 1;
    }
    elsif ($a_chrom eq 'X') {
        $rv = -1;
    }
    elsif ($b_chrom eq 'X') {
        $rv = 1;
    }
    elsif ($a_chrom eq 'Y') {
        $rv = -1;
    }
    elsif ($b_chrom eq 'Y') {
        $rv = 1;
    }
    elsif ($a_chrom eq 'MT') {
        $rv = -1;
    }
    elsif ($b_chrom eq 'MT') {
        $rv = 1;
    }
    else {
        $rv = ($a_chrom cmp $b_chrom);
    }

    return $rv;
}

1;

