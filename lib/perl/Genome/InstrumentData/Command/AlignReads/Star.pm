package Genome::InstrumentData::Command::AlignReads::Star;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command::AlignReads::Star {
    is => ['Genome::InstrumentData::Command::AlignReads'],
    has_param => [
        lsf_resource => {
            calculate => q{
                $self->bsub_rusage;
            },
            default_value => &_fallback_lsf_resource, 
        },
    ]
};

sub _fallback_lsf_resource {
    my $mem_mb = 1024 * 32;
    my $mem_kb = $mem_mb*1024;
    my $cpus   = 12;

    my $select  = "select[ncpus >= $cpus && mem >= $mem_mb] span[hosts=1]";
    my $rusage  = "rusage[mem=$mem_mb]";
    my $options = "-M $mem_kb -n $cpus";

    my $required_usage = "-R '$select $rusage' $options";

    return $required_usage;
}


1;
