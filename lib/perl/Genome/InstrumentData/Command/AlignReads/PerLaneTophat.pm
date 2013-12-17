package Genome::InstrumentData::Command::AlignReads::PerLaneTophat;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command::AlignReads::PerLaneTophat {
    is => ['Genome::InstrumentData::Command::AlignReads'],
    has_param => [
        lsf_resource => {
            calculate => q{
                $self->bsub_rusage;
            },
            default_value => &_fallback_lsf_resource, # workflow doesn't support varying this per instance
        },
    ]
};

sub _fallback_lsf_resource {
    my $mem_mb = 1024 * 28; # increased b/c we have about 16 GB available when 6 jobs run on a 96 Gb server
    my $cpus = 12;

    my $mem_kb = $mem_mb*1024;
    my $tmp_gb = 300;

    my $user = getpwuid($<);
    my $queue = $ENV{GENOME_LSF_QUEUE_ALIGNMENT_DEFAULT};

    # TODO: in-house TGI concepts like alignment_pd shouldn't be methods on the generic config module :( -ssmith
    $queue = $ENV{GENOME_LSF_QUEUE_ALIGNMENT_PROD} if (Genome::Config->can('should_use_alignment_pd') and Genome::Config->should_use_alignment_pd);

    #my $host_groups;
    #my $command = qq(bqueues -l $queue | grep ^HOSTS:);
    #$host_groups = qx($command);
    #$host_groups =~ s/\/\s+/\ /;
    #$host_groups =~ s/^HOSTS:\s+//;

    my $select  = "select[ncpus >= $cpus && mem >= $mem_mb && gtmp >= $tmp_gb] span[hosts=1]";
    my $rusage  = "rusage[mem=$mem_mb, gtmp=$tmp_gb]";
    my $options = "-M $mem_kb -n $cpus -q $queue";

    my $required_usage = "-R '$select $rusage' $options";

    #check to see if our resource requests are feasible (This uses "maxmem" to check theoretical availability)
    #factor of four is based on current six jobs per host policy this should be revisited later
    #my $select_check = "select[ncpus >= $cpus && maxmem >= " . ($mem_mb * 4) . " && gtmp >= $tmp_gb] span[hosts=1]";

    #$command = qq(bhosts -R '$select_check' $host_groups | grep ^blade);
    #$command = qq(bhosts -R '$select_check' $host_groups );
    #my @selected_blades = qx($command);

    #if (@selected_blades) {
        return $required_usage;
    #} else {
    #    die __PACKAGE__->error_message("Failed to find hosts that meet resource requirements ($required_usage). [Looked with $select_check]");
    #}
}



1;
