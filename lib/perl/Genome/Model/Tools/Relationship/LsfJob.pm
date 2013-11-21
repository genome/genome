package Genome::Model::Tools::Relationship::LsfJob;

use strict;
use warnings;
use Data::Dumper;
use Genome;           
use Genome::Info::IUB;
use POSIX;
use Cwd;
use File::Basename;
use File::Path;
use LSF::Job;

class Genome::Model::Tools::Relationship::LsfJob {
    is => 'Command',
    has_input => [
    command_line=>{
        is_optional=>0,
    },
    parents=> {
        is_optional=>1,
    },
    children=> {
        is_optional=>1,
    },
    status=> {
        default=>'Created',
    },
    queue=> {
        default => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER},
    },
    std_out=> {
    },
    std_err=>{
    },
    resource_request=>{
        is_optional=>1,
    },
    job_id=>{
        is_optional=>1,
    },
    lsf_job_object=>{
        is_optional=>1,
    },
    ],
};

sub help_brief {
}

sub help_detail {
}


sub execute {
    my $self = shift;
    my @args = ("-q " .  $self->queue);
    if($self->std_out) {
        push @args, "-o " . $self->std_out;
    }
    if($self->std_err) {
        push @args, "-e " . $self->std_err;
    }
    ##TODO: implement resource request
    my $cmd = "bsub " . join(" ", @args) . " " . $self->command_line;
    print $cmd . "\n";
    $self->job_id($self->parse_bsub_output(`$cmd`));
    if($self->job_id) {
        return $self;
    }
    else {
        return 0;
    }
}

sub bjob_l_is_job_successful {
    my $self = shift;
    my $job_id = $self->job_id;
    my @output = `bjobs -l $job_id`;
    for my $output_line (@output) {
        if ($output_line =~ m/Done successfully/) {
            return 1;
        }
    }
    return 0;
}


sub parse_bsub_output {
    my ($self, $bsub_return) = @_;
    my ($id) = ($bsub_return =~ m/<(\d+)>/);
    if($id) {
        return $id;
    }
    else {
        return 0;
    }
}

