package Genome::Model::Tools::ViromeEvent::RepeatMasker::OuterCheckResult;

use strict;
use warnings;

use Genome;
use Workflow;
use IO::File;
use File::Basename;
use File::Copy;

class Genome::Model::Tools::ViromeEvent::RepeatMasker::OuterCheckResult{
    is => 'Genome::Model::Tools::ViromeEvent',
    has_output => [
        files_to_run=> { is => 'ARRAY',  doc => 'array of files for repeat masker', is_optional => 1},
    ]
};

sub help_brief {
    return "gzhao's Repeat Masker check result";
}

sub help_detail {
    return <<"EOS"
This script will check the result of RepeatMasker which
should generate .masked file in the given directory.
EOS
}

sub execute
{
    my $self = shift;

    my $dir = $self->dir;
    my $sample_name = basename ($dir);

    my $repeat_masker_dir = $dir.'/'.$sample_name.'.fa.cdhit_out_RepeatMasker';
    unless (-d $repeat_masker_dir) {
	$self->log_event("Failed to find repeat masker dir for sample: $sample_name");
	return;
    }

    my @fastas = glob ("$repeat_masker_dir/$sample_name*.fa");
    unless (scalar @fastas > 0) {
	$self->log_event("No input fastas for repeat masker run found for sample: $sample_name");
	return;
    }

    $self->files_to_run( \@fastas );

    $self->log_event("Completed check to run repeat masker for sample: $sample_name");

    return 1;
}

1;

