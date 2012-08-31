package Genome::Model::Tools::ViromeEvent::RepeatMasker::RemoveFile;

use strict;
use warnings;

use Genome;
use Workflow;
use IO::File;
use File::Basename;

class Genome::Model::Tools::ViromeEvent::RepeatMasker::RemoveFile{
    is => 'Genome::Model::Tools::ViromeEvent',
};

sub help_brief {
    "gzhao's Repeat Masker remove files";
}

sub help_detail {
    "This script will check the result of RepeatMasker which should generate .masked file in the given directory.";
}

sub execute {
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
 
    foreach my $file (@fastas) {
	unlink $file.'.cat' if -e $file.'.cat';
	unlink $file.'.tbl' if -e $file.'.tbl';
	unlink $file.'.ref' if -e $file.'.ref';
	unlink $file.'.out' if -e $file.'.out';
    }

    $self->log_event("Remove unnecessary file for $sample_name completed");

    return 1;
}

1;

