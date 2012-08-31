package Genome::Model::Tools::ViromeEvent::CDHIT::RemoveFile;

use strict;
use warnings;

use Genome;
use Workflow;
use IO::File;
use File::Basename;

class Genome::Model::Tools::ViromeEvent::CDHIT::RemoveFile{
    is => 'Genome::Model::Tools::ViromeEvent',
};

sub help_brief {
    "Tool to remove selected files from cd-hit run directory";
}

sub help_detail {
    'Tool to remove selected files from cd-hit run directory';
}

sub execute {
    my $self = shift;
    my $dir = $self->dir;
    my $sample_name = basename($dir);
    my @files = glob ("$dir/*.bak.clstr");
    if ( not @files ) {
        $self->log_event("No cd-hit clstr files to remove for sample: $sample_name");
        return 1;
    }
    foreach ( @files ) {
	unlink $_;
    }
    $self->log_event("Removed temp files for sample: $sample_name");

    return 1;
}

1;
