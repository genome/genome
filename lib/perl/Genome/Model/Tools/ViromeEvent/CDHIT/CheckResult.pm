package Genome::Model::Tools::ViromeEvent::CDHIT::CheckResult;

use strict;
use warnings;

use Genome;
use Workflow;
use IO::File;

class Genome::Model::Tools::ViromeEvent::CDHIT::CheckResult{
    is => 'Genome::Model::Tools::ViromeEvent',
};

sub help_brief {
    "Tool to run cd-hit for virome screening pipeline"
}

sub help_detail {
    'Tool to run cd-hit for virome screening pipline';
}

sub execute {
    my $self = shift;
    my $dir = $self->dir;
    my $sample_name = File::Basename::basename ($dir);

    my $faFile = $dir.'/'.$sample_name.'.fa';
    my $cdhitReport = $faFile.'.cdhitReport';

    if (-s $cdhitReport) {
        my $result = `grep completed $cdhitReport`;
        if ($result =~ /program completed/) {
	    $self->log_event("Already completed for sample: $sample_name");
	    return 1;
        }
    }

    my $archos = `uname -a`;
    my $cd_hit_dir = ($archos =~ /64/) ? '/cd-hit-64/cd-hit-est' : '/cd-hit-32/cd-hit-est';    

    my $com = '/gscmnt/sata835/info/medseq/virome/scripts_used_by_virome' . $cd_hit_dir .
	      ' -i ' . $faFile .
	      ' -o ' . $faFile . '.cdhit_out' .
	      ' -c 0.98 -n 8 -G 0 -aS 0.98 -g 1 -r 1 -M 4000 -d 0' .
	      ' > ' . $cdhitReport;

    $self->log_event("Executing for sample: $sample_name");

    if (system($com)) { #returns 0 when successful
        $self->log_event("Failed to successfully run cd-hit for sample: $sample_name");
        return; #need to make sure this exits
    }

    $self->log_event("Cd-hit ran successfully for sample: $sample_name");

    #remove files

    return 1;
}

1;
