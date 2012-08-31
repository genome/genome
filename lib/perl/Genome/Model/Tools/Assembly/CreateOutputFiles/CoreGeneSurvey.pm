package Genome::Model::Tools::Assembly::CreateOutputFiles::CoreGeneSurvey;

use strict;
use warnings;

use Genome;
use Cwd;

class Genome::Model::Tools::Assembly::CreateOutputFiles::CoreGeneSurvey {
    is => 'Command',
    has => [
	core_gene_option => {
	    is => 'Text',
	    doc => 'Gene type, either bact or archaea',
	    valid_values => ['bact', 'archaea'],
	},
	subject_file => {
	    is => 'Text',
	    doc => 'Assembly fasta or ?? file like contigs.bases',
	},
    ],
};

sub execute {
    my $self = shift;

    #This has to run in the directory where subject file is located
    my $cwd = cwd();

    unless (-s $self->subject_file) {
	$self->error_message("Can't find subject file: ".$self->subject_file);
	return;
    }

    my $subject_file_name = File::Basename::basename($self->subject_file);
    my $dir_name = File::Basename::dirname($self->subject_file);

    chdir $dir_name;
    unlink 'Cov_30_PID_30.out.gz';

    my $cmd = "run_coregene_cov_pid_script $subject_file_name 30 0.3 -assembly -".$self->core_gene_option;
    if (system("$cmd")) {
	$self->error_message("run_coregene_cov_pid_script failed with command: $cmd");
	return;
    }

    chdir $cwd;

    return 1;
}

1;
