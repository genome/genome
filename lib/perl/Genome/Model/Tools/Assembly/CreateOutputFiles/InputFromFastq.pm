package Genome::Model::Tools::Assembly::CreateOutputFiles::InputFromFastq;

use strict;
use warnings;

use Genome;

use IO::File;

class Genome::Model::Tools::Assembly::CreateOutputFiles::InputFromFastq {
    is => 'Genome::Model::Tools::Assembly::CreateOutputFiles',
    has => [
	fastq_file => {
	    is => 'Text',
	    doc => 'Input fastq file for the assembly',
	},
	directory => {
	    is => 'Text',
	    doc => 'Assembly data directory',
	},
    ],
};

sub help_brief {
    'Tool to create input fasta and qual from fastq file for stats',
}

sub help_detail {
    "Tool to create fasta and quality files from velvet input fastq file";
}

sub execute {
    my $self = shift;

    my $root_name = File::Basename::basename($self->fastq_file);
    $root_name =~ s/\.fastq//;

    my $fasta_file = $self->directory.'/edit_dir/'.$root_name.'.fasta';
    my $qual_file = $self->directory.'/edit_dir/'.$root_name.'.fasta.qual';
    my $fq = Genome::Model::Tools::Sx->create(
        input => [ $self->fastq_file ],
        output => [ $fasta_file.':qual_file='.$qual_file ],
    );
    if ( not $fq ) {
        $self->error_message('Failed to create fastq to fasta converter');
        return;
    }
    $fq->dump_status_messages(1);
    if ( not $fq->execute ) {
        $self->error_message('Failed to execute fastq to fasta converter');
        return;
    }

    if ( not -e $fasta_file ) {
        $self->error_message('Executed fastq to fasta converter, but fasta file does not exist: '.$fasta_file);
        return;
    }
    if ( not -e $qual_file ) {
        $self->error_message('Executed fastq to fasta converter, but fasta file does not exist: '.$qual_file);
        return;
    }

    unlink $fasta_file.'.gz';
    unlink $qual_file.'.gz';

    if (system("gzip $fasta_file $qual_file")) {
        $self->error_message("Failed to zip files: $fasta_file $qual_file");
        return;
    }

    return 1;
}

1;
