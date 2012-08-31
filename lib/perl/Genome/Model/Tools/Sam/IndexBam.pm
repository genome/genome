package Genome::Model::Tools::Sam::IndexBam;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Sam::IndexBam {
    is  => 'Genome::Model::Tools::Sam',
    has_input => [
	bam_file => {
	    is  => 'Text',
	    doc => 'BAM input file to index',
	},
        bam_index_file => {
            is => 'Text',
            is_optional => 1,
            doc => 'The BAM index file to create',
        }
    ],
};

sub help_brief {
    'Tool to index BAM files';
}

sub help_detail {
    return <<EOS
    Tool to index BAM files.
EOS
}

sub execute {
    my $self = shift;

    my $bam_index_cmd = $self->samtools_path .' index '. $self->bam_file;
    if ($self->bam_index_file) {
        $bam_index_cmd .= ' '. $self->bam_index_file;
    } else {
        $self->bam_index_file($self->bam_file.'.bai');
    }
    Genome::Sys->shellcmd(
        cmd => $bam_index_cmd,
        input_files => [$self->bam_file],
        output_files => [$self->bam_index_file],
        skip_if_output_is_present => 0
    );
    return 1;
}


1;
