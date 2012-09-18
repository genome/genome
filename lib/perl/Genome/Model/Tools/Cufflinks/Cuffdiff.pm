package Genome::Model::Tools::Cufflinks::Cuffdiff;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Cufflinks::Cuffdiff {
    is => 'Genome::Model::Tools::Cufflinks',
    has_input => [
        transcript_gtf_file => {
            doc => 'A GTF format file of transcript annotation to perform differential expression tests.',
        },
        bam_file_paths => {
            doc => 'A string of bam file paths separated by white space. Supply replicate SAMs as comma separated lists for each condition: sample1_rep1.sam,sample1_rep2.sam,...sample1_repM.sam',
        },
        output_directory => {
            doc => 'The directory to write output files.',
        },
        params => {
            is_optional => 1,
            doc => 'A text string of optional parameters to pass to cuffdiff',
        },
    ],
};

sub help_synopsis {
    return <<EOS
gmt cufflinks cuffdiff --transcript-gtf-file=?  --bam-files=...
EOS
}

sub help_brief {
    "A wrapper around the cuffdiff command";
}

sub help_detail {
    return <<EOF
Cufflinks includes a program, "Cuffdiff", that you can use to find significant changes in transcript expression, splicing, and promoter use.

More information about cuffdiff can be found at http://cufflinks.cbcb.umd.edu/.
EOF
}

sub execute {
    my $self = shift;

    my $cmd = $self->cuffdiff_path .' '. ($self->params || '') .' -o '. $self->output_directory .' '. $self->transcript_gtf_file .' '. $self->bam_file_paths;
    Genome::Sys->shellcmd(
        cmd => $cmd,
    );
    return 1;
}

1;
