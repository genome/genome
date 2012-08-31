package Genome::Model::Tools::454::Sff2Fastq;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::454::Sff2Fastq {
    is => ['Command'],
    has => [
            sff_file => {
                            is => 'string',
                            doc => 'The input SFF file',
                     },
            fastq_file => {
                            is => 'string',
                            doc => 'The output fastq file to generate',
                        },
            force_untrimmed => {
                       is => 'boolean',
                       doc => 'Output untrimmed data (don\'t use 454 default trimming scheme)',
                       is_optional => 1,
		   },
        ],
};

sub help_brief {
    "constructs a single fastq from a 454 SFF file"
}

sub help_detail {
    return <<EOS
    see 'sff2fastq' usage for valid params
EOS
}

sub execute {
    my $self = shift;


    my $cmd = sprintf("sff2fastq %s -o %s %s", ($self->force_untrimmed?'-n' : ''), $self->fastq_file, $self->sff_file);

    Genome::Sys->shellcmd(
                                         cmd => $cmd,
                                         input_files => [$self->sff_file],
                                         output_files => [$self->fastq_file],
                                     );
    return 1;
}

1;


