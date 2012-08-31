package Genome::Model::Tools::Fastq::Trimq2;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Fastq::Trimq2 {
    is => 'Command',
    is_abstract  => 1,
    has_optional => [
        length_limit => {
            is  => 'Integer',
            doc => 'length limit of first Q2 position to 5 prime end (left), the reads fastq failing this limit will be thrown into xxxx.trimq2.filtered.fastq, default is 32',
            default => 32,
        },
        output_dir => {
            is  => 'Text',
            doc => 'The directory that report, filtered.fastq, trimmed_original_fastq will be written to, default is the same path as fastq_file',
        },
        trim_string => {
            is      => 'String',
            default => '#',
            doc     => 'quality trim string like #, B, i.. There are two types of fastq, sanger fastq quality = Qphred + 33, Solexa fastq quality = Qphred + 64. default is sanger format',
        },
        report_file => {
            is  => 'Text',
            doc => 'the report file name with path, default is trimq2.report in fastq_file dir',
            is_optional => 1,
        },
    ],
};


sub help_brief {
    'Tools to trim phred quality 2 sequences and qualities'
}


sub help_detail { 
    'Tools to trim phred quality 2 sequences and qualities'
}


sub create {
    my $class = shift;
    my $self  = $class->SUPER::create(@_);
    my $dir   = $self->output_dir;
    
    if ($dir) {
        unless (-d $dir) {
            $self->error_message("output dir provided does not exist: $dir");
            return;
        }
    }

    unless ($self->trim_string =~ /\S{1}/) {
        $self->error_message("quality trim string must be one character");
        return;
    }
    return $self;
}

1;


