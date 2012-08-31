package Genome::Model::Tools::BedTools::BedToBam;

use strict;
use warnings;

use Genome;

my $DEFAULT_GENOME = '/gscmnt/gc4096/info/model_data/2741951221/build101947881/genome.tsv';
my $DEFAULT_BED12 = 0;

class Genome::Model::Tools::BedTools::BedToBam {
    is => 'Genome::Model::Tools::BedTools',
    has_input => [
        input_file => {
            is => 'Text',
            doc => 'The input BED file',
        },
        output_file => {
            is => 'Text',
            doc => 'The output BAM file',
        },
        genome => {
            is => 'Text',
            doc => 'The genome reference.  See BEDTools User Manual for documentation of format.',
            is_optional => 1,
            default_value => $DEFAULT_GENOME,
        },
        bed12 => {
            is => 'Boolean',
            doc => 'The BED file is in BED12 format.  The BAM CIGAR string will reflect BED blocks.',
            default_value => $DEFAULT_BED12,
            is_optional => 1,
        },
        mapq => {
            is => 'Integer',
            doc => 'Set the mappinq quality for the BAM records. (INT) Default: 255',
            is_optional => 1,
        },
    ],
};

sub help_brief {
    'Converts BED/GFF/VCF features to BAM format.',
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
  gmt bed-tools bed-to-bam --input-file=example.bed --output-file=example.bam
EOS
}

sub help_detail {
    return <<EOS
More information about the BedTools suite of tools can be found at http://code.google.com/p/bedtools/. 
EOS
}

sub default_genome {
    return $DEFAULT_GENOME;
}

sub default_bed12 {
    return $DEFAULT_BED12;
}

sub execute {
    my $self = shift;
    my $executable = $self->bedtools_path .'/bin/bedToBam';
    unless (-e $executable) {
        die('Failed to find executable at: '. $executable);
    }
    my $params = '';
    if (defined($self->mapq)) {
        $params .= ' -mapq '. $self->mapq;
    }
    if ($self->bed12) {
        $params .= ' -bed12 ';
    }
    my $cmd = $executable .' '. $params .' -i '. $self->input_file  .' -g '. $self->genome .' > '. $self->output_file;
    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$self->input_file, $self->genome],
        output_files => [$self->output_file],
        skip_if_output_is_present => 0,
    );
    return 1;
}

1;
