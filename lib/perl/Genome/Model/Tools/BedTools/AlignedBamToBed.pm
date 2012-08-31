package Genome::Model::Tools::BedTools::AlignedBamToBed;

use strict;
use warnings;

use Genome;


class Genome::Model::Tools::BedTools::AlignedBamToBed {
    is => 'Genome::Model::Tools::BedTools',
    has_input => [
        input_file => {
            is => 'Text',
            doc => 'The input BAM file',
        },
        output_file => {
            is => 'Text',
            doc => 'The output BED file containing only aligned reads',
        },
    ],
};

sub help_brief {
    'Converts BAM files to BED format.',
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
  gmt bed-tools bam-to-bed --input-file=example.bam --output-file=example.bed
EOS
}

sub help_detail {
    return <<EOS
More information about the BedTools suite of tools can be found at http://code.google.com/p/bedtools/. 
EOS
}


sub execute {
    my $self = shift;
    my $executable = $self->bedtools_path .'/bin/bamToBed';
    unless (-e $executable) {
        die('Failed to find executable at: '. $executable);
    }
    my $input_file = $self->input_file;
    my $output_file = $self->output_file;
    #my $params = '';
    my $samtools_cmd = 'samtools view -bF 0x4 '.$input_file;
    #samtools view -bF 0x4 example.bam | /gsc/pkg/bio/bedtools/BEDTools-2.10.0-64/bin/bamToBed -i stdin > example.bed
    my $bedtools_cmd = $executable .' -i stdin '. '> ' . $output_file;
    
    my $cmd = $samtools_cmd .' | '.$bedtools_cmd;
    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$input_file],
        output_files => [$output_file],
        skip_if_output_is_present => 0,
    );
    return 1;
}

1;
