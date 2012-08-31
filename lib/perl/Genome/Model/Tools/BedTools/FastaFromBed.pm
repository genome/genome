package Genome::Model::Tools::BedTools::FastaFromBed;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::BedTools::FastaFromBed {
    is => 'Genome::Model::Tools::BedTools',
    has_input => [
        input_fasta => {
            is => 'Text',
            doc => 'Input FASTA file. DEFAULT GENOME: NCBI-human-build36.',
            default_value => '/gscmnt/gc4096/info/model_data/2741951221/build101947881/all_sequences.fa',
        },
        bed_file => {
            is => 'Text',
            doc => 'BED/GFF/VCF file of ranges to extract from input_fasta',
        },
        output_file => {
            is => 'Text',
            doc => 'Output file (can be FASTA or TAB-delimited)',
        },
        use_name => {
            is => 'Boolean',
            doc => 'Use the name column in the BED file for the FASTA headers in the output FASTA file.',
            is_optional => 1,
            default_value => 0,
        },
        output_format => {
            is => 'Text',
            doc => 'Report extract sequences in a tab-delimited format or a FASTA format.',
            is_optional => 1,
            default_value => 'fasta',
            valid_values => ['fasta','tab'],
        },
        force_strandedness => {
            is => 'Boolean',
            doc => 'If the feature occupies the antisense strand, the sequence will be reverse complemented.',
            is_optional => 1,
            default_value => 0,
        },
    ],
};

sub help_brief {
    'Extracts sequences from a FASTA file for each of the intervals defined in a BED file.',
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
  gmt bed-tools fasta-from-bed --input-bed=example.bed --input-fasta=example_in.fa --output-fasta=example_out.fa
EOS
}

sub help_detail {
    return <<EOS
More information about the BedTools suite of tools can be found at http://code.google.com/p/bedtools/. 
EOS
}

sub execute {
    my $self = shift;
    my $executable = $self->bedtools_path .'/bin/fastaFromBed';
    unless (-e $executable) {
        die('Failed to find executable at: '. $executable);
    }
    my $params;
    if ($self->use_name) {
        $params .= ' -name ';
    }
    if ($self->output_format eq 'tab') {
        $params .= ' -tab ';
    }
    if ($self->force_strandedness) {
        $params .= ' -s ';
    }
    my $cmd = $executable .' '. $params .' -fi '. $self->input_fasta  .' -bed '. $self->bed_file .' -fo '. $self->output_file;
    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$self->bed_file, $self->input_fasta],
        output_files => [$self->output_file],
        skip_if_output_is_present => 0,
    );
    return 1;
}

1;
