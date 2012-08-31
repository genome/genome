package Genome::Model::Tools::SmrtAnalysis::GffToVcf;

use strict;
use warnings;

use Genome;
use File::Basename;

class Genome::Model::Tools::SmrtAnalysis::GffToVcf {
    is => ['Genome::Model::Tools::SmrtAnalysis::Base'],
    has_input => {
        gff_file => {
            doc => 'The GFF3 format file to convert.',
            is_output => 1,
        },
    },
    has_optional_input => {
        vcf_file => {
            is => 'Text',
            doc => 'The output VCF format file.',
            is_output => 1,
        },
        global_reference => {
            is => 'Text',
            doc => 'Name of global reference to put in Meta field',
        },
    },
};

sub help_brief {
    'Utility for converting GFF3 to BED format. Currently supports regional coverage or variant .bed output.'
}

sub help_detail {
    return <<EOS
Utility for converting GFF3 to BED format. Currently supports regional coverage or variant .bed output.
EOS
}

sub execute {
    my $self = shift;

    my $cmd = $self->analysis_bin .'/gffToVcf.py';
    if (defined($self->global_reference)) {
        $cmd .= ' --globalReference="'. $self->global_reference .'"';
    }
    unless ($self->vcf_file) {
        my ($basename,$dirname,$suffix) = File::Basename::fileparse($self->gff_file,qw/\.gff/);
        $self->vcf_file($dirname .'/'. $basename .'.vcf');
    }
    $cmd .= ' '. $self->gff_file .' > '. $self->vcf_file;
    $self->shellcmd(
        cmd => $cmd,
        input_files => [$self->gff_file],
        output_files => [$self->vcf_file],
        skip_if_output_is_present => 0,
    );
    return 1;
};


1;
