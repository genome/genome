package Genome::Model::Tools::SmrtAnalysis::GffToBed;

use strict;
use warnings;

use Genome;
use File::Basename;

class Genome::Model::Tools::SmrtAnalysis::GffToBed {
    is => ['Genome::Model::Tools::SmrtAnalysis::Base'],
    has_input => {
        gff_file => {
            doc => 'The GFF3 format file to convert.',
            is_output => 1,
        },
        purpose => {
            is => 'Text',
            doc => 'The purpose of the GFF3 conversion.',
            valid_values => ['coverage','variants'],
        },
    },
    has_optional_input => {
        bed_file => {
            is => 'Text',
            doc => 'The output BED format file.',
            is_output => 1,
        },
        name => {
            is => 'Text',
            doc => 'track name to display in header',
        },
        description => {
            is => 'Text',
            doc => 'track description to display in header',
        },
        use_score => {
            is => 'Text',
            doc => 'whether or not to use score for feature display',
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

    my $cmd = $self->analysis_bin .'/gffToBed.py';
    if (defined($self->use_score)) {
        $cmd .= ' --useScore='. $self->use_score;
    }
    if (defined($self->name)) {
        $cmd .= ' --name='. $self->name;
    }
    if (defined($self->description)) {
        $cmd .= ' --description="'. $self->description .'"';
    }
    unless ($self->bed_file) {
        my ($basename,$dirname,$suffix) = File::Basename::fileparse($self->gff_file,qw/\.gff/);
        $self->bed_file($dirname .'/'. $basename .'.bed');
    }
    $cmd .= ' '. $self->purpose .' '. $self->gff_file .' > '. $self->bed_file;
    $self->shellcmd(
        cmd => $cmd,
        input_files => [$self->gff_file],
        output_files => [$self->bed_file],
        skip_if_output_is_present => 0,
    );
    return 1;
};


1;
