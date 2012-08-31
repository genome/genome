package Genome::Model::Tools::SmrtAnalysis::StsPlots;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::SmrtAnalysis::StsPlots {
    is => ['Genome::Model::Tools::SmrtAnalysis::Base'],
    has_input => {
        csv_file => {
            doc => 'The input *.sts.csv movie file.',
        },
        pdf_file => {
            doc => 'The output PDF file.',
            is_output => 1,
        },
    },
};

sub help_brief {
    'Generate heatmap based on *.sts.csv movie files';
}

sub help_detail {
    return <<EOS
This tool runs a PacBio R script for generating a summary PDF file of quality/production type metrics.
NOTE: This requires R version 2.11.1 or greater.
EOS
}

sub execute {
    my $self = shift;

    my $csv_file = $self->csv_file;
    $csv_file =~ s/ /\\ /g;
    my $pdf_file = $self->pdf_file;
    $pdf_file =~ s/ /\\ /g;
    my $cmd = 'Rscript '. $self->r_lib ."/stsPlots.R $csv_file $pdf_file";
    $self->shellcmd(
        cmd => $cmd,
        input_files => [$self->csv_file],
        output_files => [$self->pdf_file],
        skip_if_output_is_present => 0,
    );
    return 1;
};


1;
