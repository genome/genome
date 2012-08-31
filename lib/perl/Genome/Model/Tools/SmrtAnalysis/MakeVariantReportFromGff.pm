package Genome::Model::Tools::SmrtAnalysis::MakeVariantReportFromGff;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::SmrtAnalysis::MakeVariantReportFromGff {
    is  => 'Genome::Model::Tools::SmrtAnalysis::Base',
    has_input => [
        alignment_summary_gff_file => {
            is => 'Text',
            doc => 'The alignment summary GFF3 format file.',
        },
        variants_gff_file => {
            is => 'Text',
            doc => 'The variants in GFF3 format.',
        },
    ],
    has_optional_input => [
        report_xml_file => {
            is => 'Text',
            doc => 'The output xml file to write report to.',
            is_output => 1,
        },
        output_dir => {
            is => 'Text',
            doc => 'Output directory for associated files',
        },
        dots_per_inch => {
            is => 'Number',
            doc => 'dots/inch',
        },
        dumpdata => {
            is => 'Boolean',
            doc => 'Flag to dump the data backing charts',
            default_value => 0,
        },
        reference_directory => {
            is => 'Text',
            doc => 'The reference entry directorymakeCoverageReportFromGff.py',
        },
    ],
};

sub help_brief {
    'Make the reports for variants.'
}


sub help_detail {
    return <<EOS 

EOS
}

sub execute {
    my $self = shift;
    my $cmd = $self->analysis_bin .'/makeVariantReportFromGff.py';
    if ($self->dumpdata) {
        $cmd .= ' --dumpdata';
    }
    if (defined($self->dots_per_inch)) {
        $cmd .= ' --dpi='. $self->dots_per_inch;
    }
    if (defined($self->reference_directory)) {
        $cmd .= ' --refdir='. $self->reference_directory;
    }

    if (defined($self->output_dir)) {
        unless (-d $self->output_dir) {
            unless (Genome::Sys->create_directory($self->output_dir)) {
                die('Failed to create output directory: '. $self->output_dir);
            }
        }
        $cmd .= ' --output='. $self->output_dir;
        unless ($self->report_xml_file) {
            $self->report_xml_file($self->output_dir .'/coverage.xml');
        }
    }
    $cmd .= ' '. $self->alignment_summary_gff_file .' '. $self->variants_gff_file .' > '. $self->report_xml_file;
    $self->shellcmd(
        cmd => $cmd,
        input_files => [$self->alignment_summary_gff_file,$self->variants_gff_file],
        output_files => [$self->report_xml_file],
        skip_if_output_is_present => 0,
    );
    return 1;
}


1;
