package Genome::Model::Tools::SmrtAnalysis::MakeTopVariantsReport;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::SmrtAnalysis::MakeTopVariantsReport {
    is  => 'Genome::Model::Tools::SmrtAnalysis::Base',
    has_input => [
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
        how_many => {
            is => 'Number',
            doc => 'How many variants to pull out of variants.gff',
        },
        sort_batch_size => {
            is => 'Number',
            doc => 'During reading of the source gff, perform intermediate sorts of record batches of given size',
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
    if (defined($self->how_many)) {
        $cmd .= ' --howMany='. $self->how_many;
    }
    if (defined($self->sort_batch_size)) {
        $cmd .= ' --sortBatchSize='. $self->sort_batch_size;
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
    $cmd .= ' '. $self->variants_gff_file .' > '. $self->report_xml_file;
    $self->shellcmd(
        cmd => $cmd,
        input_files => [$self->variants_gff_file],
        output_files => [$self->report_xml_file],
        skip_if_output_is_present => 0,
    );
    return 1;
}


1;
