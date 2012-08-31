package Genome::Model::Tools::SmrtAnalysis::MakeControlReport;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::SmrtAnalysis::MakeControlReport {
    is  => 'Genome::Model::Tools::SmrtAnalysis::Base',
    has_input => [
        cmp_hdf5_file => {
            doc => 'The aligned reads in the format of cmp.h5',
        },
        filtered_summary_csv_file => {
            doc => 'The filter summary csv file output from FilterPlsH5',
        },
    ],
    has_optional_input => [
        report_xml_file => {
            doc => 'The output summary report in xml format.',
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
        debug => {
            is => 'Boolean',
            doc => 'Enable debug output',
            default_value => 0,
        },
    ],
};

sub help_brief {
    'Create the spike-in control report.'
}


sub help_detail {
    return <<EOS 
Create the spike-in control report.
EOS
}

sub execute {
    my $self = shift;
    my $cmd = $self->analysis_bin .'/makeControlReport.py';
    if ($self->debug) {
        $cmd .= ' --debug';
    }
    if (defined($self->dots_per_inch)) {
        $cmd .= ' --dpi='. $self->dots_per_inch;
    }
    if (defined($self->output_dir)) {
        unless (-d $self->output_dir) {
            unless (Genome::Sys->create_directory($self->output_dir)) {
                die('Failed to create output directory: '. $self->output_dir);
            }
        }
        $cmd .= ' --output='. $self->output_dir;
        unless ($self->report_xml_file) {
            $self->report_xml_file($self->output_dir .'/controlReport.xml');
        }
    }
    $cmd .= ' '. $self->filtered_summary_csv_file .' '. $self->cmp_hdf5_file .' > '. $self->report_xml_file;
    $self->shellcmd(
        cmd => $cmd,
        input_files => [$self->filtered_summary_csv_file,$self->cmp_hdf5_file],
        output_files => [$self->report_xml_file],
        skip_if_output_is_present => 0,
    );
    return 1;
}


1;
