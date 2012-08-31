package Genome::Model::Tools::SmrtAnalysis::MakeFilterStatsReport;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::SmrtAnalysis::MakeFilterStatsReport {
    is  => 'Genome::Model::Tools::SmrtAnalysis::Base',
    has_input => [
        filter_csv_file => {
            is => 'Text',
            doc => 'The csv file output from FilterPlsH5.',
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
        skip_pre_plots => {
            is => 'Boolean',
            doc => 'Flag to skip generation of pre-filter histograms.',
            default_value => 1,
        },
    ],
    has_optional_param => [
        lsf_resource => { default_value => "-g /pacbio/smrtanalysis -M 32000000 -R 'select[type==LINUX64 && mem>=32000 && tmp>=160000] rusage[mem=32000,tmp=80000]'" },
    ],
};

sub help_brief {
    'Make the reports for filtered reads/subreads.'
}


sub help_detail {
    return <<EOS 
Generates Post-Filter statistics for aligned read length, accuracy, and
optionally z-scores.     Takes as input a cmp.h5 file.
EOS
}

sub execute {
    my $self = shift;
    my $cmd = $self->analysis_bin .'/makeFilterStatsReport.py';
    if ($self->skip_pre_plots) {
        $cmd .= ' --skipPrePlots';
    }
    if ($self->dumpdata) {
        $cmd .= ' --dumpdata';
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
            $self->report_xml_file($self->output_dir .'/filterReports_filterStats.xml');
        }
    }
    $cmd .= ' '. $self->filter_csv_file .' > '. $self->report_xml_file;
    $self->shellcmd(
        cmd => $cmd,
        input_files => [$self->filter_csv_file],
        output_files => [$self->report_xml_file],
        skip_if_output_is_present => 0,
    );
    return 1;
}


1;
