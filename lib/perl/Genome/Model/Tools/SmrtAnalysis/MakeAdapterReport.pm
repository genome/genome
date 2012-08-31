package Genome::Model::Tools::SmrtAnalysis::MakeAdapterReport;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::SmrtAnalysis::MakeAdapterReport {
    is  => 'Genome::Model::Tools::SmrtAnalysis::Base',
    has_input => [
        hdf5_fofn => {
            is => 'Text',
            doc => 'A FOFN of pls.h5 or bas.h5 files.',
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
        debug => {
            is => 'Boolean',
            doc => 'Enable debug output',
            default_value => 0,
        },
    ],
    has_optional_param => [
        lsf_resource => { default_value => "-g /pacbio/smrtanalysis -M 16000000 -R 'select[type==LINUX64 && mem>=16000 && tmp>=160000] rusage[mem=16000,tmp=80000]'" },
    ],
};

sub help_brief {
    'Make the reports for adapter dimers.'
}


sub help_detail {
    return <<EOS 
Generates adapter related statistics which are helpful in quantifying the
presence of adapter dimers.
EOS
}

sub execute {
    my $self = shift;
    my $cmd = $self->analysis_bin .'/makeAdapterReport.py';
    if ($self->debug) {
        $cmd .= ' --debug';
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
            $self->report_xml_file($self->output_dir .'/filterReports_adapters.xml');
        }
    }
    $cmd .= ' '. $self->hdf5_fofn .' > '. $self->report_xml_file;
    $self->shellcmd(
        cmd => $cmd,
        input_files => [$self->hdf5_fofn],
        output_files => [$self->report_xml_file],
        skip_if_output_is_present => 0,
    );
    return 1;
}


1;
