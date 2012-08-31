package Genome::Model::Tools::SmrtAnalysis::MakeMappingStatsReport;

use strict;
use warnings;

use Genome;

my $DEFAULT_LSF_RESOURCE = "-g /pacbio/smrtanalysis -M 8000000 -R 'select[type==LINUX64 && mem>=8000 && tmp>=40000] rusage[mem=8000,tmp=20000]'";

class Genome::Model::Tools::SmrtAnalysis::MakeMappingStatsReport {
    is  => 'Genome::Model::Tools::SmrtAnalysis::Base',
    has_input => [
        cmp_hdf5_file => {
            is => 'Text',
            doc => 'The aligned reads cmp.h5 file.',
        },
        filtered_regions_fofn => {
            is => 'Text',
            doc => 'The fofn of region tables after filtering.',
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
        pie_values => { 
            is => 'Text',
            doc => 'KV pairs to plot in a pie chart',
        },              
        mode => {       
            is => 'Text',
            doc => 'How to generate report.',
            valid_values => ['internal','external'],
            default_value => 'external',
        },
    ],
    has_optional_param => [
        lsf_resource => { default_value => $DEFAULT_LSF_RESOURCE },
    ],
};

sub help_brief {
    'Generates Post-Mapping statistics for aligned read length, accuracy, and optionally z-scores.     Takes as input a cmp.h5 file.'
}


sub help_detail {
    return <<EOS 
Generates Post-Mapping statistics for aligned read length, accuracy, and optionally z-scores.     Takes as input a cmp.h5 file.
EOS
}

sub execute {
    my $self = shift;
    my $cmd = $self->analysis_bin .'/makeMappingStatsReport.py';
    if ($self->dumpdata) {
        $cmd .= ' --dumpdata';
    }
    if (defined($self->dots_per_inch)) {
        $cmd .= ' --dpi='. $self->dots_per_inch;
    }
    if (defined($self->pie_values)) {
        $cmd .= ' --pievalues="'. $self->pie_values .'"';
    }

    if (defined($self->output_dir)) {
        unless (-d $self->output_dir) {
            unless (Genome::Sys->create_directory($self->output_dir)) {
                die('Failed to create output directory: '. $self->output_dir);
            }
        }
        $cmd .= ' --output='. $self->output_dir;
        unless ($self->report_xml_file) {
            $self->report_xml_file($self->output_dir .'/quality.xml');
        }
    }
    $cmd .= ' --mode='. $self->mode .' '. $self->filtered_regions_fofn .' '. $self->cmp_hdf5_file .' > '. $self->report_xml_file;
    $self->shellcmd(
        cmd => $cmd,
        input_files => [$self->filtered_regions_fofn,$self->cmp_hdf5_file],
        output_files => [$self->report_xml_file],
        skip_if_output_is_present => 0,
    );
    return 1;
}


1;
