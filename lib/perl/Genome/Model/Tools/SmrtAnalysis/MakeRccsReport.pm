package Genome::Model::Tools::SmrtAnalysis::MakeRccsReport;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::SmrtAnalysis::MakeRccsReport {
    is  => 'Genome::Model::Tools::SmrtAnalysis::Base',
    has_input => [
        rccs_per_base_info_file => { },
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
        extra => {
            is => 'Boolean',
            doc => 'Turn on extra plots',
            default_value => 0,
        },
    ],
};

sub help_brief {
    'Create report reference circular consensus.'
}


sub help_detail {
    return <<EOS 
Create report reference circular consensus.
EOS
}

sub execute {
    my $self = shift;
    my $cmd = $self->analysis_bin .'/makeRCCSReport.py';
    if ($self->debug) {
        $cmd .= ' --debug';
    }
    if ($self->extra) {
        $cmd .= ' --extra';
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
            $self->report_xml_file($self->output_dir .'/RCCSReport.xml');
        }
    }
    $cmd .= ' '. $self->rccs_per_base_info_file .' > '. $self->report_xml_file;
    $self->shellcmd(
        cmd => $cmd,
        input_files => [$self->rccs_per_base_info_file],
        output_files => [$self->report_xml_file],
        skip_if_output_is_present => 0,
    );
    return 1;
}


1;
