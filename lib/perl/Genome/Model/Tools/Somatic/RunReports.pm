package Genome::Model::Tools::Somatic::RunReports;

use warnings;
use strict;

use Genome;
use Carp;
use IO::File;
use Genome::Sys;
use File::Basename;

class Genome::Model::Tools::Somatic::RunReports{
    is => 'Command',
    has => [
        model_id => {
            is_optional => 1,
            doc => 'the model for which to generate the variant report, must use this or build_id'
        },
        build_id => {
            is_optional => 1,
            is_input => 1,
            doc => 'the build for which to generate the variant report, must use this or model_id',
        },
        variant_report_output => {
            is => 'Text',
            is_input => 1,
            is_output => 1,
            doc => 'store the variant report HTML in the specified file'
        },
        file_summary_report_output => {
            is => 'Text',
            is_input => 1,
            is_output => 1,
            doc => 'store the file summary report HTML in the specified file'
        }
    ],
    has_optional => [
        skip_if_output_present => {
            is => 'Boolean',
            is_optional => 1,
            is_input => 1,
            default => 0,
            doc => 'enable this flag to shortcut if the report_outputs are already present. Useful for pipelines.',
        },
        
    ],
};

sub help_brief {
    "Generate reports for the somatic pipeline, currently including the variant report (cancer page) and the file summary report (line counts of all major files)",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt somatic run-reports -b 99433342 --variant ~/example_variant_report.html --file-summary ~/example_file_summary_report.html
EOS
}

sub help_detail {
    return <<EOS 
produces two HTML report listing the variants and structural variants and the line counts of all somatic output files
EOS
}

sub execute {
    my $self = shift;

    if (($self->skip_if_output_present)&&(-s $self->variant_report_output)&&(-s $self->file_summary_report_output)) {
        $self->debug_message("Skipping execution: Output is already present and skip_if_output_present is set to true");
        return 1;
    }

    my $build_id = $self->_resolve_build_id;
    $build_id or return;
    
    # Run the variant report
    my $variant_generator = Genome::Model::Somatic::Report::Variant->create(
        build_id => $build_id
    );
    $self->make_report($variant_generator, $self->variant_report_output);
    
    # Run the file summary report
    my $fs_generator = Genome::Model::Somatic::Report::FileSummary->create(
        build_id => $build_id
    );
    $self->make_report($fs_generator, $self->file_summary_report_output);

    return 1;
}

sub make_report {
    my ($self, $generator, $output) = @_;

    my $report = $generator->generate_report();

    my $xslt_file = $generator->get_xsl_file_for_html;

    my $transform = Genome::Report::XSLT->transform_report(report => $report, xslt_file => $xslt_file);
    my $html = $transform->{content};
    
    if(Genome::Sys->check_for_path_existence($output)) {
        $self->debug_message("Report html $output exists!   Moving it out of the way...");
        my $n = 1;
        my $max = 20;
        while ($n < $max and Genome::Sys->check_for_path_existence($output . '.' . $n)) {
            $n++;
        }
        if ($n == $max) {
            $self->error_message("Too many re-runs of this report!  Clean up old ones first.");
            die $self->error_message;
        }
        rename $output, "$output.$n";
        if (Genome::Sys->check_for_path_existence($output)) {
            $self->error_message("failed to move old report dir $output to $output.$n!: $!");
            die $self->error_message;
        }
    }

    my $fh = Genome::Sys->open_file_for_writing($output)
        or confess;
    $fh->print( $html );
    $fh->close;

    my $build = Genome::Model::Build->get($generator->build_id);
    if ($build->add_report($report)) {
        $self->debug_message('Saved report: '.$report);
    } else {
        $self->error_message('Error saving '.$report.'. Error: '. $build->error_message);
        die($self->error_message);
    }

    return 1;
}

sub _resolve_build_id {
    my $self = shift;
    
    my $build_id = $self->build_id;

    if ($build_id){ 
        my $build = Genome::Model::Build->get($build_id);
        
        unless($build) {
            $self->error_message("Build not found for id " . $self->build_id);
            return;
        }
    } else {
        my $model = Genome::Model->get($self->model_id);
        
        unless($model) { 
            $self->error_message("Model not found for id " . $self->model_id);
            return;
        }
        
        my $build = $model->last_succeeded_build;
        
        unless($build) {
            $self->error_message("No successful build for model");
            return;
        }
        
        $build_id = $build->build_id;
    }
    
    return $build_id;
}

1;
