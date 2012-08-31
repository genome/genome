package Genome::Model::SomaticCapture::Command::RunReports;

use warnings;
use strict;

use Genome;
use Carp;
use IO::File;
use File::Basename;

class Genome::Model::SomaticCapture::Command::RunReports {
    is => 'Genome::Command::Base',
    has => [
        build => {
            shell_args_position => 1,
            is => 'Genome::Model::Build::SomaticCapture',
            id_by => 'build_id',
            doc => 'the build on which to run the reports'
        },
        build_id => { is => 'Integer', is_input => 1, is_output => 1 },
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
        },
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
    "Generate reports for the somatic-capture pipeline, currently including the variant report (cancer page) and the file summary report (line counts of all major files)",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
genome model somatic-capture run-reports --variant ~/example_variant_report.html --file-summary ~/example_file_summary_report.html 103889438
EOS
}

sub help_detail {
    return <<EOS
produces two HTML report listing the variants and structural variants and the line counts of all somatic-capture output files
EOS
}

sub execute {
    my $self = shift;

    unless($self->build) {
        $self->error_message('Failed to resolve build.');
        return;
    }

    if (($self->skip_if_output_present)&&(-s $self->variant_report_output)&&(-s $self->file_summary_report_output)) {
        $self->status_message("Skipping execution: Output is already present and skip_if_output_present is set to true");
        return 1;
    }

    my $build_id = $self->build_id;

    # Run the variant report
    my $variant_generator = Genome::Model::SomaticCapture::Report::Variant->create(
        build_id => $build_id
    );
    $self->make_report($variant_generator, $self->variant_report_output);

    # Run the file summary report
    my $fs_generator = Genome::Model::SomaticCapture::Report::FileSummary->create(
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
        $self->status_message("Report html $output exists!   Moving it out of the way...");
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
        $self->status_message('Saved report: '.$report);
    } else {
        $self->error_message('Error saving '.$report.'. Error: '. $build->error_message);
        die($self->error_message);
    }

    return 1;
}

1;
