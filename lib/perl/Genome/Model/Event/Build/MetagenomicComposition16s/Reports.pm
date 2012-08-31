package Genome::Model::Event::Build::MetagenomicComposition16s::Reports;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';

class Genome::Model::Event::Build::MetagenomicComposition16s::Reports {
    is => 'Genome::Model::Event::Build::MetagenomicComposition16s',
};

sub execute {
    my $self = shift;

    my @report_names = (qw/ summary /);
    push @report_names, 'composition' if not $self->model->is_for_qc;

    for my $report_name ( @report_names ) {
        $self->_generate_and_save_report($report_name);
    }

    return 1;
}

sub _generate_and_save_report {
    my ($self, $name) = @_;

    $self->status_message(ucfirst($name).' report...');

    my $build = $self->build;
    my $class = 'Genome::Model::MetagenomicComposition16s::Report::'.Genome::Utility::Text::string_to_camel_case($name);
    my $generator = $class->create(
        build_id => $build->id,
    );
    unless ( $generator ) {
        $self->error_message("Could not create $name report generator for ".$build->description);
        return;
    }
    my $report = $generator->generate_report;
    unless ( $report ) {
        $self->error_message("Could not generate $name report for ".$build->description);
        return;
    }

    unless ( $build->add_report($report) ) {
        $self->error_message("Could not save $name report for ".$build->description);
    }

    my @datasets = $report->get_datasets;
    unless ( @datasets ) { # not ok
        $self->error_message("No datasets in $name report for ".$build->description);
        return;
    }
    my $file_base = sprintf(
        '%s/%s',
        $build->reports_directory,
        $report->name_to_subdirectory($report->name),
    );

    for my $dataset ( @datasets ) {
        my $dataset_name = $dataset->name;
        my $file = sprintf(
            '%s/%s.%s.tsv',
            $file_base,
            $self->model->subject_name,
            $dataset_name,
        );
        unlink $file if -e $file;
        my $fh = Genome::Sys->open_file_for_writing($file)
            or return; # not ok
        my ($svs) = $dataset->to_separated_value_string(separator => "\t");
        unless ( $svs ) { # not ok
            $self->error_message("Could not get dataset ($dataset) for $name report for ".$build->description);
            return;
        }
        $fh->print($svs);
        $fh->close;
    }

    my $xsl_file = $generator->get_xsl_file_for_html;
    if ( -e $xsl_file ) {
        my $xslt = Genome::Report::XSLT->transform_report(
            report => $report,
            xslt_file => $xsl_file,
        );
        unless ( $xslt ) {
            $self->error_message("Can't transform report to html.");
            return;
        }
        # Save
        my $html_file = $report->directory.'/report.html';
        my $fh = Genome::Sys->open_file_for_writing($html_file);
        unless ( $fh ) {
        }
        $fh->print( $xslt->{content} );
    }

    $self->status_message(ucfirst($name).' report...OK');
    
    return $report;
}

1;

#$HeadURL$
#$Id$
