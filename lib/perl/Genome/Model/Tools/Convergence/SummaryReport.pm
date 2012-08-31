package Genome::Model::Tools::Convergence::SummaryReport;

use warnings;
use strict;

use Genome;

class Genome::Model::Tools::Convergence::SummaryReport{
    is => 'Command',
    has => [
        build_id => {
            is_input => 1,
            doc => 'the build for which to generate the summary report',
        },
    ],
    has_optional => [
        _summary_report => {
            is_output => 1,
        }
    ],
};

sub help_brief {
    "generate summary report",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt convergence summary-report --build-id 101043871
EOS
}

sub help_detail {
    return <<EOS 
produces an HTML report summarizing the Convergence model
EOS
}

sub execute {
    my $self = shift;
    
    my $build_id = $self->build_id;
    $build_id or return;
    
    my $build = Genome::Model::Build->get($build_id);
    unless($build) {
        $self->error_message("Build not found for ID: " . $build_id);
        return;
    }

    my $generator = Genome::Model::Convergence::Report::Summary->create(
        build_id => $build_id
    );
    
    my $report = $generator->generate_report();
    
    $build->add_report($report);
    
    my $subdir = $build->resolve_reports_directory . '/' . $report->name_to_subdirectory($report->name);
    $self->_summary_report($subdir . '/report.html');
    
    my $xslt_file = $generator->get_xsl_file_for_html;
    
    my $transform = Genome::Report::XSLT->transform_report(report => $report, xslt_file => $xslt_file);
    my $html = $transform->{content};
    
    my $fh = Genome::Sys->open_file_for_writing($self->_summary_report);
    unless($fh){
        $self->error_message('Failed to open report file <'. $self->_summary_report . '> for writing.');
        return;
    }
    
    $fh->print( $html );
    $fh->close;
    
    return 1;
}

1;
