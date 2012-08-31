#:boberkfe the generate report detail method here is way too long, this should
#:boberkfe be broken up to divide out the work between collecting data and making graphs, etc

package Genome::Model::ReferenceAlignment::Report::DbSnpConcordance;

use strict;
use warnings;

use Genome;

use IO::String;
use Template;

class Genome::Model::ReferenceAlignment::Report::DbSnpConcordance {
    is => 'Genome::Model::Report',
    has => [
        name                        => { default_value => 'dbSNP Concordance' },
        description => {
            calculate => q|
            return "<div>Db Snp coverage for " . $self->model->name . " (build " . $self->build_id . ") as of " . UR::Time->now.'</div>';
            |,
        },

        report_templates => {
            is => 'String',
            is_many => 1,
            default_value => [
                'txt'
            ],
            doc => 'The paths of template(s) to use to format the report.  (In .tt2 format)',
        },
    ]
};

sub _add_to_report_xml 
{
    return { html => shift->generate_report_detail() }
}

sub generate_report_brief 
{
    my $self=shift;
    my $build = $self->build;
    return "<div>Db Snp coverage for " . $self->model->name . " (build " . $self->build_id . ") as of " . UR::Time->now.'</div>';
}

sub generate_report_detail 
{
    my $self = shift;
   
    my $build = $self->build;
    my $cmd = Genome::Model::ReferenceAlignment::Command::CreateMetrics::DbSnpConcordance->create(build_id => $build->id);
    $cmd->execute || die "Failed to execute DbSnpConcordance command";
    my $build_id = $build->id;
    my $model = $build->model;

    my $content;
    
    $self->status_message("Running dbSNP Concordance Report for build ".$build->id.".");
    
    my $module_path = $INC{"Genome/Model/ReferenceAlignment/Report/DbSnpConcordance.pm"};
    die 'failed to find module path!' unless $module_path;

    my $view_url = Genome::Config::base_web_uri() . "/genome/model/build/db-snp-concordance.html?id=" . $build->id;
    my $rv = <<EOF
<html>
<head>
<meta http-equiv="refresh" content="1;url=$view_url">
</head>
<body>
This report has been converted to a view.<br>
If you are not redirected automatically, click here: <a href="$view_url">$view_url</a>.
</body>
</html>
EOF
;
 
    return $rv;
}

1;
