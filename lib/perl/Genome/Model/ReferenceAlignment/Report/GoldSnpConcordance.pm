#:boberkfe shouldn't this use the template toolkit to format the report, and not
#:boberkfe html within this module?

package Genome::Model::ReferenceAlignment::Report::GoldSnpConcordance;

use strict;
use warnings;

use Genome;

use IO::String;
use Template;

class Genome::Model::ReferenceAlignment::Report::GoldSnpConcordance {
    is => 'Genome::Model::Report',
    has => [
        # the name is essentially constant
        name                        => { default_value => 'Gold_SNP_Concordance' },
        description => {
            calculate => q|
            return "<div>Gold Snp coverage for " . $self->model_name . " (build " . $self->build_id . ") as of " . UR::Time->now.'</div>';
            |,
        },
        report_templates => {
            is => 'String',
            is_many => 0,
            default_value => "",
            doc => 'The paths of template(s) to use to format the report.  (In .tt2 format)',
        },
    ]
};

sub _add_to_report_xml {
    return { html => shift->generate_report_detail() }
}

sub generate_report_brief {
    my $self=shift;
    my $build = $self->build;
    return "<div>Gold Snp coverage for " . $self->model_name . " (build " . $self->build_id . ") as of " . UR::Time->now.'</div>';
}

sub generate_report_detail {
    my $self = shift;
    my $build = $self->build;
    my $content;

    my $gold_snp_build = $build->gold_snp_build;
    if (!defined $gold_snp_build) {
        return $self->error_message("Unable to locate gold snp build for " . $build->__display_name__);
    }

    my $view_url = Genome::Config::base_web_uri() . "/genome/model/build/set/intersect-snv.html?-standard_build_id=" . $gold_snp_build->id . "&id=" . $build->id;
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
