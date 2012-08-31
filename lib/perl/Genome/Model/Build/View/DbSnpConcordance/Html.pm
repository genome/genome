package Genome::Model::Build::View::DbSnpConcordance::Html;

use strict;
use warnings;

use File::Basename;
use Data::Dumper;
use Genome;
use JSON;
use Sort::Naturally;
use Template;

my $DBSNP_METRICS_COMMAND = "genome model reference-alignment create-metrics db-snp-concordance --build";

class Genome::Model::Build::View::DbSnpConcordance::Html {
    is => 'UR::Object::View::Default::Html',
    has_constant => [
        perspective => {
            value => 'intersect',
        },
    ],
};

sub _generate_content {
    my $self = shift;

    my $content;
    eval {
        my $build = $self->subject;
        my $results_dir = $build->data_directory . "/reports";
        my $unfiltered_file = "$results_dir/dbsnp_concordance.txt";
        my $filtered_file = "$results_dir/dbsnp_concordance.filtered.txt";
        for my $f ($unfiltered_file, $filtered_file) {
            if (! -e $f) {
                die "Input file '$f' does not exist.\nTo generate the required data ".
                    "for this view, run:\n\n$DBSNP_METRICS_COMMAND " . $build->id . "\n";
            }
        }
                                 
        my $unfiltered_results = Genome::Model::Tools::Joinx::SnvConcordanceByQuality::parse_results_file($unfiltered_file);
        my $filtered_results = Genome::Model::Tools::Joinx::SnvConcordanceByQuality::parse_results_file($filtered_file);

        $content = $self->_render_view($unfiltered_results, $filtered_results);
    };
    return $self->_format_error($@) if $@;
    return $content;
}

sub _format_error {
    my ($self, $message) = @_;
    return <<EOF
<html><head><title>dbSNP Concordance - error</title></head>
<body>
<table width=100% height="95%"><tr align="center" valign="center"><td>
<pre>$message</pre>
</td></tr></table>
</body></html>
EOF
}

sub _get_results {
    my ($self, $build) = @_;
}

sub _base_path_for_templates {
    my $module = __PACKAGE__;
    $module =~ s/::/\//g;
    $module .= '.pm';
    my $module_path = $INC{$module};
    unless ($module_path) {
        die "Module " . __PACKAGE__ . " failed to find its own path!  Checked for $module in \%INC...";
    }
    return dirname($module_path);
}

sub _support_file_path {
    my ($self, $file) = @_;
    return $self->_base_path_for_templates . "/$file";
}

sub _hashref_to_json_array {
    my ($h) = @_;

    my @arr = map { [$_, $h->{$_}] } sort {$a <=> $b} keys %$h;

    return to_json(\@arr, {ascii => 1});
}

sub _render_view {
    my ($self, $unfiltered_results, $filtered_results) = @_;

    ## get CSS resources
    my $css_file = $self->_support_file_path("Intersect.css");
    my $css_fh = IO::File->new($css_file);
    unless ($css_fh) {
        die "failed to open file $css_file!"; 
    }
    my $page_css = join('',$css_fh->getlines);
    ## get javascript resources
    my $graph_script_file = $self->_support_file_path("Intersect.js");
    my $graph_script_fh = IO::File->new($graph_script_file);
    unless ($graph_script_fh) {
        die "failed to open file $graph_script_file"; 
    }
    my $graph_script = join('',$graph_script_fh->getlines);

    my $bname = $self->subject->__display_name__;
    my @vars = (
        page_title                     => "dbSNP Concordance for build '$bname'",
       
        # unfiltered results
        total_unfiltered_snps          => $unfiltered_results->{'total_snvs'},
        dbsnp_unfiltered_positions     => commify($unfiltered_results->{'total_hits'}),
        unfiltered_concordance         => sprintf("%.02f%%", $unfiltered_results->{'total_concordance'}),
        unfiltered_hit_snv_data        => _hashref_to_json_array($unfiltered_results->{'hit'}),
        unfiltered_all_snv_data        => _hashref_to_json_array($unfiltered_results->{'all'}),
        unfiltered_concordance_data    => _hashref_to_json_array($unfiltered_results->{'concordance'}),

        # filtered results
        total_filtered_snps          => $filtered_results->{'total_snvs'},
        dbsnp_filtered_positions     => commify($filtered_results->{'total_hits'}),
        filtered_concordance         => sprintf("%.02f%%", $filtered_results->{'total_concordance'}),
        filtered_hit_snv_data        => _hashref_to_json_array($filtered_results->{'hit'}),
        filtered_all_snv_data        => _hashref_to_json_array($filtered_results->{'all'}),
        filtered_concordance_data    => _hashref_to_json_array($filtered_results->{'concordance'}),

        graph_script                   => $graph_script,

        page_css                       => $page_css
    );

    my $template = $self->_support_file_path("Intersect.html.tt2");
    my $tt = Template->new({ ABSOLUTE => 1 }) || die "$Template::ERROR\n";

    $self->status_message("processing template $template");

    my $content;
    my $rv = $tt->process($template, { @vars }, \$content) || die $tt->error(), "\n";
    if ($rv != 1) {
   	    die "Bad return value from template processing for summary report generation: $rv ";
    }
    unless ($content) {
        die "No content returned from template processing!";
    }
    return $content;
}

sub commify {
    my $text = reverse $_[0];
    $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    return scalar reverse $text;
}

1;
