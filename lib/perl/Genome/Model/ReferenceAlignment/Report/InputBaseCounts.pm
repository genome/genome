#:boberkfe this ought to use the template toolkit template for formatting too
#:boberkfe not html inside here

package Genome::Model::ReferenceAlignment::Report::InputBaseCounts;

use strict;
use warnings;

use Genome;
use CGI;

use IO::String;
use Template;

my $base_template_path = __PACKAGE__->_base_path_for_templates;

class Genome::Model::ReferenceAlignment::Report::InputBaseCounts {
    is => 'Genome::Model::Report',
    has => [
        # the name is essentially constant
        name                        => { default_value => 'Input_Base_Counts' },
        description => {
            calculate => q|
            return "<div>Input Base Counts for " . $self->model_name . " (build " . $self->build_id . ") as of " . UR::Time->now.'</div>';
            |,
        },
        report_templates => {
            is => 'String',
            is_many => 0,
            default_value => "$base_template_path.html.tt2",
            doc => 'The paths of template(s) to use to format the report.  (In .tt2 format)',
        },
        test => {
            is             => 'Boolean',
            default_value  => 0,
            doc            => "Saves copies of the generated data in the pwd if they do not exist. Re-uses them on the next run(s)."
        }
    ]
};

sub _base_path_for_templates 
{
    my $module = __PACKAGE__;
    $module =~ s/::/\//g;
    $module .= '.pm';
    my $module_path = $INC{$module};
    unless ($module_path) {
        die "Module " . __PACKAGE__ . " failed to find its own path!  Checked for $module in \%INC...";
    }
    return $module_path;
}

sub _add_to_report_xml 
{
    my $self = shift;
#    return {
#        description => $self->generate_report_brief,
#        html => $self->generate_report_detail,
#    };
#  below is part of shift to templating system
    my $template = shift;

    my @templates = $self->report_templates;
    unless (@templates) {
        die "No report templates assigned! Cannot generate any content."
    }

    #my $data = { description => $self->generate_report_brief };
    my $data = {};
    
    for my $template (@templates) {
        my $content = $self->generate_report_detail($template);
        my ($format,$key);
        if ($content =~ /\<\s*HTML/i) {
            $format = 'HTML';
            $key = 'html';
        }
        if (exists $data->{$key}) {
            die "Multiple templates return content in $format format. This is not supported, sadly."
                . "  Error processing $template";
        }
        $data->{$key} = $content;
    };
    return $data;
# end of additions
}

sub generate_report_brief 
{
    my $self=shift;
#    my $build = $self->build;
    return "<div>Input Base Counts for " . $self->model_name . " (build " . $self->build_id . ") as of " . UR::Time->now.'</div>';
}

sub generate_report_detail 
{
    my $self = shift;
    my $template = shift;
    unless ($template) {
        die "Please specify which template to use for this report.";
    }

    my $build = $self->build;
    my $model = $build->model;
    
    my $module_path = $INC{"Genome/Model/ReferenceAlignment/Report/InputBaseCounts.pm"};
    die 'failed to find module path!' unless $module_path;
   
    my $r = new CGI;
    my $style = $self->get_css();

    my $content;
    my $report_content;
    my $total_kb = 0;
    my @lane_data;

    my @inputs = $self->build->model->instrument_data_inputs;
    for my $input (@inputs) {
        my $instrument_data = $input->value;
        my $lane_kb;
        if($input->filter_desc)
        {
           $lane_kb = sprintf("%d",$instrument_data->total_bases_read($input->filter_desc)/1000);
        }
        else
        {
           $lane_kb = sprintf("%d",$instrument_data->total_bases_read()/1000);
        }
        $total_kb += $lane_kb;
        push @lane_data, [$instrument_data->run_identifier, $instrument_data->subset_name, commify($lane_kb), ($instrument_data->can('read_length') ? $instrument_data->read_length : 'NA')];
    }

    my $total_gb = sprintf("%.03f", $total_kb/1000000);

    if ($model->read_trimmer_name && $model->read_trimmer_name =~ /^trimq2/) {
        my ($total_ct, $total_trim_ct) = $build->calculate_input_base_counts_after_trimq2;
        if ($total_ct and $total_trim_ct) {
            my $gb       = sprintf("%.03f", $total_ct/1000000000);
            my $trim_gb  = sprintf("%.03f", $total_trim_ct/1000000000);
            $total_gb = "$trim_gb/$gb";
        }
        else {
            $self->warning_message("Failed to get input base counts after trimq2");
        }
    }

    $report_content=<<END_CONTENT;
<h2 class="section_title">Input Base Counts Per Lane</h2>
<div class="section_content">
<table class="flowcells">
<tr><td class="summary">Total Bases</td><td class="summary">$total_gb Gb</td><td colspan="2" class="summary">&nbsp;</td></tr>
<tr><th>Flowcell</th><th>Lane</th><th>Kb</th><th>Read Length</th></tr>
END_CONTENT

    for (@lane_data) {
        my ($fc, $lane, $kb, $rl) = @$_;
        $report_content .= "<tr><td class='flowcell'>$fc</td><td>$lane</td><td>$kb</td><td>$rl</td></tr>\n";
    }

    $report_content .= "\n</table></div>";

    $build->set_metric("instrument data total kb", $total_kb);



    my @vars = (
        model_id       => $model->id,
        model_name     => $model->name,
        build_id       => $build->id,
        page_title     => "Input Base Count Report for Model " . $model->id . " (" . $model->name . "), build " .$build->id,
        style          => $style,
        report_content => $report_content
    );

    my $tt = Template->new({
        ABSOLUTE => 1,
    }) || die "$Template::ERROR\n";

    my $rv = $tt->process($template, { @vars }, \$content) || die $$tt->error(), "\n";
    if ($rv != 1) {
        die "Bad return value from template processing for summary report generation: $rv ";
    }
    unless ($report_content) {
        die "No content returned from template processing!";
    }

    my $body = IO::String->new();  
    die $! unless $body;
    $body->print($content);
    $body->seek(0, 0);
    return join('', $body->getlines);        

}

sub get_css
{
    my $module_path = $INC{"Genome/Model/ReferenceAlignment/Report/InputBaseCounts.pm"};
    die 'failed to find module path!' unless $module_path;
    
    ## get CSS resources
    my $css_file = "$module_path.html.css";
    my $css_fh = IO::File->new($css_file);
    unless ($css_fh) {
        die "failed to open file $css_file!"; 
    }
    my $page_css = join('',$css_fh->getlines);

}

sub commify {
    my $text = reverse $_[0];
    $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    return scalar reverse $text;
}


1;
