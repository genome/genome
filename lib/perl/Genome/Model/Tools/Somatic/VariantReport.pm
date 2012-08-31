package Genome::Model::Tools::Somatic::VariantReport;

use warnings;
use strict;

use Genome;
use Carp;
use IO::File;
use Genome::Sys;
use File::Basename;

class Genome::Model::Tools::Somatic::VariantReport{
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
        report_output => {
            is => 'Text',
            is_input => 1,
            is_output => 1,
            doc => 'store the report HTML in the specified file'
        }
    ],
    has_optional => [
        skip_if_output_present => {
            is => 'Boolean',
            is_optional => 1,
            is_input => 1,
            default => 0,
            doc => 'enable this flag to shortcut if the report_output is already present. Useful for pipelines.',
        },
        
    ],
};

sub help_brief {
    "make variant report",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt somatic variant-report -b 99433342 --report ~/example_variant_report.html
EOS
}

sub help_detail {
    return <<EOS 
produces an HTML report listing the variants and structural variants 
EOS
}

sub execute {
    my $self = shift;

    if (($self->skip_if_output_present)&&(-s $self->report_output)) {
        $self->status_message("Skipping execution: Output is already present and skip_if_output_present is set to true");
        return 1;
    }

    my $build_id = $self->_resolve_build_id;
    $build_id or return;

    my $generator = Genome::Model::Somatic::Report::Variant->create(
        build_id => $build_id
    );
    
    my $report = $generator->generate_report();
    
    my $xslt_file = $generator->get_xsl_file_for_html;
    
    my $transform = Genome::Report::XSLT->transform_report(report => $report, xslt_file => $xslt_file);
    my $html = $transform->{content};
    
    my $fh = Genome::Sys->open_file_for_writing($self->report_output)
        or confess;
    $fh->print( $html );
    $fh->close;
    
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
