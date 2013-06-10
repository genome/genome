package Genome::Model::Tools::Graph::MutationDiagram;

use strict;
use warnings;

use Genome;
use Genome::Model::Tools::Graph::MutationDiagram::MutationDiagram;

class Genome::Model::Tools::Graph::MutationDiagram {
    is => 'Command',
    has => [
        annotation => {
            type => 'String',
            doc => "Annotator output.  Requires --reference-transcripts option",
        },
        annotation_format => {
            type => 'String',
            doc => "annotation file format",
            valid_values => ["tgi", "vep"],
            default_value => "tgi",
        },
        reference_transcripts => {
            type => 'String',
            doc => 'name/version number of the reference transcripts set ("NCBI-human.combined-annotation/0") Defaults to "NCBI-human.combined-annotation/54_36p_v2"',
            example_values => ['NCBI-human.combined-annotation/54_36p_v2'],
            is_optional => 1,
        },
        annotation_build_id => {
            type => 'Text',
            doc => 'The id of the annotation build to use',
            is_optional => 1,
        },
        genes  => {
            type => 'String',
            doc => "comma separated list of (hugo) gene names (uppercase)--default is ALL",
            is_optional => 1
        },
        custom_domains   => {
            type => 'String',
            doc => "comma separated list of protein domains to add. Expects triplets of name,start,end.",
            is_optional => 1
        },
        output_directory => {
            type => 'Text',
            doc => 'The output directory to write .svg files in',
            default => '.',
        },
        file_prefix => {
            type => 'Text',
            doc => 'A prefix to prepend to all filenames',
            default => '',
        },
        file_suffix => {
            type => 'Text',
            doc => 'A suffix to append to all filenames (before the extension)',
            default => '',
        },
        vep_frequency_field => {
            type => 'Text',
            doc => 'For VEP annotation, the name of a field in the EXTRA column that specifies the frequency of mutations',
            default_value => 'COUNT',
        },
    ],
    has_optional => [
        max_display_frequency => {
            type => 'Number',
            doc => "The maximum number of single lollis for any one mutations. Those sites exceeding this number will be truncated. Specifying this option automatically adds the number of mutations to the labels.",
        },
        lolli_shape => {
            type => 'Text',
            valid_values => ["circle", "diamond", "square"],
            default_value => "circle",
            doc => 'shape of the lolli part of each lollipop',
        },
    ],
};

sub help_brief {
    "report mutations as a (svg) diagram"
}

sub help_synopsis {
    return <<"EOS"
gmr graph mutation-diagram  --annotation my.maf
EOS
}

sub help_detail {
    return <<"EOS"
Generates (gene) mutation diagrams from an annotation file.
EOS
}

sub execute {
    my $self = shift;
    my $anno_file = $self->annotation;
    if($anno_file) {
        my %params = (
            annotation => $anno_file,
            annotation_format => $self->annotation_format,
            annotation_build_id => $self->annotation_build_id,
            hugos => $self->genes,
            custom_domains => $self->custom_domains,
            reference_transcripts => $self->reference_transcripts,
            output_directory => $self->output_directory,
            basename => $self->file_prefix,
            suffix => $self->file_suffix,
            vep_frequency_field => $self->vep_frequency_field,
            max_display_freq => $self->max_display_frequency,
            lolli_shape => $self->lolli_shape,
        );
        if ($self->annotation_build_id) {
            $params{annotation_build_id} = $self->annotation_build_id;
        }
        elsif ($self->reference_transcripts) {
            $params{reference_transcripts} = $self->reference_transcripts;
        }
        else {
            $self->error_message("Must provide either annotation_build_id or reference_transcripts");
            return;
        }

        my $anno_obj = new Genome::Model::Tools::Graph::MutationDiagram::MutationDiagram(
            %params);
    }
    else {
        $self->error_message("Must provide annotation output format");
        return;
    }
    return 1;
}


1;

#$HeadURL: svn+ssh://svn/srv/svn/gscpan/perl_modules/trunk/Genome/Model/ReferenceAlignment/Report/MutationDiagram.pm $
#$Id: MutationDiagram.pm 53299 2009-11-20 22:45:10Z eclark $
