package Genome::Model::ImportedVariationList::Command::ImportVariants;

use warnings;
use strict;

use Genome;
use POSIX; # ceil()

class Genome::Model::ImportedVariationList::Command::ImportVariants {
    is => 'Command::V2',
    has_input => [
        input_path => {
            is => 'Text',
            doc => 'Path to the variants file to import. Should be sorted by joinx before import.',
        },
        version => {
            is => 'Text',
            doc => 'The version of the build to create',
        },
        reference_sequence_build => {
            is => 'Genome::Model::Build::ReferenceSequence',
            doc => 'The reference upon which the DBSnp build will be based'
        },
        source_name => {
            is => 'Text',
            doc => 'The name of the source of the imported variants (e.g., dbsnp, 1kg)',
        },
        variant_type => {
            is => 'Text',
            valid_values => ['snv', 'cnv', 'indel', 'sv'],
            doc => 'The type of variants in the imported file',
        },
        format => {
            is => 'Text',
            valid_values => ['vcf', 'bed'],
            doc => 'The format of the file being imported.',
        },
        description => {
            is => 'Text',
            doc => 'A description of the imported file',
        },
        model_name => {
            is => 'Text',
            doc => 'Name of the model to put the build under (overrides default model name of prefix-reference sequence model)',
            is_optional => 1,
        },
        bed_file => {
            is => 'Path',
            doc => 'Path to bed-format version of the dbsnp build',
            is_optional => 1,
        },
    ],
    has_transient_optional_output => [
        build => {
            is => 'Genome::Model::Build::ImportedVariationList',
            doc => 'Build created by this command'
        },
    ],
};

sub execute {
    my $self = shift;

    my $size_kb = int(ceil((-s $self->input_path) / 1024.0));
    if (!$size_kb) {
        $self->error_message("Input file " . $self->input_path . " does not exist or is empty!");
        return;
    }

    my $manual_result = Genome::Model::Tools::DetectVariants2::Result::Manual->create(
        reference_build_id => $self->reference_sequence_build->id,
        variant_type => $self->variant_type,
        format => $self->format,
        original_file_path => $self->input_path,
        description => $self->description,
    );

    if ($self->bed_file) {
        my $new_path = $manual_result->output_dir."/".$self->variant_type."s.hq.bed";
        $self->debug_message("Copying bed file to $new_path");
        Genome::Sys->copy_file($self->bed_file, $new_path);
        my @allocations = $manual_result->disk_allocations;
        $allocations[0]->reallocate;
    }

    my %params = (
        version => $self->version,
        prefix => $self->source_name,
        source_name => $self->source_name,
    );
    my $result_property = $self->variant_type."_result";
    $params{$result_property} = $manual_result;

    if (defined $self->model_name) {
        $params{"model_name"} = $self->model_name;
    }

    my $imported_variationlist_definition = Genome::Model::Command::Define::ImportedVariationList->create(
        %params
    );

    unless($imported_variationlist_definition->execute()){
        die($self->error_message("ImportedVariationList creation commanad failed"));
    }

    $self->build($imported_variationlist_definition->build);

   
    return 1;
}

1;


