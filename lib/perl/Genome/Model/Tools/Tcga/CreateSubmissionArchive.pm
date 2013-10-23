package Genome::Model::Tools::Tcga::CreateSubmissionArchive;

use strict;
use warnings;
use Genome;

my $NULL_CHARACTER = "->";
my $SDRF_FILE_NAME = "sdrf";

my $SDRF_MATERIAL_HEADERS = ["Extract Name", "Comment [TCGA Barcode]", "Comment [is tumor]", "Material Type", 
                    "Annotation REF", "Comment [TCGA Genome Reference]"];
my $SDRF_LIBRARY_HEADERS = ["Protocol REF", 
                    "Parameter Value [Vendor]", "Parameter Value [Catalog Number]", "Parameter Value [Produce URL]",
                    "Parameter Value [Target File URL]", "Parameter Value [Target File Format]",
                    "Parameter Value [Target File Format Version]", "Parameter Value [Probe File URL]", 
                    "Parameter Value [Probe File Format]", "Parameter Value [Probe File Format Version]",
                    "Parameter Value [Target Reference Accession]"];
my $SDRF_SEQUENCING_HEADERS = ["Protocol REF"];
my $SDRF_MAPPING_HEADERS = ["Protocol REF",
                    "Comment [Derived Data File REF]", "Comment [TCGA CGHub ID]", "Comment [TCGA CGHub metadata URL]",
                    "Comment [TCGA Include for Analysis]"], 
my $SDRF_MAPPING_DERIVED_HEADERS = ["Derived Data File", "Comment [TCGA Include for Analysis]",
                    "Comment [TCGA Data Type]", "Comment [TCGA Data Level]", "Comment [TCGA Archive Name]", 
                    "Parameter Value [Protocol Min Base Quality]", "Parameter Value [Protocol Min Map Quality]",
                    "Parameter Value [Protocol Min Tumor Coverage]", "Parameter Value [Protocol Min Normal Coverage]"];
my $SDRF_VARIANT_CALLING_HEADERS = ["Protocol REF", "Derived Data File", "Comment [TCGA Spec Version]", 
                    "Comment [TCGA Include for Analysis]", "Comment [TCGA Data Type]", "Comment [TCGA Data Level]",
                    "Comment [TCGA Archive Name]"];
my $SDRF_MAF_GENERATION_HEADERS = ["Protocol REF", "Derived Data File", "Comment [TCGA Spec Version]",
                    "Comment [TCGA Include for Analysis]", "Comment [TCGA Data Type]", "Comment [TCGA Data Level]",
                    "Comment [TCGA Archive Name]"];
my $SDRF_VALIDATION_HEADERS = ["Protocol REF", "Derived Data File", "Comment [TCGA Spec Version]",
                    "Comment [TCGA Include for Analysis]", "Comment [TCGA Data Type]", "Comment [TCGA Data Level]",
                    "Comment [TCGA Archive Name]"];

my @HEADERS = (
    ["Material", $SDRF_MATERIAL_HEADERS],
    ["Library", $SDRF_LIBRARY_HEADERS],
    ["Sequencing", $SDRF_SEQUENCING_HEADERS],
    ["Mapping", $SDRF_MAPPING_HEADERS],
    ["Mapping2", $SDRF_MAPPING_DERIVED_HEADERS],
    ["Variants", $SDRF_VARIANT_CALLING_HEADERS],
    ["Maf", $SDRF_MAF_GENERATION_HEADERS],
    ["Validation", $SDRF_VALIDATION_HEADERS],
);

class Genome::Model::Tools::Tcga::CreateSubmissionArchive {
    is => 'Command',
    has => [
        models => {
            is => "Genome::Model::ReferenceAlignment",
            is_many => 1,
        },
        output_dir => {
            is => "Text",
        },
        archive_name => {
            is => "Text",
        },
        create_archive => {
            is => "Boolean",
            default_value => 0,
        },
    ],
    has_optional => [
        maf_file => {
            is => "Text",
        },
    ],
};

sub execute {
    my $self = shift;
    my @sdrf_rows;
    my @manifest_rows;

    for my $model ($self->models) {
        my $build = $model->last_succeeded_build;
        unless($build) {
            $self->error_message("Could not resolve build from model ".$model->__display_name__);
            return;
        }
        push @sdrf_rows, $self->create_vcf_row($build);
        push @manifest_rows, $self->create_manifest_row($build);
    }

    if ($self->maf_file) {
        push @sdrf_rows, $self->create_maf_row;
    }

    $self->print_idf($self->gather_protocols(@sdrf_rows));
    $self->print_sdrf(@sdrf_rows);
    $self->print_manifest(@manifest_rows);

    return 1;
}

sub create_vcf_row {

}

sub create_maf_row {
}

sub create_manifest_row {
}

sub gather_protocols {
}

sub print_idf {
}

sub print_sdrf {
    my $self = shift;
    my @rows = @_;
    my $writer = Genome::Utility::IO::SeparatedValueWriter->create(
        output => $self->output_dir."/".$SDRF_FILE_NAME,
        separator => "\t",
        headers => $self->get_sdrf_headers,
        print_headers => 0,
    );

    for my $row (@rows) {
        $writer->write_one($self->fill_in_nulls($row));
    }

    $writer->close;

    return 1;
}

sub print_manifest {
}

sub fill_in_nulls {
    my $self = shift;
    my $row = shift;

    for my $header ($self->get_sdrf_headers) {
        unless(defined $row->{$header}) {
            $row->{$header} = $NULL_CHARACTER;
        }
    }

    return $row;
}

sub get_null_character {
    return $NULL_CHARACTER;
}

sub get_sdrf_headers {
    my $self = shift;
    my @headers;
    for my $pair (@HEADERS) {
        my ($key, $value) = @$pair;
        push @headers, map {"$key ".$_} @$value;
    }
    return @headers;
}

1;

