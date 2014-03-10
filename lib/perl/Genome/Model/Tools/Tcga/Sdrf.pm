package Genome::Model::Tools::Tcga::Sdrf;

use strict;
use warnings;
use Genome;

my $NULL_CHARACTER = "->";

my $SDRF_MATERIAL_HEADERS = ["Extract Name", "Comment [TCGA Barcode]", "Comment [is tumor]", "Material Type", 
                    "Annotation REF", "Comment [TCGA Genome Reference]"];
my $SDRF_LIBRARY_HEADERS = ["Protocol REF", 
                    "Parameter Value [Vendor]", "Parameter Value [Catalog Name]",
                    "Parameter Value [Catalog Number]", "Parameter Value [Annotation URL]", 
                    "Parameter Value [Product URL]",
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

my $CGHUB_INFO;
my $CGHUB_INFO_BY_TCGA_NAME;

class Genome::Model::Tools::Tcga::Sdrf {
    has => [
        idf => {
            is => 'Genome::Model::Tools::Tcga::Idf',
        },
        cghub_id_file => {
            is => 'File',
        },
        archive_name => {
            is => 'Text',
        },
    ],
};

sub print_sdrf {
    my $self = shift;
    my $output_file = shift;
    my @rows = @_;

    my $temp = Genome::Sys->create_temp_file_path;
    my $writer = Genome::Utility::IO::SeparatedValueWriter->create(
        output => $temp,
        separator => "\t",
        headers => [$self->get_sdrf_headers],
        print_headers => 0,
    );

    for my $row (@rows) {
        $writer->write_one($self->fill_in_nulls($row));
    }

    my $in = Genome::Sys->open_file_for_reading($temp);

    my $out = Genome::Sys->open_file_for_writing($output_file);
    $out->print(join("\t", $self->get_output_headers)."\n");
    
    while (my $line = <$in>) {
        $out->print($line);
    }
    $out->close;

    return 1;
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

sub get_output_headers {
    my $self = shift;

    my @new_headers;
    for my $header ($self->get_sdrf_headers) {
        my @fields = split(/\s/, $header);
        shift @fields;
        push @new_headers, join(" ", @fields);
    }
    return @new_headers;
}

sub create_vcf_row {
    my $self = shift;
    my $build = shift;
    my $somatic_build = shift;
    my $vcf = shift;
    my $sample_info = shift;

    my $row = $self->fill_in_common_fields($build, $somatic_build, $sample_info);

    $row->{"Variants Derived Data File"} = $vcf;
    return $row;
}

sub fill_in_common_fields {
    my $self = shift;
    my $build = shift;
    my $somatic_build = shift;
    my $sample = shift;

    my %row;
    $row{"Material Extract Name"} = $sample->{"SampleUUID"}->{content};
    unless ($row{"Material Extract Name"}) {
        die $self->error_message("No UUID in vcf from build ".$build->id);
    }
    $row{"Material Comment [TCGA Barcode]"} = $sample->{"SampleTCGABarcode"}->{content};
    unless ($row{"Material Comment [TCGA Barcode]"}) {
        die $self->error_message("No extraction label set on subject of build ".$build->id);
    }
    #remaining required fields:
    my $sample_common_name = $build->subject->common_name;
    my $is_tumor;
    if ($sample_common_name eq "normal") {
        $is_tumor = "no";
    }
    elsif ($sample_common_name eq "tumor" or $sample_common_name eq "recurrent") {
        $is_tumor = "yes";
    }
    else {
        die $self->error_message("Unrecognized sample common name ".$sample_common_name.
            " for build ".$build->id);
    }
    $row{"Material Comment [is tumor]"} = $is_tumor;
    $row{"Material Material Type"} = "DNA";
    $row{"Material Comment [TCGA Genome Reference]"} = "GRCh37-lite";
    $row{"Library Protocol REF"} = $self->idf->resolve_library_protocol();
    ($row{"Library Parameter Value [Vendor]"},
    $row{"Library Parameter Value [Catalog Name]"},
    $row{"Library Parameter Value [Catalog Number]"}) = $self->resolve_capture_reagent($build);
    $row{"Sequencing Protocol REF"} = $self->idf->resolve_sequencing_protocol();
    $row{"Mapping Protocol REF"} = $self->idf->resolve_mapping_protocol($somatic_build->processing_profile);
    $row{"Mapping Comment [Derived Data File REF]"} =  $sample->{"File"}->{content};
    $row{"Mapping Comment [TCGA CGHub ID]"} = $self->resolve_cghub_id($build);
    $row{"Mapping Comment [TCGA Include for Analysis]"} = "yes";
    $row{"Variants Protocol REF"} = $self->idf->resolve_variants_protocol($somatic_build->processing_profile);
    $row{"Variants Comment [TCGA Include for Analysis]"} = "yes";
    $row{"Variants Comment [TCGA Data Type]"} = "Mutations";
    $row{"Variants Comment [TCGA Data Level]"} = "Level 2";
    $row{"Variants Comment [TCGA Archive Name]"} = $self->archive_name;
    return \%row;
}

sub resolve_cghub_id {
    my $self = shift;
    my $build = shift;

    unless (defined $CGHUB_INFO) {
        $CGHUB_INFO = $self->load_cghub_info("BAM_path");
    }

    unless (defined $CGHUB_INFO_BY_TCGA_NAME) {
        $CGHUB_INFO_BY_TCGA_NAME = $self->load_cghub_info("TCGA_Name");
    }
    my $id = $CGHUB_INFO->{$build->whole_rmdup_bam_file};
    unless (defined $id) {
        $id = $CGHUB_INFO_BY_TCGA_NAME->{$build->subject->extraction_label};
        unless (defined $id) {
            die("CGHub id could not be resolved for build ".$build->id." with bam file ".$build->whole_rmdup_bam_file);
        }
    }
    return $id;
}

sub load_cghub_info {
    my $self = shift;
    my $load_by = shift;
    my %id_hash;

    my $reader = Genome::Utility::IO::SeparatedValueReader->create(
        input => $self->cghub_id_file,
        separator => "\t",
    );

    while (my $line = $reader->next) {
        $id_hash{$line->{$load_by}} = $line->{CGHub_ID};
    }
    return \%id_hash;
}

sub resolve_capture_reagent {
    my $self = shift;
    my $build = shift;

    my %CAPTURE_REAGENTS = (
        "SeqCap EZ Human Exome v2.0" => {
            vendor => "Nimblegen",
            name => "Nimblegen SeqCap EZ Human Exome Library v2.0",
            number => "05860504001",
        },
        "11111001 capture chip set" => {
            vendor => "Nimblegen",
            name => "Nimblegen EZ Exome v3.0",
            number => "06465692001",
        },
        "SeqCap EZ Human Exome v3.0" => {
            vendor => "Nimblegen",
            name => "Nimblegen SeqCap EZ Human Exome Library v3.0",
            number => "06465692001",
        },
    );

    my $reagent_info = $CAPTURE_REAGENTS{$build->model->target_region_set_name};
    unless (defined $reagent_info) {
        die "No reagent info for capture set name: ".$build->target_region_set_name;
    }
    return ($reagent_info->{vendor}, $reagent_info->{name}, $reagent_info->{number});
}

sub create_maf_row {
    my $self = shift;
    my $build = shift;
    my $somatic_build = shift;
    my $maf_file =shift;
    my $sample_info = shift;

    my $row = $self->fill_in_common_fields($build, $somatic_build, $sample_info);

    $row->{"Maf Protocol REF"} = $self->idf->resolve_maf_protocol;
    #Required if providing maf file:
    $row->{"Maf Derived Data File"} = $maf_file;
    $row->{"Maf Comment [TCGA Spec Version]"} = 2.3;
    $row->{"Maf Comment [TCGA Include for Analysis]"} = "yes";
    $row->{"Maf Comment [TCGA Data Type]"} = "Mutations";
    $row->{"Maf Comment [TCGA Data Level]"} = "Level 2";
    $row->{"Maf Comment [TCGA Archive Name]"} = $self->archive_name;

    return $row;
}

1;

