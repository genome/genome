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
    $row{"Library Protocol REF"} = $self->idf->resolve_library_protocol($build);
    ($row{"Library Parameter Value [Vendor]"},
    $row{"Library Parameter Value [Catalog Name]"},
    $row{"Library Parameter Value [Catalog Number]"},
    $row{"Library Parameter Value [Target File URL]"},
    $row{"Library Parameter Value [Target File Format]"},
    $row{"Library Parameter Value [Target File Format Version]"}) = $self->resolve_capture_reagent($build);
    $row{"Sequencing Protocol REF"} = $self->idf->resolve_sequencing_protocol($build);
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

sub create_maf_row {
    my $self = shift;
    my $build = shift;
    my $somatic_build = shift;
    my $maf_file = shift;
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

sub resolve_capture_reagent {
    my $self = shift;
    my $build = shift;

    if ($build->has_imported_instrument_data) {
        return (undef, undef, undef, undef, undef, undef);
    }
    unless ($build->model->target_region_set_name) {#WGS
        return (undef, undef, undef, "NA", "NA", "NA");
    }

    my %CAPTURE_REAGENTS = $self->capture_reagents;

    # Initialize to common hardcoded values
    my $total_reagent_info = {
        target_file_format => "BED",
        target_file_version => "http://genome.ucsc.edu/FAQ/FAQformat.html#format1",
    };

    # In the target_region_set_name, reagents will appear as 'reagent1 + reagent2 + reagentN'
    my @target_region_components = split /\+/, $build->model->target_region_set_name;

    for my $component (@target_region_components) {
        $component =~ s/^\s+|\s+$//g;
        my $reagent_infos = $CAPTURE_REAGENTS{$component};
        unless (defined $reagent_infos) {
            die $self->error_message("No reagent info for sub-component (%s) of capture set name (%s) from build (%s).", $component, $build->target_region_set_name, $build->id);
        }

        # Each target region component has (possibly) multiple reagents
        for my $reagent_info (@$reagent_infos) {
            for my $key ($self->mutable_capture_reagent_keys) {
                if (defined $total_reagent_info->{$key}) {
                    $total_reagent_info->{$key} = join(",", $total_reagent_info->{$key}, $reagent_info->{$key});
                } else {
                    $total_reagent_info->{$key} = $reagent_info->{$key};
                }
            }
        }
    }

    return (map { $total_reagent_info->{$_} } $self->capture_reagent_keys);
}

# The order here matters due to the way resolve_capture_reagent returns
sub capture_reagent_keys {
    my $self = shift;
    my @keys = $self->mutable_capture_reagent_keys;
    push @keys, qw(target_file_format target_file_version);
    return @keys;
}

sub mutable_capture_reagent_keys {
    return qw(reagent_vendor reagent_name catalog_number target_file);
}

# The following hash was pulled from /gsc/scripts/opt/lims/snapshots/current/lib/perl5/CGHub.pm
sub capture_reagents {
    return (
        '11111001 capture chip set' => [
            {
                reagent_vendor => 'Nimblegen',
                reagent_name   => 'Nimblegen EZ Exome v3.0',
                catalog_number => '06465692001',
                target_file    => 'http://www.nimblegen.com/downloads/annotation/ez_exome_v3/SeqCapEZ_Exome_v3.0_Design_Annotation_files.zip#SeqCap_EZ_Exome_v3_capture.bed',
            },
        ],
        '120613_HG19_EC_HPV_39235 capture oligo tube' => [
            # This is the 11111001 capture chip set Nimblegen reagent combined
            # with custom HPV IDT probes
            {
                reagent_vendor => 'Nimblegen',
                reagent_name   => 'Nimblegen EZ Exome v3.0',
                catalog_number => '06465692001',
                target_file    => 'http://www.nimblegen.com/downloads/annotation/ez_exome_v3/SeqCapEZ_Exome_v3.0_Design_Annotation_files.zip#SeqCap_EZ_Exome_v3_capture.bed',
            },
            {
                reagent_vendor => 'IDT',
                reagent_name   => '120613_HG19_EC_HPV_39235 capture oligo tube',
                catalog_number => 'NA',
                target_file    => 'ftp://genome.wustl.edu/pub/custom_capture/120613_HG19_EC_HPV_39235_capture_oligo_tube/6D44F569CD1711E1AFBE5C7646F0A7A3.bed',
            },
        ],
        '37543 capture oligo tube' => [
            {
                # Custom through Nimblegen
                reagent_vendor => 'Nimblegen',
                reagent_name   => '37543 capture oligo tube',
                catalog_number => 'NA',
                target_file    => 'ftp://genome.wustl.edu/pub/custom_capture/37543_capture_oligo_tube/7FA367D2572211E18F4237793494AFD5.bed',
            },
        ],
        '37810 capture oligo tube' => [
            {
                # Custom through Nimblegen
                reagent_vendor => 'Nimblegen',
                reagent_name   => '37810 capture oligo tube',
                catalog_number => 'NA',
                target_file    => 'ftp://genome.wustl.edu/pub/custom_capture/37810_capture_oligo_tube/9E6401595C3211E1B2416043993C62A0.bed',
            },
        ],
        'AML_KP - OID36117 capture chip set' => [
            {
                # Custom through Nimblegen
                reagent_vendor => 'Nimblegen',
                reagent_name   => 'AML_KP - OID36117 capture chip set',
                catalog_number => 'NA',
                target_file    => 'ftp://genome.wustl.edu/pub/custom_capture/AML_KP_OID36117_capture_chip_set/D64EFD6F04F011E18AB3DB401536219E.bed',
            },
        ],
        'HPV IDT all pooled probes' => [
            {
                reagent_vendor => 'IDT',
                reagent_name   => 'HPV IDT all pooled probes',
                catalog_number => 'NA',
                target_file    => 'ftp://genome.wustl.edu/pub/custom_capture/HPV_IDT_all_pooled_probes/37b1fc41fd114b64b94336ff9b4d97ae.bed',
            },
        ],
        'Nimblegen_V3_Exome_HPV_Probes_HG19 capture oligo tube' => [
            {
                reagent_vendor => 'IDT',
                reagent_name   => 'HPV_IDT_probes capture chip set',
                catalog_number => 'NA',
                target_file    => 'ftp://genome.wustl.edu/pub/custom_capture/HPV_IDT_probes_capture_chip_set/AC6217418DAB11E1BD99FC55D6BB89D5.bed',
            },
            {
                reagent_vendor => 'Nimblegen',
                reagent_name   => 'Nimblegen SeqCap EZ Human Exome Library v3.0',
                catalog_number => '06465692001',
                target_file    => 'http://www.nimblegen.com/downloads/annotation/ez_exome_v3/SeqCapEZ_Exome_v3.0_Design_Annotation_files.zip#SeqCap_EZ_Exome_v3_capture.bed',
            },
        ],
        'RT42434_combined_capture_pool' => [
            {
                reagent_vendor => 'Agilent',
                reagent_name   => 'RT42434_pool_1',
                catalog_number => 'NA',
                target_file    => 'ftp://genome.wustl.edu/pub/custom_capture/RT42434_pool_1/RT42434_pool_1.bed',
            },
            {
                reagent_vendor => 'Agilent',
                reagent_name   => 'RT42434_pool_2',
                catalog_number => 'NA',
                target_file    => 'ftp://genome.wustl.edu/pub/custom_capture/RT42434_pool_2/RT42434_pool_2.bed',
            },
            {
                reagent_vendor => 'Agilent',
                reagent_name   => 'RT42434_pool_3',
                catalog_number => 'NA',
                target_file    => 'ftp://genome.wustl.edu/pub/custom_capture/RT42434_pool_3/RT42434_pool_3.bed',
            },
            {
                reagent_vendor => 'Agilent',
                reagent_name   => 'RT42434_pool_1b',
                catalog_number => 'NA',
                target_file    => 'ftp://genome.wustl.edu/pub/custom_capture/RT42434_pool_1b/RT42434_pool_1b.bed',
            },
        ],
        'RT42434_pool_1' => [
            {
                reagent_vendor => 'Agilent',
                reagent_name   => 'RT42434_pool_1',
                catalog_number => 'NA',
                target_file    => 'ftp://genome.wustl.edu/pub/custom_capture/RT42434_pool_1/RT42434_pool_1.bed',
            },
        ],
        'RT45860 combined pool 55k (a/b) and 27k (1/2)' => [
            {
                reagent_vendor => 'Agilent',
                reagent_name   => 'RT42434_pool_1',
                catalog_number => 'NA',
                target_file    => 'ftp://genome.wustl.edu/pub/custom_capture/RT42434_pool_1/RT42434_pool_1.bed',
            },
            {
                reagent_vendor => 'Agilent',
                reagent_name   => 'RT42434_pool_2',
                catalog_number => 'NA',
                target_file    => 'ftp://genome.wustl.edu/pub/custom_capture/RT42434_pool_2/RT42434_pool_2.bed',
            },
            {
                reagent_vendor => 'Agilent',
                reagent_name   => 'RT42434_pool_3',
                catalog_number => 'NA',
                target_file    => 'ftp://genome.wustl.edu/pub/custom_capture/RT42434_pool_3/RT42434_pool_3.bed',
            },
            {
                reagent_vendor => 'Agilent',
                reagent_name   => 'RT42434_pool_1b',
                catalog_number => 'NA',
                target_file    => 'ftp://genome.wustl.edu/pub/custom_capture/RT42434_pool_1b/RT42434_pool_1b.bed',
            },
        ],
        'RT47454 TCGA OV Reorder combined pool' => [
            {
                reagent_vendor => 'Agilent',
                reagent_name   => 'TCGA OV Reorder high',
                catalog_number => 'NA',
                target_file    => 'ftp://genome.wustl.edu/pub/custom_capture/TCGA_OV_Reorder_high/TCGA_OV_Reorder_high.bed',
            },
            {
                reagent_vendor => 'Agilent',
                reagent_name   => 'TCGA OV Reorder highest',
                catalog_number => 'NA',
                target_file    => 'ftp://genome.wustl.edu/pub/custom_capture/TCGA_OV_Reorder_highest/TCGA_OV_Reorder_highest.bed',
            },
            {
                reagent_vendor => 'Agilent',
                reagent_name   => 'TCGA OV Reorder low',
                catalog_number => 'NA',
                target_file    => 'ftp://genome.wustl.edu/pub/custom_capture/TCGA_OV_Reorder_low/TCGA_OV_Reorder_low.bed',
            },
            {
                reagent_vendor => 'Agilent',
                reagent_name   => 'TCGA OV Reorder lowest',
                catalog_number => 'NA',
                target_file    => 'ftp://genome.wustl.edu/pub/custom_capture/TCGA_OV_Reorder_lowest/TCGA_OV_Reorder_lowest.bed',
            },
        ],
        'SeqCap EZ Human Exome v2.0' => [
            {
                reagent_vendor => 'Nimblegen',
                reagent_name   => 'Nimblegen SeqCap EZ Human Exome Library v2.0',
                catalog_number => '05860504001',
                target_file    => 'http://www.nimblegen.com/downloads/annotation/ez_exome_v2/SeqCapEZ_Exome_v2.0_Design_Annotation_files.zip#Design_Annotation_files/Target_Regions/SeqCap_EZ_Exome_v2.bed',
            },
        ],
        'SeqCap EZ Human Exome v3.0' => [
            {
                reagent_vendor => 'Nimblegen',
                reagent_name   => 'Nimblegen SeqCap EZ Human Exome Library v3.0',
                catalog_number => '06465692001',
                target_file    => 'http://www.nimblegen.com/downloads/annotation/ez_exome_v3/SeqCapEZ_Exome_v3.0_Design_Annotation_files.zip#SeqCap_EZ_Exome_v3_capture.bed',
            },
        ],
        'TCGA OV Massive Pool for RT48966' => [
            {
                reagent_vendor => 'Agilent',
                reagent_name   => 'RT42434_pool_1',
                catalog_number => 'NA',
                target_file    => 'ftp://genome.wustl.edu/pub/custom_capture/RT42434_pool_1/RT42434_pool_1.bed',
            },
            {
                reagent_vendor => 'Agilent',
                reagent_name   => 'RT42434_pool_2',
                catalog_number => 'NA',
                target_file    => 'ftp://genome.wustl.edu/pub/custom_capture/RT42434_pool_2/RT42434_pool_2.bed',
            },
            {
                reagent_vendor => 'Agilent',
                reagent_name   => 'RT42434_pool_3',
                catalog_number => 'NA',
                target_file    => 'ftp://genome.wustl.edu/pub/custom_capture/RT42434_pool_3/RT42434_pool_3.bed',
            },
            {
                reagent_vendor => 'Agilent',
                reagent_name   => 'RT42434_pool_1b',
                catalog_number => 'NA',
                target_file    => 'ftp://genome.wustl.edu/pub/custom_capture/RT42434_pool_1b/RT42434_pool_1b.bed',
            },
            {
                reagent_vendor => 'Agilent',
                reagent_name   => 'TCGA OV Reorder high',
                catalog_number => 'NA',
                target_file    => 'ftp://genome.wustl.edu/pub/custom_capture/TCGA_OV_Reorder_high/TCGA_OV_Reorder_high.bed',
            },
            {
                reagent_vendor => 'Agilent',
                reagent_name   => 'TCGA OV Reorder highest',
                catalog_number => 'NA',
                target_file    => 'ftp://genome.wustl.edu/pub/custom_capture/TCGA_OV_Reorder_highest/TCGA_OV_Reorder_highest.bed',
            },
            {
                reagent_vendor => 'Agilent',
                reagent_name   => 'TCGA OV Reorder low',
                catalog_number => 'NA',
                target_file    => 'ftp://genome.wustl.edu/pub/custom_capture/TCGA_OV_Reorder_low/TCGA_OV_Reorder_low.bed',
            },
            {
                reagent_vendor => 'Agilent',
                reagent_name   => 'TCGA OV Reorder lowest',
                catalog_number => 'NA',
                target_file    => 'ftp://genome.wustl.edu/pub/custom_capture/TCGA_OV_Reorder_lowest/TCGA_OV_Reorder_lowest.bed',
            },
        ],
        'TCGA OV Reorder high' => [
            {
                reagent_vendor => 'Agilent',
                reagent_name   => 'TCGA OV Reorder high',
                catalog_number => 'NA',
                target_file    => 'ftp://genome.wustl.edu/pub/custom_capture/TCGA_OV_Reorder_high/TCGA_OV_Reorder_high.bed',
            },
        ],
        'TCGA OV Reorder highest' => [
            {
                reagent_vendor => 'Agilent',
                reagent_name   => 'TCGA OV Reorder highest',
                catalog_number => 'NA',
                target_file    => 'ftp://genome.wustl.edu/pub/custom_capture/TCGA_OV_Reorder_highest/TCGA_OV_Reorder_highest.bed',
            },
        ],
        'TCGA OV Reorder low' => [
            {
                reagent_vendor => 'Agilent',
                reagent_name   => 'TCGA OV Reorder low',
                catalog_number => 'NA',
                target_file    => 'ftp://genome.wustl.edu/pub/custom_capture/TCGA_OV_Reorder_low/TCGA_OV_Reorder_low.bed',
            },
        ],
        'TCGA OV Reorder lowest' => [
            {
                reagent_vendor => 'Agilent',
                reagent_name   => 'TCGA OV Reorder lowest',
                catalog_number => 'NA',
                target_file    => 'ftp://genome.wustl.edu/pub/custom_capture/TCGA_OV_Reorder_lowest/TCGA_OV_Reorder_lowest.bed',
            },
        ],
        'WO2736953 pooled probes' => [
            {
                reagent_vendor => 'IDT',
                reagent_name   => 'WO2736953 pooled probes',
                catalog_number => 'NA',
                target_file    => 'ftp://genome.wustl.edu/pub/custom_capture/WO2736953_pooled_probes/784c240d7e6942afb8514ebdb6a950d9.bed',
            },
        ],
        'WO2768646 pooled probes' => [
            {
                reagent_vendor => 'IDT',
                reagent_name   => 'WO2768646 pooled probes',
                catalog_number => 'NA',
                target_file    => 'ftp://genome.wustl.edu/pub/custom_capture/WO2768646_pooled_probes/562e962c09834121a17bf931182c90e7.bed',
            },
        ],
        'WO2790654 pooled probes' => [
            {
                reagent_vendor => 'IDT',
                reagent_name   => 'WO2790654 pooled probes',
                catalog_number => 'NA',
                target_file    => 'ftp://genome.wustl.edu/pub/custom_capture/WO2790654_pooled_probes/1d60152280514553b6a01cd20d2b12e8.bed',
            },
        ],
        'WO2791991 pooled probes' => [
            {
                reagent_vendor => 'IDT',
                reagent_name   => 'WO2791991 pooled probes',
                catalog_number => 'NA',
                target_file    => 'ftp://genome.wustl.edu/pub/custom_capture/WO2791991_pooled_probes/ce2c7958845b4895b878a6eda8a9c521.bed',
            },
        ],
        'WO2793950 pooled probes' => [
            {
                reagent_vendor => 'IDT',
                reagent_name   => 'WO2793950 pooled probes',
                catalog_number => 'NA',
                target_file    => 'ftp://genome.wustl.edu/pub/custom_capture/WO2793950_pooled_probes/0383d24d42694f7a98b17df4f5104b4d.bed',
            },
        ],
        'WO2830729 pooled probes' => [
            {
                reagent_vendor => 'IDT',
                reagent_name   => 'WO2830729 pooled probes',
                catalog_number => 'NA',
                target_file    => 'ftp://genome.wustl.edu/pub/custom_capture/WO2830729_pooled_probes/fe74a5ca10fc4f378a733cdfd308b130.bed',
            },
        ],
        'WO2831284 pooled probes' => [
            {
                reagent_vendor => 'IDT',
                reagent_name   => 'WO2831284 pooled probes',
                catalog_number => 'NA',
                target_file    => 'ftp://genome.wustl.edu/pub/custom_capture/WO2831284_pooled_probes/b431049e36034157903ae2c75302af2d.bed',
            },
        ],
        'agilent sureselect exome version 1' => [
            {
                reagent_vendor => 'Agilent',
                reagent_name   => 'SureSelect Human All Exon 50Mb Kit',
                catalog_number => 'S02972011',
                target_file    => 'https://earray.chem.agilent.com/earray/',
            },
        ],
        'agilent sureselect exome version 2 broad' => [
            {
                reagent_vendor => 'Agilent',
                reagent_name   => 'SureSelect Human All Exon 38 Mb v2',
                catalog_number => 'S0293689',
                target_file    => 'https://earray.chem.agilent.com/earray/',
            },
        ],
        'hg18 nimblegen exome version 2' => [
            {
                reagent_vendor => 'Nimblegen',
                reagent_name   => 'hg18 nimblegen exome version 2',
                catalog_number => 'NA',
                target_file    => 'ftp://genome.wustl.edu/pub/custom_capture/hg18_nimblegen_exome_version_2/hg18_nimblegen_exome_version_2.bed',
            },
        ],
        'nimblegen exome version 1' => [
            {
                reagent_vendor => 'Nimblegen',
                reagent_name   => 'NimbleGen Sequence Capture 2.1M Human Exome Array',
                catalog_number => 'Obsolete',
                target_file    => 'http://www.nimblegen.com/downloads/annotation/seqcap_exome/2.1M_Human_Exome_Annotation.zip#2.1M_Human_Exome.bed',
            },
        ],
        'HBV_IDT_probes pooled probes' => [
            {
                reagent_vendor => 'IDT',
                reagent_name   => 'HBV_IDT_probes pooled probes',
                catalog_number => 'NA',
                target_file    => 'ftp://genome.wustl.edu/pub/custom_capture/HBV_IDT_probes_pooled_probes/40774c8461274a81b4d223161dc84936.bed',
            }
        ],
    );
}

1;
