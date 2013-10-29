package Genome::Model::Tools::Tcga::CreateSubmissionArchive;

use strict;
use warnings;
use Genome;
use Genome::File::Vcf::Reader;

my $NULL_CHARACTER = "->";
my $IDF_FILE_NAME = "idf";
my $SDRF_FILE_NAME = "sdrf";
my $MANIFEST_FILE_NAME = "MANIFEST.txt";

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

my @IDF_FIELDS = (
    "Investigation Title", "Experimental Design", "Experimental Design Term Source REF",
    "Experimental Factor Name", "Experimental Factor Type", "Person Last Name",
    "Person First Name", "Person Middle Initials", "Person Email", "Person Address",
    "Person Affiliation", "Person Roles", "PubMed ID", "Publication Author List",
    "Publication Title", "Publication Status", "Experiment Description",
    "SDRF File",
);

my %PROTOCOL_PARAMS = (
    "library preparation" => ["Vendor", "Catalog Name", "Catalog Number", "Annotation URL", "Product URL", "Target File URL", 
                            "Target File Format", "Target File Format Version", "Probe File URL", "Probe File Format",
                            "Probe File Format Version", "Target Reference Accession"],
    "sequence alignment" => ["Protocol Min Base Quality", "Protocol Min Map Quality", "Protocol Min Tumor Coverage",
                            "Protocol Min Normal Coverage"],
);

class Genome::Model::Tools::Tcga::CreateSubmissionArchive {
    is => 'Command::V2',
    has => [
        models => {
            is => 'Genome::Model::ReferenceAlignment',
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
        somatic_maf_file => {
            is => "Text",
        },
        germline_maf_file => {
            is => "Text",
        },
    ],
};

sub execute {
    my $self = shift;
    print STDERR "Executing\n";
    my @sdrf_rows;
    my %protocol_db;

    my $archive_dir = $self->output_dir."/".$self->archive_name;
    Genome::Sys->create_directory($archive_dir);

    for my $model ($self->models) {
        my $build = $model->last_succeeded_build;
        unless($build) {
            $self->error_message("Could not resolve build from model ".$model->__display_name__);
            return;
        }
        Genome::Sys->copy_file($build->data_directory."/variants/snvs.vcf.gz", $archive_dir."/".$build->id.".snvs.vcf.gz");
        Genome::Sys->copy_file($build->data_directory."/variants/indels.vcf.gz", $archive_dir."/".$build->id.".indels.vcf.gz");
        push @sdrf_rows, $self->create_snvs_vcf_row($build, $self->archive_name, \%protocol_db);
        push @sdrf_rows, $self->create_indels_vcf_row($build, $self->archive_name, \%protocol_db);
        if ($self->somatic_maf_file) {
            Genome::Sys->copy_file($self->somatic_maf_file, $archive_dir."/somatic.maf");
            push @sdrf_rows, $self->create_maf_row($build, $self->archive_name, $self->somatic_maf_file, \%protocol_db);
        }
        if ($self->germline_maf_file) {
            Genome::Sys->copy_file($self->germline_maf_file, $archive_dir."/germline.maf");
            push @sdrf_rows, $self->create_maf_row($build, $self->archive_name, $self->germline_maf_file, \%protocol_db);
        }
    }

    $self->print_idf($archive_dir."/".$IDF_FILE_NAME, \%protocol_db);
    $self->print_sdrf($archive_dir."/".$SDRF_FILE_NAME, @sdrf_rows);

    #print manifest
    my $cd_cmd = "cd $archive_dir ; md5sum * > $MANIFEST_FILE_NAME";
    Genome::Sys->shellcmd(cmd => $cd_cmd, output_files => [$archive_dir."/$MANIFEST_FILE_NAME"]);
    return 1;
}

sub create_maf_row {
    my $self = shift;
    my $build = shift;
    my $archive_name = shift;
    my $maf_file =shift;
    my $protocol_db = shift;
    my $row = $self->fill_in_common_fields($build, $archive_name, $protocol_db);

    $row->{"Maf Protocol REF"} = $self->resolve_maf_protocol($build, $protocol_db);
    #Required if providing maf file:
    $row->{"Maf Derived Data File"} = $maf_file;
    $row->{"Maf Comment [TCGA Spec Version]"} = 2.3;
    $row->{"Maf Comment [TCGA Include for Analysis]"} = "yes";
    $row->{"Maf Comment [TCGA Data Type]"} = "Mutations";
    $row->{"Maf Comment [TCGA Data Level]"} = "Level 2";
    $row->{"Maf Comment [TCGA Archive Name]"} = $archive_name;

    return $row;
}

sub resolve_maf_protocol {
    my $self = shift;
    my $build = shift;
    my $protocol_db = shift;

    unless (defined $protocol_db->{"mutation filtering annotation and curation"}) {
        $protocol_db->{"mutation filtering annotation and curation"} = [{name => "genome.wustl.edu:maf_creation:data_consolidation:01",
                                                                description => "Automatic and manual filtering and curation of variants"}];
    }
    return $protocol_db->{"mutation filtering annotation and curation"}->[0]->{name};
}

sub resolve_mapping_protocol {
    my $self = shift;
    my $build = shift;
    my $protocol_db = shift;

    my $name = "genome.wustl.edu:alignment:".$build->processing_profile->id.":01";
    my $description = $build->processing_profile->name;
    if (defined $protocol_db->{"sequence alignment"}){
        my $found = 0;
        for my $protocol (@{$protocol_db->{"variant calling"}}) {
            if ($protocol->{name} eq $name) {
                $found = 1;
                last;
            }
        }
        unless ($found) {
            push @{$protocol_db->{"variant calling"}}, {name => $name, description => $description};
        }      
    }
    else {
        $protocol_db->{"sequence alignment"} = [{name => "genome.wustl.edu:DNA_sequencing:Illumina:01",
                                                                description => "Illumina sequencing by synthesis"}];
    }
    return $name;
}

sub resolve_library_protocol {
    my $self = shift;
    my $build = shift;
    my $protocol_db = shift;

    unless (defined $protocol_db->{"library preparation"}){
        $protocol_db->{"library preparation"} = [{name => "genome.wustl.edu:DNA_extraction:Illumina_DNASeq:01",
                                                                description => "Illumina library prep"}];
    }
    return $protocol_db->{"library preparation"}->[0]->{name};
}

sub resolve_variants_protocol {
    my $self = shift;
    my $build = shift;
    my $protocol_db = shift;

    my $name = "genome.wustl.edu:variant_calling:".$build->processing_profile->id.":01";
    my $description = $build->processing_profile->name;
    if (defined $protocol_db->{"variant calling"}){
        my $found = 0;
        for my $protocol (@{$protocol_db->{"variant calling"}}) {
            if ($protocol->{name} eq $name) {
                $found = 1;
                last;
            }
        }
        unless ($found) {
            push @{$protocol_db->{"variant calling"}}, {name => $name, description => $description};
        }
    }
    else {
        $protocol_db->{"variant calling"} = [{name => $name, description => $description}];
    }
    return $name;
}

sub create_snvs_vcf_row {
    my $self = shift;
    my $build = shift;
    my $archive_name = shift;
    my $protocol_db = shift;

    my $row = $self->fill_in_common_fields($build, $archive_name, $protocol_db);

    my $snvs_vcf = $build->data_directory."/variants/snvs.vcf.gz";
    die "Couldn't find file $snvs_vcf" unless (-s $snvs_vcf);
    
    $row->{"Variants Derived Data File"} = $snvs_vcf;
    return $row;
}

sub create_indels_vcf_row {
    my $self = shift;
    my $build = shift;
    my $archive_name = shift;
    my $protocol_db = shift;

    my $row = $self->fill_in_common_fields($build, $archive_name, $protocol_db);

    my $indels_vcf = $build->data_directory."/variants/indels.vcf.gz";
    die "Couldn't find file $indels_vcf" unless (-s $indels_vcf);
    
    $row->{"Variants Derived Data File"} = $indels_vcf;
    return $row;
}

sub fill_in_common_fields {
    my $self = shift;
    my $build = shift;
    my $archive_name = shift;
    my $protocol_db = shift;

    my %row;

    my $snvs_vcf = $build->data_directory."/variants/snvs.vcf.gz";
    die "Couldn't find file $snvs_vcf" unless (-s $snvs_vcf);
    my $vcf_reader = new Genome::File::Vcf::Reader($snvs_vcf);
    my $vcf_sample_info = $vcf_reader->header->metainfo->{"SAMPLE"};
    unless (@$vcf_sample_info == 1) {
        die $self->error_message("Not exactly one SAMPLE vcf header in $snvs_vcf");
    }
    $row{"Material Extract Name"} = $vcf_sample_info->[0]->{"SampleUUID"}->{content};
    unless ($row{"Material Extract Name"}) {
        die $self->error_message("No UUID in vcf from build ".$build->id);
    }
    $row{"Material Comment [TCGA Barcode]"} = $vcf_sample_info->[0]->{"SampleTCGABarcode"}->{content};
    unless ($row{"Material Comment [TCGA Barcode]"}) {
        die $self->error_message("No extraction label set on subject of build ".$build->id);
    }
    #remaining required fields:
    my $sample_common_name = $build->subject->common_name;
    my $is_tumor;
    if ($sample_common_name eq "normal") {
        $is_tumor = "no";
    }
    elsif ($sample_common_name eq "tumor") {
        $is_tumor = "yes";
    }
    else {
        die $self->error_message("Unrecognized sample common name ".$sample_common_name.
            " for build ".$build->id);
    }
    $row{"Material Comment [is tumor]"} = $is_tumor;
    $row{"Material Material Type"} = "DNA";
    $row{"Material Comment [TCGA Genome Reference]"} = "GRCh-37lite";
    $row{"Library Protocol REF"} = $self->resolve_library_protocol($build, $protocol_db);
    #$row{"Library Parameter Value [Vendor]"} #only for non-custom
    #$row{"Library Parameter Value [Catalog Name]"} #only for non-custom
    #$row{"Library Parameter Value [Catalog Number]"} #only for non-custom
    $row{"Mapping Protocol REF"} = $self->resolve_mapping_protocol($build, $protocol_db);
    $row{"Mapping Comment [Derived Data File REF]"} =  $vcf_sample_info->[0]->{"File"}->{content};
    #$row{"Mapping Comment [TCGA CGHub ID]"}
    $row{"Mapping Comment [TCGA Include for Analysis]"} = "yes";
    $row{"Variants Protocol REF"} = $self->resolve_variants_protocol($build, $protocol_db);
    $row{"Variants Comment [TCGA Include for Analysis]"} = "yes";
    $row{"Variants Comment [TCGA Data Type]"} = "Mutations";
    $row{"Variants Comment [TCGA Data Level]"} = "Level 2";
    $row{"Variants Comment [TCGA Archive Name]"} = $archive_name;
    return \%row;
}

sub print_idf {
    my $self = shift;
    my $output_file = shift;
    my $protocol_db = shift;

    my $out = Genome::Sys->open_file_for_writing($output_file);
    my @protocol_names;
    my @protocol_types;
    my @protocol_descriptions;
    my @protocol_parameters;
    for my $protocol_type (keys %$protocol_db) {
        for my $protocol (@{$protocol_db->{$protocol_type}}) {
            push @protocol_names, $protocol->{name};
            push @protocol_types, $protocol_type;
            push @protocol_descriptions, $protocol->{description};
            push @protocol_parameters, $PROTOCOL_PARAMS{$protocol_type};
        }
    }

    $out->print(join("\t", "Protocol Name", @protocol_names)."\n");
    $out->print(join("\t", "Protocol Type", @protocol_types)."\n");
    $out->print(join("\t", "Protocol Description", @protocol_descriptions)."\n");
    $out->print(join("\t", "Protocol Parameters", map {if (defined $_){join(";", @{$_})}else {""}} @protocol_parameters)."\n");

    $out->close;
    return 1;
}

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

1;

