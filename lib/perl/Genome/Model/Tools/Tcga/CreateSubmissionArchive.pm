package Genome::Model::Tools::Tcga::CreateSubmissionArchive;

use strict;
use warnings;
use Genome;
use Genome::File::Vcf::Reader;

my $NULL_CHARACTER = "->";
my $IDF_FILE_EXTENSION = "idf";
my $SDRF_FILE_EXTENSION = "sdrf";
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

my $CGHUB_INFO;
my $CGHUB_INFO_BY_TCGA_NAME;

class Genome::Model::Tools::Tcga::CreateSubmissionArchive {
    is => 'Command::V2',
    has => [
        models => {
            is => 'Genome::Model::SomaticVariation',
            is_many => 1,
        },       
        output_dir => {
            is => "Text",
        },
        archive_name => {
            is => "Text",
        },
        archive_version => {
            is => "Text",
        },
        create_archive => {
            is => "Boolean",
            default_value => 0,
        },
        cghub_id_file => {
            is => "Text",
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
    my @sdrf_rows;
    my %protocol_db;

    my $vcf_archive_dir = $self->output_dir."/".$self->archive_name.".Level_2.".$self->archive_version;
    Genome::Sys->create_directory($vcf_archive_dir);
    my $magetab_archive_dir = $self->output_dir."/".$self->archive_name.".mage-tab.".$self->archive_version;
    Genome::Sys->create_directory($magetab_archive_dir);

    my %patient_ids;
    for my $somatic_model ($self->models) {
        my $somatic_build = $somatic_model->last_succeeded_build;
        unless($somatic_build) {
            $self->error_message("Could not resolve build from model ".$somatic_model->__display_name__);
            return;
        }
        my $normal_build = $somatic_build->normal_build;
        my $tumor_build = $somatic_build->tumor_build;

        my $patient_id = $self->resolve_patient_id($somatic_build);
        my $patient_id_counter = ++$patient_ids{$patient_id};
        
        my $snvs_vcf = $self->construct_vcf_name("snv", $patient_id, $patient_id_counter);
        my $indels_vcf = $self->construct_vcf_name("indel", $patient_id, $patient_id_counter);

        for my $variant_type (qw(snv indel)) {
            my $local_file = $somatic_build->data_directory."/variants/".$variant_type."s_tcga/".$variant_type."s_tcga.vcf";
            die "Tcga compliant $variant_type vcf not found for build ".$somatic_build->id unless(-s $local_file);
            Genome::Sys->copy_file($local_file, "$vcf_archive_dir/".$self->construct_vcf_name($variant_type, $patient_id, $patient_id_counter));
        }

        my $vcf_sample_info = $self->get_sample_info_from_vcf("$vcf_archive_dir/$snvs_vcf");
        unless (defined $vcf_sample_info) {
            die $self->error_message("No SAMPLE vcf header in $snvs_vcf for build ".$somatic_build->id);
        }

        for my $build(($normal_build, $tumor_build)) {
            my $sample_info = $self->get_info_for_sample($build->subject->extraction_label, $vcf_sample_info);
            
            for my $vcf($snvs_vcf, $indels_vcf) {
                push @sdrf_rows, $self->create_vcf_row($build, $self->archive_name.".".$self->archive_version, \%protocol_db, $self->cghub_id_file, $vcf, $sample_info);
            }

            for my $maf_type (qw(somatic germline)) {
                my $maf_accessor = $maf_type."_maf_file";
                if ($self->$maf_accessor) {
                    Genome::Sys->copy_file($self->$maf_accessor, $vcf_archive_dir."/$maf_type.maf");
                    push @sdrf_rows, $self->create_maf_row($build, $self->archive_name.".".$self->archive_version, "$maf_type.maf", \%protocol_db, $self->cghub_id_file, $sample_info);
                }
            }
        }
    }

    $self->print_idf($magetab_archive_dir."/".$self->archive_name.".".$self->archive_version.".".$IDF_FILE_EXTENSION.".txt", \%protocol_db);
    $self->print_sdrf($magetab_archive_dir."/".$self->archive_name.".".$self->archive_version.".".$SDRF_FILE_EXTENSION.".txt", @sdrf_rows);

    $self->print_manifest($vcf_archive_dir);
    $self->print_manifest($magetab_archive_dir);

    if ($self->create_archive) {
        $self->tar_and_md5_dir($vcf_archive_dir);
        $self->tar_and_md5_dir($magetab_archive_dir);
    }

    return 1;
}

sub construct_vcf_name {
    my $self = shift;
    my $variant_type = shift;
    my $patient_id = shift;
    my $patient_id_counter = shift;
    my $snvs_vcf = "genome.wustl.edu.$patient_id.$variant_type.".$patient_id_counter.".vcf";
}

sub print_manifest {
    my $self = shift;
    my $directory = shift;
    my $cd_cmd = "cd $directory ; md5sum * > $MANIFEST_FILE_NAME";
    Genome::Sys->shellcmd(cmd => $cd_cmd, output_files => [$directory."/$MANIFEST_FILE_NAME"]);
}

sub get_info_for_sample {
    my $self = shift;
    my $desired_sample = shift;
    my $sample_info_collection = shift;
    for my $sample (@$sample_info_collection) {
        if ($sample->{"ID"}->{content} eq $desired_sample) {
            return $sample;
        }
    }
    die "Info for sample $desired_sample was not available";
}

sub get_sample_info_from_vcf {
    my $self = shift;
    my $vcf_file = shift;
    my $vcf_reader = new Genome::File::Vcf::Reader($vcf_file);
    return $vcf_reader->header->metainfo->{"SAMPLE"};
}

sub tar_and_md5_dir {
    my $self = shift;
    my $dir = shift;

    my $tar_file = "$dir.tar.gz";
    my $tar_cmd = "tar -czf $tar_file $dir";
    Genome::Sys->shellcmd(cmd => $tar_cmd, output_files => [$tar_file]);
    my $md5 = Genome::Sys->md5sum($tar_file);
    my $md5_cmd = "echo $md5 > $tar_file.md5";
    Genome::Sys->shellcmd(cmd => $md5_cmd, output_files => ["$tar_file.md5"]);
    return 1;
}

sub resolve_patient_id {
    my $self = shift;
    my $build = shift;
    my $patient_id = $build->subject->source->upn;
    die "Could not resolve patient_id for build ".$build->id unless (defined $patient_id);
    return $patient_id;
}

sub create_maf_row {
    my $self = shift;
    my $build = shift;
    my $archive_name = shift;
    my $maf_file =shift;
    my $protocol_db = shift;
    my $cghub_id_file = shift;
    my $sample_info = shift;

    my $row = $self->fill_in_common_fields($build, $archive_name, $protocol_db, $cghub_id_file, $sample_info);

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

sub create_vcf_row {
    my $self = shift;
    my $build = shift;
    my $archive_name = shift;
    my $protocol_db = shift;
    my $cghub_id_file = shift;
    my $vcf = shift;
    my $sample_info = shift;

    my $row = $self->fill_in_common_fields($build, $archive_name, $protocol_db, $cghub_id_file, $sample_info);

    $row->{"Variants Derived Data File"} = $vcf;
    return $row;
}

sub fill_in_common_fields {
    my $self = shift;
    my $build = shift;
    my $archive_name = shift;
    my $protocol_db = shift;
    my $cghub_id_file = shift;
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
    $row{"Material Comment [TCGA Genome Reference]"} = "GRCh-37lite";
    $row{"Library Protocol REF"} = $self->resolve_library_protocol($build, $protocol_db);
    ($row{"Library Parameter Value [Vendor]"},
    $row{"Library Parameter Value [Catalog Name]"},
    $row{"Library Parameter Value [Catalog Number]"}) = $self->resolve_capture_reagent($build);
    $row{"Mapping Protocol REF"} = $self->resolve_mapping_protocol($build, $protocol_db);
    $row{"Mapping Comment [Derived Data File REF]"} =  $sample->{"File"}->{content};
    $row{"Mapping Comment [TCGA CGHub ID]"} = $self->resolve_cghub_id($build, $cghub_id_file);
    $row{"Mapping Comment [TCGA Include for Analysis]"} = "yes";
    $row{"Variants Protocol REF"} = $self->resolve_variants_protocol($build, $protocol_db);
    $row{"Variants Comment [TCGA Include for Analysis]"} = "yes";
    $row{"Variants Comment [TCGA Data Type]"} = "Mutations";
    $row{"Variants Comment [TCGA Data Level]"} = "Level 2";
    $row{"Variants Comment [TCGA Archive Name]"} = $archive_name;
    return \%row;
}

sub resolve_cghub_id {
    my $self = shift;
    my $build = shift;
    my $file = shift;

    unless (defined $CGHUB_INFO) {
        $CGHUB_INFO = $self->load_cghub_info($file, "BAM_path");
    }

    unless (defined $CGHUB_INFO_BY_TCGA_NAME) {
        $CGHUB_INFO_BY_TCGA_NAME = $self->load_cghub_info($file, "TCGA_Name");
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
    my $id_file = shift;
    my $load_by = shift;
    my %id_hash;

    my $reader = Genome::Utility::IO::SeparatedValueReader->create(
        input => $id_file,
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

