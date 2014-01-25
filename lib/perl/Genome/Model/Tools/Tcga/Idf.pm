package Genome::Model::Tools::Tcga::Idf;

use strict;
use warnings;
use Genome;

my %PROTOCOL_PARAMS = (
    "library preparation" => ["Vendor", "Catalog Name", "Catalog Number", "Annotation URL", "Product URL", "Target File URL", 
                            "Target File Format", "Target File Format Version", "Probe File URL", "Probe File Format",
                            "Probe File Format Version", "Target Reference Accession"],
    "sequence alignment" => ["Protocol Min Base Quality", "Protocol Min Map Quality", "Protocol Min Tumor Coverage",
                            "Protocol Min Normal Coverage"],
);

my %HARD_CODED_PROTOCOLS = (
    "mutation filtering annotation and curation" => {name => "genome.wustl.edu:maf_creation:data_consolidation:01",
                                                     description => "Automatic and manual filtering and curation of variants"},
    "nucleic acid sequencing"                    => {name => "genome.wustl.edu:DNA_sequencing:Illumina:01", 
                                                     description => "Illumina sequencing by synthesis"},
    "library preparation"                        => {name => "genome.wustl.edu:DNA_extraction:IlluminaHiSeq_DNASeq:01",
                                                     description => "Illumina library prep"},

);

my %PROCESSING_PROFILE_PROTOCOL_TYPES = (
    "sequence alignment" => "alignment",
    "variant calling" => "variant_calling",
);

my @HARD_CODED_ROW_HEADERS = (
    "Investigation Title",
    "Experimental Design",
    "Experimental Design Term Source REF",
    "Experimental Factor Name",
    "Experimental Factor Type",
    "Person Last Name",
    "Person First Name",
    "Person Email",
    "Person Affiliation",
    "Person Roles",
    "Experiment Description",
    "Protocol Term Source REF",
    "Term Source Name",
    "Term Source File",
    "Term Source Version",
);

class Genome::Model::Tools::Tcga::Idf {
    has => [
        protocols => {
        },
        sdrf_file => {
            is => "Text",
        },
    ],
};

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);

    my %p;
    $self->protocols(\%p);
    for my $key (keys %HARD_CODED_PROTOCOLS) {
        $self->_add_hard_coded_protocol($key);
    }
    return $self;
}

sub add_pp_protocols {
    my $self = shift;
    my $processing_profile = shift;

    for my $type (keys %PROCESSING_PROFILE_PROTOCOL_TYPES) {
        $self->_add_protocol_with_pp($processing_profile, $type);
    }
}

sub _add_hard_coded_protocol {
    my $self = shift;
    my $type = shift;

    unless (defined $self->protocols->{$type}) {
        $self->protocols->{$type} = [$HARD_CODED_PROTOCOLS{$type}];
    }
    return 1;
}

sub _resolve_hard_coded_protocol {
    my $self = shift;
    my $type = shift;

    return $self->protocols->{$type}->[0]->{name};
}

sub resolve_maf_protocol {
    my $self = shift;
    return $self->_resolve_hard_coded_protocol("mutation filtering annotation and curation");
}

sub resolve_sequencing_protocol {
    my $self = shift;

    return $self->_resolve_hard_coded_protocol("nucleic acid sequencing");
}

sub _resolve_protocol_with_pp {
    my $self = shift;
    my $processing_profile = shift;
    my $type = shift;

    return "genome.wustl.edu:".$PROCESSING_PROFILE_PROTOCOL_TYPES{$type}.":".$processing_profile->id.":01";
}

sub _add_protocol_with_pp {
    my $self = shift;
    my $processing_profile = shift;
    my $type = shift;

    my $name = $self->_resolve_protocol_with_pp($processing_profile, $type);
    my $description = $processing_profile->name;
    my $found = 0;
    for my $protocol (@{$self->protocols->{$type}}) {
        if ($protocol->{name} eq $name) {
            $found = 1;
            last;
        }
    }
    unless ($found) {
        push @{$self->protocols->{$type}}, {name => $name, description => $description};
    }      
    return 1;
}

sub resolve_mapping_protocol {
    my $self = shift;
    my $processing_profile = shift;
    
    return $self->_resolve_protocol_with_pp($processing_profile, "sequence alignment");
}

sub resolve_library_protocol {
    my $self = shift;

    return $self->_resolve_hard_coded_protocol("library preparation");
}

sub resolve_variants_protocol {
    my $self = shift;
    my $processing_profile = shift;

    return $self->_resolve_protocol_with_pp($processing_profile, "variant calling");
}

sub print_idf {
    my $self = shift;
    my $output_file = shift;

    my $out = Genome::Sys->open_file_for_writing($output_file);
    my @protocol_names;
    my @protocol_types;
    my @protocol_descriptions;
    my @protocol_parameters;
    for my $protocol_type (sort keys %{$self->protocols}) {
        for my $protocol (@{$self->protocols->{$protocol_type}}) {
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
    $out->print(join("\t", "SDRF Files", $self->sdrf_file)."\n");
    for my $row_header (@HARD_CODED_ROW_HEADERS) {
        $out->print("$row_header\n");
    }

    $out->close;
    return 1;
}

1;

