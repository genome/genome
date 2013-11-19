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

class Genome::Model::Tools::Tcga::Idf {
    has => [
        protocols => {
        },
    ],
};

sub create {
    my $class = shift;
    my $self = $class->SUPER::create;

    my %p;
    $self->protocols(\%p);
    return $self;
}

sub resolve_maf_protocol {
    my $self = shift;

    unless (defined $self->protocols->{"mutation filtering annotation and curation"}) {
        $self->protocols->{"mutation filtering annotation and curation"} = [{name => "genome.wustl.edu:maf_creation:data_consolidation:01",
                                                                description => "Automatic and manual filtering and curation of variants"}];
    }
    return $self->protocols->{"mutation filtering annotation and curation"}->[0]->{name};
}

sub resolve_mapping_protocol {
    my $self = shift;
    my $processing_profile = shift;

    my $name = "genome.wustl.edu:alignment:".$processing_profile->id.":01";
    my $description = $processing_profile->name;
    if (defined $self->protocols->{"sequence alignment"}){
        my $found = 0;
        for my $protocol (@{$self->protocols->{"sequence alignment"}}) {
            if ($protocol->{name} eq $name) {
                $found = 1;
                last;
            }
        }
        unless ($found) {
            push @{$self->protocols->{"sequence alignment"}}, {name => $name, description => $description};
        }      
    }
    else {
        $self->protocols->{"sequence alignment"} = [{name => "genome.wustl.edu:DNA_sequencing:Illumina:01",
                                                                description => "Illumina sequencing by synthesis"}];
    }
    return $name;
}

sub resolve_library_protocol {
    my $self = shift;

    unless (defined $self->protocols->{"library preparation"}){
        $self->protocols->{"library preparation"} = [{name => "genome.wustl.edu:DNA_extraction:Illumina_DNASeq:01",
                                                                description => "Illumina library prep"}];
    }
    return $self->protocols->{"library preparation"}->[0]->{name};
}

sub resolve_variants_protocol {
    my $self = shift;
    my $processing_profile = shift;

    my $name = "genome.wustl.edu:variant_calling:".$processing_profile->id.":01";
    my $description = $processing_profile->name;
    if (defined $self->protocols->{"variant calling"}){
        my $found = 0;
        for my $protocol (@{$self->protocols->{"variant calling"}}) {
            if ($protocol->{name} eq $name) {
                $found = 1;
                last;
            }
        }
        unless ($found) {
            push @{$self->protocols->{"variant calling"}}, {name => $name, description => $description};
        }
    }
    else {
        $self->protocols->{"variant calling"} = [{name => $name, description => $description}];
    }
    return $name;
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

    $out->close;
    return 1;
}

1;

