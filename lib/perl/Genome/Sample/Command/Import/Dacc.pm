package Genome::Sample::Command::Import::Dacc;

use strict;
use warnings;

use Genome;

require Carp;
use Data::Dumper 'Dumper';
require XML::LibXML;

class Genome::Sample::Command::Import::Dacc {
    is  => 'Genome::Sample::Command::Import::Base',
    has => [
        sra_sample_id => {
            is => 'Text',
            is_input => 1,
            shell_args_position => 1,
            doc => 'SRA id to download and import from the DACC.',
        },
        xml_files => {
            is => 'Text',
            is_many => 1,
            is_optional => 1,
            shell_args_position => 2,
            doc => 'XML files',
        },
    ],
};

sub help_brief {
    return 'import DACC samples';
}

sub execute {
    my $self = shift;

    if ( $self->xml_files ) {
        return $self->xml_import;
    }
    else {
        return $self->basic_import;
    }
}

sub basic_import {
    my $self = shift;

    $self->status_message('Import DACC Sample Basic...');

    my $ok = $self->_import(
        taxon => 'Human Metagenome',
        individual => {
            upn => 'dbGaP-'.$self->sra_sample_id,
            nomenclature => 'dbGaP',
            description => 'dbGaP individual: unknown, used SRS sample id',
        },
        sample => {
            name => $self->sra_sample_id,
            nomenclature => 'dbGaP',
            extraction_label => $self->sra_sample_id,
            extraction_type => 'genomic dna',
            cell_type => 'unknown',
            nomenclature => 'dbGaP',
        },
        library => 'extlibs',
    );
    return if not $ok;

    $self->status_message('Import...OK');

    return 1;
}

sub xml_import {
    my $self = shift;

    $self->status_message('Import DACC Sample XML...');

    my $sample_info = $self->_sample_info_from_xmls;
    return if not $sample_info;

    if ( not defined $sample_info->{gap_subject_id} ) {
        $self->error_message('No gap subject id for SRA id: '.$self->sra_sample_id);
        return;
    }

    if ( not defined $sample_info->{scientific_name} ) {
        $self->error_message('No scientific name for SRA id: '.$self->sra_sample_id);
        return;
    }

    my $ok = $self->_import(
        taxon => $sample_info->{scientific_name},
        individual => {
            upn => 'dbGaP-'.$sample_info->{gap_subject_id},
            nomenclature => 'dbGaP',
            gender => $sample_info->{sex} || 'unspecified',
            description => 'dbGaP individual: '.$sample_info->{gap_subject_id},
        },
        sample => { 
            name => $self->sra_sample_id,
            nomenclature => 'dbGaP',
            tissue_label => $sample_info->{sample_type}, 
            tissue_desc => $sample_info->{body_site}, 
            extraction_label => $sample_info->{sra_sample_id},
            extraction_type => 'genomic dna',
            extraction_desc => $sample_info->{description}, 
            cell_type => 'unknown',
            nomenclature => 'dbGaP',
        },
        library => 'extlibs',
    );
    return if not $ok;

    $self->status_message('Import...OK');

    return 1;
}

sub _sample_info_from_xmls {
    my $self = shift;

    $self->status_message('Sample info from XMLs...');

    # Get library from XMLs
    my @xml_files = $self->xml_files;
    if ( not @xml_files ) {
        $self->error_message('No XML files!');
        return;
    }

    my $libxml = XML::LibXML->new();
    my @sample_infos;
    for my $xml_file ( @xml_files ) {
        next if not -s $xml_file; # seen files w/ 0 size
        my $xml = eval { $libxml->parse_file($xml_file); };
        if ( not defined $xml ) {
            $self->error_message("Could not parse report XML from file ($xml_file): $@");
            return;
        }

        my %sample_info;

        # sample
        my ($sample_node) = grep { $_->nodeType == 1 } $xml->findnodes('RunViewer/SAMPLE');
        if ( not defined $sample_node ) {
            next; # ok...we'll try the next one
        }
        $sample_info{sra_sample_id} = $sample_node->getAttribute('accession');
        if ( not defined $sample_info{sra_sample_id} ) {
            $self->error_message('No sra sample id found in sample node in XML file: '.$xml_file);
            return;
        }
        for my $attr (qw/ description /) { 
            my ($node) = grep { $_->nodeType == 1 } $sample_node->findnodes(uc($attr));
            if ( not defined $node ) {
                $self->error_message("No library attribute ($attr) node found in XML file: $xml_file");
                return;
            }
            my ($value) = $node->to_literal;
            if ( not defined $value ) {
                $self->error_message("Got sample attribute ($attr) node, but there was not a value in it");
                return;
            }
            $sample_info{$attr} = $node->to_literal;
        }

        # sample name
        for my $attr (qw/ scientific_name /) {
            my ($node) = grep { $_->nodeType == 1 } $sample_node->findnodes('SAMPLE_NAME/'.uc($attr));
            if ( not defined $node ) {
                $self->error_message('No sample name node found for '.$attr.' in XML file: '.$xml_file);
                return;
            }
            $sample_info{ lc $node->nodeName } = $node->to_literal;
        }

        # sample attrs
        my @sample_attr_nodes = grep { $_->nodeType == 1 } $sample_node->findnodes('SAMPLE_ATTRIBUTES/*');
        if ( not @sample_attr_nodes ) {
            $self->error_message('No sample attribute nodes found in XML file: '.$xml_file);
            return;
        }
        for my $node ( @sample_attr_nodes ) {
            my ($tag_node) = $node->findnodes('TAG');
            my $attr = $tag_node->to_literal;
            #next if not grep { $attr eq $_ } @sample_attrs_in_attrs;
            my ($value_node) = $node->findnodes('VALUE');
            $sample_info{ $tag_node->to_literal } = $value_node->to_literal;
        }

        push @sample_infos, \%sample_info;
    }

    if ( not @sample_infos ) {
        $self->error_message("No sample info in XMLs");
        return;
    }

    $self->status_message('Validating sample info');

    my $main_sample_info = $sample_infos[0];
    for my $sample_info ( @sample_infos[1..$#sample_infos] ) {
        for my $attr ( keys %$sample_info ) {
            if ( not defined $main_sample_info->{$attr} ) {
                # some of these are have incorrect keys
                $main_sample_info->{$attr} = $sample_info->{$attr};
            }
            elsif ( $main_sample_info->{$attr} ne $sample_info->{$attr} ) {
                $self->error_message("Mismatching data from sample infos for $attr: '$main_sample_info->{$attr}' VS '$sample_info->{$attr}");
                print Dumper(@sample_infos);
                return;
            }
        }
    }

    $self->status_message('Sample info from XMLs...OK');

    return $main_sample_info;
}

1;

