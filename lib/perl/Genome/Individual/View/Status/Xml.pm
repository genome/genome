package Genome::Individual::View::Status::Xml;

use strict;
use warnings;
use Genome;
use Data::Dumper;
use XML::LibXML;

class Genome::Individual::View::Status::Xml {
    is => 'Genome::View::Status::Xml',
    has_constant => [
        default_aspects => {
            is => 'ARRAY',
            value => [
                'id',
                'name',
                'common_name',
                'gender',
                {
                    name => 'samples',
                    perspective => 'default',
                    toolkit => 'xml',
                    aspects => [
                        'id',
                        'name',
                        'common_name',
                        'extraction_label',
                        'patient_common_name',
                        'extraction_type',
                        'extraction_desc',
                        'cell_type',
                        'tissue_label',
                        'tissue_desc',
                        'organ_name',
                        {
                            name => 'models',
                            perspective => 'status',
                            toolkit => 'xml',
                            subject_class_name => 'Genome::Model',
                        },
                    ]
                },
                {
                    name => 'taxon',
                    perspective => 'default',
                    toolkit => 'xml',
                    aspects => [
                        'id',
                        'domain',
                        'species_name',
                        'domain',
                        'species_latin_name',
                        'strain_name',
                        'ncbi_taxon_id',
                    ]
                },
            ]
        }
    ]
};

1;

=pod

=head1 NAME

Genome::Individual::View::Status::XML - status summary for an individual in XML format

=head1 SYNOPSIS

$i = Genome::Individual->get(1234);
$v = Genome::Individual::View::Status::Xml->create(subject => $i);
$xml = $v->content;

=head1 DESCRIPTION

This view renders the summary of an individual's status in XML format.

=cut

