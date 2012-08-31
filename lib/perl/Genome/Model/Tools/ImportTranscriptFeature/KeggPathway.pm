package Genome::Model::Tools::ImportTranscriptFeature::KeggPathway;

use strict;
use warnings;
use Genome;
use SOAP::Lite;

class Genome::Model::Tools::ImportTranscriptFeature::KeggPathway {
    is => 'Genome::Model::Tools::ImportTranscriptFeature',
    has =>[
        species_name => {
            is => 'Text',
            is_optional => 1,
            default => 'hsa',
            doc => 'Kegg species abbreviation used to grab pathways, defaults to human (hsa)',
        },        
        gene_pathway_file => {
            is => 'Text',
            is_optional => 1,
            default => '-',
            doc => 'Output file for gene/pathway info, default is STDOUT',
        },
        kegg_wsdl => {
            is => 'Text',
            is_optional => 1,
            default => 'http://soap.genome.jp/KEGG.wsdl',
            doc => 'URL of Kegg WSDL used to access pathway and gene information',
        },
        pathway_description_file => {
            is => 'Text',
            is_optional => 1,
            doc => 'If set, pathway name/pathway description pairs are written to this file',
        },
    ],
    doc => "Downloads human metabolic pathway information from Kegg and stores gene id, pathway, and pathway description\n",
};

sub help_detail {
    return <<EOS
This tool connects to Kegg's pathway database and pulls information about 
metabolic pathways (name and description) and the genes that are a part of those 
pathways. For each gene that is part of a pathway, a path name/gene name line is written
to the gene pathway file. If specified, each pathway's name and description is stored in 
a separate file.

The Kegg database stores this pathway information for numerous species. This tool has been 
designed to work with all of them. To use a non-human species, look up the three character 
abbreviation at http://www.genome.jp/kegg/pathway.html and set the species name option.
EOS
}

# Connects to Kegg, grabs pathway info and entrez gene ids
sub execute {
    my $self = shift;
    my $server = SOAP::Lite->service($self->kegg_wsdl);
    die "Could not establish connection to Kegg WSDL" unless defined $server;

    my $pathways = $server->list_pathways($self->species_name);
    die "Could not get pathways from Kegg" unless defined $pathways;

    my $gene_pathway_fh = IO::File->new(">" . $self->gene_pathway_file);
    unless (defined $gene_pathway_fh) {
        print "Could not get gene pathway file handle for " . $self->gene_pathway_file . "\n";
        die $!;
    }
    
    my $pathway_description_fh;
    if (defined $self->pathway_description_file) {
        $pathway_description_fh = IO::File->new(">" . $self->pathway_description_file);
        unless (defined $pathway_description_fh) {
            print "Could not get pathway description file handle for " . $self->pathway_description_file . "\n";
            die $!;
        }
    }

    for my $path (@$pathways) {
        my $path_name = $path->{entry_id};
        print "Could not get path name for pathway, skipping\n" and next unless defined $path_name;

        my $path_description = $path->{definition};
        if (defined $pathway_description_fh and not defined $path_description) {
            print "Could not get path description for " . $path_name . ", skipping\n" and next;
        }

        $pathway_description_fh->print(join("\t", $path_name, $path_description) . "\n") if defined $pathway_description_fh;

        my $genes = $server->get_genes_by_pathway($path_name); 
        die "Could not get genes from Kegg for pathway " . $path_name unless defined $genes;

        for my $gene_name (@$genes) { 
            $gene_name = substr($gene_name, length($self->species_name) + 1); # Gene name consists of species name : gene id
                                                                              # Strip away everything except gene id
            $gene_pathway_fh->print(join("\t", $path_name, $gene_name) . "\n");
        }
    }
    
    $gene_pathway_fh->close;
    $pathway_description_fh->close if defined $pathway_description_fh;
    return 1;
}

1;

