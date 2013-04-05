package Genome::DruggableGene::Command::Import::RussLampel;

use strict;
use warnings;

use Genome;
use IO::File;

class Genome::DruggableGene::Command::Import::RussLampel {
    is => 'Genome::DruggableGene::Command::Import::Base',
    has => {
        infile => {
            is => 'Path',
            is_input => 1,
            doc => 'PATH.  tsv file corresponding to publicly available xls Russ and Lampel file',
        },
        tmp_dir => {
            is => 'Path',
            default => '/tmp',
            doc => 'Directory where temp files will be created',
        },
        genes_outfile => {
            is => 'Path',
            is_input => 1,
            example_values=> ['/gscmnt/sata132/techd/mgriffit/DruggableGenes/TSV/RussLampel_WashU_TARGETS.tsv'],
            doc => 'PATH.  Path to .tsv file for genes (targets)',
        },
        version => {
            doc => 'VERSION.  Version (date) of download/import (no proper version available)',
        },
        citation_base_url => {
            default => 'http://useast.ensembl.org/Homo_sapiens/Gene/Summary?g=', #No url available for direct linking of genes/drugs
        },
        citation_site_url => {
            default => 'http://www.ncbi.nlm.nih.gov/pubmed/16376820/',
        },
        citation_text => {
            default => "The druggable genome: an update. Russ AP, Lampel S. Drug Discov Today. 2005 Dec;10(23-24):1607-10. PMID: 16376820",
        },
    },
    doc => 'Parse a tsv file from Russ & Lampel (2005)',
};

sub _doc_license {
    my $self = shift;
    my (@y) = $self->_doc_copyright_years;  
    return <<EOS
Copyright (C) $y[0] Washington University in St. Louis.

It is released under the Lesser GNU Public License (LGPL) version 3.  See the 
associated LICENSE file in this distribution.
EOS
}

sub _doc_authors {
    return <<EOS
 Malachi Griffith, Ph.D.
 Obi Griffith, Ph.D.
 Jim Weible
EOS
}

=cut
sub _doc_credits {
    return ('','None at this time.');
}
=cut

sub _doc_see_also {
    return <<EOS
B<genome>(1)
EOS
}

sub _doc_manual_body {
    my $help = shift->help_detail;
    $help =~ s/\n+$/\n/g;
    return $help;
}

sub help_synopsis {
    return <<HELP
genome druggable-gene import russ-lampel --infile=/gscmnt/sata132/techd/mgriffit/DruggableGenes/PotentiallyDruggable/Russ_and_Lampel_2005/Andreas_Russ_Druggable_Gene_List.tsv --version="26-Jul-2011" 
HELP
}

sub help_detail {
    my $summary = <<HELP
Parse a tsv file corresponding to publicly available Russ and Lampel xls file (provided by Ron Bose).
Get gene info for each record (No categories are provided. All records will be given category: RussLampel)
HELP
}

sub execute {
    my $self = shift;
    my $high = 750000;
    UR::Context->object_cache_size_highwater($high);
    $self->input_to_tsv();
    $self->import_tsv();
    return 1;
}

sub input_to_tsv {
    my $self = shift;
    my $out_fh = IO::File->new($self->genes_outfile, 'w');

    my $headers = join("\t", 'gene_stable_id','display_id','description', 'HumanReadableName');
    $out_fh->print($headers, "\n");

    #Get the data in order
    my $infile_path = $self->infile;
    my $targets = $self->_parse_targets_file($infile_path);

    #Write data to the file
    for my $target_id (keys %{$targets}){

        #Gene Stable Id (ENSG)
        my $gene_stable_id =  $targets->{$target_id}{'gene_stable_id'};
        $gene_stable_id = 'N/A' unless $gene_stable_id;

        #Display Id
        my $display_id =  $targets->{$target_id}{'display_id'};
        $display_id = 'N/A' unless $display_id;

        #Gene Description
        my $description =  $targets->{$target_id}{'description'};
        $description = 'N/A' unless $description;

        #Term Category Human Readable Name - No category names provided, all will be set to RussLampel
        my $HumanReadableName = 'RussLampel';

        $out_fh->print(join("\t", $gene_stable_id, $display_id, $description, $HumanReadableName), "\n");
    }
    $out_fh->close;
    return 1;
}

sub _parse_targets_file {
    my $self = shift;
    my $targets_path = shift;
    my $targets = {};
    my $fh = IO::File->new($targets_path, 'r');

    while(my $line = <$fh>){
        next unless $line;
        chomp $line;
        $line =~ s/\r//g;
        if($line =~ m/^ENSG\d+/){
            my ($gene_stable_id,$display_id,$description) = split("\t", $line);
        $targets->{$gene_stable_id}{'gene_stable_id'} = $gene_stable_id;
        $targets->{$gene_stable_id}{'display_id'} = $display_id;
        $targets->{$gene_stable_id}{'description'} = $description;
        }else{
            #skip this line
            next;
        }
    }
    $fh->close;
    return ($targets);
}

sub import_tsv {
    my $self = shift;
    my $genes_outfile = $self->genes_outfile;
    my $citation = $self->_create_citation('RussLampel', $self->version, $self->citation_base_url, $self->citation_site_url, $self->citation_text, 'The druggable genome: an update (Russ & Lampel, 2005)');
    my @genes = $self->import_genes($genes_outfile, $citation);
    return 1;
}

sub import_genes {
    my $self = shift;
    my $version = $self->version;
    my $genes_outfile = shift;
    my $citation = shift;
    my @genes;
    my @headers = qw/gene_stable_id display_id description HumanReadableName/;
    my $parser = Genome::Utility::IO::SeparatedValueReader->create(
        input => $genes_outfile,
        headers => \@headers,
        separator => "\t",
        is_regex => 1,
    );
    
    $parser->next; #eat the headers
    while(my $input = $parser->next){
        my $gene_name = $self->_create_gene_name_report($input->{'gene_stable_id'}, $citation, 'RussLampel Gene Stable Id', '');
        my $human_readable_name = $input->{'HumanReadableName'};
        $human_readable_name =~ s/-/ /g;
        if ($human_readable_name eq 'RussLampel'){$human_readable_name="Druggable Genome";} #Create new generic category for such lists
        my $human_readable = $self->_create_gene_category_report($gene_name, 'Human Readable Name', uc($human_readable_name), '');
        unless ($input->{'gene_stable_id'} eq 'N/A'){
            my $ensembl_id = $self->_create_gene_alternate_name_report($gene_name, $input->{'gene_stable_id'}, 'Ensembl Gene Id', '', 'upper');
        }
        unless ($input->{'display_id'} eq 'N/A'){
            my $display_id = $self->_create_gene_alternate_name_report($gene_name, $input->{'display_id'}, 'Display Id', '', 'upper');
        }
        unless ($input->{'description'} eq 'N/A'){ 
            my $description = $self->_create_gene_alternate_name_report($gene_name, $input->{'description'}, 'Description', '', 'lower');
        }
        push @genes, $gene_name;
    }
    return @genes;
}

1;
