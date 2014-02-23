package Genome::DruggableGene::Command::Import::Tend;

use strict;
use warnings;
use Genome;

class Genome::DruggableGene::Command::Import::Tend {
  is => 'Genome::DruggableGene::Command::Import::Base',
  has => [
        infile => {
            is => 'Path',
            is_input => 1,
            doc => 'PATH.  Tab-delim text file provided by collaborators',
        },
        verbose => {
            is => 'Boolean',
            is_input => 1,
            is_optional => 1,
            default => 0,
            doc => 'Print more output while running',
        },
        interactions_outfile => {
            is => 'Path',
            is_input => 1,
            example_values => ['/gscmnt/sata132/techd/mgriffit/DruggableGenes/TSV/TEND_WashU_INTERACTIONS.tsv'],
            doc => 'PATH.  Path to .tsv file for drug gene interactions',
        },
        version => {
            is => 'Text',
            is_input => 1,
            doc => 'VERSION.  Version (date) of release of database from TEND group',
        },
        citation_base_url => {
            default => 'http://www.uniprot.org/uniprot/', #For genes, use Uniprot ID, Drugs will have to be handled as special case.
        },
        citation_site_url => {
            default => 'http://www.ncbi.nlm.nih.gov/pubmed/21804595/',
        },
        citation_text => {
            default => 'Trends in the exploitation of novel drug targets.  Rask-Andersen M, Almen MS, Schioth HB.  Nat Rev Drug Discov. 2011 Aug 1;10(8):579-90.  PMID: 21804595',
        },
    ],
    doc => 'Parse a tab-delim file from collaborators for TEND database',
};

sub _doc_copyright_years {
    (2011,2012);
}

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
 Obi Griffith, Ph.D.
 Malachi Griffith, Ph.D.
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
genome druggable-gene import tend --infile=/gscmnt/sata132/techd/mgriffit/DruggableGenes/KnownDruggable/TrendsExploitationNovelDrugTargets/nrd3478-s1.tsv --version="01-Aug-2011" --verbose
HELP
}

sub help_detail {
    my $summary = <<HELP
Parse a tab-delimited database file from the TEND publication (one time manual curation)
Get drug, interaction and gene info for each drug-gene interaction in the database
Get gene names and uniprot IDs from Entrez gene
Add official 'EntrezGene' name and id to each gene record via uniprot-to-entrez mapping
HELP
}

sub execute {
    my $self = shift;
    binmode(STDOUT, ":utf8");
    $self->input_to_tsv();
    return 1;
}

sub input_to_tsv {
    my $self = shift;
    my $infile = $self->infile;

    #Create interactions output file
    my $interactions_outfile = $self->interactions_outfile;
    my $interactions_fh = IO::File->new($interactions_outfile, 'w');
    my $interactions_header = join("\t", 'interaction_id', 'drug_name', 'indication', 'year_of_approval', 'uniprot_id', 'uniprot_accession_number', 'entrez_id', 'ensembl_id', 'gene_symbol', 'gene_description', 'target_main_class', 'target_subclass', 'number_transmembrane_helices');
    $interactions_fh->print($interactions_header, "\n");

    #Get the data in order
    my $infile_path = $infile;

    my ($interactions) = $self->_parse_interactions_file($infile_path);

    #Load UniProt to Entrez mapping information from file - This will be used to obtain Entrez IDs from UniProt accessions provided in Drugbank records
    my %UniProtMapping=%{$self->_get_uniprot_entrez_mapping()};

    #Write data to the file
    for my $interaction_id (keys %{$interactions}){

        my $interaction_id = $interactions->{$interaction_id}{'interaction_id'};
        $interaction_id = 'N/A' unless $interaction_id;

        my $drug_name = $interactions->{$interaction_id}{'drug_name'};
        $drug_name = 'N/A' unless $drug_name;

        my $indication = $interactions->{$interaction_id}{'indication'};
        $indication = 'N/A' unless $indication;
        
        my $year_of_approval = $interactions->{$interaction_id}{'year_of_approval'};
        $year_of_approval = 'N/A' unless $year_of_approval;
        
        my $uniprot_id = $interactions->{$interaction_id}{'uniprot_id'};
        $uniprot_id = 'N/A' unless $uniprot_id;

        #Retrieve Entrez/Ensembl IDs for interaction protein (if available)
        my $entrez_id = "N/A";
        my $ensembl_id = "N/A";
        if ($UniProtMapping{$uniprot_id}){
            $entrez_id = $UniProtMapping{$uniprot_id}{entrez_id};
            $ensembl_id = $UniProtMapping{$uniprot_id}{ensembl_id};
        }

        my $uniprot_accession_number = $interactions->{$interaction_id}{'uniprot_accession_number'};
        $uniprot_accession_number = 'N/A' unless $uniprot_accession_number;
        
        my $gene_symbol = $interactions->{$interaction_id}{'gene_symbol'};
        $gene_symbol = 'N/A' unless $gene_symbol;
        
        my $gene_description = $interactions->{$interaction_id}{'gene_description'};
        $gene_description = 'N/A' unless $gene_description;
        
        my $target_main_class = $interactions->{$interaction_id}{'target_main_class'};
        $target_main_class = 'N/A' unless $target_main_class;
        
        my $target_subclass = $interactions->{$interaction_id}{'target_subclass'};
        $target_subclass = 'N/A' unless $target_subclass;
        
        my $number_transmembrane_helices = $interactions->{$interaction_id}{'number_transmembrane_helices'};
        $number_transmembrane_helices = 'N/A' unless $number_transmembrane_helices;

        $interactions_fh->print(join("\t", $interaction_id, $drug_name, $indication, $year_of_approval, $uniprot_id, $uniprot_accession_number, $entrez_id, $ensembl_id, $gene_symbol, $gene_description, $target_main_class, $target_subclass, $number_transmembrane_helices), "\n");
        
    }
    $interactions_fh->close;
    return 1;
}

sub _parse_interactions_file {
    my $self = shift;
    my $targets_path = shift;
    my $interactions = {};
    my $fh = IO::File->new($targets_path, 'r');

    my $c = 0;
    while(my $line = <$fh>){
        next unless $line;
        $c++;
        my $interaction_id = "TEND".$c;
        chomp $line;
        $line =~ s/\r//g;
        my ($drug_name, $indication, $year_of_approval, $uniprot_id, $uniprot_accession_number, $gene_symbol, $gene_description, $target_main_class, $target_subclass, $number_transmembrane_helices) = split("\t", $line);

        $interactions->{$interaction_id}{'interaction_id'} = $interaction_id;
        $interactions->{$interaction_id}{'drug_name'} = $drug_name;
        $interactions->{$interaction_id}{'indication'} = $indication;
        $interactions->{$interaction_id}{'year_of_approval'} = $year_of_approval;
        $interactions->{$interaction_id}{'uniprot_id'} = $uniprot_id;
        $interactions->{$interaction_id}{'uniprot_accession_number'} = $uniprot_accession_number;
        $interactions->{$interaction_id}{'gene_symbol'} = $gene_symbol;
        $interactions->{$interaction_id}{'gene_description'} = $gene_description;
        $interactions->{$interaction_id}{'target_main_class'} = $target_main_class;
        $interactions->{$interaction_id}{'target_subclass'} = $target_subclass;
        $interactions->{$interaction_id}{'number_transmembrane_helices'} = $number_transmembrane_helices;

    }
    $fh->close;
    return ($interactions);
}

1;

