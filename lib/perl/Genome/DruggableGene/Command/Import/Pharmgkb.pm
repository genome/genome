package Genome::DruggableGene::Command::Import::Pharmgkb;

use strict;
use warnings;
use Genome;
use Text::ParseWords;

class Genome::DruggableGene::Command::Import::Pharmgkb {
  is => 'Genome::DruggableGene::Command::Import::Base',
  has => [
        relationships_file => {
            is => 'Path',
            is_input => 1,
            doc => 'PATH.  Tab-delim text relationships file obtained from PharmGKB',
        },
        genes_file => {
            is => 'Path',
            is_input => 1,
            doc => 'PATH.  Tab-delim text genes file obtained from PharmGKB',
        },
        drugs_file => {
            is => 'Path',
            is_input => 1,
            doc => 'PATH.  Tab-delim text drugs file obtained from PharmGKB',
        },
        interactions_outfile => {
            is => 'Path',
            is_input => 1,
            example_values => ['/gscmnt/sata132/techd/mgriffit/DruggableGenes/TSV/PharmGKB_WashU_INTERACTIONS.tsv'],
            doc => 'PATH.  Path to .tsv file for drug gene interactions',
        },
        version => {
            doc => 'VERSION.  Version (date) of release of data files from PharmGKB',
        },
        citation_base_url => {
            default => 'http://www.pharmgkb.org',
        },
        citation_site_url => {
            default => 'http://www.pharmgkb.org/',
        },
        citation_text => {
          default => 'From pharmacogenomic knowledge acquisition to clinical applications: the PharmGKB as a clinical pharmacogenomic biomarker resource.  McDonagh EM, Whirl-Carrillo M, Garten Y, Altman RB, Klein TE. Biomark Med. 2011 Dec;5(6):795-806. PMID: 22103613',
        },
    ],
    doc => 'Parse tab-delim files from PharmGKB database',
};

sub _doc_copyright_years {
    (2012);
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
genome druggable-gene import pharmgkb --version=\"12-Jul-2012\" --relationships-file=/gscmnt/sata132/techd/mgriffit/DruggableGenes/KnownDruggable/PharmGKB/2012-07-12/relationships/relationships.tsv --drugs-file=/gscmnt/sata132/techd/mgriffit/DruggableGenes/KnownDruggable/PharmGKB/2012-07-12/drugs/drugs.tsv --genes-file=/gscmnt/sata132/techd/mgriffit/DruggableGenes/KnownDruggable/PharmGKB/2012-07-12/genes/genes.tsv 
HELP
}

sub help_detail {
    my $summary = <<HELP
Parse tab-delimited database files from PharmGKB
Get drug, interaction and gene info for each drug-gene interaction in the database
Limited to only Drug-Gene interactions in relationships.tsv file and their corresponding drugs/genes
HELP
}

sub execute {
    my $self = shift;
    binmode(STDOUT, ":utf8");
    $self->input_to_tsv();
}

sub input_to_tsv {
    my $self = shift;
    my $relationships_file = $self->relationships_file;
    my $drugs_file = $self->drugs_file;
    my $genes_file = $self->genes_file;

    #Create interactions output file
    my $interactions_outfile = $self->interactions_outfile;
    my $interactions_fh = IO::File->new($interactions_outfile, 'w');
    my $interactions_header = join("\t", 'Interaction_id', 'Entity1_id','Entity1_type','Entity2_id','Entity2_type','Evidence','Association','PK','PD','PMIDs','Drug_Name','Generic_Names','Trade_Names','Brand_Mixtures','Drug_Type','Drug_Cross_References','SMILES','External_Vocabulary','Entrez_Id','Ensembl_Id','Gene_Name','Symbol','Alternate_Names','Alternate_Symbols','Is_VIP','Has_Variant_Annotation','Gene_Cross_References');
    $interactions_fh->print($interactions_header, "\n");

    #Get the data in order
    my $relationships = $self->_parse_relationships_file($relationships_file);
    my $drugs = $self->_parse_drugs_file($drugs_file);
    my $genes = $self->_parse_genes_file($genes_file);

    #Write data to the file
    for my $relationship_id (keys %{$relationships}){
        #Interaction Id - Created by joining Entity1_id and Entity2_id
        my $interaction_id =  $relationships->{$relationship_id}{'relationship_id'};
        $interaction_id = 'N/A' unless $interaction_id;

        #Entity1 id (Drug)
        my $Entity1_id =  $relationships->{$relationship_id}{'Entity1_id'};
        $Entity1_id = 'N/A' unless $Entity1_id;

        #Entity1 type (Drug)
        my $Entity1_type =  $relationships->{$relationship_id}{'Entity1_type'};
        $Entity1_type = 'N/A' unless $Entity1_type;

        #Entity2 id (Gene)
        my $Entity2_id =  $relationships->{$relationship_id}{'Entity2_id'};
        $Entity2_id = 'N/A' unless $Entity2_id;

        #Entity2 type (Gene)
        my $Entity2_type =  $relationships->{$relationship_id}{'Entity2_type'};
        $Entity2_type = 'N/A' unless $Entity2_type;

        #Interaction Evidence
        my $Evidence =  $relationships->{$relationship_id}{'Evidence'};
        $Evidence = 'N/A' unless $Evidence;

        #Association
        my $Association =  $relationships->{$relationship_id}{'Association'};
        $Association = 'N/A' unless $Association;

        #PK
        my $PK =  $relationships->{$relationship_id}{'PK'};
        $PK = 'N/A' unless $PK;

        #PD
        my $PD =  $relationships->{$relationship_id}{'PD'};
        $PD = 'N/A' unless $PD;

        #PMIDs
        my $PMIDs =  $relationships->{$relationship_id}{'PMIDs'};
        $PMIDs = 'N/A' unless $PMIDs;

        #Drug Name
        my $Drug_Name =  $drugs->{$Entity1_id}{'Drug_Name'};
        $Drug_Name = 'N/A' unless $Drug_Name;

        #Generic Names
        my $Generic_Names =  $drugs->{$Entity1_id}{'Generic_Names'};
        $Generic_Names = 'N/A' unless $Generic_Names;

        #Trade Names
        my $Trade_Names =  $drugs->{$Entity1_id}{'Trade_Names'};
        $Trade_Names = 'N/A' unless $Trade_Names;

        #Brand Mixtures
        my $Brand_Mixtures =  $drugs->{$Entity1_id}{'Brand_Mixtures'};
        $Brand_Mixtures = 'N/A' unless $Brand_Mixtures;

        #Drug Type
        my $Drug_Type =  $drugs->{$Entity1_id}{'Drug_Type'};
        $Drug_Type = 'N/A' unless $Drug_Type;

        #Drug Cross References
        my $Drug_Cross_References =  $drugs->{$Entity1_id}{'Drug_Cross_References'};
        $Drug_Cross_References = 'N/A' unless $Drug_Cross_References;

        #SMILES
        my $SMILES =  $drugs->{$Entity1_id}{'SMILES'};
        $SMILES = 'N/A' unless $SMILES;

        #External_Vocabulary
        my $External_Vocabulary =  $drugs->{$Entity1_id}{'External_Vocabulary'};
        $External_Vocabulary = 'N/A' unless $External_Vocabulary;

        #Entrez Id
        my $Entrez_Id =  $genes->{$Entity2_id}{'Entrez_Id'};
        $Entrez_Id = 'N/A' unless $Entrez_Id;

        #Ensembl Id
        my $Ensembl_Id =  $genes->{$Entity2_id}{'Ensembl_Id'};
        $Ensembl_Id = 'N/A' unless $Ensembl_Id;

        #Gene Name
        my $Gene_Name =  $genes->{$Entity2_id}{'Gene_Name'};
        $Gene_Name = 'N/A' unless $Gene_Name;

        #Symbol
        my $Symbol =  $genes->{$Entity2_id}{'Symbol'};
        $Symbol = 'N/A' unless $Symbol;

        #Alternate Names
        my $Alternate_Names =  $genes->{$Entity2_id}{'Alternate_Names'};
        $Alternate_Names = 'N/A' unless $Alternate_Names;

        #Alternate Symbols
        my $Alternate_Symbols =  $genes->{$Entity2_id}{'Alternate_Symbols'};
        $Alternate_Symbols = 'N/A' unless $Alternate_Symbols;

        #Is_VIP
        my $Is_VIP =  $genes->{$Entity2_id}{'Is_VIP'};
        $Is_VIP = 'N/A' unless $Is_VIP;

        #Has Variant Annotation
        my $Has_Variant_Annotation =  $genes->{$Entity2_id}{'Has_Variant_Annotation'};
        $Has_Variant_Annotation = 'N/A' unless $Has_Variant_Annotation;

        #Gene Cross References
        my $Gene_Cross_References =  $genes->{$Entity2_id}{'Gene_Cross_References'};
        $Gene_Cross_References = 'N/A' unless $Gene_Cross_References;

        $interactions_fh->print(join("\t", $interaction_id, $Entity1_id, $Entity1_type, $Entity2_id, $Entity2_type, $Evidence, $Association, $PK, $PD, $PMIDs, $Drug_Name, $Generic_Names, $Trade_Names, $Brand_Mixtures, $Drug_Type, $Drug_Cross_References, $SMILES, $External_Vocabulary, $Entrez_Id, $Ensembl_Id, $Gene_Name, $Symbol, $Alternate_Names, $Alternate_Symbols, $Is_VIP, $Has_Variant_Annotation, $Gene_Cross_References), "\n");
    }
    $interactions_fh->close;
    return 1;
}

sub _parse_relationships_file {
    my $self = shift;
    my $relationships_path = shift;
    my $relationships = {};
    my $fh = IO::File->new($relationships_path, 'r');

    while(my $line = <$fh>){
        next unless $line;
        chomp $line;
        $line =~ s/\r//g;
        if($line =~ m/\w+\s+Drug\s+\w+\s+Gene\s+/){ #Note: File contains both Drug-Gene and Gene-Drug pairs but these were determined to be duplicates
            my ($Entity1_id,$Entity1_type,$Entity2_id,$Entity2_type,$Evidence,$Association,$PK,$PD,$PMIDs) = split("\t", $line);
            if ($Association eq 'associated'){ #Require that the evidence for relationship between the entities is positive (not negative or ambiguous)
                my $relationship_id=join("_",$Entity1_id,$Entity2_id); #Assumes combination of Entity 1 and 2 for unique relationship - checked manually
                $relationships->{$relationship_id}{'relationship_id'} = $relationship_id;
                $relationships->{$relationship_id}{'Entity1_id'} = $Entity1_id;
                $relationships->{$relationship_id}{'Entity1_type'} = $Entity1_type;
                $relationships->{$relationship_id}{'Entity2_id'} = $Entity2_id;
                $relationships->{$relationship_id}{'Entity2_type'} = $Entity2_type;
                $relationships->{$relationship_id}{'Evidence'} = $Evidence;
                $relationships->{$relationship_id}{'Association'} = $Association;
                $relationships->{$relationship_id}{'PK'} = $PK;
                $relationships->{$relationship_id}{'PD'} = $PD;
                $relationships->{$relationship_id}{'PMIDs'} = $PMIDs;
            }else{
                #skip this line
                next;
            }
        }else{
            #skip this line
            next;
        }
    }
    $fh->close;
    return ($relationships);
}

sub _parse_drugs_file {
    my $self = shift;
    my $drugs_path = shift;
    my $drugs = {};
    my $fh = IO::File->new($drugs_path, 'r');

    while(my $line = <$fh>){
        next unless $line;
        chomp $line;
        $line =~ s/\r//g;
        if($line =~ m/^PA\d+/){
            my ($PharmGKB_Drug_Accession_Id,$Drug_Name,$Generic_Names,$Trade_Names,$Brand_Mixtures,$Drug_Type,$Drug_Cross_References,$SMILES,$External_Vocabulary) = split("\t", $line);
            $drugs->{$PharmGKB_Drug_Accession_Id}{'PharmGKB_Drug_Accession_Id'} = $PharmGKB_Drug_Accession_Id;
            $drugs->{$PharmGKB_Drug_Accession_Id}{'Drug_Name'} = $Drug_Name;
            $drugs->{$PharmGKB_Drug_Accession_Id}{'Generic_Names'} = $Generic_Names;
            $drugs->{$PharmGKB_Drug_Accession_Id}{'Trade_Names'} = $Trade_Names;
            $drugs->{$PharmGKB_Drug_Accession_Id}{'Brand_Mixtures'} = $Brand_Mixtures;
            $drugs->{$PharmGKB_Drug_Accession_Id}{'Drug_Type'} = $Drug_Type;
            $drugs->{$PharmGKB_Drug_Accession_Id}{'Drug_Cross_References'} = $Drug_Cross_References;
            $drugs->{$PharmGKB_Drug_Accession_Id}{'SMILES'} = $SMILES;
            $drugs->{$PharmGKB_Drug_Accession_Id}{'External_Vocabulary'} = $External_Vocabulary;
        }else{
            #skip this line
            next;
        }
    }
    $fh->close;
    return ($drugs);
}

sub _parse_genes_file {
    my $self = shift;
    my $genes_path = shift;
    my $genes = {};
    my $fh = IO::File->new($genes_path, 'r');

    while(my $line = <$fh>){
        next unless $line;
        chomp $line;
        $line =~ s/\r//g;
        if($line =~ m/^PA\d+/){
            my ($PharmGKB_Gene_Accession_Id,$Entrez_Id,$Ensembl_Id,$Gene_Name,$Symbol,$Alternate_Names,$Alternate_Symbols,$Is_VIP,$Has_Variant_Annotation,$Gene_Cross_References) = split("\t", $line);
            $genes->{$PharmGKB_Gene_Accession_Id}{'PharmGKB_Gene_Accession_Id'} = $PharmGKB_Gene_Accession_Id;
            $genes->{$PharmGKB_Gene_Accession_Id}{'Entrez_Id'} = $Entrez_Id;
            $genes->{$PharmGKB_Gene_Accession_Id}{'Ensembl_Id'} = $Ensembl_Id;
            $genes->{$PharmGKB_Gene_Accession_Id}{'Gene_Name'} = $Gene_Name;
            $genes->{$PharmGKB_Gene_Accession_Id}{'Symbol'} = $Symbol;
            $genes->{$PharmGKB_Gene_Accession_Id}{'Alternate_Names'} = $Alternate_Names;
            $genes->{$PharmGKB_Gene_Accession_Id}{'Alternate_Symbols'} = $Alternate_Symbols;
            $genes->{$PharmGKB_Gene_Accession_Id}{'Is_VIP'} = $Is_VIP;
            $genes->{$PharmGKB_Gene_Accession_Id}{'Has_Variant_Annotation'} = $Has_Variant_Annotation;
            $genes->{$PharmGKB_Gene_Accession_Id}{'Gene_Cross_References'} = $Gene_Cross_References;
        }else{
            #skip this line
            next;
        }
    }
    $fh->close;
    return ($genes);
}

1;
