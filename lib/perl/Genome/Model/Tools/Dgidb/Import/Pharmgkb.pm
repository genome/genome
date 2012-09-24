package Genome::Model::Tools::Dgidb::Import::Pharmgkb;

use strict;
use warnings;
use Genome;
use Text::ParseWords;

binmode(STDOUT, ":utf8");

class Genome::Model::Tools::Dgidb::Import::Pharmgkb {
  is => 'Genome::Model::Tools::Dgidb::Import::Base',
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
            default => '/gscmnt/sata132/techd/mgriffit/DruggableGenes/TSV/PharmGKB_WashU_INTERACTIONS.tsv',
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
B<gmt>(1)
EOS
}

sub _doc_manual_body {
    my $help = shift->help_detail;
    $help =~ s/\n+$/\n/g;
    return $help;
}

sub help_synopsis {
    return <<HELP
gmt dgidb import pharmgkb --version=\"2012-07-12\" --relationships-file=/gscmnt/sata132/techd/mgriffit/DruggableGenes/KnownDruggable/PharmGKB/2012-07-12/relationships/relationships.tsv --drugs-file=/gscmnt/sata132/techd/mgriffit/DruggableGenes/KnownDruggable/PharmGKB/2012-07-12/drugs/drugs.tsv --genes-file=/gscmnt/sata132/techd/mgriffit/DruggableGenes/KnownDruggable/PharmGKB/2012-07-12/genes/genes.tsv 
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
    $self->input_to_tsv();
    $self->import_tsv();
    $self->_destroy_and_rebuild_pubchem_and_drug_groups();
    return 1;
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
        $interaction_id = 'NA' unless $interaction_id;

        #Entity1 id (Drug)
        my $Entity1_id =  $relationships->{$relationship_id}{'Entity1_id'};
        $Entity1_id = 'NA' unless $Entity1_id;

        #Entity1 type (Drug)
        my $Entity1_type =  $relationships->{$relationship_id}{'Entity1_type'};
        $Entity1_type = 'NA' unless $Entity1_type;

        #Entity2 id (Gene)
        my $Entity2_id =  $relationships->{$relationship_id}{'Entity2_id'};
        $Entity2_id = 'NA' unless $Entity2_id;

        #Entity2 type (Gene)
        my $Entity2_type =  $relationships->{$relationship_id}{'Entity2_type'};
        $Entity2_type = 'NA' unless $Entity2_type;

        #Interaction Evidence
        my $Evidence =  $relationships->{$relationship_id}{'Evidence'};
        $Evidence = 'NA' unless $Evidence;

        #Association
        my $Association =  $relationships->{$relationship_id}{'Association'};
        $Association = 'NA' unless $Association;

        #PK
        my $PK =  $relationships->{$relationship_id}{'PK'};
        $PK = 'NA' unless $PK;

        #PD
        my $PD =  $relationships->{$relationship_id}{'PD'};
        $PD = 'NA' unless $PD;

        #PMIDs
        my $PMIDs =  $relationships->{$relationship_id}{'PMIDs'};
        $PMIDs = 'NA' unless $PMIDs;

        #Drug Name
        my $Drug_Name =  $drugs->{$Entity1_id}{'Drug_Name'};
        $Drug_Name = 'NA' unless $Drug_Name;

        #Generic Names
        my $Generic_Names =  $drugs->{$Entity1_id}{'Generic_Names'};
        $Generic_Names = 'NA' unless $Generic_Names;

        #Trade Names
        my $Trade_Names =  $drugs->{$Entity1_id}{'Trade_Names'};
        $Trade_Names = 'NA' unless $Trade_Names;

        #Brand Mixtures
        my $Brand_Mixtures =  $drugs->{$Entity1_id}{'Brand_Mixtures'};
        $Brand_Mixtures = 'NA' unless $Brand_Mixtures;

        #Drug Type
        my $Drug_Type =  $drugs->{$Entity1_id}{'Drug_Type'};
        $Drug_Type = 'NA' unless $Drug_Type;

        #Drug Cross References
        my $Drug_Cross_References =  $drugs->{$Entity1_id}{'Drug_Cross_References'};
        $Drug_Cross_References = 'NA' unless $Drug_Cross_References;

        #SMILES
        my $SMILES =  $drugs->{$Entity1_id}{'SMILES'};
        $SMILES = 'NA' unless $SMILES;

        #External_Vocabulary
        my $External_Vocabulary =  $drugs->{$Entity1_id}{'External_Vocabulary'};
        $External_Vocabulary = 'NA' unless $External_Vocabulary;

        #Entrez Id
        my $Entrez_Id =  $genes->{$Entity2_id}{'Entrez_Id'};
        $Entrez_Id = 'NA' unless $Entrez_Id;

        #Ensembl Id
        my $Ensembl_Id =  $genes->{$Entity2_id}{'Ensembl_Id'};
        $Ensembl_Id = 'NA' unless $Ensembl_Id;

        #Gene Name
        my $Gene_Name =  $genes->{$Entity2_id}{'Gene_Name'};
        $Gene_Name = 'NA' unless $Gene_Name;

        #Symbol
        my $Symbol =  $genes->{$Entity2_id}{'Symbol'};
        $Symbol = 'NA' unless $Symbol;

        #Alternate Names
        my $Alternate_Names =  $genes->{$Entity2_id}{'Alternate_Names'};
        $Alternate_Names = 'NA' unless $Alternate_Names;

        #Alternate Symbols
        my $Alternate_Symbols =  $genes->{$Entity2_id}{'Alternate_Symbols'};
        $Alternate_Symbols = 'NA' unless $Alternate_Symbols;

        #Is_VIP
        my $Is_VIP =  $genes->{$Entity2_id}{'Is_VIP'};
        $Is_VIP = 'NA' unless $Is_VIP;

        #Has Variant Annotation
        my $Has_Variant_Annotation =  $genes->{$Entity2_id}{'Has_Variant_Annotation'};
        $Has_Variant_Annotation = 'NA' unless $Has_Variant_Annotation;

        #Gene Cross References
        my $Gene_Cross_References =  $genes->{$Entity2_id}{'Gene_Cross_References'};
        $Gene_Cross_References = 'NA' unless $Gene_Cross_References;

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

sub import_tsv {
    my $self = shift;
    my $interactions_outfile = $self->interactions_outfile;
    $self->preload_objects;
    my @interactions = $self->import_interactions($interactions_outfile);
    return 1;
}

sub preload_objects {
    my $self = shift;
    my $source_db_name = 'PharmGKB';
    my $source_db_version = $self->version;

    #Let's preload anything for this database name and version so that we can avoid death by 1000 queries
    my @gene_names = Genome::DruggableGene::GeneNameReport->get(source_db_name => $source_db_name, source_db_version => $source_db_version);
    for my $gene_name (@gene_names){
        $gene_name->gene_alt_names;
        $gene_name->gene_categories;
    }
    my @drug_names = Genome::DruggableGene::DrugNameReport->get(source_db_name => $source_db_name, source_db_version => $source_db_version);
    for my $drug_name (@drug_names){
        $drug_name->drug_alt_names;
        $drug_name->drug_categories;
    }
    my @gene_ids = map($_->id, @gene_names);
    my @interactions = Genome::DruggableGene::DrugGeneInteractionReport->get(gene_id => \@gene_ids);
    for my $interaction (@interactions){
        $interaction->interaction_attributes;
    }
    return 1;
}

sub import_interactions {
    my $self = shift;
    my $interactions_outfile = shift;
    my $version = $self->version;
    my @interactions;
    my @headers = qw(Interaction_id Entity1_id Entity1_type Entity2_id Entity2_type Evidence Association PK PD PMIDs Drug_Name Generic_Names Trade_Names Brand_Mixtures Drug_Type Drug_Cross_References SMILES External_Vocabulary Entrez_Id Ensembl_Id Gene_Name Symbol Alternate_Names Alternate_Symbols Is_VIP Has_Variant_Annotation Gene_Cross_References);
    my $parser = Genome::Utility::IO::SeparatedValueReader->create(
        input => $interactions_outfile,
        headers => \@headers,
        separator => "\t",
        is_regex => 1,
    );

    my $citation = $self->_create_citation('PharmGKB', $version, $self->citation_base_url, $self->citation_site_url, $self->citation_text, 'PharmGKB - The Pharmacogenomics Knowledgebase');

    $parser->next; #eat the headers
    while(my $interaction = $parser->next){
        my $drug_name = $self->_import_drug($interaction, $citation);
        my $gene_name = $self->_import_gene($interaction, $citation);
        my $drug_gene_interaction = $self->_create_interaction_report($citation, $drug_name, $gene_name, '');
        push @interactions, $drug_gene_interaction;
        #No interaction type provided in PharmGKB relationships file
        #Set to 'n/a' to be consistent with other data sources where intereaction unknown
        my $type_attribute = $self->_create_interaction_report_attribute($drug_gene_interaction, 'Interaction Type', 'n/a');
    }
    return @interactions;
}

sub _import_drug {
    my $self = shift;
    my $interaction = shift;
    my $citation = shift;
    my $drug_accession = $self->_create_drug_name_report($interaction->{Entity1_id}, $citation, 'PharmGKB', '');
    my $primary_drug_name = $self->_create_drug_alternate_name_report($drug_accession, $interaction->{Entity1_id}, 'Primary Drug Name', '');
    my $drug_name = $self->_create_drug_alternate_name_report($drug_accession, $interaction->{Drug_Name}, 'Drug Name', '');
    my @drug_generic_names = split(",", $interaction->{Generic_Names});
    for my $drug_generic_name (@drug_generic_names){
      next if $drug_generic_name eq 'NA';
      my $synonym_association = $self->_create_drug_alternate_name_report($drug_accession, $drug_generic_name, 'Drug Generic Name', '');
    }
    my @drug_tradenames = split(",", $interaction->{Trade_Names});
    for my $drug_tradename (@drug_tradenames){
      next if $drug_tradename eq 'NA';
      my $tradename_association = $self->_create_drug_alternate_name_report($drug_accession, $drug_tradename, 'Drug Trade Name', '');
    }
    my @drug_cross_references = split(",", $interaction->{Drug_Cross_References});
    for my $drug_cross_reference (@drug_cross_references){
      next if $drug_cross_reference eq 'NA';
      my @data_pair = split(":", $drug_cross_reference);
      my $cross_ref_type=$data_pair[0];
      my $cross_ref_value=$data_pair[1];
      my $cross_reference_association = $self->_create_drug_alternate_name_report($drug_accession, $cross_ref_value, $cross_ref_type, '');
    }
    unless($interaction->{SMILES} eq 'NA'){
        my $SMILES_association = $self->_create_drug_alternate_name_report($drug_accession, $interaction->{SMILES}, 'SMILES', '');
    }

    my @external_vocabs = split(",", $interaction->{External_Vocabulary});
    for my $external_vocab (@external_vocabs){
      next if $external_vocab eq 'NA';
      my @data_pair = split(":", $external_vocab);
      my $external_vocab_type=$data_pair[0];
      my $external_vocab_value=$data_pair[1];
      my $external_vocab_association = $self->_create_drug_category_report($drug_accession, $external_vocab_type, $external_vocab_value,'');
    }
    unless($interaction->{Drug_Type} eq 'NA'){
        my $drug_type_association = $self->_create_drug_category_report($drug_accession, 'Drug Type', $interaction->{Drug_Type}, '');
    }
    return $drug_accession;
}

sub _import_gene {
    my $self = shift;
    my $interaction = shift;
    my $citation = shift;
    my $gene_accession = $self->_create_gene_name_report($interaction->{Entity2_id}, $citation, 'PharmGKB Gene Accession', '');
    my $Entrez_Id_association = $self->_create_gene_alternate_name_report($gene_accession, $interaction->{Entrez_Id}, 'Entrez Gene Id', '');
    my $Ensembl_Id_association = $self->_create_gene_alternate_name_report($gene_accession, $interaction->{Ensembl_Id}, 'Ensembl Gene Id', '');
    my $Gene_Name_association = $self->_create_gene_alternate_name_report($gene_accession, $interaction->{Gene_Name}, 'Gene Name', '');
    my $Gene_Symbol_association = $self->_create_gene_alternate_name_report($gene_accession, $interaction->{Symbol}, 'Gene Symbol', '');
    my @Alternate_Names = quotewords(',', 0, $interaction->{Alternate_Names});
    for my $Alternate_Name (@Alternate_Names){
        next if $Alternate_Name eq 'NA';
        my $alt_name_association = $self->_create_gene_alternate_name_report($gene_accession, $Alternate_Name, 'Alternate Name','');
    }
    my @Alternate_Symbols = quotewords(',', 0, $interaction->{Alternate_Symbols});
        for my $Alternate_Symbol (@Alternate_Symbols){      
        next if $Alternate_Symbol eq 'NA';
        my $alt_symbol_association = $self->_create_gene_alternate_name_report($gene_accession, $Alternate_Symbol, 'Alternate Symbol','');
    }
#   my @gene_cross_references = split(",", $interaction->{Gene_Cross_References});
#    for my $gene_cross_reference (@gene_cross_references){
#        next if $gene_cross_reference eq 'NA';
#        my @data_pair = split(":", $gene_cross_reference);
#        my $cross_ref_type=join("_", "PharmGKB", $data_pair[0]);
#        my $cross_ref_value=$data_pair[1];
#        my $cross_reference_association = $self->_create_gene_alternate_name_report($gene_accession, $cross_ref_value, $cross_ref_type, '');
#    }
    unless($interaction->{Is_VIP} eq 'NA'){
    my $is_vip_association = $self->_create_gene_category_report($gene_accession, 'Is VIP', $interaction->{Is_VIP},'');
    }
    unless($interaction->{Has_Variant_Annotation} eq 'NA'){
    my $has_var_annot_association = $self->_create_gene_category_report($gene_accession, 'Has Variant Annotation', $interaction->{Has_Variant_Annotation},'');
    return $gene_accession;
    }
}

1;
