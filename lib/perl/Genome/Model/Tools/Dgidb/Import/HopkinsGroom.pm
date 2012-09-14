package Genome::Model::Tools::Dgidb::Import::HopkinsGroom;

use strict;
use warnings;

use Genome;
use IO::File;

my $high = 750000;
UR::Context->object_cache_size_highwater($high);

class Genome::Model::Tools::Dgidb::Import::HopkinsGroom {
    is => 'Genome::Model::Tools::Dgidb::Import::Base',
    has => {
        infile => {
            is => 'Path',
            is_input => 1,
            doc => 'PATH.  text file obtained from publication and curated/updated',
        },
        hopkins_term_file => {
            is => 'Path',
            is_input => 1,
            doc => 'Path to HopkinsGroom term file (provides DGIDB Human Readable Names for HopkinsGroom terms)', 
        },
        tmp_dir => {
            is => 'Path',
            default => '/tmp',
            doc => 'Directory where temp files will be created',
        },
        genes_outfile => {
            is => 'Path',
            is_input => 1,
            default => '/gscmnt/sata132/techd/mgriffit/DruggableGenes/TSV/HopkinsGroom_WashU_TARGETS.tsv',
            doc => 'PATH.  Path to .tsv file for genes (targets)',
        },
        version => {
            doc => 'VERSION.  Version (date) of last update of HopkinsGroom mapping to current Interpro/protein/gene IDs',
        },
        citation_base_url => {
            default => 'http://www.ncbi.nlm.nih.gov/pubmed/12209152', 
        },
        citation_site_url => {
            default => 'http://www.ncbi.nlm.nih.gov/pubmed/12209152', 
        },
        citation_text => {
            default => "The druggable genome. Hopkins AL, Groom CR. Nat Rev Drug Discov. 2002 Sep;1(9):727-30.",
        },
    },
    doc => 'Parse a text file from Hopkins & Groom publication',
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
gmt dgidb import hopkins-groom --infile=/gscuser/ogriffit/Projects/DruggableGenes/PotentiallyDruggable/Hopkins_and_Groom_2002/HopkinsGroomGenes.tsv --hopkins-term-file=/gscuser/ogriffit/Projects/DruggableGenes/PotentiallyDruggable/Hopkins_and_Groom_2002/HopkinsGroomTerms2DGIDB.txt --version=11Sep2012
HELP
}

sub help_detail {
    my $summary = <<HELP
Parse a text file from the Hopkins and Groom publication
Get gene and category info for each record
HELP
}

sub execute {
    my $self = shift;
    $self->input_to_tsv();
    $self->import_tsv();
    return 1;
}

sub input_to_tsv {
    my $self = shift;
    my $out_fh = IO::File->new($self->genes_outfile, 'w');

    my $headers = join("\t", 'relationship_id', 'Interpro_Acc','Interpro_Name','Interpro_Short_Name','Interpro_Type','DGIDB_Human_Readable','Uniprot_Acc','Uniprot_Id','Uniprot_Protein_Name','Uniprot_Gene_Name','Uniprot_Evidence','Uniprot_Status','Entrez_Id','Ensembl_Id');
    $out_fh->print($headers, "\n");

    #Get the data in order
    my $infile_path = $self->infile;
    my $terms_file_path = $self->hopkins_term_file;
    my $targets = $self->_parse_targets_file($infile_path);
    my $terms = $self->_parse_terms_file($terms_file_path);

    #Write data to the file
    for my $target_id (keys %{$targets}){
        #relationship id
        my $relationship_id = $targets->{$target_id}{'relationship_id'};
        $relationship_id = 'NA' unless $relationship_id;

        #Interpro Accession
        my $Interpro_Acc = $targets->{$target_id}{'Interpro_Acc'};
        $Interpro_Acc = 'NA' unless $Interpro_Acc;

        #Uniprot Accession
        my $Uniprot_Acc = $targets->{$target_id}{'Uniprot_Acc'};
        $Uniprot_Acc = 'NA' unless $Uniprot_Acc;

        #Uniprot Id
        my $Uniprot_Id = $targets->{$target_id}{'Uniprot_Id'};
        $Uniprot_Id = 'NA' unless $Uniprot_Id;

        #Uniprot Protein Name
        my $Uniprot_Protein_Name = $targets->{$target_id}{'Uniprot_Protein_Name'};
        $Uniprot_Protein_Name = 'NA' unless $Uniprot_Protein_Name;

        #Uniprot Gene Name
        my $Uniprot_Gene_Name = $targets->{$target_id}{'Uniprot_Gene_Name'};
        $Uniprot_Gene_Name = 'NA' unless $Uniprot_Gene_Name;

        #Uniprot Evidence
        my $Uniprot_Evidence = $targets->{$target_id}{'Uniprot_Evidence'};
        $Uniprot_Evidence = 'NA' unless $Uniprot_Evidence;

        #Uniprot Status
        my $Uniprot_Status = $targets->{$target_id}{'Uniprot_Status'};
        $Uniprot_Status = 'NA' unless $Uniprot_Status;

        #Entrez Id
        my $Entrez_Id = $targets->{$target_id}{'Entrez_Id'};
        $Entrez_Id = 'NA' unless $Entrez_Id;

        #Ensembl Id
        my $Ensembl_Id = $targets->{$target_id}{'Ensembl_Id'};
        $Ensembl_Id = 'NA' unless $Ensembl_Id;

        #Term Category Name
        my $Name = $terms->{$Interpro_Acc}{'Name'};
        $Name = 'NA' unless $Name;

        #Term Category Short Name
        my $Short_Name = $terms->{$Interpro_Acc}{'Short_Name'};
        $Short_Name = 'NA' unless $Short_Name;

        #Term Category Type
        my $Type = $terms->{$Interpro_Acc}{'Type'};
        $Type = 'NA' unless $Type;

        #DGIDB Human Readable
        my $DGIDB_Human_Readable = $terms->{$Interpro_Acc}{'DGIDB_Human_Readable'};
        $DGIDB_Human_Readable = 'NA' unless $DGIDB_Human_Readable;

        $out_fh->print(join("\t", $relationship_id,$Interpro_Acc,$Name,$Short_Name,$Type,$DGIDB_Human_Readable,$Uniprot_Acc,$Uniprot_Id,$Uniprot_Protein_Name,$Uniprot_Gene_Name,$Uniprot_Evidence,$Uniprot_Status,$Entrez_Id,$Ensembl_Id), "\n");
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
        if($line =~ m/^\w+\s+\w+_HUMAN/){
            my ($Uniprot_Acc, $Uniprot_Id, $Uniprot_Protein_Name, $Uniprot_Gene_Name, $Uniprot_Evidence, $Uniprot_Status, $Entrez_Id, $Ensembl_Id, $Interpro_Acc) = split("\t", $line);
            my $relationship_id=join("_",$Uniprot_Acc,$Interpro_Acc);
            $targets->{$relationship_id}{'relationship_id'} = $relationship_id;
            $targets->{$relationship_id}{'Uniprot_Acc'} = $Uniprot_Acc;
            $targets->{$relationship_id}{'Uniprot_Id'} = $Uniprot_Id;
            $targets->{$relationship_id}{'Uniprot_Protein_Name'} = $Uniprot_Protein_Name;
            $targets->{$relationship_id}{'Uniprot_Gene_Name'} = $Uniprot_Gene_Name;
            $targets->{$relationship_id}{'Uniprot_Evidence'} = $Uniprot_Evidence;
            $targets->{$relationship_id}{'Uniprot_Status'} = $Uniprot_Status;
            $targets->{$relationship_id}{'Entrez_Id'} = $Entrez_Id;
            $targets->{$relationship_id}{'Ensembl_Id'} = $Ensembl_Id;
            $targets->{$relationship_id}{'Interpro_Acc'} = $Interpro_Acc;
        }else{
            #skip this line
            next;
        }
    }
    $fh->close;
    return ($targets);
}

sub _parse_terms_file {
    my $self = shift;
    my $terms_path = shift;
    my $terms = {};
    my $fh = IO::File->new($terms_path, 'r');
    my $header = <$fh>;
    while(my $line = <$fh>){
        next unless $line;
        chomp $line;
        $line =~ s/\r//g;
        if($line =~ m/^IPR\d+/){
          my ($Interpro_Acc,$Name,$Short_Name,$Type,$Count,$DGIDB_Human_Readable) = split("\t", $line);
          $terms->{$Interpro_Acc}{'Interpro_Acc'} = $Interpro_Acc;
          $terms->{$Interpro_Acc}{'Name'} = $Name;
          $terms->{$Interpro_Acc}{'Short_Name'} = $Short_Name;
          $terms->{$Interpro_Acc}{'Type'} = $Type;
          $terms->{$Interpro_Acc}{'DGIDB_Human_Readable'} = $DGIDB_Human_Readable;
        }else{
            #skip this line
            next;
        }
    }
    $fh->close;
    return ($terms);
}

sub import_tsv {
    my $self = shift;
    my $genes_outfile = $self->genes_outfile;
    my $citation = $self->_create_citation('HopkinsGroom', $self->version, $self->citation_base_url, $self->citation_site_url, $self->citation_text);
    my @genes = $self->import_genes($genes_outfile, $citation);
    return 1;
}

sub import_genes {
    my $self = shift;
    my $version = $self->version;
    my $genes_outfile = shift;
    my $citation = shift;
    my @genes;
    my @headers = qw/relationship_id Interpro_Acc Interpro_Name Interpro_Short_Name Interpro_Type DGIDB_Human_Readable Uniprot_Acc Uniprot_Id Uniprot_Protein_Name Uniprot_Gene_Name Uniprot_Evidence Uniprot_Status Entrez_Id Ensembl_Id/;
    my $parser = Genome::Utility::IO::SeparatedValueReader->create(
        input => $genes_outfile,
        headers => \@headers,
        separator => "\t",
        is_regex => 1,
    );
    
    $parser->next; #eat the headers
    while(my $hopkins_input = $parser->next){
        #Create gene record with all alternate names
        my $gene_name = $self->_create_gene_name_report($hopkins_input->{'Uniprot_Acc'}, $citation, 'HopkinsGroom_gene_name', '');
        my $uniprot_id = $self->_create_gene_alternate_name_report($gene_name, $hopkins_input->{'Uniprot_Id'}, 'Uniprot_Id', '');
        my $uniprot_protein_name = $self->_create_gene_alternate_name_report($gene_name, $hopkins_input->{'Uniprot_Protein_Name'}, 'Uniprot_Protein_Name', '');
        unless ($hopkins_input->{'Uniprot_Gene_Name'} eq 'NA'){
          my $uniprot_gene_name = $self->_create_gene_alternate_name_report($gene_name, $hopkins_input->{'Uniprot_Gene_Name'}, 'Uniprot_Gene_Name', '');
        }
        unless ($hopkins_input->{'Entrez_Id'} eq 'NA'){
          my $entrez_id = $self->_create_gene_alternate_name_report($gene_name, $hopkins_input->{'Entrez_Id'}, 'Entrez_Id', '');
        }
        unless ($hopkins_input->{'Ensembl_Id'} eq 'NA'){
          my $ensembl_string=$hopkins_input->{'Ensembl_Id'};
          $ensembl_string=~s/\s//;
          my @ensembl_ids=split(";", $ensembl_string);
          foreach my $ensembl_id (@ensembl_ids){
            my $ensembl_id_entry = $self->_create_gene_alternate_name_report($gene_name, $ensembl_id, 'Ensembl_Id', '');
          }
        }
        #Put all genes in HopkinsGroom category as well as any others
        my $human_readable_all = $self->_create_gene_category_report($gene_name, 'human_readable_name', 'HOPKINSGROOM', ''); 
        my $human_readable_name = $hopkins_input->{'DGIDB_Human_Readable'};
        $human_readable_name =~ s/-/ /g;
        $human_readable_name =~ s/\// /g;
        unless ($human_readable_name eq 'NA'){
          my $human_readable = $self->_create_gene_category_report($gene_name, 'human_readable_name', uc($human_readable_name), '');
        }
        #Add additional category details
        my $Interpro_Acc = $self->_create_gene_category_report($gene_name, 'Interpro_Acc', $hopkins_input->{'Interpro_Acc'}, '');
        my $Uniprot_Evidence = $self->_create_gene_category_report($gene_name, 'Uniprot_Evidence', $hopkins_input->{'Uniprot_Evidence'}, '');
        my $Uniprot_Status = $self->_create_gene_category_report($gene_name, 'Uniprot_Status', $hopkins_input->{'Uniprot_Status'}, '');
        my $Interpro_Name = $self->_create_gene_category_report($gene_name, 'Interpro_Name', $hopkins_input->{'Interpro_Name'}, '');
        my $Interpro_Short_Name = $self->_create_gene_category_report($gene_name, 'Interpro_Short_Name', $hopkins_input->{'Interpro_Short_Name'}, '');
        my $Interpro_Type = $self->_create_gene_category_report($gene_name, 'Interpro_Type', $hopkins_input->{'Interpro_Type'}, '');
        push @genes, $gene_name;
    }
    return @genes;
}

1;
