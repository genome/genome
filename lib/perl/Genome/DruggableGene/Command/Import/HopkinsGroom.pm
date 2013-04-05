package Genome::DruggableGene::Command::Import::HopkinsGroom;

use strict;
use warnings;

use Genome;
use IO::File;

class Genome::DruggableGene::Command::Import::HopkinsGroom {
    is => 'Genome::DruggableGene::Command::Import::Base',
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
            example_values => ['/gscmnt/sata132/techd/mgriffit/DruggableGenes/TSV/HopkinsGroom_WashU_TARGETS.tsv'],
            doc => 'PATH.  Path to .tsv file for genes (targets)',
        },
        version => {
            doc => 'VERSION.  Version (date) of last update of HopkinsGroom mapping to current Interpro/protein/gene IDs',
        },
        citation_base_url => {
            default => 'http://www.uniprot.org/uniprot/', #HopkinsGroom genes are uniprot ids - link to there for lack of better place
        },
        citation_site_url => {
            default => 'http://www.ncbi.nlm.nih.gov/pubmed/12209152/', 
        },
        citation_text => {
            default => "The druggable genome. Hopkins AL, Groom CR. Nat Rev Drug Discov. 2002 Sep;1(9):727-30. PMID: 12209152",
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
genome druggable-gene import hopkins-groom --infile=/gscuser/ogriffit/Projects/DruggableGenes/PotentiallyDruggable/Hopkins_and_Groom_2002/HopkinsGroomGenes.tsv --hopkins-term-file=/gscuser/ogriffit/Projects/DruggableGenes/PotentiallyDruggable/Hopkins_and_Groom_2002/HopkinsGroomTerms2DGIDB.txt --version="11-Sep-2012"
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
    $self->status_message(".tsv file created for use in 'rake dgidb:import:hopkins_groom <.tsv_file_path> <source_db_version>'");

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
        $relationship_id = 'N/A' unless $relationship_id;

        #Interpro Accession
        my $Interpro_Acc = $targets->{$target_id}{'Interpro_Acc'};
        $Interpro_Acc = 'N/A' unless $Interpro_Acc;

        #Uniprot Accession
        my $Uniprot_Acc = $targets->{$target_id}{'Uniprot_Acc'};
        $Uniprot_Acc = 'N/A' unless $Uniprot_Acc;

        #Uniprot Id
        my $Uniprot_Id = $targets->{$target_id}{'Uniprot_Id'};
        $Uniprot_Id = 'N/A' unless $Uniprot_Id;

        #Uniprot Protein Name
        my $Uniprot_Protein_Name = $targets->{$target_id}{'Uniprot_Protein_Name'};
        $Uniprot_Protein_Name = 'N/A' unless $Uniprot_Protein_Name;

        #Uniprot Gene Name
        my $Uniprot_Gene_Name = $targets->{$target_id}{'Uniprot_Gene_Name'};
        $Uniprot_Gene_Name = 'N/A' unless $Uniprot_Gene_Name;

        #Uniprot Evidence
        my $Uniprot_Evidence = $targets->{$target_id}{'Uniprot_Evidence'};
        $Uniprot_Evidence = 'N/A' unless $Uniprot_Evidence;

        #Uniprot Status
        my $Uniprot_Status = $targets->{$target_id}{'Uniprot_Status'};
        $Uniprot_Status = 'N/A' unless $Uniprot_Status;

        #Entrez Id
        my $Entrez_Id = $targets->{$target_id}{'Entrez_Id'};
        $Entrez_Id = 'N/A' unless $Entrez_Id;

        #Ensembl Id
        my $Ensembl_Id = $targets->{$target_id}{'Ensembl_Id'};
        $Ensembl_Id = 'N/A' unless $Ensembl_Id;

        #Term Category Name
        my $Name = $terms->{$Interpro_Acc}{'Name'};
        $Name = 'N/A' unless $Name;

        #Term Category Short Name
        my $Short_Name = $terms->{$Interpro_Acc}{'Short_Name'};
        $Short_Name = 'N/A' unless $Short_Name;

        #Term Category Type
        my $Type = $terms->{$Interpro_Acc}{'Type'};
        $Type = 'N/A' unless $Type;

        #DGIDB Human Readable
        my $DGIDB_Human_Readable = $terms->{$Interpro_Acc}{'DGIDB_Human_Readable'};
        $DGIDB_Human_Readable = 'N/A' unless $DGIDB_Human_Readable;

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

1;
