package Genome::DruggableGene::Command::Import::Dgene;

use strict;
use warnings;

use Genome;
use IO::File;

my $high = 750000;
UR::Context->object_cache_size_highwater($high);

class Genome::DruggableGene::Command::Import::Dgene {
    is => 'Genome::DruggableGene::Command::Import::Base',
    has => {
        infile => {
            is => 'Path',
            is_input => 1,
            doc => 'PATH.  csv file provided by collaborators',
        },
        dgene_term_file => {
            is => 'Path',
            is_input => 1,
            doc => 'Path to dGene term file (provides DGIDB Human Readable Names for dGene terms)', 
        },
        tmp_dir => {
            is => 'Path',
            default => '/tmp',
            doc => 'Directory where temp files will be created',
        },
        genes_outfile => {
            is => 'Path',
            is_input => 1,
            default => '/gscmnt/sata132/techd/mgriffit/DruggableGenes/TSV/dGene_WashU_TARGETS.tsv',
            doc => 'PATH.  Path to .tsv file for genes (targets)',
        },
        version => {
            doc => 'VERSION.  Version (date) of release of database from dGene group',
        },
        citation_base_url => {
            default => 'http://www.ncbi.nlm.nih.gov/gene?term=', #Since dGene is based on Entrez genes, dGene records can be linked out to there
        },
        citation_site_url => {
            default => 'http://hematology.wustl.edu/faculty/Bose/BoseBio.html', #To be updated upon publication
        },
        citation_text => {
            default => "The Druggable Gene List, dGENE, provides a Rapid Filter for Cancer Genome Sequencing Data. Kumar R, Chang L, Ellis MJ, Bose R. Manuscript in preparation.",
        },
    },
    doc => 'Parse a csv file from dGene (Ron Bose)',
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
genome druggable-gene import dgene --infile=/gscmnt/sata132/techd/mgriffit/DruggableGenes/PotentiallyDruggable/RonBose/02Aug2012/DRUGGABLE_GENOME_08022012_2257.txt --dgene-term-file=/gscmnt/sata132/techd/mgriffit/DruggableGenes/PotentiallyDruggable/RonBose/02Aug2012/TargetDgeneTerms.tsv --version="02-Aug-2012" 
HELP
}

sub help_detail {
    my $summary = <<HELP
Parse a csv file from the dGene group (produced for a publication by Ron Bose et al)
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

    my $headers = join("\t", 'dgene_class', 'dgene_short_name', 'dgene_full_name', 'human_readable_name', 'dGeneID', 'tax_id', 'GeneID', 'Symbol', 'LocusTag', 'Synonyms', 'dbXrefs', 'chromosome', 'map_location', 'description', 'gene_type', 'symbol_from_authority', 'full_name', 'nomenclature_status', 'other_designations', 'modification_date');
    $out_fh->print($headers, "\n");

    #Get the data in order
    my $infile_path = $self->infile;
    my $terms_file_path = $self->dgene_term_file;
    my $targets = $self->_parse_targets_file($infile_path);
    my $terms = $self->_parse_terms_file($terms_file_path);


    #Write data to the file
    for my $target_id (keys %{$targets}){
        #dGene Id
        my $dGeneID =  $targets->{$target_id}{'dGeneID'};
        $dGeneID = 'N/A' unless $dGeneID;

        #Target Gene Taxonomy ID
        my $tax_id =  $targets->{$target_id}{'tax_id'};
        $tax_id = 'N/A' unless $tax_id;

        #Target Gene Entrez Id
        my $GeneID =  $targets->{$target_id}{'GeneID'};
        $GeneID = 'N/A' unless $GeneID;

        #Target Gene Symbol
        my $Symbol = $targets->{$target_id}{'Symbol'};
        $Symbol = 'N/A' unless $Symbol;

        #Target Gene Locus Tag
        my $LocusTag = $targets->{$target_id}{'LocusTag'};
        $LocusTag = 'N/A' unless $LocusTag;

        #Target Gene Synonyms
        my $Synonyms = $targets->{$target_id}{'Synonyms'};
        $Synonyms = 'N/A' unless $Synonyms;

        #Target Gene dbXrefs
        my $dbXrefs = $targets->{$target_id}{'dbXrefs'};
        $dbXrefs = 'N/A' unless $dbXrefs;

        #Target Gene Chromosome
        my $chromosome = $targets->{$target_id}{'chromosome'};
        $chromosome = 'N/A' unless $chromosome;

        #Target Gene Map Location
        my $map_location = $targets->{$target_id}{'map_location'};
        $map_location = 'N/A' unless $map_location;

        #Target Gene Description
        my $description = $targets->{$target_id}{'description'};
        $description = 'N/A' unless $description;

        #Target Gene Type
        my $gene_type = $targets->{$target_id}{'gene_type'};
        $gene_type = 'N/A' unless $gene_type;

        #Target Gene Symbol From Authority
        my $symbol_from_authority = $targets->{$target_id}{'symbol_from_authority'};
        $symbol_from_authority = 'N/A' unless $symbol_from_authority;

        #Target Gene Full Name
        my $full_name = $targets->{$target_id}{'full_name'};
        $full_name = 'N/A' unless $full_name;

        #Target Gene Nomenclature Status
        my $nomenclature_status = $targets->{$target_id}{'nomenclature_status'};
        $nomenclature_status = 'N/A' unless $nomenclature_status;

        #Target Gene Other Designations
        my $other_designations = $targets->{$target_id}{'other_designations'};
        $other_designations = 'N/A' unless $other_designations;

        #Target Gene Modification Date
        my $modification_date = $targets->{$target_id}{'modification_date'};
        $modification_date = 'N/A' unless $modification_date;

        #Target Gene Druggable Class (Term Category)
        my $class = $targets->{$target_id}{'class'};
        $class = 'N/A' unless $class;

        #Term Category Short Name
        my $ShortName = $terms->{$class}{'ShortName'};
        $ShortName = 'N/A' unless $ShortName;

        #Term Category Full Name
        my $FullName = $terms->{$class}{'FullName'};
        $FullName = 'N/A' unless $FullName;

        #Term Category Human Readable Name
        my $HumanReadableName = $terms->{$class}{'HumanReadableName'};
        $HumanReadableName = 'N/A' unless $HumanReadableName;

        $out_fh->print(join("\t", $class,$ShortName,$FullName,$HumanReadableName,$dGeneID, $tax_id,$GeneID,$Symbol,$LocusTag,$Synonyms,$dbXrefs,$chromosome,$map_location,$description,$gene_type,$symbol_from_authority,$full_name,$nomenclature_status,$other_designations,$modification_date), "\n");
        
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
        if($line =~ m/^dGene\d+/){
            my ($dGeneID, $tax_id,$GeneID,$Symbol,$LocusTag,$Synonyms,$dbXrefs,$chromosome,$map_location,$description,$gene_type,$symbol_from_authority,$full_name,$nomenclature_status,$other_designations,$modification_date,$class) = split("\t", $line);
        $targets->{$dGeneID}{'dGeneID'} = $dGeneID;
        $targets->{$dGeneID}{'tax_id'} = $tax_id;
        $targets->{$dGeneID}{'GeneID'} = $GeneID;
        $targets->{$dGeneID}{'Symbol'} = $Symbol;
        $targets->{$dGeneID}{'LocusTag'} = $LocusTag;
        $targets->{$dGeneID}{'Synonyms'} = $Synonyms;
        $targets->{$dGeneID}{'dbXrefs'} = $dbXrefs;
        $targets->{$dGeneID}{'chromosome'} = $chromosome;
        $targets->{$dGeneID}{'map_location'} = $map_location;
        $targets->{$dGeneID}{'description'} = $description;
        $targets->{$dGeneID}{'gene_type'} = $gene_type;
        $targets->{$dGeneID}{'symbol_from_authority'} = $symbol_from_authority;
        $targets->{$dGeneID}{'full_name'} = $full_name;
        $targets->{$dGeneID}{'nomenclature_status'} = $nomenclature_status;
        $targets->{$dGeneID}{'other_designations'} = $other_designations;
        $targets->{$dGeneID}{'modification_date'} = $modification_date;
        $targets->{$dGeneID}{'class'} = $class;
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
        if($line =~ m/\S+\t\S+\t.+\t\S+/){
            my ($ShortName,$FullName,$HumanReadableName,$DgeneClass) = split("\t", $line);
        $terms->{$DgeneClass}{'DgeneClass'} = $DgeneClass;
        $terms->{$DgeneClass}{'ShortName'} = $ShortName;
        $terms->{$DgeneClass}{'FullName'} = $FullName;
        $terms->{$DgeneClass}{'HumanReadableName'} = $HumanReadableName;
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
    my $citation = $self->_create_citation('dGene', $self->version, $self->citation_base_url, $self->citation_site_url, $self->citation_text, 'dGENE - The Druggable Gene List');
    my @genes = $self->import_genes($genes_outfile, $citation);
    return 1;
}

sub import_genes {
    my $self = shift;
    my $version = $self->version;
    my $genes_outfile = shift;
    my $citation = shift;
    my @genes;
    my @headers = qw/dgene_class dgene_short_name dgene_full_name human_readable_name dGeneID tax_id GeneID Symbol LocusTag Synonyms dbXrefs chromosome map_location description gene_type symbol_from_authority full_name nomenclature_status other_designations modification_date/;
    my $parser = Genome::Utility::IO::SeparatedValueReader->create(
        input => $genes_outfile,
        headers => \@headers,
        separator => "\t",
        is_regex => 1,
    );
    
    $parser->next; #eat the headers
    while(my $dgene_input = $parser->next){
        my $gene_name = $self->_create_gene_name_report($dgene_input->{'GeneID'}, $citation, 'dGene Gene Id', '');
        my $human_readable_name = $dgene_input->{'human_readable_name'};
        $human_readable_name =~ s/-/ /g;
        my $human_readable = $self->_create_gene_category_report($gene_name, 'Human Readable Name', uc($human_readable_name), '');
        my $symbol = $self->_create_gene_alternate_name_report($gene_name, $dgene_input->{'Symbol'}, 'Gene Symbol', '', 'upper');
        my $entrez_id = $self->_create_gene_alternate_name_report($gene_name, $dgene_input->{'GeneID'}, 'Entrez Gene Id', '', 'upper');
        push @genes, $gene_name;
    }
    return @genes;
}

1;
