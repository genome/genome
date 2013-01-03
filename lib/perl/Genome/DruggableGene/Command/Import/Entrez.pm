package Genome::DruggableGene::Command::Import::Entrez;

use strict;
use warnings;

use Genome;
use Term::ANSIColor qw(:constants);
use XML::Simple;


class Genome::DruggableGene::Command::Import::Entrez {
    is => 'Genome::DruggableGene::Command::Import::Base',
    has => [
        genes_outfile => {
            is => 'Path',
            is_input => 1,
            default => '/gscmnt/sata132/techd/mgriffit/DruggableGenes/TSV/Entrez_WashU_TARGETS.tsv',
            doc => 'PATH.  Path to .tsv file for genes (targets)',
        },
        gene_info_file => {
            is => 'Path',
            is_input => 1,
            doc => 'PATH. #ftp://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz',
        },
        gene2accession_file => {
            is => 'Path',
            is_input => 1,
            doc => 'PATH.  #ftp://ftp.ncbi.nih.gov/gene/DATA/gene2accession.gz',
        },
        citation_base_url => {
            default => 'http://www.ncbi.nlm.nih.gov/gene?term=',
        },
        citation_site_url => {
            default => 'http://www.ncbi.nlm.nih.gov/gene',
        },
        citation_text => {
            default => "Entrez Gene: gene-centered information at NCBI. Maglott D, Ostell J, Pruitt KD, Tatusova T. Nucleic Acids Res. 2011 Jan;39(Database issue)52-7. Epub 2010 Nov 28. PMID: 21115458.",
        },
    ],
    doc => 'Import genes from Entrez gene info file',
};

sub _doc_copyright_years {
    (2011);
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
 Malachi Griffith, Ph.D.
EOS
}

=cut
sub _doc_credits {
    return ('','None at this time.');
}
=cut

sub _doc_see_also {
    return <<EOS
B<go>(1)
EOS
}

sub _doc_manual_body {
    my $help = shift->help_detail;
    $help =~ s/\n+$/\n/g;
    return $help;
}

sub help_synopsis {
    return <<HELP
genome druggable-gene import entrez --gene-info-file=/gscmnt/sata132/techd/mgriffit/DruggableGenes/EntrezGene/gene_info.human --gene2accession-file=/gscmnt/sata132/techd/mgriffit/DruggableGenes/EntrezGene/gene2accession.human  --version="17-Sep-2012"
HELP
}

sub help_detail {
    my $summary = <<HELP
The Entrez gene database forms the preliminary basis for DGIDB gene groups. All human Entrez gene IDs are imported along with their corresponding symbols, synonyms and Ensembl Gene Id where available. Gene entities from other datasources will be grouped into these Entrez based gene groups by matching their identifiers via Entrez ID, Ensembl ID, gene symbol or synonyms. 
HELP
}

sub execute {
    my $self = shift;
    binmode(STDOUT, ":utf8");
    my $high = 750000;
    UR::Context->object_cache_size_highwater($high);
    $self->input_to_tsv();
    $self->import_tsv();
    return 1;
}

sub import_tsv {
    my $self = shift;
    my $targets_outfile = $self->genes_outfile;
    $self->preload_objects;
    my @gene_name_reports = $self->import_genes($targets_outfile);
    return 1;
}

sub import_genes {
    my $self = shift;
    my $version = $self->version;
    my $gene_outfile = shift;
    my @gene_name_reports;
    my @headers = qw/ entrez_id entrez_gene_symbol entrez_gene_synonyms ensembl_ids description/;
    my $parser = Genome::Utility::IO::SeparatedValueReader->create(
        input => $gene_outfile,
        headers => \@headers,
        separator => "\t",
        is_regex=> 1,
    );
    my $citation = $self->_create_citation('Entrez', $version, $self->citation_base_url, $self->citation_site_url, $self->citation_text, 'NCBI Entrez Gene');

    $parser->next; #eat the headers
    while(my $gene = $parser->next){
        my $gene_name = $gene->{entrez_id};
        my $gene_name_report = $self->_create_gene_name_report($gene_name, $citation, 'Entrez Gene Id', '');
        push @gene_name_reports, $gene_name_report;
        my $description = $gene->{description};
        my $desc_alt = $self->_create_gene_alternate_name_report($gene_name_report, $gene->{description}, 'Gene Description', '', 'lower');
        my $gene_name_alt = $self->_create_gene_alternate_name_report($gene_name_report, $gene->{entrez_id}, 'Entrez Gene Id', '', 'upper');
        my $gene_symbol_association = $self->_create_gene_alternate_name_report($gene_name_report, $gene->{entrez_gene_symbol}, 'Gene Symbol', '', 'upper');
        my @entrez_gene_synonyms = split(',', $gene->{entrez_gene_synonyms});
        for my $entrez_gene_synonym (@entrez_gene_synonyms){
            if ($entrez_gene_synonym and $entrez_gene_synonym ne 'N/A'){
                my $gene_alternate_name_report = $self->_create_gene_alternate_name_report($gene_name_report, $entrez_gene_synonym, 'Gene Synonym', '', 'upper');
            }
        }
        my @ensembl_gene_ids = split(',', $gene->{ensembl_ids});
        for my $ensembl_gene_id (@ensembl_gene_ids){
            if ($ensembl_gene_id and $ensembl_gene_id ne 'N/A'){
                my $ensembl_alternate_name_report = $self->_create_gene_alternate_name_report($gene_name_report, $ensembl_gene_id, 'Ensembl Gene Id', '', 'upper');
            }
        }
    }
    return @gene_name_reports;
}

sub preload_objects {
    my $self = shift;
    my $source_db_name = 'Entrez';
    my $source_db_version = $self->version;

    #Let's preload anything for this database name and version so that we can avoid death by 1000 queries
    my @gene_name_reports = Genome::DruggableGene::GeneNameReport->get(source_db_name => $source_db_name, source_db_version => $source_db_version);
    for my $gene_name_report (@gene_name_reports){
        $gene_name_report->gene_alt_names;
        $gene_name_report->gene_categories;
    }
    return 1;
}

sub help_usage_complete_text {
    my $self = shift;
    my $usage = $self->SUPER::help_usage_complete_text(@_);
    return GREEN . $usage . RESET;
}

sub input_to_tsv {
    my $self = shift;

    #Parse Entrez flatfiles
    my $entrez_data = $self->loadEntrezData();
    $self->{entrez_data}=$entrez_data;

    my $targets_outfile = $self->genes_outfile;
    open(TARGETS, ">$targets_outfile") || die "\n\nCould not open outfile: $targets_outfile\n\n";
    binmode(TARGETS, ":utf8");

    my $targets_header = join("\t", 'entrez_id', 'entrez_gene_symbol', 'entrez_gene_synonyms', 'ensembl_ids', 'description');
    print TARGETS "$targets_header\n";

    my %entrez_ids = %{$entrez_data->{'entrez_ids'}};
    for my $entrez_id (keys %entrez_ids){
        my %entrez_id_names = %{$entrez_ids{$entrez_id}};
        my $synonyms = $entrez_id_names{synonyms_array};
        $synonyms = join(",", @$synonyms);
        my $ensembl_ids = $entrez_id_names{ensembl_array};
        $ensembl_ids = join(",", @$ensembl_ids);
        my $symbol = $entrez_id_names{symbol};
        my $description = $entrez_id_names{description};
        print TARGETS join("\t", $entrez_id, $symbol, $synonyms, $ensembl_ids, $description), "\n";
    }

    close(TARGETS);

    print "\n\n";
    return 1;
}

#######################################################################################################################################################################
#Load Entrez Data from flatfiles                                                                                                                                      #
#######################################################################################################################################################################
sub loadEntrezData {
    my $self = shift;
    my %edata;

    my %entrez_map;      #Entrez_id -> symbol, synonyms, Ensembl gene id
    my %symbols_map;     #Symbols   -> entrez_id(s)
    my %synonyms_map;    #Synonyms  -> entrez_id(s)
    my %p_acc_map;       #Protein accessions -> entrez_id

    my $gene2accession_file = $self->gene2accession_file;
    my $gene_info_file = $self->gene_info_file;
    open (GENE, "$gene_info_file") || die "\n\nCould not open gene_info file: $gene_info_file\n\n";
    while(<GENE>){
        chomp($_);
        if ($_ =~ /^\#/){
            next();
        }
        my @line = split("\t", $_);
        my $tax_id = $line[0];
        #Skip all non-human records
        unless ($tax_id eq "9606"){
            next();
        }
        my $entrez_id = $line[1];
        my $symbol = $line[2];
        my $synonyms = $line[4];
        my $xref = $line[5]; #we will ignore all cross references except Ensembl for now
        my $description = $line[8]; #Be aware this string contains commas, semicolons, and doubtless many other potentially problematic characters 

        #my @synonyms_array = grep{$_ ne '-'} (split("\\|", $synonyms), map{$_ =~ s/Ensembl://; $_} grep{$_ =~ /ENSG/ } split(/\|/, $xref));
        my @synonyms_array = grep{$_ ne '-'} split("\\|", $synonyms);
        my @ensembl_array = grep{$_ ne '-'} map{$_ =~ s/Ensembl://; $_} grep{$_ =~ /ENSG/ } split(/\|/, $xref);
        $synonyms = join("|", @synonyms_array);
        if ($synonyms eq ''){
            $synonyms = 'N/A';
        }
        my $ensembl_string = join("|", @ensembl_array);
        if ($ensembl_string eq ''){
          $ensembl_string = 'N/A';
        }
        my %synonyms_hash;
        foreach my $syn (@synonyms_array){
            $synonyms_hash{$syn} = 1;
        }
        my %ensembl_hash;
        foreach my $ensembl (@ensembl_array){
          $ensembl_hash{$ensembl} = 1;
        }

        #Store entrez info keyed on entrez id
        $entrez_map{$entrez_id}{symbol} = $symbol;
        $entrez_map{$entrez_id}{description} = $description;
        $entrez_map{$entrez_id}{synonyms_string} = $synonyms;
        $entrez_map{$entrez_id}{synonyms_array} = \@synonyms_array;
        $entrez_map{$entrez_id}{synonyms_hash} = \%synonyms_hash;
        $entrez_map{$entrez_id}{ensembl_string} = $ensembl_string;
        $entrez_map{$entrez_id}{ensembl_array} = \@ensembl_array;
        $entrez_map{$entrez_id}{ensembl_hash} = \%ensembl_hash;

        #Store entrez info keyed on symbol
        if ($symbols_map{$symbol}){
            my $ids = $symbols_map{$symbol}{entrez_ids};
            $ids->{$entrez_id} = 1;
        }else{
            my %tmp;
            $tmp{$entrez_id} = 1;
            $symbols_map{$symbol}{entrez_ids} = \%tmp;
        }

        #Store synonym to entrez_id mappings
        foreach my $syn (@synonyms_array){
            if ($synonyms_map{$syn}){
                my $ids = $synonyms_map{$syn}{entrez_ids};
                $ids->{$entrez_id} = 1;
            }else{
                my %tmp;
                $tmp{$entrez_id} = 1;
                $synonyms_map{$syn}{entrez_ids} = \%tmp;
            }
        }
    }
    close (GENE);

    open (ACC, "$gene2accession_file") || die "\n\nCould not open gene2accession file: $gene2accession_file\n\n";
    while(<ACC>){
        chomp($_);
        if ($_ =~ /^\#/){
            next();
        }
        my @line = split("\t", $_);
        my $tax_id = $line[0];
        #Skip all non-human records
        unless ($tax_id eq "9606"){
            next();
        }
        my $entrez_id = $line[1];
        my $prot_id = $line[5];
        #If the prot is not defined, skip
        if ($prot_id eq "-"){next();}
        #Clip the version number
        if ($prot_id =~ /(\w+)\.\d+/){
            $prot_id = $1;
        }
        if ($p_acc_map{$prot_id}){
            my $ids = $p_acc_map{$prot_id}{entrez_ids};
            $ids->{$entrez_id} = 1;
        }else{
            my %tmp;
            $tmp{$entrez_id} = 1;
            $p_acc_map{$prot_id}{entrez_ids} = \%tmp;
        }
    }
    close (ACC);

    $edata{'entrez_ids'} = \%entrez_map;
    $edata{'symbols'} = \%symbols_map;
    $edata{'synonyms'} = \%synonyms_map;
    $edata{'protein_accessions'} = \%p_acc_map;

    return(\%edata);
}

1;
