package Genome::Model::Tools::Dgidb::Import::Entrez;

use strict;
use warnings;

use Genome;
use Term::ANSIColor qw(:constants);
use XML::Simple;

binmode(STDOUT, ":utf8");

my $high = 750000;
UR::Context->object_cache_size_highwater($high);

class Genome::Model::Tools::Dgidb::Import::Entrez {
    is => 'Genome::Model::Tools::Dgidb::Import::Base',
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
gmt dgidb import entrez --gene_info_file=/gscmnt/sata132/techd/mgriffit/DruggableGenes/EntrezGene/gene_info --gene2accession_file=/gscmnt/sata132/techd/mgriffit/DruggableGenes/EntrezGene/gene2accession.gz  --version 3
HELP
}

sub help_detail {
#TODO: make this accurate
    my $summary = <<HELP
WRITE ME
HELP
}

sub execute {
    my $self = shift;
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
    my @headers = qw/ entrez_id entrez_gene_symbol entrez_gene_synonyms /;
    my $parser = Genome::Utility::IO::SeparatedValueReader->create(
        input => $gene_outfile,
        headers => \@headers,
        separator => "\t",
        is_regex=> 1,
    );
    my $citation = $self->_create_citation('Entrez', $version, $self->citation_base_url, $self->citation_site_url, $self->citation_text);

    $parser->next; #eat the headers
    while(my $gene = $parser->next){
        my $gene_name = $gene->{entrez_id};
        my $gene_name_report = $self->_create_gene_name_report($gene_name, $citation, 'entrez_id', '');
        push @gene_name_reports, $gene_name_report;
        my $gene_symbol_association = $self->_create_gene_alternate_name_report($gene_name_report, $gene->{entrez_gene_symbol}, 'entrez_gene_symbol', '');
        my @entrez_gene_synonyms = split(',', $gene->{entrez_gene_synonyms});
        for my $entrez_gene_synonym (@entrez_gene_synonyms){
            if ($entrez_gene_synonym and $entrez_gene_synonym ne 'na'){
                my $gene_alternate_name_report = $self->_create_gene_alternate_name_report($gene_name_report, $entrez_gene_synonym, 'entrez_gene_synonym', '');
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

    my $targets_header = join("\t", 'entrez_id', 'entrez_gene_symbol', 'entrez_gene_synonyms');
    print TARGETS "$targets_header\n";

    my %entrez_ids = %{$entrez_data->{'entrez_ids'}};
    for my $entrez_id (keys %entrez_ids){
        my %entrez_id_names = %{$entrez_ids{$entrez_id}};
        my $synonyms = $entrez_id_names{synonyms_array};
        $synonyms = join(",", @$synonyms);
        my $symbol = $entrez_id_names{symbol};
        print TARGETS join("\t", $entrez_id, $symbol, $synonyms), "\n";
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

    my %entrez_map;      #Entrez_id -> symbol, synonyms
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

        my @synonyms_array = grep{$_ ne '-'} (split("\\|", $synonyms), map{$_ =~ s/Ensembl://; $_} grep{$_ =~ /ENSG/ } split(/\|/, $xref));
        $synonyms = join("|", @synonyms_array);
        if ($synonyms eq ''){
            $synonyms = 'na';
        }

        my %synonyms_hash;
        foreach my $syn (@synonyms_array){
            $synonyms_hash{$syn} = 1;
        }

        #Store entrez info keyed on entrez id
        $entrez_map{$entrez_id}{symbol} = $symbol;
        $entrez_map{$entrez_id}{synonyms_string} = $synonyms;
        $entrez_map{$entrez_id}{synonyms_array} = \@synonyms_array;
        $entrez_map{$entrez_id}{synonyms_hash} = \%synonyms_hash;

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
