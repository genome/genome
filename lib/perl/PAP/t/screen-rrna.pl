#!/usr/bin/env genome-perl

use strict;
use warnings;

use PAP::Command::RRnaScreen;
use Getopt::Long;
use File::Slurp;
use Carp;
use BAP::DB::GeneTag;
use BAP::DB::Tag;
use BAP::DB::CodingGene;
use BAP::DB::RnaGene;
use BAP::DB::TrnaGene;

my $fastafile = undef;
my $devflag = 0;
my $blastdb = "/gscmnt/278/analysis/HGMI/rRNA_testing/16s_23srnadb";
my $outfile = undef;

GetOptions(
           "fasta=s" => \$fastafile,
           "dev"     => \$devflag,
           "blast-db=s" => \$blastdb,
           "output=s" => \$outfile,
);

if($devflag)
{
    $BAP::DB::DBI::db_env = 'dev';
}

unless(defined($fastafile) )
{
    print STDERR "usage: $0 --fasta <fastafile> --output <outfile>\n";
    print STDERR "   [--dev] [--blast-db <rrna blastdb>]\n";
    croak "need to provide fasta file and output arguments";
}

my $rrs = PAP::Command::RRnaScreen->create(
                                           'fasta_file' => $fastafile,
                                           'dev_flag' => $devflag,
                                           'blast_db' => $blastdb,
                                           'upload_to_biosql' => 0,
) or croak "can't create rrna screen : $!";

$rrs->execute() or croak "can't execute rrna screen : $!";
my $dead_genes = $rrs->dead_genes();

my @lines = ();
foreach my $dead (@{$dead_genes})
{
    push(@lines,"Sequence ".$dead."\n");
    push(@lines,"Dead\n");
    push(@lines,"\n");
    # mark as dead in database
    tag_gene($dead);
    # should do a commit here or something.
}

BAP::DB::DBI->dbi_commit();
if(defined($outfile))
{
    write_file($outfile,@lines);
}


# subroutines
sub tag_gene
{
    my ($gene_name) = @_;
    my $gene_iterator = BAP::DB::CodingGene->search({'gene_name' => $gene_name});
    if(defined($gene_iterator))
    {
        # tag the gene;
        my $gene = $gene_iterator->next;
        my ($tag) = BAP::DB::Tag->search({'tag_name'  => 'Dead',
                                          'tag_value' => 'rrna fragment hit'});
        my ($genetag) = BAP::DB::GeneTag->insert({'gene_id' => $gene->gene_id,
                                                  'tag_id'  => $tag->tag_id
                                                 });
    }
    else
    {
        carp "can't find $gene_name in mgap database!";
    }
    return;
}

__END__

# $Id$
