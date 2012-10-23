package Genome::Model::Tools::Hgmi::DumpProtein;

use strict;
use warnings;

use Genome;
use Command;
use Carp;

use Bio::SeqIO;
use Bio::Seq;
use Bio::DB::BioDB;
use Bio::DB::Query::BioQuery;

use File::Slurp;
use File::Temp qw/ tempfile tempdir /;
use DateTime;
use List::MoreUtils qw/ uniq /;
use IPC::Run;
use Workflow::Simple;
use Data::Dumper;


UR::Object::Type->define(
    class_name => __PACKAGE__,
    is         => 'Command',
    has        => [
        'locus_tag' => {
            is  => 'String',
            doc => "taxonomy name"
        },

        'pep_file' => {
            is          => 'String',
            doc         => "fasta file of gene proteins",

        },
        'dev' => {
            is => 'Boolean',
            doc => "use development databases",
            is_optional => 1,
        },
    ]
);

sub help_brief
{
    "Dumps to file the protein sequences.";
}

sub help_synopsis
{
    return <<"EOS"
Bridges between HGMI tools and PAP.  This tool pulls data from biosql,
then initializes and runs the PAP workflow.
EOS
}

sub help_detail
{
    return <<"EOS"
Bridges between HGMI tools and PAP.  This tool loads predictions from mgap to
biosql, pulls data from biosql, then initializes and runs the PAP workflow.
EOS
}

sub execute
{
    my $self = shift;

    $self->get_gene_peps();
    
    return 1;
}

sub get_gene_peps
{
    my $self = shift;
    # this needs to handle switching to
    # either dev or prod
    my $dbadp;

    if($self->dev())
    {
        $dbadp = Bio::DB::BioDB->new(
            -database => 'biosql',
            -user     => 'sg_user',
            -pass     => 'sgus3r',
            -dbname   => 'DWDEV',
            -driver   => 'Oracle'
        );
    }
    else
    {
        $dbadp = Bio::DB::BioDB->new(
            -database => 'biosql',
            -user     => 'sg_user',
            -pass     => 'sg_us3r',
            -dbname   => 'DWRAC',
            -driver   => 'Oracle'
        );
    }

    my $file = $self->pep_file();
    my $seqout = new Bio::SeqIO(
        -file   => ">$file",
        -format => "fasta"
    );

    my $adp = $dbadp->get_object_adaptor("Bio::SeqI");

    $adp->dbh->{'LongTruncOk'} = 0;
    $adp->dbh->{'LongReadLen'} = 1000000000;

    my $query = Bio::DB::Query::BioQuery->new();
    $query->datacollections( [ "Bio::PrimarySeqI s", ] );

    my $locus_tag = $self->locus_tag;
    $query->where( ["s.display_id like '$locus_tag%'"] );
    my $res = $adp->find_by_query($query);
    while ( my $seq = $res->next_object() )
    {
        my $gene_name = $seq->display_name();

        #print $gene_name, "\n";
        my @feat = $seq->get_SeqFeatures();
GENE:        foreach my $f (@feat)
        {
            my $display_name = $f->display_name();
            #print STDERR $display_name," ", $f->primary_tag,"\n";
            next GENE if $f->primary_tag ne 'gene';
            next if $f->has_tag('Dead');
            my $ss;
            $ss = $seq->subseq( $f->start, $f->end );

            unless(defined($ss)) { 
                die "failed to fetch sequence for '$display_name' from BioSQL";
            }
            
            my $pep = $self->make_into_pep( $f->strand,
                $ss,
                #$seq->subseq($f->start,$f->end),
                $display_name );
            $seqout->write_seq($pep);
            #print STDERR "sequence should be written out\n";
        }
    }
    if(! -f $file )
    {
        print STDERR "the fasta file $file doesn't exist! SendToPap.pm\n";
        return 0;
    }
    unless( -s $file > 0 )
    {
        print STDERR "the fasta file $file still empty! SendToPap.pm\n";
    }

    return 1;
}

sub make_into_pep
{
    my ( $self, $strand, $subseq, $name ) = @_;

    my $seq = new Bio::Seq(
        -display_id => $name,
        -seq        => $subseq
    );
    my $newseq;
    if ( $strand < 0 )
    {
        $newseq = $seq->revcom->translate->seq;
    }
    else
    {
        $newseq = $seq->translate->seq;
    }
    my $seqobj = new Bio::Seq(
        -display_id => $name,
        -seq        => $newseq
    );
    return $seqobj;
}


1;

# $Id$
