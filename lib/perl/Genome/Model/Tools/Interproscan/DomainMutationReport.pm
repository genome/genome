package Genome::Model::Tools::Interproscan::DomainMutationReport;

use strict;
use warnings;

# this probably need a bit more documentation

use Genome;
use Command;
use Carp;
use IO::File;
use File::Slurp;
use Text::CSV_XS;
use File::Temp qw/ tempfile /;
use List::MoreUtils qw/ uniq /;
use SnpDom;

UR::Object::Type->define(
    class_name => __PACKAGE__,
    is         => 'Command',
    has        => [
        'maf' => {
            is       => 'String',
            doc      => "mutation annotation file",
        },
        'filter' => {
            is          => 'String',
            doc         => "filter expression for domains",
            is_optional => 1
        },
        'output' => {
            is       => 'String',
            doc      => "output file name",
        },
        'gcol' => {
            is          => 'Integer',
            doc         => "hugo gene name column",
            is_optional => 1,
            default     => 0
        },
        'tcol' => {
            is          => 'Integer',
            doc         => "transcript name column",
            is_optional => 1,
            default     => 3
        },
        'acol' => {
            is          => 'Integer',
            doc         => "Amino acid change column",
            is_optional => 1,
            default     => 12
        },
        'sepchar' => {
            is          => 'String',
            doc         => "separator character in matrix",
            is_optional => 1,
            default     => "\t"
        },
    ]
);

sub help_brief
{
"tool for generating reports of how many mutations fall into specific protein domains";
}

# there's much similarity with this and AppendDomain.pm
sub execute
{
    my $self = shift;

    # preprocess a temp file for the 'snp file'.

    my $tmpfile = $self->extract_snp_data2tmpfile();

    # do SnpDom stuff.
    my $snpdomains = new SnpDom( '-inc-ts' => 1, );
    my $muthash;
    my %pfamlen;
    my $mutcounts;
    my $flag = 1;    # probably need to get rid of this...

    # read in mutations
    $snpdomains->read_muts( $tmpfile, "1,2" );
    $muthash = $snpdomains->get_all_mutations();
    $snpdomains->mutation_in_dom( \%pfamlen, $self->filter );
    $mutcounts = $snpdomains->get_mutation_counts();

    my $counts  = $self->gene_counts( $muthash, $flag );
    my $peplens = undef;

    $peplens = $snpdomains->get_lengths_from_db();

    my @outputlines;
    foreach my $g ( sort keys %$mutcounts )
    {
        foreach my $dom ( sort keys %{ $mutcounts->{$g} } )
        {
            if ( !exists( $counts->{$g} ) )
            {
                print STDERR "NO MUTATION COUNTS FOR ", $g, "\n";
                $counts->{$g} = 0;
            }
            if ( !exists( $peplens->{$g} ) )
            {
                print STDERR "NO PEPLENGTHS FOR ", $g, "\n";
                $peplens->{$g} = 0;
            }
            my $result               = $peplens->{$g} - $pfamlen{$dom}{$g};
            my $mutations_not_in_dom = $counts->{$g} - $mutcounts->{$g}->{$dom};
            if ( $mutations_not_in_dom < 0 )
            {
                warn "negative mut count!";
            }
            my ( $ts_name, $hugo ) = split( /,/x, $g );
            my $line = join( "\t",
                $hugo, $ts_name, $dom, $mutcounts->{$g}->{$dom},
                $mutations_not_in_dom, $pfamlen{$dom}{$g}, $result );

            #my $line = join("\t",$g,$dom, $mutcounts->{$g}->{$dom},
            #                    $mutations_not_in_dom,
            #                    $pfamlen{$dom}{$g},
            #                    $result);

            $line .= "\n";

            push( @outputlines, $line );
        }
    }

    write_file( $self->output, @outputlines );
    unlink($tmpfile);

    return 1;
}

sub extract_snp_data2tmpfile
{
    my $self = shift;
    my ( $fh, $tmpfile ) =
      tempfile( "gt-dom-mut-repXXXXXX", SUFFIX => '.dat' );
    my $c = new Text::CSV_XS( { sep_char => $self->sepchar } );
    my @lines = read_file( $self->maf );
    my @tmp;
    foreach my $l (@lines)
    {
        $c->parse($l);
        my @f = $c->fields();

        my $nr =
            $f[ $self->gcol ] . "\t"
          . $f[ $self->tcol ] . ","
          . $f[ $self->gcol ] . "\t"
          . $f[ $self->acol ] . "\n";
        push( @tmp, $nr );
    }
    write_file( $tmpfile, @tmp );
    return $tmpfile;
}

sub gene_counts
{
    my $self = shift;
    my ( $href, $_flag ) = @_;
    my $cref;
    foreach my $g ( sort keys %{$href} )
    {
        if ( !defined($_flag) )
        {
            my ( $transcript, $gene ) = split( /,/x, $g );
            if ( !exists( $cref->{$gene} ) )
            {
                $cref->{$gene} = $#{ $href->{$g} } + 1;
            }
            else
            {
                $cref->{$gene} = $cref->{$gene} + $#{ $href->{$g} } + 1;
            }
        }
        else
        {

            #my ($transcript,$gene) = split(/,/x,$g);
            if ( !exists( $cref->{$g} ) )
            {
                $cref->{$g} = $#{ $href->{$g} } + 1;
            }
            else
            {
                $cref->{$g} = $cref->{$g} + $#{ $href->{$g} } + 1;
            }
        }

    }
    return $cref;
}

1;

# $Id$
