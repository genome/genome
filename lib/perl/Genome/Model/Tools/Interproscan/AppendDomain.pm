package Genome::Model::Tools::Interproscan::AppendDomain;

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
            is  => 'String',
            doc => "mutation annotation file",
        },
        'filter' => {
            is  => 'String',
            doc => "filter expression for domains",
        },
        'output' => {
            is  => 'String',
            doc => "output file name",
        },
        'gcol' => {
            is          => 'Integer',
            doc         => "hugo gene name column",
            default     => 0,
            is_optional => 1
        },
        'tcol' => {
            is          => 'Integer',
            doc         => "transcript name column",
            default     => 3,
            is_optional => 1
        },
        'acol' => {
            is          => 'Integer',
            doc         => "Amino acid change column",
            default     => 12,
            is_optional => 1
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
"tool for appending interproscan domain results to mutation annotation file";
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt append-domain --maf maf-file --output outfile \
    [--sepchar separator] [--tcol transcriptcolumn] \
    [--acol aminoacidchangecolumn] [--gcol genenamecolumn] \
    [--filter filterregex]
EOS
}

sub help_detail {                           
    return <<EOS 
This appends inerproscan domain results to mutations that fall into protein domains
EOS
}


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

    # read in mutations
    $snpdomains->read_muts( $tmpfile, "1,2" );
    $muthash = $snpdomains->get_all_mutations();
    $snpdomains->mutation_in_dom( \%pfamlen, $self->filter );
    $mutcounts = $snpdomains->get_mutation_counts();

    my @mutrecs = read_file( $self->maf )
      or croak "can't read " . $self->maf . " : $!";
    chomp(@mutrecs);

    my $c = new Text::CSV_XS( { sep_char => $self->sepchar } );
    my %hash;

    # start appending to the records.
    foreach my $line (@mutrecs)
    {
        $c->parse($line);
        my @f        = $c->fields();
        my $gene     = $f[ $self->gcol ];
        my $tname    = $f[ $self->tcol ];
        my $aachange = $f[ $self->acol ];
        push( @{ $hash{$gene}{$tname} }, $aachange );
    }

    my @new;

    foreach my $l (@mutrecs)
    {
        $c->parse($l);
        my @f        = $c->fields();
        my $gene     = $f[ $self->gcol ];
        my $tname    = $f[ $self->tcol ];
        my $aachange = $f[ $self->acol ];
        my $obj      = $snpdomains->get_mut_obj( $tname . "," . $gene );
        if ( defined($obj) )
        {
            my $doms = $obj->get_domain($aachange);
            if ( defined($doms) )
            {
                $f[ $#f + 1 ] = join( ":", uniq @$doms );
            }
        }
        $c->combine(@f);
        push( @new, $c->string() . "\n" );
    }

    write_file( $self->output, @new );
    unlink($tmpfile);

    return 1;
}

sub extract_snp_data2tmpfile
{
    my $self = shift;
    my ( $fh, $tmpfile ) =
      tempfile( "gt-append-domainXXXXXX", SUFFIX => '.dat' );
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

1;

# $Id$
