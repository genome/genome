#!/usr/bin/env genome-perl

use strict;
use warnings;

use Data::Dumper;
use Finfo::CommandLineOptions;
use Finishing::Assembly::Commands::AnalyzeJoin;
use Finishing::Assembly::Commands::CreateAssembly;
use Finishing::Assembly::Commands::CreateOrganism;
use Finishing::Assembly::Commands::ProjectCheckout;
use Finishing::Assembly::Commands::DeploySchema;
use Finishing::Assembly::Commands::ExportChromosomeAgp;
use Finishing::Assembly::Commands::ImportChromosomeAgp;
use Finishing::Assembly::Commands::JoinContigs;
use Finishing::Assembly::Commands::RemoveAndReplaceAce;

use Finfo::Logging 'fatal_msg';

require Term::ANSIColor;
require Text::Wrap;

my $function = shift @ARGV;
my @functions = _functions();
unless ( $function and $function !~ /^\-/ )
{
    $Text::Wrap::columns = 100;

    print Text::Wrap::wrap
    (
        " ", # 1st line indent, 
        "  ", # all other lines indent, 
        Term::ANSIColor::colored("Usage\n", 'bold'),
        "\tfin-assem-admin <function> <options>\n\n\n",
    ),
    Text::Wrap::wrap
    (
        " ",
        "  ",
        Term::ANSIColor::colored("For specific function help\n", 'bold'),
        "\tfin-assem-admin <function> --help\n\n\n",
    ),
    Text::Wrap::wrap
    (
        " ",
        "  ",
        Term::ANSIColor::colored("Valid functions\n", 'bold'),
        join(', ', @functions),
        "\n",
    );

    exit 1;
}

Finfo::Validate->validate
(
    attr => 'function',
    value => $function,
    isa => [ 'in_list', @functions ],
    msg => 'fatal',
);

my $class = _function_to_class($function);
my $clo = Finfo::CommandLineOptions->new
(
    classes => [ $class ],
    add_q => 1,
);
my $opts = $clo->get_options;
my $command = $class->new( %{ $opts->{$class} } );
$command->execute;

exit 0;

###########################################################

sub _classes
{
    return 
    (qw/
        AnalyzeJoin
        ProjectCheckout
        CreateOrganism
        CreateAssembly
        DeploySchema
        ExportChromosomeAgp
        ImportChromosomeAgp
        JoinContigs
        RemoveAndReplaceAce
        /);
}

sub _function_to_class
{
    my ($function) = @_;

    return 'Finishing::Assembly::Commands::' . join('', map { ucfirst } split(/\-/, $function));
}

sub _functions
{
    my @functions;
    foreach my $class ( _classes() )
    {
        my @words = $class =~ /([A-Z](?:[A-Z]*(?=$|[A-Z][a-z])|[a-z]*))/g;
        push @functions, join('-', map { lc } @words);
    }

    return @functions;
}

sub search
{
    my ($factory, $opts) = @_;

    if ( my $query = $opts->{search}->{query} )
    {
        main->fatal_msg("Not a search query: $query") if $query =~ /delete|drop|insert/i;
        my $dbh = $factory->storage->dbh;
        my $sth = $dbh->prepare($query);
        main->fatal_msg( $DBI::errstr ) unless $sth;
        $sth->execute
            or main->fatal_msg( $DBI::errstr );

        my $columns = $sth->{NAME_uc};
        my $format = _get_column_format($columns);

        print sprintf($format, @$columns);

        foreach my $aryref ( @{ $sth->fetchall_arrayref } )
        {
            print sprintf($format, map { $_ || 'NULL' } @$aryref);
        }
    }
    elsif ( my $table = delete $opts->{search}->{table} )
    {
        $table = join('', map { ucfirst } split(/_/, $table));
        my @search_params;
        my $rs = $factory->resultset(ucfirst $table)->search
        (
        );
        my @columns = $rs->result_source->columns;
        my $format = _get_column_format(\@columns);
        print sprintf($format, @columns);
        while ( my $obj = $rs->next )
        {
            print sprintf($format, map { $obj->$_ || 'NULL' } @columns);
        }
    }

    return 1;
}

sub _get_column_format
{
    my ($columns) = @_;

    mai->fatal_msg("Need columns") unless $columns;
    
    my $format;
    foreach my $name ( @$columns )
    {
        $format .= sprintf('%%-%ds', ( $name =~ /id/i and length($name) < 8 ) ? 7 : 20);
    }
    $format .= "\n";

    return $format;
}

