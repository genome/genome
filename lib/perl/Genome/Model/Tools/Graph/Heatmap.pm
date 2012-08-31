
# Rename the final word in the full class name <---
package Genome::Model::Tools::Graph::Heatmap;

use strict;
use warnings;

use Genome;
use Command;
use Carp;
use IO::File;
use File::Slurp;
use Text::CSV_XS;

UR::Object::Type->define(
    class_name => __PACKAGE__,
    is         => 'Command',
    has        => [
        'matrix' => {
            is  => 'String',
            doc => "mutation matrix file",
        },
        'image' => {
            is  => 'String',
            doc => "image output file (png)",
        },
        'columns' => {
            is  => 'Integer',
            doc => "number of columns to use from matrix",
        },
        'sepchar' => {
            is          => 'String',
            doc         => "separator character in matrix",
            is_optional => 1,
            default     => "\t"
        },

    ],
);

sub help_brief
{
    "tool to generate heatmaps of mutation matrix files via R";
}

sub execute
{
    my $self = shift;

    #    unless($self->sepchar)
    #    {
    #        $self->sepchar("\t");
    #    }
    my $maxmutations = $self->validate_matrix( $self->matrix );
    unless ($maxmutations)
    {
        croak "matrix did not validate, check format";
    }

    # this needs to be nailed down better.
    my @colors = qw/blue pink red yellow orange purple/;

    my $color_vector =
      $self->array2vector( \@colors, { limit => $maxmutations } );
    my $matrix  = $self->matrix;
    my $image   = $self->image;
    my $columns = $self->columns;
    my $rcmd    = [
        qq(
                    x <- read.table( "$matrix" ,header=T )
                    mt <- as.matrix(x[,1:$columns])
                    png("$image", width=1600, height=1200)
                    tmt <- t(mt)
                    image(c(1:length(rownames(tmt))),
                          c(1:length(colnames(tmt))),
                          tmt,axes=FALSE, col = $color_vector,
                          xlab="",ylab="")
                    axis(side=1,at=c(1:length(rownames(tmt))),
                         labels=rownames(t(mt)),
                         las=2)
                    dev.off()
                    ),
    ];

    my ( $fh, $tmpfile ) = Genome::Sys->create_temp_file('gmt_graph_heatmap.R');

    write_file( $tmpfile, @$rcmd );

    system("/usr/bin/R --no-save -q < $tmpfile >/dev/null 2>&1");

    if ($@)
    {
        print STDERR "problem with R?";
        print STDERR $@, "\n";
        exit 2;
    }

    unlink($tmpfile);

    return 1;
}

sub validate_matrix
{
    my $self  = shift;
    my @lines = read_file( $self->matrix );
    chomp @lines;
    my $c = new Text::CSV_XS( { sep_char => $self->sepchar } );

    my $max         = 0;
    my $header_line = undef;

    foreach my $line (@lines)
    {
        $c->parse($line);
        my @f = $c->fields();

        if ( $f[0] )
        {
            foreach my $idx ( 1 .. $#f )
            {
                if ( $f[$idx] > $max )
                {
                    $max = $f[$idx];
                }
            }
        }
        else
        {

            # header line
            $header_line = $line;
        }
    }

    unless ( defined($header_line) )
    {
        return undef;
    }

    # what to return, $max?
    return $max;
}

sub array2vector
{
    my $self = shift;
    my ( $aref, $optref ) = @_;
    my $str = "c(";

    foreach my $i ( 0 .. $optref->{limit} )
    {
        $str = $str . qq("$aref->[$i]",);
    }
    $str =~ s/,$//;
    $str = $str . ")";
    return $str;
}

1;

# $Id$
