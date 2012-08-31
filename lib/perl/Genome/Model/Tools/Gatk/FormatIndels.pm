package Genome::Model::Tools::Gatk::FormatIndels;

#####################################################################################################################################
# FormatIndels - Convert GATK indel calls to a format that can be used with the WU annotator
#
#   AUTHOR:   Dan Koboldt (dkoboldt@genome.wustl.edu)
#
#   CREATED:  10/23/2009 by D.K.
#   MODIFIED: 06/05/2012 by C.K.
#
#####################################################################################################################################

use strict;
use warnings;
use IO::File;
use Genome;

class Genome::Model::Tools::Gatk::FormatIndels {
    is => 'Command',
    has => [
        variants_file => { is => 'Text', doc => "File of indel predictions", is_optional => 0, is_input => 1 },
        output_file => { is => 'Text', doc => "Output file to receive formatted lines", is_optional => 0, is_input => 1, is_output => 1 },
    ],
};

sub help_brief {
    "Formats indels for the annotation pipeline"
}

sub help_synopsis {
    return <<EOS
 gmt gatk format-indel-bed --variants-file GATK.out --output-file GATK.out.formatted
EOS
}

sub execute
{
    my $self = shift;

    # Grab arguments
    my $variants_file = $self->variants_file;
    my $output_file = $self->output_file;
    my @formatted = (); # This is a buffer to hold the output till we're ready to print to file

    # Parse the input file and convert each of them to an annotator-friendly format
    my $inFh = IO::File->new( $variants_file ) or die "Input file not found: $variants_file. $!";
    while( my $line = $inFh->getline ) {
        chomp( $line );
        my @cols = split( /\t/, $line );
        my $num_cols = scalar( @cols );
        my ( $chrom, $start, $stop, $indel, @rest ) = @cols;

        # Skip headers and any lines that are unusually short. There should be 10 columns for germline calls, 17 for somatic calls
        next if( $num_cols < 10 || $chrom =~ m/^(CHROM|REF)/ );

        # Fix unsupported chromosome names if necessary
        $chrom =~ s/chr//;
        $chrom = "MT" if( $chrom eq "M" );

        my $ref = my $var = "-";
        if( substr( $indel, 0, 1 ) eq '+' ) { # This is an insertion
            $stop = $start + 1;
            ( $var ) = $indel =~ m/^[+-](\w+)/;
        }
        else { # This is a deletion
            $start++;
            ( $ref ) = $indel =~ m/^[+-](\w+)/;
        }

        # Add the reformatted line to an array for printing to a file later
        push( @formatted, "$chrom\t$start\t$stop\t$ref\t$var\t" . join( "\t", @rest ) . "\n" );
    }
    $inFh->close;

    # Print the reformatted indels to the user-specified output file
    my $outFh = IO::File->new( $output_file, ">" ) or die "Cannot open $output_file. $!";
    $outFh->print( @formatted );
    $outFh->close;

    # Sort the formatted indels by loci
    my $cmd_obj = Genome::Model::Tools::Joinx::Sort->create(
        input_files => [ $output_file ],
        output_file => "$output_file.sorted",
    );
    $cmd_obj->execute;
    system( "mv -f $output_file.sorted $output_file" );

    return 1;
}

1;
