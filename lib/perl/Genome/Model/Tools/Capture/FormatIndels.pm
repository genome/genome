package Genome::Model::Tools::Capture::FormatIndels;

#####################################################################################################################################
# FormatIndels - Convert VarScan indel calls to a format that can be used with the WU annotator
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

class Genome::Model::Tools::Capture::FormatIndels {
    is => 'Command',
    has => [
        variants_file => { is => 'Text', doc => "File of indel predictions", is_optional => 0, is_input => 1 },
        output_file => { is => 'Text', doc => "Output file to receive formatted lines", is_optional => 0, is_input => 1, is_output => 1 },
        append_line => { is => 'Boolean', doc => "If set to 1, appends the additional columns to output file", is_optional => 1, is_input => 1, default => 1 },
    ],
};

sub help_brief {
    "Formats indels for the annotation pipeline"
}

sub help_synopsis {
    return <<EOS
 gmt capture format-indels --variants-file varScan.out --output-file varScan.out.formatted
EOS
}

sub execute {
    my $self = shift;

    # Grab arguments
    my $variants_file = $self->variants_file;
    my $output_file = $self->output_file;
    my @formatted = (); # This is a buffer to hold the output till we're ready to print to file

    # Parse the indels
    my $inFh = IO::File->new( $variants_file ) or die "Input file not found: $variants_file. $!";
    while( my $line = $inFh->getline ) {
        chomp( $line );
        my ( $chrom, $start, $refbase, $indel, @rest ) = split( /\t/, $line );

        next if( $chrom =~ m/^(CHROM|REF)/ );
        # Fix unsupported chromosome names if necessary
        $chrom =~ s/^chr//;
        $chrom = "MT" if( $chrom eq "M" );

        my $allele1 = my $allele2 = "";
        my $stop = 0;
        my $ref = my $var = "";
        my $indel_type = my $indel_size = "";

        if( $refbase =~ /[0-9]/ ) {
            $stop = $refbase;
            if($indel =~ '/') {
                my @tempArray = split(/\//, $indel);
                $var = $tempArray[0];
                $var = $tempArray[1] if($tempArray[1] ne '*');

                if(substr($var, 0, 1) eq '+') {
                    $ref = "-";
                    $var =~ s/[^ACGTN]//g;
                }
                elsif(substr($var, 0, 1) eq '-') {
                    $ref = $var;
                    $ref =~ s/[^ACGTN]//g;
                    $var = "-";
                }
                elsif($tempArray[0] eq '*') {
                    $ref = "-";
                    $var = $tempArray[1];
                }
                elsif($tempArray[1] eq '*') {
                    $ref = $tempArray[0];
                    $var = "-";
                }
            }
            else {
                $ref = $indel;
                $var = $rest[0];
            }
        }
        else {
            if( $indel =~ '/' ) {
                my @tempArray = split(/\//, $indel);
                $var = $tempArray[0];
                $var = $tempArray[1] if($tempArray[1] ne '*');

                if(substr($var, 0, 1) eq '+') {
                    $ref = "-";
                    $var =~ s/[^ACGTN]//g;
                }
                elsif(substr($var, 0, 1) eq '-') {
                    $ref = $var;
                    $ref =~ s/[^ACGTN]//g;
                    $var = "-";
                }
            }
            elsif( $indel =~ 'INS' || $indel =~ 'DEL' ) {
                my @indelContents = split(/\-/, $indel);
                if( $indel =~ 'INS' ) {
                    $ref = "-";
                    $var = $indelContents[2];
                }
                else {
                    $ref = $indelContents[2];
                    $var = "-";
                }
            }
            else {
                $ref = $refbase;
                $var = $indel;
            }
        }

        # Correct alleles
        if( $ref eq '-' || $var eq '-' ) {
            $allele1 = $ref;
            $allele2 = $var;

            if($ref eq '-') {
                $indel_type = "INSERTION";
                $indel_size = length($var);
            }
            else {
                $indel_type = "DELETION";
                $indel_size = length($ref);
            }
        }
        elsif(substr($var, 0, 1) eq '+') {
            $allele1 = "-";
            $allele2 = uc($var);
            $allele2 =~ s/[^ACGTN]//g;
            $indel_type = "INSERTION";
            $indel_size = length($allele2);
        }
        elsif(substr($var, 0, 1) eq '-') {
            $allele2 = "-";
            $allele1 = uc($var);
            $allele1 =~ s/[^ACGTN]//g;
            $indel_type = "DELETION";
            $indel_size = length($allele1);
        }
        else {
            warn "Unable to format $line\n";
        }

        # If no chr stop, calculate it
        if(!$stop) {
            if($indel_type eq "INSERTION") {
                $stop = $start + 1;
            }
            else {
                $start++;
                $stop = $start + $indel_size - 1;
            }
        }

		# Ensure that insertions are zero-based
	if($indel_type eq "INSERTION")
	{
		$start-- if($start == $stop);		
	}


        # Construct the reformatted line and buffer it up for output
        $line = "$chrom\t$start\t$stop\t$allele1\t$allele2";
        $line .= "\t" . join( "\t", @rest ) if( $self->append_line );
        $line .= "\n";
        push( @formatted, $line );
    }
    $inFh->close;

    # Print the reformatted indels to the user-specified output file
    my $outFh = IO::File->new( $output_file, ">" ) or die "Cannot open $output_file. $!";
    $outFh->print( @formatted );
    $outFh->close;

    # Sort the re-formatted indels by loci
    my $cmd_obj = Genome::Model::Tools::Joinx::Sort->create(
        input_files => [ $output_file ],
        output_file => "$output_file.sorted",
    );
    $cmd_obj->execute;
    system( "mv -f $output_file.sorted $output_file" );

    return 1;
}

1;
