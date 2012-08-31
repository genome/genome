package Genome::Model::Tools::Capture::FormatVcfSnps;

#####################################################################################################################################
# FormatVcfSnps - Convert VarScan indel calls to a format that can be used with the WU annotator
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

class Genome::Model::Tools::Capture::FormatVcfSnps {
    is => 'Command',
    has => [
        vcf_file => { is => 'Text', doc => "VCF file of indel predictions", is_optional => 0, is_input => 1 },
        output_file => { is => 'Text', doc => "Output file to receive formatted lines", is_optional => 0, is_input => 1, is_output => 1 },
        append_line => { is => 'Boolean', doc => "If set to 1, appends the additional columns to output file", is_optional => 1, is_input => 1, default => 1 },
    ],
};

sub help_brief {
    "Formats SNPs in a VCF file for the annotation pipeline"
}

sub help_synopsis {
    return <<EOS
 gmt capture format-vcf-snps --vcf-file varScan.out --output-file varScan.out.formatted
EOS
}

sub execute {
    my $self = shift;

    # Grab arguments
    my $vcf_file = $self->vcf_file;
    my $output_file = $self->output_file;
    my @formatted = (); # This is a buffer to hold the output till we're ready to print to file

    # Parse the indels
    my $inFh = IO::File->new( $vcf_file ) or die "Input file not found: $vcf_file. $!";
    while( my $line = $inFh->getline ) {
        chomp( $line );

	if(substr($line, 0, 1) eq '#')
	{
		## Skip VCF header ##
	}
	else
	{
		my ( $chrom, $position, $id, $ref, $alt, @rest ) = split( /\t/, $line );
	
		next if( $chrom =~ m/^(CHROM|REF)/ );
		# Fix unsupported chromosome names if necessary
		$chrom =~ s/^chr//;
		$chrom = "MT" if( $chrom eq "M" );
	
		my $chr_start = my $chr_stop = 0;
		
		my @vars = split(/\,/, $alt);
		foreach my $var (@vars)
		{
			# Construct the reformatted line and buffer it up for output
			$line = "$chrom\t$position\t$position\t$ref\t$var\n";
			push( @formatted, $line );					
		}


	}

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
