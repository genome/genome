package Genome::Model::Tools::EpitopePrediction::AddCoords;

use strict;
use warnings;
use Data::Dumper;
use Genome;
use File::Basename qw/fileparse/;

class Genome::Model::Tools::EpitopePrediction::AddCoords {
    is        => ['Genome::Model::Tools::EpitopePrediction::Base'],
    has_input => [
        parsed_file => {
            is  => 'Text',
            doc => 'Parsed output file from Netmhc',
        },
	somatic_variation_build => {
            is => 'Genome::Model::Build::SomaticVariation',
            is_optional => 1,
            doc => 'The somatic variation build to use for analysis',
        },
        input_tsv_file => {
            is => 'Text',
            is_optional => 1,
            doc => 'The custom input tsv file to use for analysis if no somatic variation build is used',
        },

],
    has_output => [
        coords_file => {
            is => 'Text',
            doc => 'File to write the parsed output',
		}    
	],
};

sub help_brief {
    "Takes in the \"Parsed output file from NetMHC as well as the SNVs file and output a file with genomic coordinates for all epitopes."
}

sub execute {
    my $self = shift;
    if (!defined($self->somatic_variation_build) && !defined($self->input_tsv_file)) {
        die $self->error_message("Either somatic variation build or input tsv file needs to be provided");
    }

    if (defined($self->somatic_variation_build)) {
        if (defined($self->input_tsv_file)) {
            die $self->error_message("Custom tsv file cannot be used in combination with somatic variation build");
        }
        else {
            my $tsv_file = File::Spec->join(
                $self->somatic_variation_build->data_directory,
                'effects',
                'snvs.hq.tier1.v1.annotated.top.header'
            );
            $self->status_message("Somatic variation build given. Setting input_tsv_file to $tsv_file");
            $self->input_tsv_file($tsv_file);
        }
		}


	my $parsed_fh  = Genome::Sys->open_file_for_reading($self->parsed_file);
	my $var_fh  = Genome::Sys->open_file_for_reading($self->input_tsv_file);
	my $coords_fh =  Genome::Sys->open_file_for_writing($self->coords_file);

	my %lookup;

	while (my $var_line = $var_fh->getline) {
        	chomp $var_line;
    		my @f = split( "\t", $var_line );
    		if ($f[13] eq 'missense') {
        	my $key = $f[15];
        	$key =~ s/p\.//g;
        	$lookup{$key}->{chr}      = $f[0];
        	$lookup{$key}->{start}    = $f[1];
        	$lookup{$key}->{stop}     = $f[2];
        	$lookup{$key}->{gene}     = $f[6];
        	$lookup{$key}->{mutation} = $f[15];
        	$lookup{$key}->{type}     = $f[13];
    		$lookup{$key}->{transcript} = $f[7];
		$lookup{$key}->{reference} = $f[3];
		$lookup{$key}->{var}=	$f[4];
		$lookup{$key}->{ensgene}= $f[23];

			}
		}
	 while (my $candidates_line = $parsed_fh->getline) {
		chomp $candidates_line;
		if ($candidates_line !~ /^Mode/)
		{
    			my @f   = split( "\t", $candidates_line );
    			my $key = $f[5];
    			if ($lookup{$key}) {
			print $candidates_line ."\n";	        
		print $coords_fh join(
                   	"\t",
        	           $lookup{$key}->{chr},
                	   $lookup{$key}->{start},
				$lookup{$key}->{stop},
		   $lookup{$key}->{reference},
		   $lookup{$key}->{var},
		   $candidates_line,
		   $lookup{$key}->{transcript},
		 $lookup{$key}->{ensgene}
			) . "\n";
   			 }
    				else {
       					 print "[could not resolve: $key]\n";
   				 }
		}
	}
			
close($parsed_fh);
close ($var_fh);
return 1;

}

1;
