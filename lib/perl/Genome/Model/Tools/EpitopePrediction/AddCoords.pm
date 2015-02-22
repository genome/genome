package Genome::Model::Tools::EpitopePrediction::AddCoords;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::EpitopePrediction::AddCoords {
    is        => ['Genome::Model::Tools::EpitopePrediction::Base'],
    has_input => [
        epitope_file => {
            is  => 'FilePath',
            doc => 'Consolidated and filtered (based on binding affinities) epitope file from all lengths & alleles',
        },
	somatic_variation_build => {
            is => 'Genome::Model::Build::SomaticVariation',
            is_optional => 1,
            doc => 'The somatic variation build to use for analysis',
        },
        input_tsv_file => {
            is => 'FilePath',
            is_optional => 1,
            doc => 'The custom input tsv file to use for analysis if no somatic variation build is used',
        },

],
    has_output => [
        coords_file => {
            is => 'FilePath',
            doc => 'File to write the output with genomic coordinates',
		}    
	],
};

sub help_brief {
    "Takes in the candidate epitope file as well as the SNVs file and outputs a file with genomic coordinates for all epitopes."
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
			$self->warning_message("Somatic variation build given. Setting input_tsv_file to $tsv_file");
			$self->input_tsv_file($tsv_file);
		}
	}


	my $epitope_fh  = Genome::Sys->open_file_for_reading($self->epitope_file);
	my $var_fh  = Genome::Sys->open_file_for_reading($self->input_tsv_file);
	my $coords_fh =  Genome::Sys->open_file_for_writing($self->coords_file);

	my %lookup;

	while (my $var_line = $var_fh->getline) {
		chomp $var_line;
		my @f = split( "\t", $var_line );
		if ($f[13] eq 'missense') {
			my $key = $f[6].'_'.$f[15];
			$key =~ s/p\.//g;

			my %index = (
					chr        =>  0,
					start      =>  1,
					stop       =>  2,
					reference  =>  3,
					var        =>  4,
					gene       =>  6,
					transcript =>  7,
					type       => 13,
					mutation   => 15,
					ensgene    => 23,
				    );
			for my $col (keys %index) {
				my $i = $index{$col};
				$lookup{$key}->{$col} = $f[$i];
			}
		}
	}


	$self->print_header($coords_fh);

	while (my $candidates_line = $epitope_fh->getline) {
		chomp $candidates_line;
		if ($candidates_line !~ /^Mode/)
		{
			my @f   = split( "\t", $candidates_line );
			my $key = $f[3].'_'.$f[5];
			if ($lookup{$key}) {
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

				$self->warning_message('[could not resolve: %s]', $key)
			}
		}
	}

	close($epitope_fh);
	close ($var_fh);
	return 1;

}
sub print_header {
	my $self      = shift;
	my $output_fh = shift;

	print $output_fh join("\t",
			"Chr",
			"Start",
			"Stop",
			"Ref",
			"Var",
			"Mode",
			"Sample",
			"Length",
			"Gene Name",
			"Point Mutation",
			"Sub-peptide Position",
			"MT score",
			"WT score",
			"MT epitope seq",
			"WT epitope seq",
			"Fold change",
			"Transcript Name",
			"Ensembl Gene ID"
			) . "\n";
}
1;
