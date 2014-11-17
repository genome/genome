package Genome::Model::Tools::EpitopePrediction::BindingFilter;
# Based on myParseReportsTopGenesTopAlleleStringent_4.pl by TWYLIE (Feb 2012)

use strict;
use warnings;
use Genome;
use File::Basename;

class Genome::Model::Tools::EpitopePrediction::BindingFilter {
    is        => ['Genome::Model::Tools::EpitopePrediction::Base'],
    has_input => [
        fof_file => {
            is  => 'FilePath',
            doc => 'FOF containing list of parsed epitope files for different allele-length combinations (same sample)',
        },

],
    has_output => [
        output_file => {
            is => 'FilePath',
            doc => 'Output .xls file containing list of filtered epitopes based on binding affinity for each allele-length combination',
                }
        ],
};

sub help_brief {
    "Takes in a FOF with parsed NetMHC files for different allele-length combinations and outputs best candidates based on binding affinity."
}

sub execute {
        my $self = shift;
       

	my $fof_fh  = Genome::Sys->open_file_for_reading($self->fof_file);
	my $out_fh =  Genome::Sys->open_file_for_writing($self->output_file);

	my %prediction;
	my $threshold = 500 ;
	while (my $file = $fof_fh->getline) {
# [0] Gene Name
# [1] Point Mutation
# [2] Sub-peptide Position
# [3] MT score
# [4] WT score
# [5] MT epitope seq
# [6] WT epitope seq
# [7] Fold change
		chomp $file ;
		my $basename = basename( $file);
		my @f      = split( /\./, $basename );
		my $sample = $f[0];
		$sample =~ s/_netmhc//g;
		my $allele = $f[1];
		my $length = $f[2];

		my $mode = 'filtered';
		if ($file =~ /nofilter/i) { $mode = 'not filtered' }

		my $i = 0;

		my $parsed_fh  = Genome::Sys->open_file_for_reading($file);

		while (my $line = $parsed_fh->getline) {
			chomp $line;
			$i++;
			next if ($i == 1);  # skips headers....
				my (
						$gene_name,
						$point_mutation,
						$sub_peptide_mutation,
						$mt_score,
						$wt_score,
						$mt_epitope_seq,
						$wt_epitope_seq,
						$fold_change,
				   ) = split( "\t", $line );


			$gene_name = {
				gene_name            => $gene_name,
				allele               => $allele,
				point_mutation       => $point_mutation,
				sub_peptide_mutation => $sub_peptide_mutation,
				mt_score             => $mt_score,
				wt_score             => $wt_score,
				mt_epitope_seq       => $mt_epitope_seq,
				wt_epitope_seq       => $wt_epitope_seq,
				fold_change          => $fold_change,
			};
			push( @{ $prediction{$mode}->{$sample}->{$length}->{genes} }, \$gene_name );
		}
		close ($parsed_fh);

	}
	close ($fof_fh);

	my %best;
	foreach my $mode (sort keys %prediction) {
		if ($mode eq 'not filtered') { next ;}
		foreach my $sample (sort keys %{ $prediction{$mode} }) {
			foreach my $length (sort keys %{ $prediction{$mode}->{$sample} }) {
				foreach my $gene (sort @{ $prediction{$mode}->{$sample}->{$length}->{genes} }) {
# BEST
					unless( !$best{$sample}->{$$gene->{gene_name}}->{SCORE} ) {
						if ($$gene->{mt_score} < $best{$sample}->{$$gene->{gene_name}}->{SCORE}) {
							$best{$sample}->{$$gene->{gene_name}}->{SCORE} = $$gene->{mt_score};
							$best{$sample}->{$$gene->{gene_name}}->{GENES} = [];
							$$gene->{sample} = $sample;
							$$gene->{length} = $length;
							$$gene->{mode}   = $mode;
							push( @{ $best{$sample}->{$$gene->{gene_name}}->{GENES} }, $gene );
						}
						elsif ($$gene->{mt_score} == $best{$sample}->{$$gene->{gene_name}}->{SCORE}) {
							$best{$sample}->{$$gene->{gene_name}}->{SCORE} = $$gene->{mt_score};
							$$gene->{sample} = $sample;
							$$gene->{length} = $length;
							$$gene->{mode}   = $mode;
							push( @{ $best{$sample}->{$$gene->{gene_name}}->{GENES} }, $gene );
						}
					}
					else {
						$best{$sample}->{$$gene->{gene_name}}->{SCORE} = $$gene->{mt_score};
						$$gene->{sample} = $sample;
						$$gene->{length} = $length;
						$$gene->{mode}   = $mode;
						push( @{ $best{$sample}->{$$gene->{gene_name}}->{GENES} }, $gene );
					}

				}
			}
		}
	}

# REPORTING
	foreach my $sample (sort keys %best) {
		print $out_fh join(
				"\t",
				'Mode',
				'Sample',
				'Length',
				'Gene Name',
				'Allele',
				'Point Mutation',
				'Sub Peptide Position',
				'MT Score',
				'WT Score',
				'MT Epitope Seq',
				'WT Epitope Seq',
				'Fold Change',
				) . "\n";
		foreach my $gene (sort keys %{ $best{$sample} }) {
			foreach my $entry (@{ $best{$sample}->{$gene}->{GENES} }) {
				print $out_fh join(
						"\t",
						$$entry->{mode},
						$$entry->{sample},
						$$entry->{length},
						$$entry->{gene_name},
						$$entry->{allele},
						$$entry->{point_mutation},
						$$entry->{sub_peptide_mutation},
						$$entry->{mt_score},
						$$entry->{wt_score},
						$$entry->{mt_epitope_seq},
						$$entry->{wt_epitope_seq},
						$$entry->{fold_change},
						) . "\n" if ($$entry->{mt_score} < $threshold);
			}
		}
		close( $out_fh );
	}
	return 1;
}
__END__
