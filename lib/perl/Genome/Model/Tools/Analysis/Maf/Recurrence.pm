package Genome::Model::Tools::Analysis::Maf::Recurrence;

use strict;
use warnings;

use FileHandle;

use Genome;
use Genome::Model::Tools::Analysis::Maf::Helpers qw/getMutationSamples byGeneTranscript load_aa_changes load_annotation/;

my %stats = ();
my $max_proximity = 0;

class Genome::Model::Tools::Analysis::Maf::Recurrence {
	is => 'Command',                       
	
	has => [
		maf_file	=> { is => 'Text', doc => "Original MAF file" },
		annotation_file	=> { is => 'Text', doc => "Full annotation for variants in MAF file" },
		output_file	=> { is => 'Text', doc => "Output file for recurrence report", is_optional => 1 },
		output_maf	=> { is => 'Text', doc => "Output file of complete MAF with recurrent sites only", is_optional => 1 },		
		max_proximity	=> { is => 'Text', doc => "Maximum aa distance between mutations [0]", is_optional => 1 },
		verbose		=> { is => 'Text', doc => "Print verbose output", is_optional => 1 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {
    "Performs a recurrency analysis on mutations in a MAF file"                 
}

sub help_synopsis {
    return <<EOS
This command performs a proximity analysis on mutations in a MAF file
EXAMPLE:	gt analysis maf proximity --maf-file original.maf --output-file proximity-genes.tsv
EOS
}

sub help_detail {
    return <<EOS 

EOS
}

sub execute {
	my $self = shift;

	## Get required parameters ##
	my $maf_file = $self->maf_file;
	$max_proximity = $self->max_proximity if(defined($self->max_proximity));

	if(!(-e $maf_file))
	{
		die "Error: MAF file not found!\n";
	}

	if($self->output_file)
	{
		open(OUTFILE, ">" . $self->output_file) or die "Can't open outfile: $!\n";
	}

	if($self->output_maf)
	{
		open(OUTMAF, ">" . $self->output_maf) or die "Can't open outfile: $!\n";
	}

	
	my %mutated_aa_positions = load_aa_changes($maf_file, $self->annotation_file, \%stats);

	$stats{'num_aa_positions'} = $stats{'num_recurrent_genes'} = $stats{'num_recurrent_positions'} = 0;
	
	my %gene_counted = ();
	
	my %gene_aa_pos_recurrent_list = ();
	

	foreach my $aa_key (sort byGeneTranscript keys %mutated_aa_positions)
	{
		$stats{'num_aa_positions'}++;
		
		my @aa_lines = split(/\n/, $mutated_aa_positions{$aa_key});
		my $num_lines = @aa_lines;
		
		if($num_lines > 1)	
		{
			(my $gene, my $transcript_name, my $aa_position) = split(/\t/, $aa_key);


			## Parse the lines ##
			my %sample_counted = ();
			my $num_unique_samples = 0;
			
			my $sample_list = getMutationSamples($mutated_aa_positions{$aa_key});
			my @sample_list = split(/\n/, $sample_list);
			$num_unique_samples = @sample_list;
			
			my $sample_lines = "";
			my $sample_maf_lines = "";

			my $this_recurrent_list = "";

			foreach my $maf_line (@aa_lines)
			{
				my @lineContents = split(/\t/, $maf_line);
				my $tumor_sample = $lineContents[15];
				
				if(!$sample_counted{$tumor_sample})
				{
					$this_recurrent_list .= $maf_line . "\n";

#					$num_unique_samples++;
					my $ref = $lineContents[10];
					my $var = $lineContents[11];
					$var = $lineContents[12] if($var eq $ref);
					my $sample_line = "";

					$sample_line = join("\t", $gene, $transcript_name, $aa_position, $lineContents[2], $lineContents[4], $lineContents[5], $lineContents[6], $ref, $var, $lineContents[15]);						

					$sample_maf_lines .= join("\t", $gene, $transcript_name, $aa_position) . "\t" . $maf_line . "\n";
					$sample_lines .= $sample_line . "\n";
					$sample_counted{$tumor_sample} = 1;
				}
			}

			if($num_unique_samples > 1)
			{
				## Check to make sure that this exact recurrent list wasn't already printed ##
				
				if($gene_aa_pos_recurrent_list{"$gene\t$aa_position"} && $gene_aa_pos_recurrent_list{"$gene\t$aa_position"} eq $this_recurrent_list)
				{
					warn "NOT reporting duplicate find for $gene $aa_position\n" if($self->verbose);
				}
				else
				{
					$stats{'num_recurrent_positions'}++;
					print "$gene\t$transcript_name\t$aa_position\t$num_unique_samples\n" if($self->verbose);
					
					if($self->output_file)
					{
						print OUTFILE "$sample_lines";
					}
	
					if($self->output_maf)
					{
						print OUTMAF "$sample_maf_lines";
					}
					
					if(!$gene_counted{$gene})
					{
						$stats{'num_recurrent_genes'}++;
						$gene_counted{$gene} = 1;
					}
					
					$gene_aa_pos_recurrent_list{"$gene\t$aa_position"} = $this_recurrent_list;
				}
				

			}

		}
	}


	close(OUTFILE) if($self->output_file);
	close(OUTMAF) if($self->output_maf);

	print $stats{'num_aa_positions'} . " different amino acid positions\n";
	print $stats{'num_recurrent_positions'} . " recurrent positions\n";
	print $stats{'num_recurrent_genes'} . " recurrent genes\n";	
}

1;
