
package Genome::Model::Tools::Analysis::Indels::BuildContigs;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# SearchRuns - Search the database for runs
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	04/01/2009 by D.K.
#	MODIFIED:	04/01/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;
use POSIX;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Analysis::Indels::BuildContigs {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		variant_file	=> { is => 'Text', doc => "Indels in annotation format" },
		contig_size	=> { is => 'Text', doc => "Size of reference/variant contigs to generate", is_optional => 1, default => 150 },
		reference	=> { is => 'Text', doc => "Size of reference/variant contigs to generate", is_optional => 1, default => "/gscmnt/839/info/medseq/reference_sequences/NCBI-human-build36/all_sequences.fasta" },
		output_file	=> { is => 'Text', doc => "Output of reference/variant contig FASTAs", is_optional => 1 },
        output_reference => { is => 'Boolean', doc => 'Whether or not to output the reference contigs', is_optional => 1, default => 1, },
        samtools_compatible => { is => 'Boolean', doc => 'Whether or not to output the contig names with underscores to make them compatible with samtools', is_optional => 1, default => 0 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Build reference/variant contigs for indel validation"                 
}

sub help_synopsis {
    return <<EOS
This command builds reference/variant contigs for indel validation
EXAMPLE:	gmt analysis indels build-contigs --variant-file [indels.formatted.tsv] --output-file [indels.contigs.fasta]
EOS
}

sub help_detail {                           # this is what the user will see with the longer version of help. <---
    return <<EOS 

EOS
}


################################################################################################
# Execute - the main program logic
#
################################################################################################

sub execute {                               # replace with real execution logic.
	my $self = shift;

    # Check that we're on a 64-bit system and can run with the deployed samtools
    unless (POSIX::uname =~ /64/) {
        $self->error_message("Must run on a 64 bit machine");
        die;
    }

	## Get required parameters ##
	my $variant_file = $self->variant_file;
	my $contig_size = $self->contig_size;
	my $reference = $self->reference;
	my $output_file = $self->output_file;
    my $samtools_compatible = $self->samtools_compatible;

    my %indel_names;

	my %stats = ();
	$stats{'num_indels'} = 0;
    $stats{'num_dups'} = 0;

	if($output_file)
	{
		open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";
	}
	
	my $input = new FileHandle ($variant_file);
	my $lineCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		

		(my $chrom, my $chr_start, my $chr_stop, my $ref, my $var) = split(/\t/, $line);
		$chrom =~ s/[^0-9XYMNT\_]//g;
		$chr_start =~ s/[^0-9]//g if($chr_start);
		$chr_stop =~ s/[^0-9]//g if($chr_stop);
	
		if($chrom && $chr_start && $chr_stop)
		{
			$stats{'num_indels'}++;
			my $indel_type = my $indel_size = my $allele = "";
			
			if($ref eq "0" || $ref eq "-")
			{
				$indel_type = "Ins";
				$indel_size = length($var);
                $allele = uc($var);

				## Build indel name ##			
				my $indel_name = "$chrom:$chr_start-$chr_stop:$indel_type:$indel_size:$allele";
                if($samtools_compatible) {
                    $indel_name =~ s/[:-]/_/g;  #replace dashes and colons with underscores
                }

                if(exists($indel_names{$indel_name})) {
                    $self->error_message("Skipping indel with duplicate $indel_name");
                    $stats{'num_dups'}++;
                    next;
                }
                else {
                    $indel_names{$indel_name} = 1;
                }
				
				my $flank_size = ($contig_size) / 2;
				$flank_size = sprintf("%d", $flank_size);
				my $flank5_start = $chr_start - $flank_size + 1;
				my $flank5_stop = $chr_start;
				my $flank3_start = $chr_stop;
				my $flank3_stop = $chr_stop + $flank_size - 1;
						
				## Get upstream flank ##

				my $flank5_seq = `samtools faidx $reference $chrom:$flank5_start-$flank5_stop | grep -v \">\"`;
				$flank5_seq =~ s/[^A-Za-z]//g if($flank5_seq);

				## Get downstream flank ##

				my $flank3_seq = `samtools faidx $reference $chrom:$flank3_start-$flank3_stop | grep -v \">\"`;
				$flank3_seq =~ s/[^A-Za-z]//g if($flank3_seq);

				## Here's a sequence check of the un-interrupted reference, if needed ##
				my $check_seq = `samtools faidx $reference $chrom:$flank5_start-$flank3_stop | grep -v \">\"`;
				$check_seq =~ s/[^A-Za-z]//g if($check_seq);

				## Build Reference Contig ##
				my $reference_contig = $flank5_seq . $flank3_seq;

				## Build Variant Contig ##
				my $variant_contig = $flank5_seq . $var . $flank3_seq;

				## Trim variant contig length back to size of indel contig ##
				
				my $end = 0;
				while(length($variant_contig) > length($reference_contig))
				{
					## While reference contig is longer, trim from start and then from end ##
					if($end)
					{
						$variant_contig = substr($variant_contig, 0, length($variant_contig) - 1);
						$end = 0;
					}
					else
					{
						$variant_contig = substr($variant_contig, 1, length($variant_contig) - 1);
						$end = 1;
					}
				}
                if($self->output_reference) {
    				print OUTFILE ">" . $indel_name . "_ref\n";
	    			print OUTFILE $reference_contig . "\n";
                }

				print OUTFILE ">" . $indel_name . "_var\n";
				print OUTFILE $variant_contig . "\n";

#				print "$chrom\t$chr_start\t$chr_stop\t$ref\t$var\t$indel_type\t$flank_size\n";

			}
			else
			{
				$indel_type = "Del";
				$indel_size = length($ref);
                $allele = uc($ref);
				## Build indel name ##			
				my $indel_name = "$chrom:$chr_start-$chr_stop:$indel_type:$indel_size:$allele";

                if($samtools_compatible) {
                    $indel_name =~ s/[:-]/_/g;  #replace dashes and colons with underscores
                }

                if(exists($indel_names{$indel_name})) {
                    $self->error_message("Skipping indel with duplicate $indel_name");
                    $stats{'num_dups'}++;
                    next;
                }
                else {
                    $indel_names{$indel_name}=1;
                }
				
				## Fix chromosome stop position, which sometimes is the base after the deletion stops ##
				$chr_stop = $chr_start + $indel_size - 1;
				
#				my $flank_size = ($contig_size - $indel_size) / 2;
				my $flank_size = ($contig_size) / 2;
				$flank_size = sprintf("%d", $flank_size);
				my $flank5_start = $chr_start - $flank_size;
				my $flank5_stop = $chr_start - 1;
				my $flank3_start = $chr_stop + 1;
				my $flank3_stop = $chr_stop + $flank_size;
						
				## Get upstream flank ##

				my $flank5_seq = `samtools faidx $reference $chrom:$flank5_start-$flank5_stop | grep -v \">\"`;
				$flank5_seq =~ s/[^A-Za-z]//g if($flank5_seq);

				## Get downstream flank ##

				my $flank3_seq = `samtools faidx $reference $chrom:$flank3_start-$flank3_stop | grep -v \">\"`;
				$flank3_seq =~ s/[^A-Za-z]//g if($flank3_seq);


				## Here's a sequence check of the un-interrupted reference, if needed ##
#				my $check_seq = `samtools faidx $reference $chrom:$flank5_start-$flank3_stop | grep -v \">\"`;
#				$check_seq =~ s/[^A-Za-z]//g if($check_seq);

				## Build Reference Contig ##
				my $reference_contig = $flank5_seq . $ref . $flank3_seq;

				## Build Variant Contig ##
				
				my $variant_contig = $flank5_seq . $flank3_seq;

				## Trim reference contig length back to size of indel contig ##
				
				my $end = 0;
				while(length($reference_contig) > length($variant_contig))
				{
					## While reference contig is longer, trim from start and then from end ##
					if($end)
					{
						$reference_contig = substr($reference_contig, 0, length($reference_contig) - 1);
						$end = 0;
					}
					else
					{
						$reference_contig = substr($reference_contig, 1, length($reference_contig) - 1);
						$end = 1;
					}
				}

				## Print contigs to outfile
                if($self->output_reference) {
    				print OUTFILE ">" . $indel_name . "_ref\n";
	    			print OUTFILE $reference_contig . "\n";
                }

				print OUTFILE ">" . $indel_name . "_var\n";
				print OUTFILE $variant_contig . "\n";

			}
			
			
			print "$chrom\t$chr_start\t$chr_stop\t$ref\t$var\t$indel_type-$indel_size\n";
			
		}
	}
	
	close($input);

	print $stats{'num_indels'} . " indels in file\n";
	print $stats{'num_dups'} . " duplicated (and skipped) indels in file\n";
	print "Contigs printed to $output_file\n";
	close(OUTFILE) if($output_file);

	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}




1;

