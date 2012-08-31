
package Genome::Model::Tools::Capture::WiggleFromPileup;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# WiggleFromQ20 - Converts a pileup file to a three column file of chromosome, position, and bases with q>20.
#					
#	AUTHOR:		Will Schierding
#
#	CREATED:	06/24/2010 by W.S.
#	MODIFIED:	06/24/2010 by W.S.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;
use FileHandle;
use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Capture::WiggleFromPileup {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		normal_bam	=> { is => 'Text', doc => "Normal Bam File", is_optional => 0, is_input => 1 },
		tumor_bam	=> { is => 'Text', doc => "Tumor Bam File", is_optional => 0, is_input => 1 },
		reference_fasta	=> { is => 'Text', doc => "Reference Genome", is_optional => 1, is_input => 1, default =>"/gscmnt/839/info/medseq/reference_sequences/NCBI-human-build36/all_sequences.fa"},
		regions_file	=> { is => 'Text', doc => "Tab-delimited list of regions", is_optional => 0, is_input => 1 },
		min_depth_normal	=> { is => 'Text', doc => "Minimum Q20 depth for Normal [6]", is_optional => 1, is_input => 1, default => 6 },
		min_depth_tumor	=> { is => 'Text', doc => "Minimum Q20 depth for Tumor [8]", is_optional => 1, is_input => 1, default => 8 },
		output_dir     => { is => 'Text', doc => "Output directory to receive pileup, q20, and wiggle files", is_optional => 1, is_input => 1, is_output => 1 },
		output_file     => { is => 'Text', doc => "Output wiggle file (per-base qual>min coverage)", is_optional => 1, is_input => 1, is_output => 1, default =>"wigglefile.out" },
		gzip_after	=> { is => 'Text', doc => "If set to 1, compress the file after building", is_optional => 1, is_input => 1, default => 0 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Input bam files, builds a wiggle coverage file for a list of ROIs"
}

sub help_synopsis {
    return <<EOS
Builds a wiggle coverage file for a list of ROIs - Takes in two bam files, calls samtools pileup, creates Q20, then converts to wiggle
EXAMPLE:	gmt capture wiggle-from-pileup --normal-bam [normal bam file] --tumor-bam [tumor bam file] --regions-file [targets.tsv] --output-file [patient.wiggle]
TEST: gmt capture wiggle-from-pileup --normal-bam /gscmnt/sata870/info/model_data/2857976750/build102718369/alignments/102718369_merged_rmdup.bam --tumor-bam /gscmnt/sata860/info/model_data/2857976747/build102718323/alignments/102718323_merged_rmdup.bam --regions-file /gscmnt/sata836/info/medseq/TCGA-OV-Exome-Capture/targets/broad-exome.refseq.nochr.txt --output-file /gscmnt/sata836/info/medseq/TCGA-OV-Exome-Capture/TopUpComplete/wiggle_files/test.wiggle
EOS
}

sub help_detail {                           # this is what the user will see with the longer version of help. <---
    return <<EOS 
    Takes in two bam files, calls samtools pileup, creates Q20, then converts to wiggle coverage file for a list of ROIs                 
EOS
}


################################################################################################
# Execute - the main program logic
#
################################################################################################

sub execute {                               # replace with real execution logic.
	my $self = shift;
	#bam files
	my $normal_bam = $self->normal_bam;
	my $tumor_bam = $self->tumor_bam;
	# reference genome for samtools pileup
	my $bam_ref = $self->reference_fasta;
	my $output_dir = $self->output_dir;
	my $regions_file = $self->regions_file;
	my $output_file = $self->output_file;
	my $min_depth_normal = $self->min_depth_normal;
	my $min_depth_tumor = $self->min_depth_tumor;

	#set coverage at q20
	my $min_base_qual = 20;
	my $min_coverage = 0;

	unless (-e $bam_ref) {
		die "No reference fasta file\n";
	    }

#	my $q20_normal_file = $output_dir.'/q20_coverage_normal.txt';
#	my $q20_tumor_file = $output_dir.'/q20_coverage_tumor.txt';
#	# Open Output
#	unless (open(NORMAL_Q20,">$q20_normal_file")) {
#		die "Could not open output file '$q20_normal_file' for writing";
#	  }
#	# Open Output
#	unless (open(TUMOR_Q20,">$q20_tumor_file")) {
#		die "Could not open output file '$q20_tumor_file' for writing";
#	  }
	## Open outfile ##
	my $outfile = $output_dir."/".$output_file;
	open(OUTFILE, ">$outfile") or die "Can't open outfile: $!\n";

	print "Checking for bam index file...\n";

	my $normal_bam_index = $normal_bam.'.bai';
	my $tumor_bam_index = $tumor_bam.'.bai';
	unless (-e $normal_bam_index) {
		print "No normal bam index file\n";
		my $cmd_bam_bai = "samtools index $normal_bam";
		system($cmd_bam_bai);
	    }
	unless (-e $tumor_bam_index) {
		print "No tumor bam index file\n";
		my $cmd_bam_bai = "samtools index $tumor_bam";
		system($cmd_bam_bai);
	    }

	print "Bam indexing complete, checking bam file...\n";
	#get chromosome lengths
	my $tempdir = "/tmp/$output_file/";
	mkdir($tempdir);
	my $tmp_name = "$tempdir/tempbamheader.txt";
	my @command = `samtools view -H $normal_bam > $tmp_name`;
	my %chr_len;
	my $header = new FileHandle ($tmp_name);
	while (my $headline = <$header>) {
		$headline =~ s/\s+/\t/g;
		my @lineContents = split(/\t/, $headline);
		my $linetype = $lineContents[0];
		my $chrom = $lineContents[1];
		my $chr_length = $lineContents[2];
		unless ($linetype eq '@SQ') { next;};
		$chrom =~ s/SN://g;
		$chr_length =~ s/LN://g;
		$chr_len{$chrom} = $chr_length;
	}
	close($header);
	my $cmd = "rm $tmp_name";
	system($cmd);

	my $lineCounter_normal = 0;
	my $lineCounter_tumor = 0;
	my $out_normal = "$tempdir/temppileup_normal.txt";
	my $out_tumor = "$tempdir/temppileup_tumor.txt";
	my $pileupfh;
#	my %stats = ();
#	$stats{'bases'} = $stats{'covered'} = $stats{'not_covered'} = 0;
	my @chromes = (sort keys %chr_len);
	print "Starting pileup and wiggle file creation...\n";
	for my $chr (@chromes) {
		my $length = $chr_len{$chr};
		my $num_div = ($length / 25000000);
		my @reg_divides;
		my $count = 1;
		my $plusone = (int($num_div)+1);

		foreach my $math (1..$plusone) {
			my $end = ($math * 25000000);
			my $range;
			if ($math <= int($num_div)) {
				$range = "$count-$end";
			}
			elsif ($math == $plusone) {
				$range = "$count";
			}
			push (@reg_divides, $range);
			$count = ($end + 1);
		}
#		foreach my $poop (@reg_divides) {
#			print "$poop\n";
#		}
		for my $region (@reg_divides) {
			my %normal_coverage = ();
			my %tumor_coverage = ();

			print "Loading normal coverage...\n";
			#Call Samtools pileup for normal
			my $samcmd = "samtools view -uh $normal_bam $chr:$region | samtools pileup -v -f $bam_ref - >$out_normal";
			print "$samcmd\n";
			system ($samcmd);
			$pileupfh = new FileHandle ($out_normal);
		
			while (<$pileupfh>) {
				chomp;
				my $line = $_;
				$lineCounter_normal++;		
			
				my @lineContents = split(/\t/, $line);			
				my $chrom = $lineContents[0];
				my $position = $lineContents[1];
				my $ref_base = $lineContents[2];
				my $depth = $lineContents[3];
				my $qualities = $lineContents[5];
		
				## Go through each quality ##
				
				my @qualities = split(//, $qualities);
				my $num_quals = 0;
				my $qual_coverage = 0;
				
				foreach my $code (@qualities) {
					my $qual_score = ord($code) - 33;
					$num_quals++;
		
					if($qual_score >= $min_base_qual)
					{
						$qual_coverage++;
					}
				}
			
				if($qual_coverage >= $min_coverage) {
#					print NORMAL_Q20 "$chrom\t$position\t$qual_coverage\n";			
					my $wiggle_chr = $chrom;
					$wiggle_chr =~ s/chr//;
					$wiggle_chr = "MT" if($chrom eq "M");
					my $key = "$wiggle_chr\t$position";
					$normal_coverage{$key} = $qual_coverage;
				}
			}
			close($pileupfh);

			my $cmd = "rm $out_normal";
			system($cmd);
		
			print "Loading tumor coverage...\n";
		
			#Call Samtools pileup for tumor
			$samcmd = "samtools view -uh $tumor_bam $chr:$region | samtools pileup -v -f $bam_ref - >$out_tumor";
			print "$samcmd\n";
			system ($samcmd);
			$pileupfh = new FileHandle ($out_tumor);
			while (<$pileupfh>) {
				chomp;
				my $line = $_;
				$lineCounter_tumor++;		
			
				my @lineContents = split(/\t/, $line);			
				my $chrom = $lineContents[0];
				my $position = $lineContents[1];
				my $ref_base = $lineContents[2];
				my $depth = $lineContents[3];
				my $qualities = $lineContents[5];
	
				## Go through each quality ##
				
				my @qualities = split(//, $qualities);
				my $num_quals = 0;
				my $qual_coverage = 0;
				
				foreach my $code (@qualities) {
					my $qual_score = ord($code) - 33;
					$num_quals++;
		
					if($qual_score >= $min_base_qual)
					{
						$qual_coverage++;
					}
				}
				
				if($qual_coverage >= $min_coverage) {
#					print TUMOR_Q20 "$chrom\t$position\t$qual_coverage\n";			
					my $wiggle_chr = $chrom;
					$wiggle_chr =~ s/chr//;
					$wiggle_chr = "MT" if($chrom eq "M");
					my $key = "$wiggle_chr\t$position";
					$tumor_coverage{$key} = $qual_coverage;
				}
			}
			close($pileupfh);
	
			$cmd = "rm $out_tumor";
			system($cmd);
		
			print "Parsing regions file...\n";	
			my $input = new FileHandle ($regions_file);
			my $lineCounter = 0;
			while (<$input>)
			{
				chomp;
				my $line = $_;
				$lineCounter++;		
		
				(my $chrom, my $chr_start, my $chr_stop) = split(/\t/, $line);

				#set and test coverage for only coverage within the current iteration loop ==> $region (@reg_divides)
				if ($chrom ne $chr) {
					next;
				}

				my $rangestart;
				my $rangeend;
				if ($region =~ m/(\d+)-(\d+)/) {
					$rangestart = $1;
					$rangeend = $2;
				}
				elsif ($region =~ m/(\d)/) {
					$rangestart = $region;
					$rangeend = $chr_len{$chrom};
				}
				else {
					die "crap, parsed region wrong";
				}
				if ($chr_start > $rangeend) {
					next;
				}
				if ($chr_stop < $rangestart) {
					next;
				}
				if ($chr_stop > $rangeend) {
					$chr_stop = $rangeend;
				}

				if ($chr_start < $rangestart) {
					$chr_start = $rangestart;
				}
		
				## Print wiggle header ##
				print OUTFILE "fixedStep chrom=chr$chrom start=$chr_start step=1\n";			
		
				for(my $position = $chr_start; $position <= $chr_stop; $position++)
				{
					my $key = "$chrom\t$position";
#					$stats{'bases'}++;
		
					## Determine if coverage is met ##
					
					if($normal_coverage{$key} && $tumor_coverage{$key} && $normal_coverage{$key} >= $min_depth_normal && $tumor_coverage{$key} >= $min_depth_tumor)
					{				
						print OUTFILE "1\n";
#						$stats{'covered'}++;
					}
					else
					{
						print OUTFILE "0\n";
#						$stats{'not_covered'}++;
					}
				}
			}
			close($input);
		}
	}

	close(OUTFILE);
	rmdir($tempdir);
#	print $stats{'bases'} . " bases in ROI\n";
#	print $stats{'covered'} . " bases covered >= " . $min_depth_normal . "x in normal and >= " . $min_depth_tumor . "x in tumor\n";
#	print $stats{'not_covered'} . " bases NOT covered >= " . $min_depth_normal . "x in normal and >= " . $min_depth_tumor . "x in tumor\n";

	if($self->gzip_after)
	{
		print "Compressing $outfile...\n";
		system("gzip $outfile");
#		print "Compressing $q20_normal_file...\n";
#		system("gzip $q20_normal_file"); 
#		print "Compressing $q20_tumor_file...\n";
#		system("gzip $q20_tumor_file");  
	}
	
	return 1;
}




1;

