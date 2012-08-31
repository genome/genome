
package Genome::Model::Tools::Capture::GermlineMafStatistics;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# GermlineMafStatistics - Take Maf File and Generate Standard Statistics
#					
#	AUTHOR:		Will Schierding (wschierd@genome.wustl.edu)
#
#	CREATED:	02/28/2011 by W.S.
#	MODIFIED:	02/28/2011 by W.S.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;
use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it
use POSIX;

## Declare global statistics hash ##

my %stats = ();
my %stats2 = ();

class Genome::Model::Tools::Capture::GermlineMafStatistics {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		maf_file	=> { is => 'Text', doc => "Maf File To Read (or space separated list of Maf files to merge then read" , is_optional => 0},
		regions_file	=> { is => 'Text', doc => "Regions File -- Used to Get Targeted Region Basic Stats" , is_optional => 1},
		output_file	=> { is => 'Text', doc => "Output File With Statistics" , is_optional => 0},
		r_script_output_file     => { is => 'Text', doc => "R script built and run by this module (R not run if this option is not set)", is_optional => 1, is_input => 1},
		include_filterfailed     => { is => 'Text', doc => "Include stats output of filterfailed stats", is_optional => 1, is_input => 1, default => 0},
		plot_output_file_basename     => { is => 'Text', doc => "PDF output file basenames (include full path, will have suffixes appended automatically)", is_optional => 1, is_input => 1, is_output => 1 },
		skip_if_output_is_present     => { is => 'Text', doc => "Skip if Output is Present", is_optional => 1, is_input => 1, default => 0},
		genelist	=> { is => 'Text', doc => "One gene per line -- Used to Get Stats on these particular genes" , is_optional => 1},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Take Maf File and Generate Standard Statistics -- for GERMLINE projects"                 
}

sub help_synopsis {
    return <<EOS
Generate MAF File, Get dbsnp output, and strandfilter -- for GERMLINE events
EXAMPLE:	gmt capture germline-maf-statistics
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

	my $output_file = $self->output_file;
	my $maf_list = $self->maf_file;
	my $r_script_output_file = $self->r_script_output_file;
	my $skip_if_output_is_present = $self->skip_if_output_is_present;
	my $genelist;
	if ($self->genelist) {
		$genelist = $self->genelist;
	}
	if($self->r_script_output_file) {
		unless($self->plot_output_file_basename){
			die "Can't define r script without defining output files\n";
		}
	}
#	my %model_hash;
#	my $model_input = new FileHandle ($model_list_file);
#	while (my $line = <$model_input>) {
#		chomp($line);
#		$line =~ s/\s+/\t/g;
#		my ($model_id, $sample_name, $build_id, $build_status, $builddir) = split(/\t/, $line);
#		$model_hash{$sample_name} = "$model_id\t$build_id\t$builddir";
#	}
#	my $model_count = 0;
#	foreach (sort keys %model_hash) {
#		$model_count++;
#	}
#	print "Model List Loaded, $model_count Models in List\n";
	my %genehash;
	my %genehash_indel;
	my @excluded_genes;
	if ($self->genelist) {
		my $gene_file = new FileHandle ($genelist);
		while (my $line = <$gene_file>) {
			chomp($line);
			if ($line =~ m/exclude/) {
				$line =~ s/exclude//;
				$line =~ s/\(//;
				$line =~ s/\)//;
				@excluded_genes = split(/ /, $line);
			}
			else {
				$genehash{'Strandfilter_Passed'}{'novel'}{$line} = 0;
				$genehash{'Strandfilter_Passed'}{'dbsnp'}{$line} = 0;
				$genehash{'Strandfilter_Failed'}{'novel'}{$line} = 0;
				$genehash{'Strandfilter_Failed'}{'dbsnp'}{$line} = 0;
				$genehash_indel{'Strandfilter_Passed'}{'INS'}{$line} = 0;
				$genehash_indel{'Strandfilter_Passed'}{'DEL'}{$line} = 0;
				$genehash_indel{'Strandfilter_Failed'}{'INS'}{$line} = 0;
				$genehash_indel{'Strandfilter_Failed'}{'DEL'}{$line} = 0;
			}
		}
	}

	my $total_targeted_region = 0;
	if ($self->regions_file) {
		my $regions_file = $self->regions_file;
		my $region_input = new FileHandle ($regions_file);
#		my $header = <$region_input>;
		while (my $line = <$region_input>) {
			chomp($line);
			my ($region_chr,$region_start,$region_stop,$region_name) = split(/\t/, $line);
			my $area = $region_stop - $region_start + 1;
			$total_targeted_region += $area;
		}
	}

	## Open the outfile ##
	open(OUTFILE, ">$output_file") or die "Can't open output file: $!\n";
	if ($self->include_filterfailed) {
		my $output_file_failed = $output_file.".failed";
		open(OUTFILE2, ">$output_file_failed") or die "Can't open output file: $!\n";
	}
	my %variant_status;
	my @maffiles = split(/\s/, $maf_list);
	my %variants;
	foreach my $maf_file (@maffiles) {
		my %samples;
		my $input = new FileHandle ($maf_file);
		my $header = <$input>;
		while (my $line = <$input>) {
			chomp($line);
			my ($gene_name,$gene_id,$center,$ref_build,$chromosome,$chr_start,$chr_stop,$strand,$mutation_type,$variant_type,$ref,$tumor_gt_allele1,$tumor_gt_allele2,$dbsnp_rs,$dbsnp_status,$sample_name,$sample_name2,$ref2,$ref3,$val1,$val2,$val3,$val4,$strandfilter_status,$val_status,$mut_status,$val_method,$sequence_phase,$sequence_source,$score,$bam_file,$sequencer, @annotation) = split(/\t/, $line);
			unless ($self->include_filterfailed) {
				if ($strandfilter_status eq 'Strandfilter_Failed') {
					next;
				}
			}
			$samples{$sample_name}++;
			$stats{'TotalVariants'}++;
			$stats{$strandfilter_status}{$variant_type}++;
			$stats{$strandfilter_status}{'total'}++;
			my $variant = "$chromosome\t$chr_start\t$chr_stop\t$ref\t$tumor_gt_allele2";
			my $variant_pos = "$chromosome\t$chr_start\t$chr_stop";
			$variant_status{$variant_type}++;
			if ($variant_type eq 'SNP') {
				if ($dbsnp_rs =~ m/novel/i) {
					my $dbsnp = 'novel';
					$stats{$dbsnp}++;
					$stats{$strandfilter_status}{$dbsnp}++;
					$stats2{$strandfilter_status}{$dbsnp}{$variant_type}++;
					$variants{$strandfilter_status}{$dbsnp}{'variant'}{$variant}++;
					$variants{$strandfilter_status}{$dbsnp}{'position'}{$variant_pos}++;

					if (defined($genehash{$strandfilter_status}{$dbsnp}{$gene_name})) {
						$genehash{$strandfilter_status}{$dbsnp}{$gene_name}++;
					}
					else {
						my $skip =0;
						foreach my $genelist (@excluded_genes) {
							if ($gene_name =~ m/^$genelist/) {
								$skip = 1;
							}
						}
						unless ($skip) {
							my @genes = sort keys %{$genehash{'Strandfilter_Passed'}{'novel'}};
							foreach my $genelist (@genes) {
								#add alternatives
								$genelist =~ s/@//;
								$genelist =~ s/\d+$//;
								$genelist =~ s/-$//;
								$genelist =~ s/\d+$//;
								if ($gene_name =~ m/^$genelist/) {
									$genehash{$strandfilter_status}{$dbsnp}{$gene_name}++;
									unless (defined($genehash{$strandfilter_status}{'dbsnp'}{$gene_name})) {
										$genehash{$strandfilter_status}{'dbsnp'}{$gene_name} = 0;
									}
								}
							}
						}
					}
				}
				else {
					my $dbsnp = 'dbsnp';
					$stats{$dbsnp}++;
					$stats{$strandfilter_status}{$dbsnp}++;
					$stats2{$strandfilter_status}{$dbsnp}{$variant_type}++;
					$variants{$strandfilter_status}{$dbsnp}{'variant'}{$variant}++;
					$variants{$strandfilter_status}{$dbsnp}{'position'}{$variant_pos}++;
					if (defined($genehash{$strandfilter_status}{$dbsnp}{$gene_name})) {
						$genehash{$strandfilter_status}{$dbsnp}{$gene_name}++;
					}
					else {
						my $skip =0;
						foreach my $genelist (@excluded_genes) {
							if ($gene_name =~ m/^$genelist/) {
								$skip = 1;
							}
						}
						unless ($skip) {
							my @genes = sort keys %{$genehash{'Strandfilter_Passed'}{'novel'}};
							foreach my $genelist (@genes) {
								#add alternatives
								$genelist =~ s/@//;
								$genelist =~ s/\d+$//;
								$genelist =~ s/-$//;
								$genelist =~ s/\d+$//;
								if ($gene_name =~ m/^$genelist/) {
									$genehash{$strandfilter_status}{$dbsnp}{$gene_name}++;
									unless (defined($genehash{$strandfilter_status}{'novel'}{$gene_name})) {
										$genehash{$strandfilter_status}{'novel'}{$gene_name} = 0;
									}
								}
							}
						}
					}
				}
			}
			else { #indels, only for gene list looks
				my $dbsnp = $variant_type;
				if ($dbsnp ne 'INS' && $dbsnp ne 'DEL') {
					print "$dbsnp\n";
				}
				if (defined($genehash_indel{$strandfilter_status}{$dbsnp}{$gene_name})) {
					$genehash_indel{$strandfilter_status}{$dbsnp}{$gene_name}++;
				}
				else {
					my $skip =0;
					foreach my $genelist (@excluded_genes) {
						if ($gene_name =~ m/^$genelist/) {
							$skip = 1;
						}
					}
					unless ($skip) {
						my @genes = sort keys %{$genehash{'Strandfilter_Passed'}{'novel'}};
						foreach my $genelist (@genes) {
							#add alternatives
							$genelist =~ s/@//;
							$genelist =~ s/\d+$//;
							$genelist =~ s/-$//;
							$genelist =~ s/\d+$//;
							if ($gene_name =~ m/^$genelist/) {
								$genehash_indel{$strandfilter_status}{$dbsnp}{$gene_name}++;
							}
						}
					}
				}
			}
		}
		foreach my $sample_name (sort keys %samples) {
			$stats{'samples'}++;
		}
	}

#	print Data::Dumper::Dumper(%genehash_indel);

	my $sum = 0;
	my $sum2 = 0;
	foreach my $strandfilter_status (sort keys %variants) {
		foreach my $dbsnp (sort keys %{$variants{$strandfilter_status}}) {
			foreach my $variant (sort keys %{$variants{$strandfilter_status}{$dbsnp}{'variant'}}) {
				if (defined($variants{$strandfilter_status}{$dbsnp}{'variant'}{$variant}) && $variants{$strandfilter_status}{$dbsnp}{'variant'}{$variant} == 1) {
					$stats2{$strandfilter_status}{$dbsnp}{'singletons'}++;
					$sum += $variants{$strandfilter_status}{$dbsnp}{'variant'}{$variant};
				}
				elsif (defined($variants{$strandfilter_status}{$dbsnp}{'variant'}{$variant}) && $variants{$strandfilter_status}{$dbsnp}{'variant'}{$variant} > 1) {
					$stats2{$strandfilter_status}{$dbsnp}{'nonsingletons'}++;
					$sum2 += $variants{$strandfilter_status}{$dbsnp}{'variant'}{$variant};
				}
			}
			foreach my $variant (sort keys %{$variants{$strandfilter_status}{$dbsnp}{'position'}}) {
				$stats2{$strandfilter_status}{$dbsnp}{'TotalVariantsSNV'}++;
			}
		}
	}

#	print Data::Dumper::Dumper(%stats);
	print OUTFILE "Category\tTotal\tPerSample\n";
	print OUTFILE "Total Samples:\t$stats{'samples'}\t\n";
	if ($self->regions_file) {
		print OUTFILE "Total Size of Targeted Regions:\t$total_targeted_region\t\n";
	}
	else {
		print OUTFILE "Total Size of Targeted Regions:\tTarget Region Not Specified\t\n";
	}

	print OUTFILE "Total Unfiltered Variants:\t$stats{'TotalVariants'}"."\t". sprintf("%.2f", ($stats{'TotalVariants'} / $stats{'samples'}))."\n";
	print OUTFILE "Total Variants Passing Strand/Indel Filter:\t$stats{'Strandfilter_Passed'}{'total'}"."\t".sprintf("%.2f", ($stats{'Strandfilter_Passed'}{'total'} / $stats{'samples'}))."\n";

	foreach my $variant_type (sort keys %variant_status) {
		print OUTFILE "Total FilterPassed $variant_type"."s:\t$stats{'Strandfilter_Passed'}{$variant_type}"."\t".sprintf("%.2f", ($stats{'Strandfilter_Passed'}{$variant_type} / $stats{'samples'}))."\n";
#		if ($variant_type eq "SNP") {
#			print OUTFILE "Total FilterPassed Novel $variant_type"."s:\t$stats2{'Strandfilter_Passed'}{'novel'}{$variant_type}"."\t".sprintf("%.2f", ($stats2{'Strandfilter_Passed'}{'novel'}{$variant_type} / $stats{'samples'}))."\n";
#			print OUTFILE "Total FilterPassed $variant_type"."s in dbSNP:\t$stats2{'Strandfilter_Passed'}{'dbsnp'}{$variant_type}"."\t".sprintf("%.2f", ($stats2{'Strandfilter_Passed'}{'dbsnp'}{$variant_type} / $stats{'samples'}))."\n";
#		}
	}
	print OUTFILE "Total FilterPassed Novel SNVs:\t$stats{'Strandfilter_Passed'}{'novel'}"."\t".sprintf("%.2f", ($stats{'Strandfilter_Passed'}{'novel'} / $stats{'samples'}))."\n";
	print OUTFILE "Total FilterPassed SNVs in dbSNP:\t$stats{'Strandfilter_Passed'}{'dbsnp'}"."\t".sprintf("%.2f", ($stats{'Strandfilter_Passed'}{'dbsnp'} / $stats{'samples'}))."\n";
#	print OUTFILE "Total Novel SNVs:\t$stats{'novel'}"."\t".sprintf("%.2f", ($stats{'novel'} / $stats{'samples'}))."\n";
#	print OUTFILE "Total SNVs in dbSNP:\t$stats{'dbsnp'}"."\t".sprintf("%.2f", ($stats{'dbsnp'} / $stats{'samples'}))."\n";
#	print OUTFILE "\n";
	print OUTFILE "Total Number of SNVs That are Singletons:\t$sum"."\t". sprintf("%.2f", ($sum / $stats{'samples'}))."\n";
	print OUTFILE "Total Number of SNVs That are Not Singletons:\t$sum2"."\t". sprintf("%.2f", ($sum2 / $stats{'samples'}))."\n";
	my $singleton_sites = 0;
	my $nonsingleton_sites = 0;
	foreach my $strandfilter_status (sort keys %variants) {
		if ($strandfilter_status eq 'Strandfilter_Failed') {
			unless ($self->include_filterfailed) {
				next;
			}
		}
		foreach my $dbsnp (sort keys %{$variants{$strandfilter_status}}) {
			my $var1 = "Total $dbsnp $strandfilter_status SNV Sites:\t$stats2{$strandfilter_status}{$dbsnp}{'TotalVariantsSNV'}"."\t". sprintf("%.2f", ($stats2{$strandfilter_status}{$dbsnp}{'TotalVariantsSNV'} / $stats{'samples'}))."\n";
			my $var2 = "Total $strandfilter_status $dbsnp Singletons:\t$stats2{$strandfilter_status}{$dbsnp}{'singletons'}"."\t".sprintf("%.2f", ($stats2{$strandfilter_status}{$dbsnp}{'singletons'} / $stats{'samples'}))."\n";
			my $var3 = "Total $strandfilter_status $dbsnp nonSingletons:\t$stats2{$strandfilter_status}{$dbsnp}{'nonsingletons'}"."\t".sprintf("%.2f", ($stats2{$strandfilter_status}{$dbsnp}{'nonsingletons'} / $stats{'samples'}))."\n";
			unless ($strandfilter_status eq 'Strandfilter_Failed') {
				$singleton_sites += $stats2{$strandfilter_status}{$dbsnp}{'singletons'};
				$nonsingleton_sites += $stats2{$strandfilter_status}{$dbsnp}{'nonsingletons'};
			}
			if ($self->include_filterfailed) {
				print OUTFILE2 $var1;
				print OUTFILE2 $var2;
				print OUTFILE2 $var3;
			}
			else {
				print OUTFILE $var1;
				print OUTFILE $var2;
				print OUTFILE $var3;
			}
		}
	}
	print OUTFILE "Total Number of Sites That are Singletons:\t$singleton_sites"."\t".sprintf("%.2f", ($singleton_sites / $stats{'samples'}))."\n";
	print OUTFILE "Total Number of Sites That are Not Singletons:\t$nonsingleton_sites"."\t".sprintf("%.2f", ($nonsingleton_sites / $stats{'samples'}))."\n";

	if ($self->genelist) {
		foreach my $strandfilter_status (sort keys %genehash) {
			foreach my $gene (sort keys %{$genehash{$strandfilter_status}{'novel'}}) {
				foreach my $dbsnp (sort keys %{$genehash{$strandfilter_status}}) {
					if ($self->include_filterfailed && $strandfilter_status eq 'Strandfilter_Failed') {
						print OUTFILE2 "$dbsnp $gene variants found:\t".$genehash{$strandfilter_status}{$dbsnp}{$gene}."\t".sprintf("%.2f", ($genehash{$strandfilter_status}{$dbsnp}{$gene} / $stats{'samples'}))."\n";
					}
					elsif ($strandfilter_status eq 'Strandfilter_Passed') {
						print OUTFILE "$dbsnp $gene variants found:\t".$genehash{$strandfilter_status}{$dbsnp}{$gene}."\t".sprintf("%.2f", ($genehash{$strandfilter_status}{$dbsnp}{$gene} / $stats{'samples'}))."\n";
					}
				}
#				foreach my $dbsnp (sort keys %{$genehash{$strandfilter_status}}) {
#					if ($self->include_filterfailed && $strandfilter_status eq 'Strandfilter_Failed') {
#						print OUTFILE2 "$dbsnp $gene variants found:\t".$genehash_indel{$strandfilter_status}{$dbsnp}{$gene}."\t".sprintf("%.2f", ($genehash_indel{$strandfilter_status}{$dbsnp}{$gene} / $stats{'samples'}))."\n";
#					}
#					elsif ($strandfilter_status eq 'Strandfilter_Passed') {
#						print OUTFILE "$dbsnp $gene variants found:\t".$genehash_indel{$strandfilter_status}{$dbsnp}{$gene}."\t".sprintf("%.2f", ($genehash_indel{$strandfilter_status}{$dbsnp}{$gene} / $stats{'samples'}))."\n";
#					}
#				}
			}
		}
	}

	if ($self->include_filterfailed) {
		print OUTFILE2 "Total Variants Failing Strand/Indel Filter:\t$stats{'Strandfilter_Failed'}{'total'}"."\t".sprintf("%.2f", ($stats{'Strandfilter_Failed'}{'total'} / $stats{'samples'}))."\n";
	}
	if ($self->include_filterfailed) {
		foreach my $variant_type (sort keys %variant_status) {
			print OUTFILE2 "Total FilterFailed $variant_type"."s:\t$stats{'Strandfilter_Failed'}{$variant_type}"."\t".sprintf("%.2f", ($stats{'Strandfilter_Failed'}{$variant_type} / $stats{'samples'}))."\n";
#			if ($variant_type eq "SNP") {
#				print OUTFILE2 "Total FilterFailed Novel $variant_type"."s:\t$stats2{'Strandfilter_Failed'}{'novel'}{$variant_type}"."\t".sprintf("%.2f", ($stats2{'Strandfilter_Failed'}{'novel'}{$variant_type} / $stats{'samples'}))."\n";
#				print OUTFILE2 "Total FilterFailed $variant_type"."s in dbSNP:\t$stats2{'Strandfilter_Failed'}{'dbsnp'}{$variant_type}"."\t".sprintf("%.2f", ($stats2{'Strandfilter_Failed'}{'dbsnp'}{$variant_type} / $stats{'samples'}))."\n";
#			}
		}
	}
	if ($self->include_filterfailed) {
		print OUTFILE2 "Total FilterFailed Novel SNVs:\t$stats{'Strandfilter_Failed'}{'novel'}"."\t".sprintf("%.2f", ($stats{'Strandfilter_Failed'}{'novel'} / $stats{'samples'}))."\n";
		print OUTFILE2 "Total FilterFailed SNVs in dbSNP:\t$stats{'Strandfilter_Failed'}{'dbsnp'}"."\t".sprintf("%.2f", ($stats{'Strandfilter_Failed'}{'dbsnp'} / $stats{'samples'}))."\n";
	}

	close(OUTFILE);
	if ($self->include_filterfailed) {
		close(OUTFILE2);
	}


	if($self->r_script_output_file) {
		my $outfile_base = $self->plot_output_file_basename;
		# Open Output
		unless (open(R_COMMANDS,">$r_script_output_file")) {
		    die "Could not open output file '$r_script_output_file' for writing";
		  }
#		print R_COMMANDS 'options(echo = FALSE)'."\n";#suppress output to stdout
		print R_COMMANDS 'sink("/dev/null");'."\n";
		print R_COMMANDS "x <- read.table(file=\'$output_file\',sep=\"\\t\",row.names=1,blank.lines.skip = TRUE,header=TRUE);"."\n";
		print R_COMMANDS 'samples <- x$Total[1];'."\n";
		print R_COMMANDS 'colors <- c("red","blue");'."\n";
		my $piechart = $outfile_base.".piecharts.pdf";
		print R_COMMANDS "pdf(file=\"$piechart\",width=10,height=7.5);"."\n";
		print R_COMMANDS 'par(mfrow=c(2,2));'."\n";

		print R_COMMANDS 'pie(c(x$PerSample[8],x$PerSample[9]),labels=c((x[8,2]),(x[9,2])),col=colors,main="dbSNP versus Novel Variants (per sample)");'."\n";
		print R_COMMANDS 'dbsnp <- x$PerSample[9] / x$PerSample[7];'."\n";
		print R_COMMANDS 'mtext(paste(round(dbsnp,digits=4)*100,"% dbSNPs",sep=""), padj=-0.5);'."\n";
		print R_COMMANDS 'mtext(paste(samples,"Samples",sep=" "), padj=-0.7, side=1);'."\n";
		print R_COMMANDS 'legend(x="topright", title = "Variant Type", c("Novel SNVs", "dbSNP SNVs"),col=colors,cex=0.8,fill=colors);'."\n";

		print R_COMMANDS 'colors <- c("darkgreen","greenyellow");'."\n";
		print R_COMMANDS 'pie(c(x$PerSample[5],x$PerSample[6]),labels=c((x[5,2]),(x[6,2])),col=colors,main="Indels (per sample)");'."\n";
		print R_COMMANDS 'mtext(paste(samples,"Samples",sep=" "), padj=-0.7, side=1);'."\n";
		print R_COMMANDS 'DEL <- x$PerSample[5] / (x$PerSample[5]+x$PerSample[6]);'."\n";
		print R_COMMANDS 'INS <- x$PerSample[6] / (x$PerSample[5]+x$PerSample[6]);'."\n";
		print R_COMMANDS 'dellab <- paste(round(DEL,digits=4)*100,"% Deletions",sep="");'."\n";
		print R_COMMANDS 'inslab <- paste(round(INS,digits=4)*100,"% Insertions",sep="");'."\n";
		print R_COMMANDS 'mtext(paste(inslab,dellab), padj=-0.5);'."\n";
		print R_COMMANDS 'legend(x="topright", title = "Variant Type", c("Deletions", "Insertions"),col=colors,cex=0.8,fill=colors);'."\n";

		print R_COMMANDS 'colors <- c("darkred","darksalmon");'."\n";
		print R_COMMANDS 'pie(c(x$PerSample[16],x$PerSample[17]),labels=c((x[16,2]),(x[17,2])),col=colors,main="Novel Single versus Recurrent Variant Sites (per sample)");'."\n";
		print R_COMMANDS 'mtext(paste(samples,"Samples",sep=" "), padj=-0.7, side=1);'."\n";
		print R_COMMANDS 'singleton <- x$PerSample[16] / x$PerSample[15];'."\n";
		print R_COMMANDS 'mtext(paste(round(singleton,digits=4)*100,"% Novel are Singletons",sep=""), padj=-0.5);'."\n";
		print R_COMMANDS 'legend(x="topright", title = "Variant Type", c("Novel Singletons", "Novel Recurrent"),col=colors,cex=0.8,fill=colors);'."\n";

		print R_COMMANDS 'colors <- c("darkblue","steelblue1");'."\n";
		print R_COMMANDS 'pie(c(x$PerSample[13],x$PerSample[14]),labels=c((x[13,2]),(x[14,2])),col=colors,main="dbSNP Single versus Recurrent Variant Sites (per sample)");'."\n";
		print R_COMMANDS 'mtext(paste(samples,"Samples",sep=" "), padj=-0.7, side=1);'."\n";
		print R_COMMANDS 'singleton <- x$PerSample[13] / x$PerSample[12];'."\n";
		print R_COMMANDS 'mtext(paste(round(singleton,digits=4)*100,"% dbSNP are Singletons",sep=""), padj=-0.5);'."\n";
		print R_COMMANDS 'legend(x="topright", title = "Variant Type", c("dbSNP Singletons", "dbSNP Recurrent"),col=colors,cex=0.8,fill=colors);'."\n";

#		print R_COMMANDS 'mtext("Outer Margin Area", side=2, line=2, cex=2, col="blue", outer=TRUE);'."\n";

		if ($self->genelist) {
			my @genes = sort keys %{$genehash{'Strandfilter_Passed'}{'novel'}};
			my $size = @genes;

			my $start = my $init = 20;
			my $stop = $start + (2 * $size) - 1;
			print R_COMMANDS 'par(mfrow=c(1,1));'."\n";
			print R_COMMANDS 'counts = 0;'."\n";
			print R_COMMANDS 'genes = 0;'."\n";
			print R_COMMANDS 'countstotal = 0;'."\n";
			print R_COMMANDS 'genestotal = 0;'."\n";
			print R_COMMANDS 'countsindel = 0;'."\n";
			print R_COMMANDS 'genesindel = 0;'."\n";
#			print R_COMMANDS 'counts1 = 0;'."\n";
#			print R_COMMANDS 'counts2 = 0;'."\n";
#			print R_COMMANDS 'genes1 = 0;'."\n";
#			print R_COMMANDS 'genes2 = 0;'."\n";

			my $count = 1;
			while ($start <= $stop) {
				my $gene = $genes[(($start - $init) + 0)/2];
#				print R_COMMANDS "counts1[$count] <- x\$PerSample[$start+0]"."\n";
#				print R_COMMANDS "counts2[$count] <- x\$PerSample[$start+1]"."\n";
#				print R_COMMANDS "genes1[$count] <- \"$gene dbsnp\""."\n";
#				print R_COMMANDS "genes2[$count] <- \"$gene novel\""."\n";
				print R_COMMANDS "counts[$start+0-$init+1] <- x\$PerSample[$start+0]"."\n";
				print R_COMMANDS "counts[$start+1-$init+1] <- x\$PerSample[$start+1]"."\n";
				print R_COMMANDS "genes[$start+0-$init+1] <- \"$gene\""."\n";
				print R_COMMANDS "genes[$start+1-$init+1] <- \"\""."\n";
				print R_COMMANDS "countstotal[$count] <- (x\$PerSample[$start+0] + x\$PerSample[$start+1])"."\n";
				print R_COMMANDS "genestotal[$count] <- \"$gene\""."\n";

				$start = ($start + 2);
				$count++;
			}

			$count = 1;
			my @legend;
			my $strandfilter_status = 'Strandfilter_Passed';
			my @indel_types = sort keys %{$genehash_indel{$strandfilter_status}};
			my $number_of_indel_types = @indel_types;
			foreach my $gene (@genes) {
				foreach my $dbsnp (sort keys %{$genehash_indel{$strandfilter_status}}) {
					if ($gene eq $genes[0]) {
						push(@legend, "\"$dbsnp\"");
					}
					unless (defined ($genehash_indel{$strandfilter_status}{$dbsnp}{$gene})) {
						$genehash_indel{$strandfilter_status}{$dbsnp}{$gene} = 0;
					}
					my $avg = ($genehash_indel{$strandfilter_status}{$dbsnp}{$gene} / $stats{'samples'});
					print R_COMMANDS "countsindel[$count] <- $avg"."\n";
					if ($count % 2) {
						print R_COMMANDS "genesindel[$count] <- \"$gene\""."\n";
					}
					else {
						print R_COMMANDS "genesindel[$count] <- \"\""."\n";
					}
					$count++;
				}
			}
			my $types = join(",",@legend);

#			print R_COMMANDS 'counts <- cbind(counts1, counts2)'."\n";
#			print R_COMMANDS 'genes <- cbind(genes1, genes2)'."\n";
			my $scale = (1 / $size);
			$scale = 0.5;
			my $at = ($size * 2);
			print R_COMMANDS 'colors <- c("darkblue","darkred");'."\n";
			print R_COMMANDS "barplot(counts,main=\"Gene Variant Distribution (per sample)\",names.arg=genes,beside=TRUE,legend=c(\"dbsnp\",\"novel\"),col=colors,cex.names=$scale, las=2)"."\n";
			print R_COMMANDS 'mtext(paste(samples,"Samples",sep=" "), padj=-0.5);'."\n";
			print R_COMMANDS 'colors <- c("purple");'."\n";
			print R_COMMANDS "barplot(countstotal,main=\"Gene Variant Distribution (per sample)\",names.arg=genestotal,xlab=\"Genes\",beside=TRUE,cex.names=$scale, las=2, col=colors)"."\n";
			print R_COMMANDS 'mtext(paste(samples,"Samples",sep=" "), padj=-0.5);'."\n";
			print R_COMMANDS "colors <- topo.colors($number_of_indel_types);"."\n";
			print R_COMMANDS "barplot(countsindel,main=\"Gene Variant Distribution (per sample)\",names.arg=genesindel,beside=TRUE,legend=c($types),col=colors,cex.names=$scale, las=2)"."\n";
			print R_COMMANDS 'mtext(paste(samples,"Samples",sep=" "), padj=-0.5);'."\n";
#			print R_COMMANDS "axis(side=1, at=1:$at, labels=genes, las=2,cex.axis=$scale)"."\n";

			$start = $init;
			$stop = $start + (2 * $size) - 1;
			while ($start <= $stop) {
				my $gene1 = $genes[(($start - $init) + 0)/2];
				my $gene2 = $genes[(($start - $init) + 2)/2];
				my $gene3 = $genes[(($start - $init) + 4)/2];
				my $gene4 = $genes[(($start - $init) + 6)/2];

				print R_COMMANDS 'par(mfrow=c(2,2));'."\n";
				print R_COMMANDS 'colors <- c("blue","cornflowerblue");'."\n";

				if ($genehash{'Strandfilter_Passed'}{'novel'}{$gene1} == 0 && $genehash{'Strandfilter_Passed'}{'dbsnp'}{$gene1} == 0) {
					print R_COMMANDS 'frame();'."\n";
					print R_COMMANDS 'box("plot", col="black");'."\n";
					print R_COMMANDS "text(.5, .5, \"No Variants for $gene1\");"."\n";
				}
				else {
					print R_COMMANDS "pie(c(x\$PerSample[$start+0],x\$PerSample[$start+1]),labels=c((x[$start+0,2]),(x[$start+1,2])),col=colors,main=\"$gene1 dbSNP versus Novel Variants (per sample)\",sub=paste(samples,\"Samples\",sep=\" \"));"."\n";
					print R_COMMANDS "dbsnp <- x\$PerSample[$start+0] / (x\$PerSample[$start+0] + x\$PerSample[$start+1]);"."\n";
					print R_COMMANDS 'mtext(paste(round(dbsnp,digits=4)*100,"% dbSNPs",sep=""), padj=-0.5);'."\n";
					print R_COMMANDS 'legend(x="topright", title = "Variant Type", c("dbSNP SNVs","Novel SNVs"),col=colors,cex=0.8,fill=colors);'."\n";
				}

				if ($genehash{'Strandfilter_Passed'}{'novel'}{$gene2} == 0 && $genehash{'Strandfilter_Passed'}{'dbsnp'}{$gene2} == 0) {
					print R_COMMANDS 'frame();'."\n";
					print R_COMMANDS 'box("plot", col="black");'."\n";
					print R_COMMANDS "text(.5, .5, \"No Variants for $gene2\");"."\n";
				}
				else {
					print R_COMMANDS "pie(c(x\$PerSample[$start+2],x\$PerSample[$start+3]),labels=c((x[$start+2,2]),(x[$start+3,2])),col=colors,main=\"$gene2 dbSNP versus Novel Variants (per sample)\",sub=paste(samples,\"Samples\",sep=\" \"));"."\n";
					print R_COMMANDS "dbsnp <- x\$PerSample[$start+2] / (x\$PerSample[$start+2] + x\$PerSample[$start+3]);"."\n";
					print R_COMMANDS 'mtext(paste(round(dbsnp,digits=4)*100,"% dbSNPs",sep=""), padj=-0.5);'."\n";
					print R_COMMANDS 'legend(x="topright", title = "Variant Type", c("dbSNP SNVs","Novel SNVs"),col=colors,cex=0.8,fill=colors);'."\n";
				}

				if ($genehash{'Strandfilter_Passed'}{'novel'}{$gene3} == 0 && $genehash{'Strandfilter_Passed'}{'dbsnp'}{$gene3} == 0) {
					print R_COMMANDS 'frame();'."\n";
					print R_COMMANDS 'box("plot", col="black");'."\n";
					print R_COMMANDS "text(.5, .5, \"No Variants for $gene3\");"."\n";
				}
				else {
					print R_COMMANDS "pie(c(x\$PerSample[$start+4],x\$PerSample[$start+5]),labels=c((x[$start+4,2]),(x[$start+5,2])),col=colors,main=\"$gene3 dbSNP versus Novel Variants (per sample)\",sub=paste(samples,\"Samples\",sep=\" \"));"."\n";
					print R_COMMANDS "dbsnp <- x\$PerSample[$start+4] / (x\$PerSample[$start+4] + x\$PerSample[$start+5]);"."\n";
					print R_COMMANDS 'mtext(paste(round(dbsnp,digits=4)*100,"% dbSNPs",sep=""), padj=-0.5);'."\n";
					print R_COMMANDS 'legend(x="topright", title = "Variant Type", c("dbSNP SNVs","Novel SNVs"),col=colors,cex=0.8,fill=colors);'."\n";
				}

				if ($genehash{'Strandfilter_Passed'}{'novel'}{$gene4} == 0 && $genehash{'Strandfilter_Passed'}{'dbsnp'}{$gene4} == 0) {
					print R_COMMANDS 'frame();'."\n";
					print R_COMMANDS 'box("plot", col="black");'."\n";
					print R_COMMANDS "text(.5, .5, \"No Variants for $gene4\");"."\n";
				}
				else {
					print R_COMMANDS "pie(c(x\$PerSample[$start+6],x\$PerSample[$start+7]),labels=c((x[$start+6,2]),(x[$start+7,2])),col=colors,main=\"$gene4 dbSNP versus Novel Variants (per sample)\",sub=paste(samples,\"Samples\",sep=\" \"));"."\n";
					print R_COMMANDS "dbsnp <- x\$PerSample[$start+6] / (x\$PerSample[$start+6] + x\$PerSample[$start+7]);"."\n";
					print R_COMMANDS 'mtext(paste(round(dbsnp,digits=4)*100,"% dbSNPs",sep=""), padj=-0.5);'."\n";
					print R_COMMANDS 'legend(x="topright", title = "Variant Type", c("dbSNP SNVs","Novel SNVs"),col=colors,cex=0.8,fill=colors);'."\n";
				}

				$start = ($start + 8);
			}
		}


	 	print R_COMMANDS 'devoff <- dev.off();'."\n";
		print R_COMMANDS "q()\n";

		my $cmd = "R --vanilla --slave \< $r_script_output_file";
		my $return = Genome::Sys->shellcmd(
	           cmd => "$cmd",
	           output_files => [$piechart],
	           skip_if_output_is_present => $skip_if_output_is_present,
		   );
		unless($return) { 
			$self->error_message("Failed to execute: Returned $return");
			die $self->error_message;
		}
		return $return;
	}
}


=cut
Number of snps, ins, del
Avg Per Sample
Total Number of Sites with Variant
#unique, #dbsnp
	Of each, #snv, #ins, #del, #other
	average per sample
average allele frequency

given roi file, what is target region...
How many singletons are there per sample? 
%dbsnp




=cut

#my ($Sample, $SNPsCalled, $WithGenotype, $MetMinDepth, $Reference, $RefMatch, $RefWasHet, $RefWasHom, $Variant, $VarMatch, $HomWasHet, $HetWasHom, $VarMismatch, $VarConcord, $RareHomConcord, $OverallConcord) = split(/\t/, $qc_line);

################################################################################################
# SUBS
#
################################################################################################


#############################################################
# IUPAC to base - convert IUPAC code to variant base
#
#############################################################

sub iupac_to_base
{
	(my $allele1, my $allele2) = @_;
	
	return($allele2) if($allele2 eq "A" || $allele2 eq "C" || $allele2 eq "G" || $allele2 eq "T");
	
	if($allele2 eq "M")
	{
		return("C") if($allele1 eq "A");
		return("A") if($allele1 eq "C");
	}
	elsif($allele2 eq "R")
	{
		return("G") if($allele1 eq "A");
		return("A") if($allele1 eq "G");		
	}
	elsif($allele2 eq "W")
	{
		return("T") if($allele1 eq "A");
		return("A") if($allele1 eq "T");		
	}
	elsif($allele2 eq "S")
	{
		return("C") if($allele1 eq "G");
		return("G") if($allele1 eq "C");		
	}
	elsif($allele2 eq "Y")
	{
		return("C") if($allele1 eq "T");
		return("T") if($allele1 eq "C");		
	}
	elsif($allele2 eq "K")
	{
		return("G") if($allele1 eq "T");
		return("T") if($allele1 eq "G");				
	}	
	
	return($allele2);
}

#############################################################
# ParseBlocks - takes input file and parses it
#
#############################################################

sub trv_to_mutation_type
{
	my $trv_type = shift(@_);
	
	return("Missense_Mutation") if($trv_type eq "missense");	
	return("Nonsense_Mutation") if($trv_type eq "nonsense" || $trv_type eq "nonstop");	
	return("Silent") if($trv_type eq "silent");		
	return("Splice_Site_SNP") if($trv_type eq "splice_site");
	return("Splice_Site_Indel") if($trv_type eq "splice_site_del");		
	return("Splice_Site_Indel") if($trv_type eq "splice_site_ins");		
	return("Frame_Shift_Del") if($trv_type eq "frame_shift_del");		
	return("Frame_Shift_Ins") if($trv_type eq "frame_shift_ins");		
	return("In_Frame_Del") if($trv_type eq "in_frame_del");		
	return("In_Frame_Ins") if($trv_type eq "in_frame_ins");		
	return("RNA") if($trv_type eq "rna");		

	warn "Unknown mutation type $trv_type\n";
	return("Unknown");
}


################################################################################################
# Execute - the main program logic
#
################################################################################################

sub byChrPos
{
	my ($chrom_a, $pos_a) = split(/\t/, $a);
	my ($chrom_b, $pos_b) = split(/\t/, $b);
	
	$chrom_a =~ s/X/23/;
	$chrom_a =~ s/Y/24/;
	$chrom_a =~ s/MT/25/;
	$chrom_a =~ s/M/25/;
	$chrom_a =~ s/[^0-9]//g;

	$chrom_b =~ s/X/23/;
	$chrom_b =~ s/Y/24/;
	$chrom_b =~ s/MT/25/;
	$chrom_b =~ s/M/25/;
	$chrom_b =~ s/[^0-9]//g;

	$chrom_a <=> $chrom_a
	or
	$pos_a <=> $pos_b;
}

1;
