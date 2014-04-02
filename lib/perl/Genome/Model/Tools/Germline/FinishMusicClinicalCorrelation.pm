
package Genome::Model::Tools::Germline::FinishMusicClinicalCorrelation;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# FinishMusicClinicalCorrelation - Generate Tables and Statistics from Music Clinical Correlation Output
#					
#	AUTHOR:		Will Schierding (wschierd@genome.wustl.edu)
#
#	CREATED:	09/29/2010 by W.S.
#	MODIFIED:	09/29/2010 by W.S.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;
use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

## Declare global statistics hash ##

my %stats = ();

class Genome::Model::Tools::Germline::FinishMusicClinicalCorrelation {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
        input_file => { 
            is => 'Text',
            doc => "File of Results of clinical-correlation tool",
            is_input => 1,
            is_optional => 0,
        },
	clinical_data_file => { 
            is => 'Text',
            doc => "Same input as in clinical-correlation tool",
            is_input => 1,
            is_optional => 0,
        },
        maf_file => { 
            is => 'Text',
            doc => "List of mutations in MAF format",
            is_input => 1,
            is_optional => 0,
        },
        fdr_cutoff => { 
            is => 'Text',
            doc => "What FDR Cutoff to use?",
            is_input => 1,
            is_optional => 0,
	    default => 0.05,
        },
        project_name => { 
            is => 'Text',
            doc => "Project Name",
            is_input => 1,
            is_optional => 0,
	    default => "No Project Name",
        },
        output_file => {
            is => 'Text',
            is_output => 1,
            doc => "Results of clinical-correlation tool",
        },
        output_pdf_image_file => {
            is => 'Text',
            is_output => 1,
            doc => "Results of clinical-correlation tool",
        },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Generate Tables and Statistics from Music Clinical Correlation Output -- for GERMLINE projects"                 
}

sub help_synopsis {
    return <<EOS
Generate Tables and Statistics from Music Clinical Correlation Output -- for GERMLINE projects
EXAMPLE:	gmt germline finish-music-clinical-correlation
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
	my $input_file = $self->input_file;
	my $clinical_data_file = $self->clinical_data_file;
	my $maf_file = $self->maf_file;
	my $output_file = $self->output_file;
	my $output_pdf_image_file = $self->output_pdf_image_file;
	my $project_name = $self->project_name;
	my %stats;
	my %significant_hash;
	my %gene_hash;
	my %clinical_hash;
	my %log10position;

	my %maf_stuff;
	my $maf_input = new FileHandle ($maf_file);
	my $header = <$maf_input>;
	while (my $line = <$maf_input>) {
		chomp($line);
		my ($Hugo_Symbol, $Entrez_Gene_Id, $GSC_Center, $NCBI_Build, $Chromosome, $Start_position, $End_position, $Strand, $Variant_Classification, $Variant_Type, $Reference_Allele, $Variant_Allele1, $Variant_Allele2, $dbSNP_RS, $dbSNP_Val_Status, $Sample_Barcode1, $Sample_Barcode2, $Match_Norm_Seq_Allele1, $Match_Norm_Seq_Allele2, $Validation_Allele1, $Validation_Allele2, $Match_Norm_Validation_Allele1, $Match_Norm_Validation_Allele2, $Verification_Status, $Validation_Status, $Mutation_Status, $Validation_Method, $Sequencing_Phase, $Sequence_Source, $Score, $BAM_file, $Sequencer, $chromosome_name, $start, $stop, $reference, $variant, $type, $gene_name, $transcript_name, $transcript_species, $transcript_source, $transcript_version, $strand, $transcript_status, $trv_type, $c_position, $amino_acid_change, $ucsc_cons, $domain, $all_domains, $deletion_substructures, $transcript_error) = split(/\t/, $line);
		if ($Verification_Status eq 'Strandfilter_Failed') { next;}
		unless ($Variant_Type eq 'SNP') { next;}
		(my $sample) = $Sample_Barcode1 =~ m/(H_HY-\d\d\d\d\d)/;
		my $genevariant = "$Hugo_Symbol"."_"."$Chromosome"."_"."$Start_position"."_"."$End_position"."_"."$Reference_Allele"."_"."$Variant_Allele2";
		my $genevariant2 = "$Chromosome"."_"."$Start_position"."_"."$End_position";
		$maf_stuff{$genevariant} = $dbSNP_RS;
		$maf_stuff{$genevariant2} = $dbSNP_RS;
	}
	close($maf_input);

	my %fdr_cutoff;
	my $signif_value = $self->fdr_cutoff; ###WHAT FDR CUTOFF TO USE FOR SIGNIFICANCE TEST
	my $file_input = new FileHandle ($input_file);
	$header = <$file_input>;
	while (my $line = <$file_input>) {
		chomp($line);
		my ($geneinfo, $clinical, $method, $sample_number, $score, $pvalue, $fdr, $bonferonni) = split(/,/, $line);
		if ($pvalue =~ m/^0.000000$/) {
			warn "hit p-value cutoff, rounding up to 0.000001\n";
			$pvalue = 0.000001;
		}
		my ($gene, $chr, $start, $stop, $ref, $var) = split(/_/, $geneinfo);
		$gene_hash{'total'}{$gene}++;
		$clinical_hash{'total'}{$clinical}++;
		if ($ref eq '.' || $var eq '.') {
			$stats{'indel'}++;
		}
		else {
			$stats{'snv'}++;
		}
		if ($pvalue ne 'NA') {
			if ($pvalue == 0.000000) {
				$pvalue = 0.000001;
			}
			my $log = -log10($pvalue);
			$log10position{$clinical}{$chr}{$start}{$stop} = $log;
		}
		if ($fdr ne 'NA' && $fdr <= $signif_value) {
			$fdr_cutoff{$clinical}{$pvalue} = $fdr;
			$stats{'significant'}++;
			$gene_hash{'significant'}{$gene}++;
			$clinical_hash{'significant'}{$clinical}++;
			if ($ref eq '.' || $var eq '.') {
				$stats{'sigindel'}++;
			}
			else {
				$stats{'sigsnv'}++;
			}
			$significant_hash{$clinical}{$gene}++;
		}
	}

	my $total_genes = keys %{$gene_hash{'total'}};
	my $total_clins = keys %{$clinical_hash{'total'}};

	my $count = my $total = 0;
	foreach my $gene (sort keys %{$gene_hash{'significant'}}) {
		$count++;
		$total += $gene_hash{'significant'}{$gene};
	}
    if ($count) {
    	$stats{'avgvar'} = ($total / $count);
    } else {
        $stats{'avgvar'} = 0;
    }
	$stats{'sigvar'} = $count;
	$stats{'totalvar'} = $total_genes;

	$count = $total = 0;
	foreach my $clinical (sort keys %{$clinical_hash{'significant'}}) {
		$count++;
		$total += $clinical_hash{'significant'}{$clinical};
	}

    if ($count) {
    	$stats{'avgclin'} = ($total / $count);
    } else {
    	$stats{'avgclin'} = 0;
    }
	$stats{'sigclin'} = $count;
	$stats{'totalclin'} = $total_clins;

	## Open the outfile ##
	open(OUTFILE, ">$output_file") or die "Can't open output file: $!\n";

	print OUTFILE "#Total Genes\t$stats{'totalvar'}\n";
	print OUTFILE "#Total Phenotypes\t$stats{'totalclin'}\n";
	print OUTFILE "#Total SNVs\t$stats{'snv'}\n";
	print OUTFILE "#Total Indels\t$stats{'indel'}\n" if (defined $stats{'indel'});
	print OUTFILE "#Significant Variants\t$stats{'significant'}\n";
	print OUTFILE "#Significant Genes\t$stats{'sigvar'}\n";
	print OUTFILE "#Significant SNVs\t$stats{'sigsnv'}\n";
	print OUTFILE "#Significant Indels\t$stats{'sigindel'}\n" if (defined $stats{'sigindel'});
	print OUTFILE "#Significant Phenotypes\t$stats{'sigclin'}\n";
	print OUTFILE "#Average Significant Variants Per Significant Gene\t$stats{'avgvar'}\n";
	print OUTFILE "#Average Significant Variants Per Significant Phenotype\t$stats{'avgclin'}\n";

	print OUTFILE "#Clinical Variable";
	foreach my $gene (sort keys %{$gene_hash{'total'}}) {
		print OUTFILE "\t$gene";
	}
	print OUTFILE "\n";

	foreach my $clinical (sort keys %{$clinical_hash{'total'}}) {
		print OUTFILE "$clinical";
		foreach my $gene (sort keys %{$gene_hash{'total'}}) {
			if (defined $significant_hash{$clinical}{$gene}) {
				print OUTFILE "\t$significant_hash{$clinical}{$gene}";
			}
			else {
				print OUTFILE "\t0";
			}
		}
		print OUTFILE "\n";
	}
	foreach my $clinical (sort keys %log10position) {
		print OUTFILE "#Clinical Phenotype: $clinical\n";
		print OUTFILE "#Chromosome\tStart\tStop\tLog10(pvalue)\n";
		foreach my $chr (sort { $a <=> $b } keys %{$log10position{$clinical}}) {
			foreach my $start (sort { $a <=> $b } keys %{$log10position{$clinical}{$chr}}) {
				foreach my $stop (sort { $a <=> $b } keys %{$log10position{$clinical}{$chr}{$start}}) {
					print OUTFILE "$chr\t$start\t$stop\t$log10position{$clinical}{$chr}{$start}{$stop}\n";
				}
			}
		}
	}

	my ($tfh_R,$temp_path_R) = Genome::Sys->create_temp_file;
	unless($tfh_R) {
		$self->error_message("Unable to create temporary file $!");
		die;
	}
	$temp_path_R =~ s/\:/\\\:/g;

#	print $tfh_R 'sink("/dev/null")'."\n";
	print $tfh_R "genome=\"$project_name\";"."\n";
	print $tfh_R 'library(fpc);'."\n";
	print $tfh_R 'library(scatterplot3d);'."\n";
	print $tfh_R "pdf(file=\"$output_pdf_image_file\",width=10,height=7.5);"."\n";
	print $tfh_R "clindatatable <- read.table(\"$clinical_data_file\", row.names = 1, header = TRUE, sep = \"\\t\");"."\n";
	print $tfh_R 'par(mfrow=c(1,1));'."\n";

	my @clinicals = sort keys %{$clinical_hash{'total'}};
	my $num_clins = scalar(@clinicals);
	my @genes = sort keys %{$gene_hash{'total'}};
	my $num_genes = scalar(@genes);

	my %connections;
	my @correlate = @clinicals;
	for (my $i = 1; $i <= ($num_genes - 1); $i++) {
		my $gene1 = pop(@correlate);
		foreach my $nextgene (@correlate) {
			$connections{$gene1}{$nextgene}++;
		}
	}
	my $factors;
	$factors += $_ foreach 1..($num_clins - 1);
	my $splits = int($factors / 5)+1;
#	print $tfh_R "par(mfrow=c(5,$splits));"."\n";
        print $tfh_R 'corr = 0;'."\n";
	my @connect_names;
	my @connect_values;
	foreach my $clinical1 (sort keys %connections) {
		foreach my $clinical2 (sort keys %{$connections{$clinical1}}) {
			my $name = "$clinical1"."_"."$clinical2";
			if ($clinical1 eq 'hdlres') {
				print $tfh_R "corr\$$name <- cor.test(-clindatatable\$$clinical1,clindatatable\$$clinical2,na.rm=TRUE);"."\n";
				push(@connect_names,$name);
				push(@connect_values,"corr\$$name\$estimate");
			}
			elsif ($clinical2 eq 'hdlres') {
				print $tfh_R "corr\$$name <- cor.test(clindatatable\$$clinical1,-clindatatable\$$clinical2,na.rm=TRUE);"."\n";
				push(@connect_names,$name);
				push(@connect_values,"corr\$$name\$estimate");
			}
			else {
				print $tfh_R "corr\$$name <- cor.test(clindatatable\$$clinical1,clindatatable\$$clinical2,na.rm=TRUE);"."\n";
				push(@connect_names,$name);
				push(@connect_values,"corr\$$name\$estimate");
			}
		}
	}
	my $clins = "stuff <- cbind(".join(",",@connect_values).")";
	my $clins_names = "stuff2 <- cbind(\"".join("\",\"",@connect_names)."\")";
	print $tfh_R "$clins;"."\n";
	print $tfh_R "$clins_names;"."\n";
	print $tfh_R "b <- barplot(height=stuff,main=\"Gene Distribution (per clinical category)\",col=\"black\",las=2,xaxt=\"n\",ylim=c(-1,1));"."\n";
	print $tfh_R "axis(1,at=b,labels=stuff2, las=2, pos=-.75);"."\n";
	print $tfh_R "text(b, colMeans(stuff), labels=signif(stuff,3), cex=1, pos=4, col=\"red\",offset=0,srt=90);"."\n";
	my $number = 0;
	while ($number <= ($num_clins - 1)) {
		my ($tfh_step,$tpath_step) = Genome::Sys->create_temp_file;
		unless($tfh_step) {
			$self->error_message("Unable to create temporary file $!");
			die;
		}
		$tpath_step =~ s/\:/\\\:/g;

		if ($number == 0) {
			print $tfh_step "Clinical Variable";
			foreach my $gene (sort keys %{$gene_hash{'total'}}) {
				print $tfh_step "\t$gene";
			}
			print $tfh_step "\n";
			foreach my $clinical (@clinicals) {
				print $tfh_step "$clinical";
				foreach my $gene (sort keys %{$gene_hash{'total'}}) {
					if (defined $significant_hash{$clinical}{$gene}) {
						print $tfh_step "\t$significant_hash{$clinical}{$gene}";
					}
					else {
						print $tfh_step "\t0";
					}
				}
				print $tfh_step "\n";
			}

			my $genelist = "c(\"".join("\",\"",@genes)."\")";
			print $tfh_R "siggenetable <- read.table(\"$tpath_step\", row.names = 1, header = TRUE, sep = \"\\t\");"."\n";
#			print $tfh_R "maxcount <- max(siggenetable);"."\n";
			my @sums;
			foreach (my $i = 1; $i <= ($num_clins - 1); $i++) {
				push(@sums,"sum(siggenetable[$i,])");
			}
			push(@sums,"sum(siggenetable[$num_clins,])");
			print $tfh_R "maxcount <- max(c(".join(",",@sums)."));"."\n";
			print $tfh_R "maxcount2 <- max(siggenetable);"."\n";
#			print $tfh_R 'table_labels <- labels(siggenetable);'."\n";
#			print $tfh_R 'genes <- table_labels[2];'."\n";
#			print $tfh_R 'categories <- table_labels[1];'."\n";
#			foreach (my $i; $i <= ($num_clins - 1); $i++) {
#				print $tfh_R "siggenes$i <- c(genes\$names[$i],\"\",\"\",\"\",\"\",\"\",\"\",\"\",\"\");"."\n";
#			}
#			print $tfh_R "siggenes_total <- c(siggenes,siggenes1,siggenes2,siggenes3,siggenes4,siggenes5,siggenes6,siggenes7,siggenes8);"."\n";
#			print $tfh_R 'colors <- c("darkred","red","orange","yellow","green","darkgreen","blue","darkblue","violet");'."\n";
			my $scale = 0.5;
			print $tfh_R "barplot(t(siggenetable),main=\"Gene Distribution (per clinical category)\",col=rainbow($num_genes),las=2,ylim=c(0,maxcount*1.25));"."\n";
			my $cex = (1 / $num_genes);
			if ($cex < 0.1) {
				$cex = 0.1;
#				print $tfh_R "legend(x=\"top\", title = \"Genes\", legend=$genelist,col=rainbow($num_genes),cex=$cex,pch=19,ncol=100);"."\n";
				print $tfh_R "legend(x=\"topright\", title = \"Genes\", legend=\"Individual Genes Too Many to Show\",col=\"Red\",cex=1,pch=19);"."\n";
			}
			elsif ($cex > 0.5) {
				$cex = 0.5;
				print $tfh_R "legend(x=\"top\", title = \"Genes\", legend=$genelist,col=rainbow($num_genes),cex=$cex,pch=19,ncol=10);"."\n";
			}
			else {
				print $tfh_R "legend(x=\"top\", title = \"Genes\", legend=$genelist,col=rainbow($num_genes),cex=$cex,pch=19,ncol=10);"."\n";
			}


			print $tfh_R "barplot(as.matrix(siggenetable),main=\"Gene Distribution per clinical category\",legend=TRUE,cex.names=$cex,col=rainbow($num_clins),border=rainbow($num_clins),las=2,beside=TRUE,ylim=c(0,maxcount2*1.5));"."\n";
			close($tfh_step);
		}
		($tfh_step,$tpath_step) = Genome::Sys->create_temp_file;
		unless($tfh_step) {
			$self->error_message("Unable to create temporary file $!");
			die;
		}
#if ($number == 0) {
#$tpath_step = '/gscuser/wschierd/Deleteme/poop.txt';
#open($tfh_step, ">$tpath_step") or die "Can't open output file: $!\n";
#}

		$tpath_step =~ s/\:/\\\:/g;
		my $clinical = $clinicals[$number];
		my $column = ($number +1);
		print $tfh_R "clinmean <- signif(mean(clindatatable\$$clinical,na.rm=TRUE),3);"."\n";
		print $tfh_R "clinsd <- signif(sd(clindatatable\$$clinical,na.rm=TRUE),3);"."\n";
		print $tfh_R "clinmin <- signif(min(clindatatable\$$clinical,na.rm=TRUE),3);"."\n";
		print $tfh_R "clinmax <- signif(max(clindatatable\$$clinical,na.rm=TRUE),3);"."\n";
		print $tfh_R "xmin <- (clinmean - (4 * clinsd));"."\n";
		print $tfh_R "xmax <- (clinmean + (4 * clinsd));"."\n";
		print $tfh_R "hist(clindatatable\$$clinical,main=\"$clinical Distribution (per clinical category)\",xlab=\"$clinical\",xlim=c(xmin,xmax));"."\n";
		print $tfh_R "mtext(paste(\"Mean:\", clinmean, \"Std Dev:\", clinsd, \"Min:\", clinmin, \"Max:\", clinmax, sep=\" \"), padj=-0.5);"."\n";
#		print $tfh_R "mtext(paste(\"Mean:\", clinmean, \"Std Dev:\", clinsd, sep=\"\"), padj=-0.5);"."\n";
		my $line = 1;
		print $tfh_step "#Clinical Phenotype: $clinical\n";
		print $tfh_step "Line\tChromosome\tStart\tStop\tLog10\tdbsnp\n";
		my @positions;
		my @tickmarks;
		my $last_chr=0;
		foreach my $chr (sort { $a <=> $b } keys %{$log10position{$clinical}}) {
			foreach my $start (sort { $a <=> $b } keys %{$log10position{$clinical}{$chr}}) {
				foreach my $stop (sort { $a <=> $b } keys %{$log10position{$clinical}{$chr}{$start}}) {
					my $var_search = "$chr"."_"."$start"."_"."$stop";
					print $tfh_step "$line\t$chr\t$start\t$stop\t$log10position{$clinical}{$chr}{$start}{$stop}\t$maf_stuff{$var_search}\n";
					if ($chr ne $last_chr) {
						push(@positions,$chr);
						push(@tickmarks,$line);
						$last_chr = $chr;
					}
					$line++;
				}
			}
		}
		my $position_string = "\"".join("\", \"",@positions)."\"";
		my $tickmarks_string = join(", ",@tickmarks);

#			my $pcutoff = 0;
#			my $fdrcutoff = 0;
#line graph, x is line number, y is log10
#		foreach my $current_p (sort { $a <=> $b } keys %{$fdr_cutoff{$clinical}}) {
#			if ($fdr_cutoff{$clinical}{$current_p} eq 'NA' || $fdr_cutoff{$clinical}{$current_p} >= $signif_value) {
#				next;
#			}
#			else {
#				$fdrcutoff = $fdr_cutoff{$clinical}{$current_p};
#				$pcutoff = $current_p;
#			}
#		}

		my @pvalues = (sort { $a <=> $b } keys %{$fdr_cutoff{$clinical}});
		my $log;
		if (scalar(@pvalues) > 0) {
			my @pvalues2 = @pvalues;
			my $pcutoff = pop(@pvalues2);
			my $fdrcutoff = $fdr_cutoff{$clinical}{$pcutoff};
			my $multiplier = ($signif_value / $fdrcutoff);
			my $pcutoff_adj = ($pcutoff * $multiplier);
			$pcutoff_adj = (($pcutoff_adj + $pcutoff + $pcutoff + $pcutoff) / 4);
#			my $multiplier2 = (0.10 / $fdrcutoff);
#			my $pcutoff_adj2 = ($pcutoff * $multiplier2);
##			$log = -log10($pcutoff_adj);
			$log = -log10($pcutoff_adj);
#			my $log2 = -log10($pcutoff_adj2);
#			if ($fdrcutoff == 0.000500) {
#				print "\n$pcutoff\n$pcutoff_adj\n$pcutoff_adj2\n$fdrcutoff\n$log\n$log2\n";
#				print "\n$pcutoff\n$fdrcutoff\n$log\n";
#			}
		}

		print $tfh_R "clinical <- read.table(\"$tpath_step\", header = TRUE, sep = \"\\t\");"."\n"; #row.names = 1, 
		print $tfh_R "plot.default(x=clinical\$Line,y=clinical\$Log10,xlab=\"Relative Chromosomal Position\",ylab=\"-log10(pvalue)\", main=paste(genome, \" \", \"$clinical\"),lty=1,lwd = 1,cex=0.4, col=\"#000000FF\",type=\"l\",xaxt=\"n\");"."\n";
		print $tfh_R "axis(1, at=c($tickmarks_string), labels=c($position_string), las=2, cex=0.1);"."\n";
		if (scalar(@pvalues) > 0) {
			print $tfh_R "z1=subset(clinical, clinical\$Log10 > $log - .001);"."\n";
#z2 <-switch(y, fruit = "banana", vegetable = "broccoli", "Neither")
#			print $tfh_R "z2=subset(clinical, clinical\$Log10 > $log & clinical\$dbsnp != \"novel\");"."\n";
			print $tfh_R "text(z1\$Line, z1\$Log10 - 0.3, labels=paste(signif(z1\$Log10,3),\"\\n\",\"dbsnp:\\n\",z1\$dbsnp,\"\\n\",\"Chr: \",z1\$Chromosome,\"\\n\",\"Position:\",\"\\n\",z1\$Start,sep=\"\"), cex=0.6, pos=4, col=\"red\",offset=0.1);"."\n";
			print $tfh_R "abline($log, 0, col=\"red\");"."\n";
			print $tfh_R "fdrlinevalue=signif($log,4);"."\n";
			print $tfh_R "text(-40, $log-.1, labels=paste(\"FDR\<$signif_value @ ~\",fdrlinevalue), cex=0.8, pos=4, col=\"red\");"."\n";
		}
#fdr=p.adjust(tt[,"p"],method="fdr");
		close $tfh_step;
#		}
		$number++;

	}
	print $tfh_R 'devoff <- dev.off();'."\n";
	print $tfh_R "q()\n";
	close $tfh_R;

	my $cmd = "R --vanilla --slave \< $temp_path_R";
	my $return = Genome::Sys->shellcmd(
           cmd => "$cmd",
           output_files => [$output_pdf_image_file],
           skip_if_output_is_present => 0,
	   );
	unless($return) { 
		$self->error_message("Failed to execute: Returned $return");
		die $self->error_message;
	}

	return $return;
}

sub log10 {
	my $n = shift;
	return log($n)/log(10);
}

1;
