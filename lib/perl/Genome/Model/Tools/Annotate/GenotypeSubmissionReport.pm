package Genome::Model::Tools::Annotate::GenotypeSubmissionReport;

use strict;
use warnings;
use Genome;
use IPC::Run;

class Genome::Model::Tools::Annotate::GenotypeSubmissionReport {

    is => 'Command',                       
    has => [ 
	gts => {
	    type  =>  'String',
	    doc   =>  "provide the genotype submission file, Alleles assumed to be oriented on the positive strand and coordinates are on the genomic scale",
	},
	organism => {
	    type  =>  'String',
	    doc   =>  "provide the organism either mouse or human; default is human",
	    is_optional  => 1,
	    default => 'human',
	},
	version => {
	    type  =>  'String',
	    doc   =>  "provide the imported annotation version; default for human is 54_36p_v2 and for mouse is 54_37g_v2",
	    is_optional  => 1,
	    default => '54_36p_v2',
	},
	order => {
	    type  =>  'String',
	    doc   =>  "provide the column number from 1 to 1+n for chromosome,coordinate,sample,allele1,allele2 if different than the default gts format;2,4,6,7,8",
	    is_optional  => 1,
	    default => '2,4,6,7,8',
	},
	delimiter => {
	    type  =>  'String',
	    doc   =>  "provide the column seperator if other than a tab",
	    is_optional  => 1,
	    default => "\t",
	},
	ref_id => {
	    type  =>  'String',
	    doc   =>  "provide the ref sample name or part of the name only contained with in the ref sample. Default will select as ref sample any sample ending in .c1",
	    is_optional  => 1,
	    default => '.c1',
	},
	ref_dir => {
	    type  =>  'String',
	    doc   =>  "provide the full path to were the genomic reference sequence sits. This will be used to verified the reference bases in the genotype submission file",
	    is_optional  => 1,
	    default => '/gscmnt/sata180/info/medseq/biodb/shared/Hs_build36_mask1c/',
	},
	root => {
	    type  =>  'String',
	    doc   =>  "provide the root name for the output default root name from gts file",
	    is_optional  => 1,
	},

	nosift => {
   	    type  =>  'Boolean',
	    doc   => 'use this option if you do not want to run sift and polyphen',
	    is_optional  => 1,
        },
	nodbsnp => {
    	    type  =>  'Boolean',
	     doc => 'use this option if you do not want to look for variant related dbsnp info',
	    is_optional  => 1,
        },
	polyphred => {
 	    type  =>  'String',
            doc  => 'enter the polyphred file if you would like to see read genotypes and scores from polyphred in the summary.txt',
	},
	polyscan => {
  	    type  =>  'String',
	    doc  => 'enter the polyscan file if you would like to see read genotypes and scores from polyscan in the summary.txt',
	},
	indels => {
   	    type  =>  'String',
            doc => 'enter the indel.txt file if you would like to see recorded read genotypes and scores from  the indel.txt',
	},
	refseq => {
      	    type  =>  'String',
	    doc  => 'enter the refseq file to convert coordinates for parsing polyphred and polyscan files',
	},
	ace => {
       	    type  =>  'String',
	    doc  => 'enter the consed ace file for the collection of assebmly stats for the summary.txt',
	},
	pretty_source_1 => {
            type  =>  'String',
            doc  => 'pretty-source-1 to set the sample start pos',
	    is_optional  => 1,
	    default => 1,
        },
        pretty_source_2 => {
            type  =>  'String',
            doc  => 'pretty-source-2 to set the sample stop pos default is 10 this value defines the end of the sample name in the population file',
	    is_optional  => 1,
	    default => 20,
        },
    ], 
};


sub help_brief {
    return <<EOS
	gmt annotate genotype-submission-report -gts
EOS
}

sub help_synopsis {
    return <<EOS
	gmt annotate genotype-submission-report -gts -ace -indels -polyphred -polyscan -refseq

EOS
}

sub help_detail { 
    return <<EOS 
	This tool was design to read in a genotype submmision file, identify all the variants, annotate them, check for dbsnps and screen the nonsynonymous snps with sift and polyphen. However, it will work on any list of variants with this info chromosome,coordinate,sample,allele1,allele2. allele1 and allele2 may be substituted with ref_allele and variant allele but they should be in the positive orientation with respect to the reference. If the order of these variants is different from the default order based on the genotype submission file format then use the -order option also see the delimiter option if your file is not tab delimited. ref_id is not a manditory requirement. 
EOS
}

sub execute {

    my $self = shift;

    my $gts = $self->gts;
    unless (-f $gts) {$self->error_message("unable to see the genotype submission file $gts");return;}

    my $root = $self->root;
    unless ($root) {
	if ($gts =~ /\.genotype/) {
	    $root = (split(/\.genotype/,$gts))[0];
	} else {
	    $root = $gts;
	}
    }

#######################################################
##Step 1  Extracting variants from the genotype submission file
   
    my ($genotypes,$gt_counts) = &parse_gts($self);
    
#######################################################
##step 2 generate a list of variants for annotation
    $self->status_message("annotate variants from the genotype submission file");
    my ($annotation,$dbsnp,$SIFT_Prediction,$PPH_Prediction) = &annotate_variant_alleles($self,$genotypes);
    
#######################################################
##step 3 generate the summary

    my ($report_prettybase_snp_counts,$report_prettybase_snp_summary) = &report_prettybase($self,$genotypes);

    my ($summary_txt,$annotation_out,$nav) = &prep_summary($self,$genotypes,$annotation,$dbsnp,$report_prettybase_snp_summary,$PPH_Prediction,$SIFT_Prediction);

    &write_annotation_out($self,$root,$annotation_out);

    my ($pretty_source_1,$pretty_source_2) = &pretty_source($self);

    my $indel_file = $self->indels;
    my $indels;
    if (-f $indel_file) {
	$self->status_message("parsing the indel.txt file");
	($indels) = &parse_indel_file($self,$indel_file,$pretty_source_1,$pretty_source_2);
    } else {
	$self->error_message("unable to see the indel.txt file $indel_file");
    }

    my $polyphred_file = $self->polyphred;
    my $polyphred_reads;
    if (-f $polyphred_file) {
	$self->status_message("parsing the polyphred file");
	($polyphred_reads) = &parse_polyphred($self,$polyphred_file,$pretty_source_1,$pretty_source_2);
    } else {
	$self->error_message("unable to see the polyphred file $polyphred_file");
    }

    my $polyscan_file = $self->polyscan;
    my $polyscan_reads;
    if (-f $polyscan_file) {
	$self->status_message("parsing the polyscan file");
	($polyscan_reads) = &parse_polyscan($self,$polyscan_file,$pretty_source_1,$pretty_source_2);
    } else {
	$self->error_message("unable to see the polyscan file $polyscan_file");
    }

    my $refseq_file = $self->refseq;
    my $refseq_header;
    if (-f $refseq_file) {
	$self->status_message("parsing the refseq fasta file");
	my ($refseq_info) = Genome::Model::Tools::RefSeq::Fasta->execute(refseq_fasta => $refseq_file, no_stdout => 1);
	$refseq_header = $refseq_info->result;
    } else {
	$self->error_message("unable to see the refseq file $refseq_file");
    }


    my $ace_file = $self->ace;
    my $ace_reference;
    if (-f $refseq_file) {
	$self->status_message("parsing the Ace file");
	my $ace_reference_info = Genome::Model::Tools::Consed::AceReference->execute(ace_file => $ace_file, assembly_stats => "1" , derive_sample_name => "1" , no_stdout => "1");
	$ace_reference = $ace_reference_info->result;
	unless ($ace_reference) {$self->error_message("unable to get the ace reference");}
    } else {
	$self->error_message("unable to see the Ace file $ace_file");
    }

    my ($autoreport) = &mk_auto_report($self,$root,$genotypes,$report_prettybase_snp_counts,$refseq_header,$ace_reference,$summary_txt,$annotation,$dbsnp,$SIFT_Prediction,$PPH_Prediction,$indels,$polyphred_reads,$polyscan_reads);

    my $summary_output_file = "$root.summary.txt";
    my $sum_fh;

    open($sum_fh,">$summary_output_file") || $self->error_message("unable to open the file $summary_output_file for writing will exit the program now") && return;
    $self->status_message("Preparing the summary.txt");

    my $ref_head = $refseq_header->{header};
    my $title_line = "This is a summary for the analysis of $ref_head\n";
    print $sum_fh qq($title_line\n);
    
    #get_assembly_info replaced by ace_reference
    
    my $reseqid = $ace_reference->{reseqid};
    my $total_reads_analyzed = $ace_reference->{assembly_read_count};
    
    print $sum_fh qq(1. Assembly Information\n);
    print $sum_fh join(" ", "Refseq ID", $reseqid, "\n");
    print $sum_fh join(" ",	"Total Number of Reads Analyzed ", $total_reads_analyzed, "\n");
    
    print $sum_fh "\n2. Amplicon Statistics\n";
    print $sum_fh "Amplicon_Coordinates, Amplicon_ID, Reads_Analyzed\n";
    
    foreach my $amplicon_name (sort keys %{$ace_reference->{amplicon}}) {
	my $amplicon_coordinates = $ace_reference->{amplicon}->{$amplicon_name}->{amplicon_coordinates};
	my $number_of_reads_in_amp = $ace_reference->{amplicon}->{$amplicon_name}->{reads_in_amp};
	
	print $sum_fh "$amplicon_coordinates, $amplicon_name, $number_of_reads_in_amp\n";
	
    }
    #print $sum_fh "Other $non_pcr_reads\n";
	
    print $sum_fh "\n3. Sample Statistics\n";
    print $sum_fh "Sample_ID, Reads_Analyzed\n";
    foreach my $sample (sort keys %{$ace_reference->{sample}}) {
	my $number_of_reads_in_assembly = $ace_reference->{sample}->{$sample};
	print $sum_fh "$sample, $number_of_reads_in_assembly\n";
    }
    
    print $sum_fh "\n\n\n########################################################################################################\n\n\n";
    
    print $sum_fh "\n\n\n4. SNP Information (Non-synonymous and Splice Site Substitutions)\n";
    print $sum_fh "---------------------------------------------------------------------------------------";
    print $sum_fh "\nRefseq Pos, Genomic Pos, Refseq Allele, dbSNP ID\n\tTranscript, Region, Alleles, Amino_Acid, Sift Prediction, PPH Prediction, Type\n\t\tSample Genotype\n\t\t\tRead Name, Polyscan Score, Polyphred Score, Read Position\n";
    print $sum_fh "---------------------------------------------------------------------------------------\n";
    
    my $type = "snp";
    #&get_read_summary($self,$sum_fh,$type,$refseq_header,$summary_txt,$annotation,$dbsnp,$SIFT_Prediction,$PPH_Prediction,$indels,$polyphred_reads,$polyscan_reads);
    &get_read_summary($self,$sum_fh,$type,$genotypes,$refseq_header,$summary_txt,$annotation,$dbsnp,$SIFT_Prediction,$PPH_Prediction,$indels,$polyphred_reads,$polyscan_reads);
    
    print $sum_fh "\n\n\n########################################################################################################\n\n\n";
    print $sum_fh "\n\n\n5. INDEL Information (Non-synonymous and Splice Site Substitutions)\n";
    print $sum_fh "---------------------------------------------------------------------------------------";
    print $sum_fh "\nRefseq Pos, Genomic Pos, Refseq Allele, dbSNP ID\n\tTranscript, Region, Alleles, Amino_Acid, Type\n\t\tSample Genotype\n\t\t\tRead Name, Polyscan Score, Polyphred Score, Read Position\n";
    print $sum_fh "---------------------------------------------------------------------------------------\n";
    
    $type = "indel";
    #&get_read_summary($self,$sum_fh,$type,$refseq_header,$summary_txt,$annotation,$dbsnp,$SIFT_Prediction,$PPH_Prediction,$indels,$polyphred_reads,$polyscan_reads);
    &get_read_summary($self,$sum_fh,$type,$genotypes,$refseq_header,$summary_txt,$annotation,$dbsnp,$SIFT_Prediction,$PPH_Prediction,$indels,$polyphred_reads,$polyscan_reads);
    print $sum_fh "\n\n\n########################################################################################################\n\n\n";
    
    print $sum_fh "\n6. Genotype Frequencies\n";
    
    open (SS, "$report_prettybase_snp_summary") || $self->error_message("couldn't open the snp summary");
    while (<SS>) {
	chomp;
	my $ss_line=$_;
	print $sum_fh "$ss_line\n";
    }
    close(SS);

    print $sum_fh "\n\n\n########################################################################################################\n\n\n";
   
    close $sum_fh;
    
    $self->status_message("The summary.txt has been written to $summary_output_file");
    
    my @command = ["rm" , "snp_summary.txt" , "snp_counts.txt"];
    &ipc_run(@command);
    
    $self->status_message("done");

    return 1;
    
}

sub prep_summary {
    my ($self,$genotypes,$annotation,$dbsnp,$report_prettybase_snp_summary,$PPH_Prediction,$SIFT_Prediction) = @_;

    my $summary_txt;
    my $annotation_out;
    my $nav;

    foreach my $chr (sort keys %{$genotypes}) {
	foreach my $pos (sort {$a<=>$b} keys %{$genotypes->{$chr}}) {
	    
	    my $ref_sample = $genotypes->{$chr}->{$pos}->{ref_sample};
	    unless ($ref_sample) {$ref_sample = "reference";$genotypes->{$chr}->{$pos}->{ref_sample}=$ref_sample;}
	    my $ref_base = $genotypes->{$chr}->{$pos}->{ref_base};
	    
	    foreach my $sample (sort keys %{$genotypes->{$chr}->{$pos}->{GT}}) {
		my $genotype = $genotypes->{$chr}->{$pos}->{GT}->{$sample};
		my $v_allele = $genotypes->{$chr}->{$pos}->{variant_sample}->{$sample};
		
		if ($v_allele) {
		    my @vas = split(/\:/,$v_allele);

		    use List::MoreUtils qw/ uniq /;
		    @vas = uniq @vas;
		    for my $variant_allele (@vas) {
			my $variant_indel_allele;
			my $end;
			if ($variant_allele =~ /\+/ || $variant_allele =~ /\-/) {
			    my ($ref,$var);
			    ($end,$ref,$var) = &indel_alleles($self,$variant_allele,$pos,$chr);
			    $ref_base = $ref;
			    $variant_indel_allele = $variant_allele;
			    $variant_allele = $var;
			}
			
			foreach my $transcript (sort keys %{$annotation}) {
			    my $aa = $annotation->{$transcript}->{$pos}->{$variant_allele}->{aa};
			    my ($polyphen,$sift);
			    if ($aa) {
				my ($prot) = $aa =~ /p\.(\S+)/;
				if ($prot) {
				    $polyphen = $PPH_Prediction->{$transcript}->{$prot};
				    $sift = $SIFT_Prediction->{$transcript}->{$prot};
				}
			    } else {
				$aa ="-";
			    }
			    
			    unless ($polyphen) {$polyphen="-";}
			    unless ($sift) {$sift="-";}
			    
			    my $chromosome = $annotation->{$transcript}->{$pos}->{$variant_allele}->{chr};
			    my $stop = $annotation->{$transcript}->{$pos}->{$variant_allele}->{stop};
			    my $ref = $annotation->{$transcript}->{$pos}->{$variant_allele}->{ref};
			    my $variant_type = $annotation->{$transcript}->{$pos}->{$variant_allele}->{variant_type};
			    my $gene = $annotation->{$transcript}->{$pos}->{$variant_allele}->{gene};
			    my $source = $annotation->{$transcript}->{$pos}->{$variant_allele}->{source};
			    my $tv = $annotation->{$transcript}->{$pos}->{$variant_allele}->{tv};
			    my $strand = $annotation->{$transcript}->{$pos}->{$variant_allele}->{strand};
			    my $Traans_stat = $annotation->{$transcript}->{$pos}->{$variant_allele}->{Traans_stat};
			    my $trv_type = $annotation->{$transcript}->{$pos}->{$variant_allele}->{trv_type};
			    my $c_pos = $annotation->{$transcript}->{$pos}->{$variant_allele}->{c_pos};
			    my $cons_score = $annotation->{$transcript}->{$pos}->{$variant_allele}->{cons_score};
			    my $domain = $annotation->{$transcript}->{$pos}->{$variant_allele}->{domain};
			    my $dbsnp_hit = $dbsnp->{$chr}->{$pos}->{$variant_allele};

			    my ($rs_id,$dbsnp_submittor,$dbsnp_alleles,$dbsnp_allele_match);
			    if ($dbsnp_hit) {
				if ($dbsnp_hit eq "no_hit") {
				    $rs_id="-";$dbsnp_submittor="-";$dbsnp_alleles="-";$dbsnp_allele_match="-";
				} else {
				    my $match;
				    ($rs_id,$dbsnp_submittor,$match) = split(/\,/,$dbsnp_hit);
				    ($dbsnp_alleles,$dbsnp_allele_match) = split(/\:/,$match);
				}
			    } else {
				$rs_id="-";$dbsnp_submittor="-";$dbsnp_alleles="-";$dbsnp_allele_match="-";
			    }

			    if ($chromosome && $stop && $ref && $variant_type && $gene && $source && $tv && $strand && $Traans_stat && $trv_type && $c_pos && $cons_score && $domain && $dbsnp_hit) {

				my $out = qq($chromosome\t$pos\t$stop\t$ref\t$variant_allele\t$sample\t$variant_type\t$genotype\t$source\t$gene\t$transcript\t$strand\t$Traans_stat\t$trv_type\t$c_pos\t$aa\t$polyphen\t$sift\t$cons_score\t$domain\t$rs_id\t$dbsnp_submittor\t$dbsnp_alleles\t$dbsnp_allele_match);
				
				$annotation_out->{$chr}->{$pos}->{$sample}->{$out}=1; ## To write annotation.out
				
				if ($trv_type =~ /nonsense/ ||
				    $trv_type =~ /missense/ ||
				    $trv_type =~ /splice_site/ ||
				    $trv_type =~ /splice_region/ ||
				    $trv_type =~ /nonstop/ ||
				    $trv_type =~ /cryptic_splice_site/ ||
				    $trv_type =~ /splice_site_del/ ||
				    $trv_type =~ /splice_site_ins/ ||
				    $trv_type =~ /splice_region_del/ ||
				    $trv_type =~ /splice_region_ins/ ||
				    $trv_type =~ /frame_shift_del/ ||
				    $trv_type =~ /in_frame_del/ ||
				    $trv_type =~ /frame_shift_ins/ ||
				    $trv_type =~ /in_frame_ins/)
				{
				    
####################################################
### Integrate method to produce a hash of pairs with an option input pairs file then we can get_somatic_status
### my ($tumor,$normal) = &get_paired($sample);  ####make new method to produce a hash of pairs with an option input pairs file
### ($somatic_status,$tgt,$ngt) = &get_somatic_status($ref,$t_gt,$n_gt);
### then generate the navigator to key off tumor only
####################################################

				    my $info = "$variant_type:$trv_type:$aa:$ref_base:$v_allele:$genotype";
				    $nav->{$chr}->{$pos}->{$sample}=$info;
				    
				    if ($variant_type =~ /INS/ || $variant_type =~ /DEL/) {
					    $summary_txt->{indel}->{$chr}->{$pos}->{$variant_allele}->{$transcript}->{$sample}=$out;
				    } else {
					$summary_txt->{snp}->{$chr}->{$pos}->{$variant_allele}->{$transcript}->{$sample}=$out;
				    }
				}
			    } else {
				my ($stop,$variant_type);
				if ($end) {
				    $stop = $end;
				    if ($variant_allele =~ /\-/) {
					$variant_type = "DEL";
				    } else {
					$variant_type = "INS";
				    }
				} else {
				    $stop = $pos;
				    $variant_type = "SNP";
				}
				
				my $out = qq($chr\t$pos\t$stop\t$ref_base\t$variant_allele\t$sample\t$variant_type\t$genotype\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t$rs_id\t$dbsnp_submittor\t$dbsnp_alleles\t$dbsnp_allele_match);
				$annotation_out->{$chr}->{$pos}->{$sample}->{$out}=1;
			    }
			}
		    }
		}
	    }
	}
    }
    
    return ($summary_txt,$annotation_out,$nav);
    
}


sub mk_nav {

    my ($self,$nav) = @_;
    open(NAV,">snpsAindels_ss_ns.navlist") || $self->error_message("couldn't write to the nav list") && return;
    $self->status_message("writing a snpsAindels_ss_ns.navlist");
    foreach my $chr (sort keys %{$nav}) {
	foreach my $pos (sort {$a<=>$b} keys %{$nav->{$chr}}) {
	    foreach my $tumor (sort keys %{$nav->{$chr}->{$pos}}) {
		my $info = $nav->{$chr}->{$pos}->{$tumor};
		my ($normal,$comment) = split(/\,/,$info);
		
		print NAV qq($tumor,$pos,$normal,$comment\n);  ###Need to modify this for none paired samples
	    }
	}
    }
    close NAV;
    return 1;
}

sub write_annotation_out {

    my $header = "chromosome\tpos\tstop\tref_allele\tvariant_allele\tsample_name\tvariant_type\tgenotype\tsource\tgene\ttranscript\tstrand\tTraans_stat\ttrv_type\tc_pos\taa\tpolyphen\tsift\tcons_score\tdomain\trs_id\tdbsnp_submittor\tdbsnp_alleles\tdbsnp_allele_match";

    my ($self,$root,$annotation_out) = @_;
    open(OUT,">$root.annotation.out") || $self->error_message("couldn't write to $root.annotation.out") && return;

    print OUT qq($header\n);
    foreach my $chr (sort keys %{$annotation_out}) {
	foreach my $pos (sort {$a<=>$b} keys %{$annotation_out->{$chr}}) {
	    foreach my $sample (sort keys %{$annotation_out->{$chr}->{$pos}}) {
		foreach my $out (sort keys %{$annotation_out->{$chr}->{$pos}->{$sample}}) {
		    print OUT qq($out\n);
		}
	    }
	}
    }
    close OUT;
    $self->status_message("your annotation results are ready in the file $root.annotation.out");
    return 1;
}



sub mk_auto_report {

    my ($self,$root,$genotypes,$report_prettybase_snp_counts,$refseq_header,$ace_reference,$summary_txt,$annotation,$dbsnp,$SIFT_Prediction,$PPH_Prediction,$indels,$polyphred_reads,$polyscan_reads) = @_;
    my $snp_counts;
    
    $self->status_message("collecting snp counts");

    if (-f  $report_prettybase_snp_counts) {
	open(SC,$report_prettybase_snp_counts) ||  $self->error_message("couldn't open the report_prettybase_snp_counts file") && return;;
	my (@coords,@AA,@AC,@AG,@AT,@CC,@CG,@CT,@GG,@GT,@TT,@star,@het_indel,@hom_indel,@NN,@XX,@Total);
	my (@CA,@GA,@GC,@TA,@TC,@TG);
	
	while (<SC>) {
	    chomp;
	    my $line = $_;
	    my @array = split(/\t/,$line);
	    
	    my $title = $array[0];
	    unless ($title =~ /\S/) {$title = "coords";}
	    if ($title eq "***") {$title = "star";}
	    
	    if ($title eq "coords") {@coords = @array;}
	    if ($title eq "AA") {@AA = @array;}
	    if ($title eq "AC") {@AC = @array;}
	    if ($title eq "AG") {@AG = @array;}
	    if ($title eq "AT") {@AT = @array;}
	    if ($title eq "CC") {@CC = @array;}
	    if ($title eq "CA") {@CA = @array;}
	    if ($title eq "CG") {@CG = @array;}
	    if ($title eq "CT") {@CT = @array;}
	    if ($title eq "GG") {@GG = @array;}
	    if ($title eq "GA") {@GA = @array;}
	    if ($title eq "GC") {@GC = @array;}
	    if ($title eq "GT") {@GT = @array;}
	    if ($title eq "TT") {@TT = @array;}
	    if ($title eq "TA") {@TA = @array;}
	    if ($title eq "TC") {@TC = @array;}
	    if ($title eq "TG") {@TG = @array;}
	    if ($title eq "star") {@star = @array;}
	    if ($title eq "het_indel") {@het_indel = @array;}
	    if ($title eq "hom_indel") {@hom_indel = @array;}
	    if ($title eq "NN") {@NN = @array;}
	    if ($title eq "XX") {@XX = @array;}
	    if ($title eq "Total") {@Total = @array;}
	    
	}

	my $n = @coords;
	my $y = $n - 1;
	for my $x (1..$y) {
	    my $coord = $coords[$x];$snp_counts->{$coord}->{coords}=$coord;
	    if (@AA) {my $AA = $AA[$x];$snp_counts->{$coord}->{AA}=$AA;}
	    if (@AC) {my $AC = $AC[$x];$snp_counts->{$coord}->{AC}=$AC;}
	    if (@AG) {my $AG = $AG[$x];$snp_counts->{$coord}->{AG}=$AG;}
	    if (@AT) {my $AT = $AT[$x];$snp_counts->{$coord}->{AT}=$AT;}
	    if (@CC) {my $CC = $CC[$x];$snp_counts->{$coord}->{CC}=$CC;}
	    if (@CG) {my $CG = $CG[$x];$snp_counts->{$coord}->{CG}=$CG;}
	    if (@GT) {my $GT = $GT[$x];$snp_counts->{$coord}->{GT}=$GT;}
	    if (@TT) {my $TT = $TT[$x];$snp_counts->{$coord}->{TT}=$TT;}
	    if (@CA) {my $CA = $CA[$x];$snp_counts->{$coord}->{CA}=$CA;}
	    if (@CT) {my $CT = $CT[$x];$snp_counts->{$coord}->{CT}=$CT;}
	    if (@CC) {my $GG = $GG[$x];$snp_counts->{$coord}->{GG}=$GG;}
	    if (@GA) {my $GA = $GA[$x];$snp_counts->{$coord}->{GA}=$GA;}
	    if (@GC) {my $GC = $GC[$x];$snp_counts->{$coord}->{GC}=$GC;}
	    if (@GT) {my $GT = $GT[$x];$snp_counts->{$coord}->{GT}=$GT;}
	    if (@TA) {my $TA = $TA[$x];$snp_counts->{$coord}->{TA}=$TA;}
	    if (@TC) {my $TC = $TC[$x];$snp_counts->{$coord}->{TC}=$TC;}
	    if (@TG) {my $TG = $TG[$x];$snp_counts->{$coord}->{TG}=$TG;}
	    if (@star) {my $star = $star[$x];$snp_counts->{$coord}->{star}=$star;}
	    if (@het_indel) {my $het_indel = $het_indel[$x];$snp_counts->{$coord}->{het_indel}=$het_indel;}
	    if (@hom_indel) {my $hom_indel = $hom_indel[$x];$snp_counts->{$coord}->{hom_indel}=$hom_indel;}
	    if (@NN) {my $NN = $NN[$x];$snp_counts->{$coord}->{NN}=$NN;}
	    if (@XX) {my $XX = $XX[$x];$snp_counts->{$coord}->{XX}=$XX;}
	    if (@Total) {my $Total = $Total[$x];$snp_counts->{$coord}->{Total}=$Total;}
	}
    } else {
	$self->error_message("no snp counts");
    }

   
    open(AR,">$root.autoreport.tsv") || $self->error_message("couldn't write to $root.autoreport.tsv") && return;

    $self->status_message("preparing the auto report");
    
    my $arrays;
    my $sample_genotypes;
    
    foreach my $transcript (sort keys %{$annotation}) {
	my (@ref_range,@range_a,@ref_a,@variant_allele_a,@variant_type_a,@gene_a,@source_a,@tv_a,@strand_a,@Traans_stat_a,@trv_type_a,@c_pos_a,@aa_a,@cons_score_a,@domain_a,@dbsnp_a,@sift_a,@pph_a);
	
	foreach my $chr (sort keys %{$genotypes}) {
	    foreach my $pos (sort {$a<=>$b} keys %{$genotypes->{$chr}}) {
		my $start = $pos;
		my $range = $pos;
		
		foreach my $variant_allele (sort keys %{$genotypes->{$chr}->{$pos}->{variant_allele}}) {
		    if ($variant_allele =~ /\+/ || $variant_allele =~ /\-/) {
			my ($stop,$ref,$var) = &indel_alleles($self,$variant_allele,$pos,$chr);
			$variant_allele = $var;
		    }
		    my $var = $variant_allele;

		    my $stop = $annotation->{$transcript}->{$start}->{$var}->{stop};
		    
		    unless ($stop) {$stop = $pos;}
		    unless ($stop == $pos) {
			$range = "$start-$stop";
		    }
		    my $ref_startpos = &ref_coord($refseq_header,$start);
		    my $ref_range = $ref_startpos;
		    my $ref_stoppos = &ref_coord($refseq_header,$stop);
		    unless ($ref_startpos == $ref_stoppos) {
			$ref_range = "$ref_startpos-$ref_stoppos";
		    }

		    push (@ref_range,$ref_range);
		    push (@range_a,$range);
		    
		    my $ref = $annotation->{$transcript}->{$start}->{$var}->{ref};
		    unless($ref){$ref=$genotypes->{$chr}->{$pos}->{ref_base};}
		    push (@ref_a,$ref);
		    
		    push(@variant_allele_a,$variant_allele);
		    
		    my $variant_type = $annotation->{$transcript}->{$start}->{$var}->{variant_type};
		    unless($variant_type){$variant_type="-";}
		    push(@variant_type_a,$variant_type);
		    
		    my $gene = $annotation->{$transcript}->{$start}->{$var}->{gene};
		    unless($gene){$gene="-";}
		    push(@gene_a,$gene);
		    
		    my $source = $annotation->{$transcript}->{$start}->{$var}->{source};
		    unless($source){$source="-";}
		    push(@source_a,$source);
		    
		    my $tv = $annotation->{$transcript}->{$start}->{$var}->{tv};
		    unless($tv){$tv="-";}
		    push(@tv_a,$tv);
		    
		    my $strand = $annotation->{$transcript}->{$start}->{$var}->{strand};
		    unless($strand){$strand="-";}
		    push(@strand_a,$strand);
		    
		    my $Traans_stat = $annotation->{$transcript}->{$start}->{$var}->{Traans_stat};
		    unless($Traans_stat){$Traans_stat="-";}
		    push(@Traans_stat_a,$Traans_stat);
		    
		    my $trv_type = $annotation->{$transcript}->{$start}->{$var}->{trv_type};
		    unless($trv_type){$trv_type="-";}
		    push(@trv_type_a,$trv_type);
		    
		    my $c_pos = $annotation->{$transcript}->{$start}->{$var}->{c_pos};
		    unless($c_pos){$c_pos="-";}
		    push(@c_pos_a,$c_pos);
		    
		    my $aa = $annotation->{$transcript}->{$start}->{$var}->{aa};
		    unless($aa){$aa="-";}
		    push(@aa_a,$aa);
		    
		    $aa =~ s/p\.//;
		    my $sift = $SIFT_Prediction->{$transcript}->{$aa};
		    unless($sift) {$sift="-";}
		    push(@sift_a,$sift);
		    
		    my $pph = $PPH_Prediction->{$transcript}->{$aa};
		    unless($pph) {$pph="-";}
		    push(@pph_a,$pph);
		    
		    my $cons_score = $annotation->{$transcript}->{$start}->{$var}->{cons_score};
		    unless($cons_score){$cons_score="-";}
		    push(@cons_score_a,$cons_score);
		    
		    my $domain = $annotation->{$transcript}->{$start}->{$var}->{domain};
		    unless($domain){$domain="-";}
		    push(@domain_a,$domain);
		    
		    my $dbsnp = $dbsnp->{$chr}->{$start}->{$var}; #$annotation->{$transcript}->{$start}->{$var}->{dbsnp};
		    unless ($dbsnp) {$dbsnp = "-";}
		    push(@dbsnp_a,$dbsnp);

		    $sample_genotypes->{ref}->{$pos}->{$variant_allele}=$ref;

		    foreach my $sample (sort keys %{$genotypes->{$chr}->{$pos}->{GT}}) {
			my $gt = $genotypes->{$chr}->{$pos}->{GT}->{$sample};
			$sample_genotypes->{$sample}->{$pos}->{$variant_allele}=$gt;
			
		    }
		}
	    }
	}
	
	my $ref_range = join ("\t", @ref_range);
	$arrays->{$transcript}->{ref_range} = $ref_range;
	my $range = join ("\t", @range_a);
	$arrays->{$transcript}->{range} = $range;
	my $ref = join ("\t", @ref_a);
	$arrays->{$transcript}->{ref} = $ref;
	my $variant_allele = join ("\t", @variant_allele_a);
	$arrays->{$transcript}->{variant_allele} = $variant_allele;
	my $variant_type = join ("\t", @variant_type_a);
	$arrays->{$transcript}->{variant_type}=$variant_type;
	my $gene = join ("\t", @gene_a);
	$arrays->{$transcript}->{gene}=$gene;
	my $source = join ("\t", @source_a);
	$arrays->{$transcript}->{source}=$source;
	my $tv = join ("\t", @tv_a);
	$arrays->{$transcript}->{tv}=$tv;
	my $strand = join ("\t", @strand_a);
	$arrays->{$transcript}->{strand}=$strand;
	my $Traans_stat = join ("\t", @Traans_stat_a);
	$arrays->{$transcript}->{Traans_stat}=$Traans_stat;
	my $trv_type = join ("\t", @trv_type_a);
	$arrays->{$transcript}->{trv_type}=$trv_type;
	my $c_pos = join ("\t", @c_pos_a);
	$arrays->{$transcript}->{c_pos}=$c_pos;
	my $aa = join ("\t", @aa_a);
	$arrays->{$transcript}->{aa}=$aa;
	my $cons_score = join ("\t", @cons_score_a);
	$arrays->{$transcript}->{cons_score}=$cons_score;
	my $domain = join ("\t", @domain_a);
	$arrays->{$transcript}->{domain}=$domain;
	my $dbsnp = join ("\t", @dbsnp_a);
	$arrays->{$transcript}->{dbsnp}=$dbsnp;
	
	my $sift = join ("\t", @sift_a);
	$arrays->{$transcript}->{sift}=$sift;
	my $pph = join ("\t", @pph_a);
	$arrays->{$transcript}->{pph}=$pph;
	
    }

    foreach my $transcript (sort keys %{$arrays}) {

	my $ref_range = $arrays->{$transcript}->{ref_range};
	my $range = $arrays->{$transcript}->{range};
	my $ref = $arrays->{$transcript}->{ref};
	my $variant_allele = $arrays->{$transcript}->{variant_allele};
	my $variant_type = $arrays->{$transcript}->{variant_type};
	my $gene = $arrays->{$transcript}->{gene};
	my $source = $arrays->{$transcript}->{source};
	my $tv = $arrays->{$transcript}->{tv};
	my $strand = $arrays->{$transcript}->{strand};
	my $Traans_stat = $arrays->{$transcript}->{Traans_stat};
	my $trv_type = $arrays->{$transcript}->{trv_type};
	my $c_pos = $arrays->{$transcript}->{c_pos};
	my $aa = $arrays->{$transcript}->{aa};
	my $cons_score = $arrays->{$transcript}->{cons_score};
	my $domain = $arrays->{$transcript}->{domain};
	my $dbsnp = $arrays->{$transcript}->{dbsnp};
	my $sift = $arrays->{$transcript}->{sift};
	my $pph = $arrays->{$transcript}->{pph};
	print AR qq(\n\n$transcript\n);
	print AR qq(ref_range\t$ref_range\nrange\t$range\nref_allele\t$ref\nvariant_allele\t$variant_allele\nvariant_type\t$variant_type\n);
	print AR qq(gene\t$gene\nsource\t$source\ntv\t$tv\nstrand\t$strand\nTraans_stat\t$Traans_stat\n);
	print AR qq(trv_type\t$trv_type\nc_pos\t$c_pos\naa\t$aa\nsift\t$sift\npolyphen\t$pph\ncons_score\t$cons_score\ndomain\t$domain\ndbsnp\t$dbsnp\n);
	
    }
    
    print AR qq(\n);
    print AR qq(\n);
    
    my (@AA,@AC,@AG,@AT,@CC,@CG,@CT,@GG,@GT,@TT,@star,@het_indel,@hom_indel,@NN,@XX,@Total);
    my (@CA,@GA,@GC,@TA,@TC,@TG);

    my $reseqid = $ace_reference->{reseqid};
    print AR qq($reseqid);

    foreach my $pos (sort {$a<=>$b} keys %{$sample_genotypes->{ref}}) {
	foreach my $variant_allele (sort keys %{$sample_genotypes->{ref}->{$pos}}) {
	    
	    my $AA = $snp_counts->{$pos}->{AA};unless($AA){$AA=0;}push(@AA,$AA);
	    my $AC = $snp_counts->{$pos}->{AC};unless($AC){$AC=0;}push(@AC,$AC);
	    my $AG = $snp_counts->{$pos}->{AG};unless($AG){$AG=0;}push(@AG,$AG);
	    my $AT = $snp_counts->{$pos}->{AT};unless($AT){$AT=0;}push(@AT,$AT);
	    my $CC = $snp_counts->{$pos}->{CC};unless($CC){$CC=0;}push(@CC,$CC);
	    my $CA = $snp_counts->{$pos}->{CA};unless($CA){$CA=0;}push(@CA,$CA);
	    my $CG = $snp_counts->{$pos}->{CG};unless($CG){$CG=0;}push(@CG,$CG);
	    my $CT = $snp_counts->{$pos}->{CT};unless($CT){$CT=0;}push(@CT,$CT);
	    my $GG = $snp_counts->{$pos}->{GG};unless($GG){$GG=0;}push(@GG,$GG);
	    my $GA = $snp_counts->{$pos}->{GA};unless($GA){$GA=0;}push(@GA,$GA);
	    my $GC = $snp_counts->{$pos}->{GC};unless($GC){$GC=0;}push(@GC,$GC);
	    my $GT = $snp_counts->{$pos}->{GT};unless($GT){$GT=0;}push(@GT,$GT);
	    my $TT = $snp_counts->{$pos}->{TT};unless($TT){$TT=0;}push(@TT,$TT);
	    my $TA = $snp_counts->{$pos}->{TA};unless($TA){$TA=0;}push(@TA,$TA);
	    my $TC = $snp_counts->{$pos}->{TC};unless($TC){$TC=0;}push(@TC,$TC);
	    my $TG = $snp_counts->{$pos}->{TG};unless($TG){$TG=0;}push(@TG,$TG);
	    my $star = $snp_counts->{$pos}->{star};unless($star){$star=0;}push(@star,$star);
	    my $het_indel = $snp_counts->{$pos}->{het_indel};unless($het_indel){$het_indel=0;}push(@het_indel,$het_indel);
	    my $hom_indel = $snp_counts->{$pos}->{hom_indel};unless($hom_indel){$hom_indel=0;}push(@hom_indel,$hom_indel);
	    my $NN = $snp_counts->{$pos}->{NN};unless($NN){$NN=0;}push(@NN,$NN);
	    my $XX = $snp_counts->{$pos}->{XX};unless($XX){$XX=0;}push(@XX,$XX);
	    my $Total = $snp_counts->{$pos}->{Total};unless($Total){$Total=0;}push(@Total,$Total);
	    
    
	    my ($ref_allele) = split(/\s/,$sample_genotypes->{ref}->{$pos}->{$variant_allele});
	    print AR qq(\t$ref_allele);
	}
    }

    print AR qq(\n);
    
    foreach my $sample (sort keys %{$sample_genotypes}) {
	unless ($sample  eq "ref") {
	    print AR qq($sample);
	    foreach my $pos (sort {$a<=>$b} keys %{$sample_genotypes->{$sample}}) {
		foreach my $variant_allele (sort keys %{$sample_genotypes->{$sample}->{$pos}}) {
		    my $genotype = $sample_genotypes->{$sample}->{$pos}->{$variant_allele};
		    print AR qq(\t$genotype);
		}
	    }
	    print AR qq(\n);
	}
    }

    print AR qq(\n);
    
    my $AA = join ("\t", @AA);print AR qq(AA\t$AA\n);
    my $AC = join ("\t", @AC);print AR qq(AC\t$AC\n);
    my $AG = join ("\t", @AG);print AR qq(AG\t$AG\n);
    my $AT = join ("\t", @AT);print AR qq(AT\t$AT\n);
    
    my $CC = join ("\t", @CC);print AR qq(CC\t$CC\n);
    my $CA = join ("\t", @CA);print AR qq(CA\t$CA\n);
    my $CT = join ("\t", @CT);print AR qq(CT\t$CT\n);
    my $CG = join ("\t", @CG);print AR qq(CG\t$CG\n);
    
    my $GG = join ("\t", @GG);print AR qq(GG\t$GG\n);
    my $GA = join ("\t", @GA);print AR qq(GA\t$GA\n);
    my $GC = join ("\t", @GC);print AR qq(GC\t$GC\n);
    my $GT = join ("\t", @GT);print AR qq(GT\t$GT\n);
    
    my $TT = join ("\t", @TT);print AR qq(TT\t$TT\n);
    my $TA = join ("\t", @TA);print AR qq(TA\t$TA\n);
    my $TC = join ("\t", @TC);print AR qq(TC\t$TC\n);
    my $TG = join ("\t", @TG);print AR qq(TG\t$TG\n);
    
    my $star = join ("\t", @star);print AR qq(***\t$star\n);
    my $het_indel = join ("\t", @het_indel);print AR qq(het_indel\t$het_indel\n);
    my $hom_indel = join ("\t", @hom_indel);print AR qq(hom_indel\t$hom_indel\n);
    my $NN = join ("\t", @NN);print AR qq(NN\t$NN\n);
    my $XX = join ("\t", @XX);print AR qq(XX\t$XX\n);
    my $Total = join ("\t", @Total);print AR qq(Total\t$Total\n);
    print AR qq(\n);
    
    close(AR);
    
    my $unique_tutto = @Total;
    my $raw_no_ss = ($unique_tutto / "249");
    my $no_ss = ((split(/\./,$raw_no_ss))[0] + 1);
    my $report_number = 1;
    my $start_array = 1;
    my $end_array = 248;
    my $last_end_array=$unique_tutto;
    if ($no_ss > 1) {
	while ($no_ss > 0) {
	    open(NEWREPORT, ">$root.autoreport_$report_number.tsv");
	    open(AR,"$root.autoreport.tsv") || die "couldn't open the autoreport\n)";
	    while (<AR>) {
		chomp;
		my $line = $_;
		my $info = (split(/\t/, $line))[0];
		
		if ($info) {
		} else {
		    $info = ' ';
		}
		
		my @array = (split(/\t/,$line))[$start_array..$end_array];
		my $newreport =join("\t" ,@array);
		
		print NEWREPORT ("$info\t$newreport\n");
		
	    }
	    close(NEWREPORT);
	    close(AR);
	    $start_array=($end_array+1);
	    my $next_end_array=($start_array+248);
	    if ($next_end_array < $last_end_array) {
		$end_array=$next_end_array;
	    } else {
		$end_array=$last_end_array;
	    }
	    $no_ss--;
	    $report_number++;
	}
    }
    
    if (-f "$root.autoreport_1.tsv") {
	my $tsv_s = `ls $root.autoreport*.tsv`;
	chomp $tsv_s;
	$self->status_message("View results in $tsv_s");
    } else {
	$self->status_message("View results in $root.autoreport.tsv");
    }

    return unless (-f "$root.autoreport.tsv");
    return 1;   
}

###$genotypes->{$chr}->{$pos}->{ref_base}
###$genotypes->{$chr}->{$pos}->{ref_sample}
###$genotypes->{$chr}->{$pos}->{variant_allele}}->{$v}=1;
###$genotypes->{$chr}->{$pos}->{variant_sample}->{$sample}=$av;$av = join':',@variant;
###$genotypes->{$chr}->{$pos}->{GT}->{$sample}=$gt;$gt = "hom_indel";$gt = "het_indel";my $gt = join( ':',sort ($allele1,$allele2) );
###$genotypes->{$chr}->{$pos}->{line}->{$sample}=$line;
   
#######clean up
    
sub get_paired {

    my ($sample) = @_;
    my ($prefix,$id) = $sample =~ /(\S+)\-(\S+)/;

    my ($tumor,$normal);
    if ($id =~ /(\S+)t$/) {
	$tumor = $sample;
	my $pair = $prefix . "-" . "$1" . "n";
	if ($id =~ /^\D(\d\d\d)\D/) {
	    $pair = $prefix . "-" . "0" . "$1" . "n";
	}
	$normal = $pair;
	
    } elsif ($id =~ /(\S+)n$/) {
	$normal = $sample;
	my $pair = $prefix . "-" . "$1" . "t";
	$tumor = $pair;
    }
    
    unless ($tumor) {$tumor = $sample; $normal = $sample;}
    return ($tumor,$normal);
}

sub get_somatic_status {

    my ($ref,$t_gt,$n_gt) = @_;

    my ($theho,$tgt);
    if ($t_gt) {
	($theho,$tgt) = &get_heho($ref,$t_gt);
    }

    my ($nheho,$ngt);
    if ($n_gt) {
	($nheho,$ngt) = &get_heho($ref,$n_gt);
    }

    unless($tgt) {$tgt="NTG";}
    unless($ngt) {$ngt ="NNG";}
    unless ($theho) {$theho="UND";}
    unless ($nheho) {$nheho="UND";}
    
    my $nxx;
    if ($theho eq "XX") {$theho="UND";}
    if ($theho eq "NN") {$theho="UND"; }
    if ($nheho eq "XX") {$nheho="UND";$nxx="xx"; }
    if ($nheho eq "NN") {$nheho="UND"; }
    
    
    #my ($theho,$nheho) = @_;
    
    my $tri_a;
    if ($theho eq "Tri_allelic") { $tri_a=1; }
    
    unless (($theho eq "WT") || ($theho eq "HR") || ($theho eq "HET")) {$theho="UND";}
    unless (($nheho eq "WT") || ($nheho eq "HR") || ($nheho eq "HET")) {$nheho="UND";}
    
    my $somatic_status;
    if ($theho eq $nheho) {
	if ($theho eq "UND") {
	    $somatic_status = "UK";
	}
	if ($theho eq "WT") {
	    $somatic_status = "WT";
	}
	if ($theho eq "HET") {
	    $somatic_status = "G";
	    unless ($tgt eq $ngt) {
		$somatic_status = "O";
	    }
	}
	if ($theho eq "HR") {
	    $somatic_status = "G";
	    unless ($tgt eq $ngt) {
		$somatic_status = "O";
	    }
	}
    } else {
	if ($nheho eq "UND") {
	    $somatic_status = "V";
	    if ($theho eq "WT") {
		$somatic_status = "UK";
	    }
	}
	
	if ($theho eq "UND") {
	    if ($nheho eq "WT") {
		$somatic_status = "UK";
	    } else {
		$somatic_status = "V";
	    }
	}
	
	unless (($theho eq "UND") || ($nheho eq "UND")) {
	    
	    if (($nheho eq "WT") && ($theho eq "HET")) {
		$somatic_status = "S";
	    }
	    if (($nheho eq "WT") && ($theho eq "HR")) {
		$somatic_status = "S";
	    }
	    if (($nheho eq "HET") && ($theho eq "WT")) {
		$somatic_status = "LOH";
	    }
	    if (($nheho eq "HET") && ($theho eq "HR")) {
		$somatic_status = "LOH";
	    }
	    if (($nheho eq "HR") && ($theho eq "WT")) {
		$somatic_status = "O";
	    }
	    if (($nheho eq "HR") && ($theho eq "HET")) {
		$somatic_status = "O";
	    }
	}
    }
    unless ($somatic_status) {$somatic_status= "UP";}
    
    
    #if (($theho eq "WT") && $nxx) { $filter_nxxtwt->{$pos}->{$sample}=1; }
    
    
    if ($tri_a) {
	my $status = "$somatic_status" . "_Tri_allelic";
	
	return ($status,$tgt,$ngt);
	
    } else {
	return ($somatic_status,$tgt,$ngt);
    }
}

sub get_heho {
    
    my ($ref,$gt) = @_;
    my $Reference_Allele = $ref;
    
    my ($a1,$a2) = $gt =~ /^(\S+)\:(\S+)$/;
    my $heho;
    
    unless ($a1 && $a2) {
	$heho = "No_gt";
	return ($heho,$gt);
    }
    
    my $sgt = join( '',sort ($a1,$a2) );
    
    if ($a1 eq $a2) {
	if ($a1 eq $Reference_Allele) {
	    $heho = "WT";
	} else {
	    $heho = "HR";
	}
	if ($sgt =~ /XX/) {
	    $heho = "XX";
	}
	elsif ($sgt =~ /NN/) {
	    $heho = "NN";
	}
	
    } else {
	unless (($a1 eq $Reference_Allele) || ($a2 eq $Reference_Allele)) { $heho = "Tri_allelic"; return($heho,$sgt);}
	$heho = "HET";
    }
    
    return($heho,$sgt);
    
}


sub ipc_run {

    my (@command) = @_;

    my ($in, $out, $err);
    my ($obj) = IPC::Run::run(@command, \$in, \$out, \$err);
    if ($err) {
#	print qq($err\n);
    }
    if ($out) {
	return ($out);
	#print qq($out\n);
    }
}

sub parse_indel_file { 

    my $indels;
    my ($self,$indel_file,$pretty_source_1,$pretty_source_2) = @_;

    unless (-f $indel_file) {return;}
    open(ID,$indel_file) || $self->error_message("couldn't read the indel_file") && return;
    while (<ID>) {
	chomp;
	my $line = $_;
	my ($contig,$trace,$range,$indel,$hetho,$size,$bases,$score) = (split(/\s+/,$line))[0,1,2,3,4,5,6,7];
	my $sample = substr($trace,$pretty_source_1,$pretty_source_2);
	my ($end,$start);
	if ($indel eq "ins") {
	    ($end,$start) = (split(/\-/,$range))[0,1];
	} else {
	    ($start,$end) = (split(/\-/,$range))[0,1];
	}
	$indels->{$sample}->{$trace}->{$start}->{hetho} = $hetho;
	$indels->{$sample}->{$trace}->{$start}->{indel} = $indel;
	$indels->{$sample}->{$trace}->{$start}->{range} = $range;
	$indels->{$sample}->{$trace}->{$start}->{end} = $end;
	$indels->{$sample}->{$trace}->{$start}->{size} = $size;
	$indels->{$sample}->{$trace}->{$start}->{bases} = $bases;
	$indels->{$sample}->{$trace}->{$start}->{score} = $score;
    }
    close ID;
    return ($indels);
}
    
sub parse_polyphred {

    my $polyphred_reads;

    my ($self,$polyphred_file,$pretty_source_1,$pretty_source_2) = @_;
    unless (-f $polyphred_file) { return; } 
    my $parse_snps = "off";
    my $parse_pfindels = "off";
    my $parse_mgt = "off";
    
    open (PF, "$polyphred_file") || $self->error_message("couldn't read the polyphred_file") && return;
    while (<PF>) {
	chomp;
	my $line=$_;
	if (($line =~ /BEGIN_GENOTYPE/) || ($line =~ /BEGIN_COLUMNGENOTYPE/)) { #BEGIN_INDELGENOTYPE BEGIN_MANUALGENOTYPE
	    $parse_snps = "on";
	}
	if ($line =~ /BEGIN_INDELGENOTYPE/) {
	    $parse_pfindels = "on";
	}
	if ($line =~ /BEGIN_MANUALGENOTYPE/) {
	    $parse_mgt = "on";
	}
	if ($line =~ /END/) {
	    $parse_snps = "off";
	    $parse_mgt = "off";
	    $parse_pfindels = "off";
	}
	
	if ($parse_snps eq "on") {
	    unless ($line =~ /BEGIN/) {
		my ($pos, $trace_pos, $read, $allele_1, $allele_2, $pf_score, $ps_score) = (split(/\s+/,$line))[0, 2, 3, 4, 5, 6, 7];
		unless ($ps_score) {$ps_score = 0;}
		my $gt = "$allele_1\:$allele_2";
		my $sample = substr($read,$pretty_source_1,$pretty_source_2);
		$polyphred_reads->{$pos}->{$sample}->{$read}->{ps_score}=$ps_score;
		$polyphred_reads->{$pos}->{$sample}->{$read}->{pf_score}=$pf_score;
		$polyphred_reads->{$pos}->{$sample}->{$read}->{gt}=$gt;
		$polyphred_reads->{$pos}->{$sample}->{$read}->{trace_pos}=$trace_pos;
	    }
	}
	
	if ($parse_pfindels eq "on") {
	    unless ($line =~ /BEGIN/) {
		my ($pos,$read) = (split(/[\s]+/,$line))[0,3];
		my $sample = substr($read,$pretty_source_1,$pretty_source_2);
		$polyphred_reads->{$pos}->{$sample}->{$read}->{pf_indel}=$line;
	    }
	}
	
	if ($parse_mgt eq "on") {
	    unless ($line =~ /BEGIN/) {
		my ($pos,$trace_pos,$read,$gt) = (split(/[\s]+/,$line))[0,2,3,4];
		my $sample = substr($read,$pretty_source_1,$pretty_source_2);
		$polyphred_reads->{$pos}->{$sample}->{$read}->{mgt}=$gt;
	    }
	}
    }
    close (PF);
    return ($polyphred_reads);
}


sub parse_polyscan {

    my $polyscan_reads;

    my ($self,$polyscan_file,$pretty_source_1,$pretty_source_2) = @_;
    unless (-f $polyscan_file) { return; } 

    my $snp_section = "off";
    my $indel_section = "off";
    
    open (PS, "$polyscan_file") || $self->error_message("couldn't read the polyscan_file") && return;;
    while (<PS>) {
	chomp;
	my $line=$_;
	if ( $line =~ /BEGIN_SNP/ ) {
	    $snp_section = "on";
	} 
	if ( $line =~ /BEGIN_INDEL/ ) {
	    $indel_section = "on";
	} 
	if ( $line =~ /END/ ) {
	    $indel_section = "off";
	    $snp_section = "off";
	}
	
	if ($snp_section eq "on") {
	    unless ( $line =~ /BEGIN/) {
		my ($pos,$trace_pos,$read,$gt,$ps_score) = split(/[\s]+/,$line);
		my $sample = substr($read,$pretty_source_1,$pretty_source_2);
		$polyscan_reads->{$pos}->{$sample}->{$read}->{ps_score}=$ps_score;
		$polyscan_reads->{$pos}->{$sample}->{$read}->{gt}=$gt;
		$polyscan_reads->{$pos}->{$sample}->{$read}->{trace_pos}=$trace_pos;
	    }
	}
	
	if ($indel_section eq "on") {
	    unless ( $line =~ /BEGIN/) {
		my ($pos,$trace_pos,$read,$size,$type,$gt,$ps_score,$peak_diff) = split(/[\s]+/,$line);
		my $sample = substr($read,$pretty_source_1,$pretty_source_2);
		$polyscan_reads->{$pos}->{$sample}->{$read}->{ps_score}=$ps_score;
		$polyscan_reads->{$pos}->{$sample}->{$read}->{gt}=$gt;
		$polyscan_reads->{$pos}->{$sample}->{$read}->{trace_pos}=$trace_pos;
		
	    }
	}
    }
    close(PS);
    return ($polyscan_reads);
}

sub ref_coord {
    
    my ($ref_head,$pos) = @_;

    my $orientation = $ref_head->{orientation};
    my $genomic_coord = $ref_head->{genomic_coord};

    my $ref_pos;
    if ($orientation eq "plus") {
	$ref_pos = $pos - $genomic_coord;
    } elsif ($orientation eq "minus") {
	$ref_pos = $genomic_coord - $pos;
    }
    return $ref_pos;
}

sub pretty_source {

    my ($self) = @_;
    
    my $pretty_source_1 = $self->pretty_source_1;
    my $pretty_source_2 = $self->pretty_source_2;

    if ($pretty_source_1) {
	my $d = $pretty_source_1;
	$pretty_source_1 = ($d - 1);
    } else {
	$pretty_source_1 = 0;
    }
    if ($pretty_source_2) {
	if ($pretty_source_1 != 0) {
	    my $d = $pretty_source_2;
	    $pretty_source_2 = ( $d - $pretty_source_1 );
	}
    } else {
	$pretty_source_2 = 10;
	if ($pretty_source_1 != 1) {
	    my $d = $pretty_source_2;
	    $pretty_source_2 = $d - $pretty_source_1;
	}
    }
    return ($pretty_source_1,$pretty_source_2);   
}


#####################################################################################################################


sub get_read_summary {

    #my ($type) = @_;
    my ($self,$sum_fh,$type,$genotypes,$refseq_header,$summary_txt,$annotation,$dbsnp,$SIFT_Prediction,$PPH_Prediction,$indels,$polyphred_reads,$polyscan_reads) = @_;

#->{snp}
#->{$type}
#$type is either snp or indel;
    
    foreach my $chr (sort keys %{$summary_txt->{$type}}) {
	print $sum_fh qq($type\'s on Chromosome $chr\n);
	foreach my $pos (sort {$a<=>$b} keys %{$summary_txt->{$type}->{$chr}}) {
	    my ($ref_pos) = &ref_coord($refseq_header,$pos);
	    
	    my $ref = $genotypes->{$chr}->{$pos}->{ref_base};
	    foreach my $variant_allele (sort keys %{$summary_txt->{$type}->{$chr}->{$pos}}) {
		my $dbsnp_info = $dbsnp->{$chr}->{$pos}->{$variant_allele};
		if ($dbsnp_info) {
		    print $sum_fh qq($ref_pos $pos $ref $dbsnp_info\n);
		} else {
		    print $sum_fh qq($ref_pos $pos $ref\n);
		}
		
		my @ta;
		foreach my $transcript (sort keys %{$summary_txt->{$type}->{$chr}->{$pos}->{$variant_allele}}) {
		    push(@ta,$transcript);
		    #my $chr = $annotation->{$transcript}->{$pos}->{$variant_allele}->{chr}; #=$chr;
		    my $stop = $annotation->{$transcript}->{$pos}->{$variant_allele}->{stop}; #=$stop;
		    my $dba_ref = $annotation->{$transcript}->{$pos}->{$variant_allele}->{ref}; #=$ref;
		    my $variant_type = $annotation->{$transcript}->{$pos}->{$variant_allele}->{variant_type}; #=$variant_type;
		    my $gene = $annotation->{$transcript}->{$pos}->{$variant_allele}->{gene}; #=$gene;
		    my $source = $annotation->{$transcript}->{$pos}->{$variant_allele}->{source}; #=$source;
		    my $tv = $annotation->{$transcript}->{$pos}->{$variant_allele}->{tv}; #=$tv;
		    my $strand = $annotation->{$transcript}->{$pos}->{$variant_allele}->{strand}; #=$strand;
		    my $Traans_stat = $annotation->{$transcript}->{$pos}->{$variant_allele}->{Traans_stat}; #=$Traans_stat;
		    my $trv_type = $annotation->{$transcript}->{$pos}->{$variant_allele}->{trv_type}; #=$trv_type;
		    my $c_pos = $annotation->{$transcript}->{$pos}->{$variant_allele}->{c_pos}; #=$c_pos;
		    my $aa = $annotation->{$transcript}->{$pos}->{$variant_allele}->{aa}; #=$aa;
		    my $cons_score = $annotation->{$transcript}->{$pos}->{$variant_allele}->{cons_score}; #=$cons_score;
		    my $domain = $annotation->{$transcript}->{$pos}->{$variant_allele}->{domain}; #=$domain;
		    
		    if ($aa =~ /p\./) {
			my ($prot) = $aa =~ /p\.(\S+)/;
			my $polyphen = $PPH_Prediction->{$transcript}->{$prot};
			my $sift = $SIFT_Prediction->{$transcript}->{$prot};
			if ($polyphen || $sift) {
			    unless ($sift) {$sift = "-";}
			    unless ($polyphen) {$polyphen = "-";}
			    print $sum_fh qq(\t$transcript $gene $trv_type $aa \($polyphen\) \($sift\) $cons_score $domain\n);
			} else {
			    print $sum_fh qq(\t$transcript $gene $trv_type $cons_score $domain\n);
			}
		    } else {
			print $sum_fh qq(\t$transcript $gene $trv_type $cons_score\n);
		    }
		}

		my $printed_sample;
		for my $trans (@ta) {

		    foreach my $sample (sort keys %{$summary_txt->{$type}->{$chr}->{$pos}->{$variant_allele}->{$trans}}) {
			my $out = $summary_txt->{$type}->{$chr}->{$pos}->{$variant_allele}->{$trans}->{$sample};
			if ($out) {
			    unless ($printed_sample->{$pos}->{$sample}) {
				
				my $gt = $genotypes->{$chr}->{$pos}->{GT}->{$sample};
				print $sum_fh qq(\t\t$sample $gt\n);
				$printed_sample->{$pos}->{$sample} = 1;
				
				foreach my $read (sort keys %{$polyphred_reads->{$ref_pos}->{$sample}}) {
				    
				    my $ps_score = $polyphred_reads->{$ref_pos}->{$sample}->{$read}->{ps_score};
				    unless ($ps_score) {$ps_score = $polyscan_reads->{$ref_pos}->{$sample}->{$read}->{ps_score};}
				    unless ($ps_score) {$ps_score = "-";}
				    my $pf_score = $polyphred_reads->{$ref_pos}->{$sample}->{$read}->{pf_score};
				    unless ($pf_score) {$pf_score = "-";}
				    my $gt = $polyphred_reads->{$ref_pos}->{$sample}->{$read}->{gt};
				    unless ($gt) {$gt = "-";}
				    my $psgt = $polyscan_reads->{$ref_pos}->{$sample}->{$read}->{gt};
				    unless ($psgt) {$psgt = "-";}
				    my $trace_pos = $polyphred_reads->{$ref_pos}->{$sample}->{$read}->{trace_pos};
				    unless ($trace_pos) {$trace_pos = "-";}
				    my $pf_indel = $polyphred_reads->{$ref_pos}->{$sample}->{$read}->{pf_indel};
				    unless ($pf_indel) {$pf_indel = "-";}
				    my $mgt = $polyphred_reads->{$ref_pos}->{$sample}->{$read}->{mgt};
				    unless ($mgt) {$mgt = "-";}

				    if ($type eq "indel") {
					my ($idps_gt,$idps_score) = &check_indel_score($sample,$read,$ref_pos,$ref,$indels);
					my $pf_indel = $polyphred_reads->{$pos}->{$sample}->{$read}->{pf_indel};
					#$polyphred_reads->{$pos}->{$sample}->{$read}->{mgt}
					if ($pf_indel) { $gt = "indel"; }
					if ($idps_gt && $idps_score) {
					    $psgt = $idps_gt;
					    $ps_score = $idps_score;
					}
				    }
				    #Read Name, Polyscan Score, Polyphred Score, Read Position
				    #H_GP-0001nPCR0007157_046a.b1 (AA 99) (AC 80)  119
				    
				    print $sum_fh qq(\t\t\t$read \($psgt $ps_score\) \($gt $pf_score\) $trace_pos\n);
				    
				}
				foreach my $read (sort keys %{$polyscan_reads->{$ref_pos}->{$sample}}) {
				    unless ($polyphred_reads->{$ref_pos}->{$sample}->{$read}) {
					my $ps_score = $polyscan_reads->{$ref_pos}->{$sample}->{$read}->{ps_score};
					unless ($ps_score) {$ps_score = "-";}
					my $gt = $polyscan_reads->{$ref_pos}->{$sample}->{$read}->{gt};
					unless ($gt) {$gt = "-";}
					my $trace_pos = $polyscan_reads->{$ref_pos}->{$sample}->{$read}->{trace_pos};
					unless ($trace_pos) {$trace_pos = "-";}
					my $pfgt = "-";
					my $pf_score = "-";
					
					my ($idps_gt,$idps_score);
					if ($type eq "indel") {
					    my ($idps_gt,$idps_score) = &check_indel_score($sample,$read,$ref_pos,$ref,$indels);
					    my $pf_indel = $polyphred_reads->{$pos}->{$sample}->{$read}->{pf_indel};
					    #$polyphred_reads->{$pos}->{$sample}->{$read}->{mgt}
					    if ($pf_indel) { $pfgt = "indel"; }
					    if ($idps_gt && $idps_score) {
						$gt = $idps_gt;
						$ps_score = $idps_score;
					    }
					}
					print $sum_fh qq(\t\t\t$read \($gt $ps_score\) \($pfgt $pf_score\) $trace_pos\n);
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
    return 1;
}

sub check_indel_score {

    my ($sample,$read,$ref_pos,$ref,$indels) = @_;

    my ($idps_gt);
    my $hetho = $indels->{$sample}->{$read}->{$ref_pos}->{hetho};# = $hetho;
    my $indel = $indels->{$sample}->{$read}->{$ref_pos}->{indel};#  = $indel;
    #$indels->{$sample}->{$trace}->{$start}->{range} = $range;
    #$indels->{$sample}->{$trace}->{$start}->{end} = $end;
    #$indels->{$sample}->{$trace}->{$start}->{size} = $size;
    my $bases = $indels->{$sample}->{$read}->{$ref_pos}->{bases};#  = $bases;
    my $idps_score = $indels->{$sample}->{$read}->{$ref_pos}->{score};#  = $score;
    #$ref
    if ($hetho) {
	if ($hetho eq "het") {
	    if ($indel eq "del") {
		$idps_gt = "$bases" . ":" . "-" . $bases;
	    } else {
		$idps_gt = "$ref" . ":" . "+" . $bases;
	    }
	} else {
	    if ($indel eq "del") {
		$idps_gt = "-" . $bases . ":" . "-" . $bases;
	    } else {
		$idps_gt = "+" . $bases . ":" . "+" . $bases;
	    }
	}
    }
    
    #my $pf_indel = $polyphred_reads->{$ref_pos}->{$sample}->{$read}->{pf_indel};
    #$polyphred_reads->{$pos}->{$sample}->{$read}->{mgt}
    #if ($pf_indel) { $pfgt = "indel"; }
    #if ($idps_gt && $idps_score) {
#	$gt = $idps_gt;
#	$ps_score = $idps_score;
   # }
    return ($idps_gt,$idps_score);   
}

#################################

sub parse_gts {

    my ($self) = @_;

    my $gts = $self->gts;
    my $delimiter = $self->delimiter;
    my $order = $self->order;
    my ($chr_n,$pos_n,$sample_n,$allele1_n,$allele2_n) = split(/\,/,$order);

    my $organism = $self->organism;
    my $ref_id = $self->ref_id;

    my ($genotypes,$gt_counts);
    open(GTS,$gts) || $self->error_message("could't open the genotype submission file") && return; 
    $self->status_message("Extracting variants from the genotype submission file");

    while (<GTS>) {
	chomp;
	my $line = $_;
	my ($chr,$pos,$sample,$allele1,$allele2) = (split(/$delimiter/,$line))[$chr_n - 1,$pos_n - 1,$sample_n - 1,$allele1_n - 1,$allele2_n - 1];

	$chr =~ s/C//;

	my $ref_base = $genotypes->{$chr}->{$pos}->{ref_base};
	unless ($ref_base) {
	    $ref_base = &get_ref_base($self,$pos,$pos,$chr);
	    $genotypes->{$chr}->{$pos}->{ref_base} = $ref_base;
	}

	if ($ref_id eq ".c1" && $sample =~ /.c1$/) {
	    $genotypes->{$chr}->{$pos}->{ref_sample}=$sample;
	} elsif ($sample =~ /$ref_id/) {
	    $genotypes->{$chr}->{$pos}->{ref_sample}=$sample;
	} else {
	    
	    $allele1 =~ s/\'//;
	    $allele2 =~ s/\'//;
	    
	    my @variant;
	    my ($va1,$va2);
	    if ($allele1 =~ /[ACGTX\-\+]/) {
		unless ($allele1 eq $ref_base) {
		    if ($allele1 =~ /[ACGTX]/) {
			push(@variant,$allele1);
			$va1 = $allele1;
		    }
		}
	    }
	    if ($allele2 =~ /[ACGTX\-\+]/) {
		unless ($allele2 eq $ref_base) {
		    if ($allele2 =~ /[ACGTX]/) {
			push (@variant,$allele2);
			$va2 = $allele2;
		    }
		}
	    }
	    $genotypes->{$chr}->{$pos}->{line}->{$sample}=$line;
	    my $gt = join( ':',sort ($allele1,$allele2) );
	    $genotypes->{$chr}->{$pos}->{GT}->{$sample}=$gt;

	    if ($va1) {
		if ($va1 =~ /\-/ || $va1 =~ /\+/) {
		    $gt = "het_indel";
		}
	    }
	    if ($va2) {
		if ($va2 =~ /\-/ || $va2 =~ /\+/) {
		    if ($gt eq "het_indel") {
			$gt = "hom_indel";
		    } else {
			$gt = "het_indel";
		    }
		}
	    }
	    $gt_counts->{$gt}->{$chr}->{$pos}++;
	    if (@variant) {
		my $av = join':',@variant;
		for my $v (@variant) {
		    $genotypes->{$chr}->{$pos}->{variant_allele}->{$v}=1;
		    $genotypes->{$chr}->{$pos}->{variant_sample}->{$sample}=$av;
		}
	    }
	}
    }
    close (GTS);
    return ($genotypes,$gt_counts);
}

###$genotypes->{$chr}->{$pos}->{ref_base}
###$genotypes->{$chr}->{$pos}->{ref_sample}
###$genotypes->{$chr}->{$pos}->{variant_allele}}->{$v}=1;
###$genotypes->{$chr}->{$pos}->{variant_sample}->{$sample}=$av;$av = join':',@variant;
###$genotypes->{$chr}->{$pos}->{GT}->{$sample}=$gt;$gt = "hom_indel";$gt = "het_indel";my $gt = join( ':',sort ($allele1,$allele2) );
###$genotypes->{$chr}->{$pos}->{line}->{$sample}=$line;

sub annotate_variant_alleles {
    
    my ($self,$genotypes) = @_;
    
    my $version = $self->version;
    my $organism = $self->organism;
    
    if ($organism eq "mouse" && $version eq "54_36p_v2") { $version = "54_37g_v2"; }
    
    my $list = "GTS.annotation.list";
    open(LIST,">$list") || $self->error_message("couldn't write the annotation list") && return;
    foreach my $chr (sort keys %{$genotypes}) {
	foreach my $pos (sort {$a<=>$b} keys %{$genotypes->{$chr}}) {

	    my $ref_base = $genotypes->{$chr}->{$pos}->{ref_base};
	    foreach my $variant_allele (sort keys %{$genotypes->{$chr}->{$pos}->{variant_allele}}) {
		
		if ($variant_allele =~ /\+/ || $variant_allele =~ /\-/) {
		    my ($stop,$ref,$var) = &indel_alleles($self,$variant_allele,$pos,$chr);
		    
		    unless ($chr && $pos && $stop && $ref && $var) {print qq($chr $pos $ref_base allele =>$variant_allele\n);}
		    
		    print LIST qq($chr\t$pos\t$stop\t$ref\t$var\n);
		    
		} else {

		    unless ($chr && $pos && $ref_base && $variant_allele) {
			print qq($chr $pos $ref_base allele =>$variant_allele\n);
		    }
		    print LIST qq($chr\t$pos\t$pos\t$ref_base\t$variant_allele\n);
		}
	    }
	}
    }
    close (LIST);
    my $annotated_list = "GTS.annotated.list";
    
    #print qq(retrieving annotation from the dw for the variant list\n);
    $self->status_message("running annotate transcript-variants");

    my @command = ["gmt" , "annotate" , "transcript-variants" , "--variant-file" , "$list" , "--output-file" , "$annotated_list" , "--reference-transcripts" , "NCBI-$organism.combined-annotation/$version", "--annotation-filter" , "none"];

    &ipc_run(@command);

    my $dbsnp_out = "GTS.annotated.list.dbsnp";
    $self->status_message("running annotate lookup-variants");
    @command = ["gmt" , "annotate" , "lookup-variants" , "-variant-file" , "$list" , "--output-file" , "$dbsnp_out", "--report-mode" , "full"];

    &ipc_run(@command);

    my ($annotation,$dbsnp,$sift_vars) = &parse_annotation($self,$annotated_list,$dbsnp_out);
    my ($SIFT_Prediction,$PPH_Prediction);
    if ($sift_vars) {
	($SIFT_Prediction,$PPH_Prediction) = &run_sift_and_poylphen($self,$sift_vars);
    }

    @command = ["rm" , "GTS.annotation.list" , "GTS.annotated.list" , "GTS.annotated.list.dbsnp" , "SIFT_LIST" , "POLYPHEN_LIST"];
    &ipc_run(@command);

    return ($annotation,$dbsnp,$SIFT_Prediction,$PPH_Prediction);
}

sub indel_alleles {

    my ($self,$allele,$pos,$chr) = @_;

    my ($ref,$var,$end);
    if ($allele =~ /\+([\S]+)/) { ##INS
	$var = $1;
	$ref = "-";
	$end = $pos + 1;
    } elsif ($allele =~ /\-([\S]+)/) { ##DEL
	my $bases = $1;
	my $length = length($bases); 
	$end = ($pos - 1) + $length;
	$ref = &get_ref_base($self,$pos,$end,$chr);
	$var = "-";
    }
    return ($end,$ref,$var);
}

sub get_ref_base {

    my ($self,$chr_start,$chr_stop,$chr_name) = @_;
    my $organism = $self->organism;
    my $RefDir = $self->ref_dir;

    use Bio::DB::Fasta;

    if ($organism eq "mouse" && $RefDir eq "/gscmnt/sata180/info/medseq/biodb/shared/Hs_build36_mask1c/") {
	$RefDir = "/gscmnt/sata147/info/medseq/rmeyer/resources/MouseB37/";
    }

    my $refdb = Bio::DB::Fasta->new($RefDir);

    my $seq = $refdb->seq($chr_name, $chr_start => $chr_stop);
    $seq =~ s/([\S]+)/\U$1/;

    return $seq;
}


sub parse_annotation {
    
    my ($self,$annotation_file,$dbsnp_out) = @_;
    my $dbsnp;
    
    if ($dbsnp_out && -e $dbsnp_out) {
	#print qq(reading the dbsnp hit results\n);
	open(DBSNP,"$dbsnp_out") || $self->error_message("couldn't open the dbsnp results");
	while (<DBSNP>) {
	    chomp;
	    my $line = $_;
	    my ($chr,$start,$stop,$ref_a,$var_a,$info) = split(/[\s]+/,$line);
	    
	    my $pi = $dbsnp->{$chr}->{$start}->{$var_a};
	    if ($pi) {
		unless ($pi =~ /$info/) {
		    $dbsnp->{$chr}->{$start}->{$var_a}="$pi\:\:$info";
		}
	    } else {
		$dbsnp->{$chr}->{$start}->{$var_a}=$info;
	    }
	}
	close (DBSNP);
    } else {
	$self->error_message("no dbsnps query result");
    }

    my $sift_vars;
    my $annotation;
    open(ANO,"$annotation_file") || $self->error_message("couldn't open the annotated variants list");
    while (<ANO>) {
	chomp;
	my $line = $_;
	my ($chromosome,$start,$stop,$ref,$var,$variant_type,$gene,$transcript,$transcript_species,$source,$tv,$strand,$Traans_stat,$trv_type,$c_pos,$aa,$cons_score,$domain) = split(/[\s]+/,$line); ##get extended annotation
#chromosome_name	start	stop	reference	variant	type	gene_name	transcript_name	transcript_species	transcript_source	transcript_version	strand	transcript_status	trv_type	c_position	amino_acid_change	ucsc_cons	domain	all_domains	flank_annotation_distance_to_transcript	intron_annotation_substructure_ordinal	intron_annotation_substructure_size	intron_annotation_substructure_position
	
	unless ($chromosome eq "chromosome_name") {
	    
	    if ($transcript =~ /\d+/) {
		
		$annotation->{$transcript}->{$start}->{$var}->{chr}=$chromosome;
		$annotation->{$transcript}->{$start}->{$var}->{stop}=$stop;
		$annotation->{$transcript}->{$start}->{$var}->{ref}=$ref;
		$annotation->{$transcript}->{$start}->{$var}->{variant_type}=$variant_type;
		$annotation->{$transcript}->{$start}->{$var}->{gene}=$gene;
		$annotation->{$transcript}->{$start}->{$var}->{source}=$source;
		$annotation->{$transcript}->{$start}->{$var}->{tv}=$tv;
		$annotation->{$transcript}->{$start}->{$var}->{strand}=$strand;
		$annotation->{$transcript}->{$start}->{$var}->{Traans_stat}=$Traans_stat;
		$annotation->{$transcript}->{$start}->{$var}->{trv_type}=$trv_type;
		$annotation->{$transcript}->{$start}->{$var}->{c_pos}=$c_pos;
		$annotation->{$transcript}->{$start}->{$var}->{aa}=$aa;
		$annotation->{$transcript}->{$start}->{$var}->{cons_score}=$cons_score;
		$annotation->{$transcript}->{$start}->{$var}->{domain}=$domain;
		
		
		if (($variant_type eq "SNP") && ($aa =~ /p\.[A-Z][\d]+[A-Z]/)) {
		    my ($var) = $aa =~ /p\.([A-Z][\d]+[A-Z])/;
		    my ($n) = $aa =~ /p\.[A-Z]([\d]+)[A-Z]/;
		    
		    unless ($var =~ /[B|J|O|U|X|Z]/ ) {
			$sift_vars->{$transcript}->{$n}->{$var}="$chromosome,$aa";
		    }
		}
	    }
	}
    } close (ANO);
    
    return ($annotation,$dbsnp,$sift_vars);
}



sub run_sift_and_poylphen {

    my ($self,$sift_vars) = @_;

    my $version = $self->version;
    my $organism = $self->organism;

    if ($organism eq "mouse" && $version eq "54_36p_v2") { $version = "54_37g_v2"; }

    my $SIFT_Prediction;
    my $PPH_Prediction;
    
    foreach my $transcript (sort keys %{$sift_vars}) {
	open(OUT,">$transcript.nonsynonymous.snps.list") || $self->error_message("couldn't write to $transcript.nonsynonymous.snps.list") && next;
	
	my $chr;
	foreach my $n (sort {$a<=>$b} keys %{$sift_vars->{$transcript}}) {
	    foreach my $var (sort keys %{$sift_vars->{$transcript}->{$n}}) {
		print OUT qq($var\n);
	    }
	}
	close OUT;
	
	next unless (-f "$transcript.nonsynonymous.snps.list");

	my @command = ["rm" , "$transcript.protein.time.error"];
	&ipc_run(@command);

	@command = ["gmt" , "snp" , "screen-nonsynonymous-snps" , "--list" , "$transcript.nonsynonymous.snps.list" , "--transcript" , "$transcript" , "--organism" , "$organism" , "--version" , "$version"];
	
	$self->status_message("will now screen nonsynonymous snps in $transcript");

	&ipc_run(@command);
	
	my $sift_results_file = "$transcript.protein.SIFTprediction";
	
	if (-f $sift_results_file) {
	    open(SIFT,$sift_results_file) || $self->error_message("couldn't read $sift_results_file");
	    
	    #print qq(SIFT results are in $sift_resutls_file\n);
	    while (<SIFT>) {
		chomp;
		my $line = $_;
		my ($prot,$td,$n1,$n2,$n3,$n4) = split(/\s+/,$line);
		
		$SIFT_Prediction->{$transcript}->{$prot} = "$td $n1 $n2 $n3 $n4";
	    }
	    close SIFT;
	} else {
	    $self->error_message("no SIFT results were found for $transcript");
	}
	
	
	my $polyphen_results_file = "$transcript.protein.PPHprediction";
	if (-f $polyphen_results_file) {
	    
	    #print qq(polypheb results are in $polyphen_results_file\n);
	    
	    open(PPH,$polyphen_results_file) || $self->error_message("couldn't read $polyphen_results_file");
	    while (<PPH>) {
		chomp;
		my $line = $_;
		my ($snp_id,$status,$protein_id,$position,$aa1,$aa2,$prediction,$basis,$effect,$site,$region,$phat,$score_delta,$num_observ,$num_struct_init,$num_struct_filt,$pdb_id,$res_num,$chain_id,$ali_ide,$ali_len,$acc_normed,$sec_str,$map_region,$delta_volume,$delta_prop,$b_fact,$num_h_bonds,$het_cont_ave_num,$het_cont_min_dist,$inter_cont_ave_num,$inter_cont_min_dist,$sites_cont_ave_num,$sites_cont_min_dist) = split(/\t/,$line);
		my $prot = "$aa1$position$aa2";
		$PPH_Prediction->{$transcript}->{$prot} = "$prediction $score_delta";
	    }
	    close PPH;
	} else {
	    $self->error_message("no polyphen results were found for $transcript");
	}
	
	#SIFT_Prediction "R132H DEoLETERIOUS 0.00 3.07 210 281"
	#PPH_Prediction "R132H probably damaging ? "
	
	@command = ["rm" , "$transcript.nonsynonymous.snps.list" , "$transcript.fasta" , "SIFT_$transcript.protein.fasta.out" , "$transcript.protein.PPHprediction" , "$transcript.protein.SIFTprediction" , "$transcript.protein.SIFTprediction.error" , "$transcript.protein.alignedfasta" , "$transcript.protein.clumped" , "$transcript.protein.clumped.consensuskey" , "$transcript.protein.clumped.error" , "$transcript.protein.fasta" , "$transcript.protein.fasta.out" , "$transcript.protein.fasta.query" , "$transcript.protein.fasta.query.globalX" , "$transcript.protein.fasta.query.globalX.error" , "$transcript.protein.fasta.query.out" , "$transcript.protein.fasta.query.selectedclumped" , "$transcript.protein.fasta.query.selectedclumped.error" , "$transcript.protein.fasta.query.selectedclumped.log" , "$transcript.protein.fasta.query.unfiltered" , "$transcript.protein.time.error"];
	
	&ipc_run(@command);
	
    }
    
    return($SIFT_Prediction,$PPH_Prediction);
}

sub report_prettybase {

    my ($self,$genotypes) = @_;

###$genotypes->{$chr}->{$pos}->{ref_base}
###$genotypes->{$chr}->{$pos}->{ref_sample}
###$genotypes->{$chr}->{$pos}->{variant_allele}}->{$v}=1;
###$genotypes->{$chr}->{$pos}->{variant_sample}->{$sample}=$av;$av = join':',@variant;
###$genotypes->{$chr}->{$pos}->{GT}->{$sample}=$gt;$gt = "hom_indel";$gt = "het_indel";my $gt = join( ':',sort ($allele1,$allele2) );
###$genotypes->{$chr}->{$pos}->{line}->{$sample}=$line;

    open(GPB,">genomic_prettybase.txt") || $self->error_message("couldn't write the genomic_prettybase.txt file") && return;
    $self->status_message("writting the genomic_prettybase");
    foreach my $chr (sort keys %{$genotypes}) {
	foreach my $pos (sort {$a<=>$b} keys %{$genotypes->{$chr}}) {
	    my $variant_location;
	    foreach my $variant_allele (sort keys %{$genotypes->{$chr}->{$pos}->{variant_allele}}) {
		$variant_location = 1;
	    }
	    next unless $variant_location;

	    my $ref_base = $genotypes->{$chr}->{$pos}->{ref_base};
	    my $ref_sample = $genotypes->{$chr}->{$pos}->{ref_sample};

	    print GPB qq($pos\t$ref_sample\t$ref_base\t$ref_base\n);

	    foreach my $sample (sort keys %{$genotypes->{$chr}->{$pos}->{GT}}) {
		unless ($sample eq $ref_sample) {
		    my $gt = $genotypes->{$chr}->{$pos}->{GT}->{$sample};

		    my ($a1,$a2) = split(/\:/,$gt);

		    if ($a1) {
			if ($a1 =~ /\-/ || $a1 =~ /\+/) {
			    $gt = "het_indel";
			}
		    }
		    if ($a2) {
			if ($a2 =~ /\-/ || $a2 =~ /\+/) {
			    if ($gt eq "het_indel") {
				$gt = "hom_indel";
			    } else {
				$gt = "het_indel";
			    }
			}
		    }

		    #my ($a1,$a2);
		    if ($gt eq "hom_indel") {
			$a1 = "-";
			$a2 = "-";
		    } elsif ($gt eq "het_indel") {
			$a1 = "$ref_base";
			$a2 = "-";
		    } #else {
			#($a1,$a2) = split(/\:/,$gt);
		    #}
		    print GPB qq($pos\t$sample\t$a1\t$a2\n);
		}
	    }
	}
    }
    close GPB;

    $self->status_message("getting the prettybase counts");
    my @command = ["report_msa_prettybase" , "genomic_prettybase.txt"];
    &ipc_run(@command);

    my $report_prettybase_snp_counts = "snp_counts.txt";
    my $report_prettybase_snp_summary = "snp_summary.txt";

    @command = ["rm" , "genomic_prettybase.txt" , "snp_report.txt" , "snp_ratios.txt" , "snp_macro.bas"];
    &ipc_run(@command);

    unless (-f $report_prettybase_snp_counts && -f $report_prettybase_snp_summary) {$self->error_message("couldn't find the output from report prettybase");}
    return($report_prettybase_snp_counts,$report_prettybase_snp_summary);

}

sub get_polyscan_positions {
    
    my $polyscan_positions;

    my ($self,$polyscanfile) = @_;
    
    use MG::IO::Polyscan;
    use MG::IO::Polyscan::Contig;
    
    my $polyscan=new MG::IO::Polyscan();
    $polyscan->readIn($polyscanfile);  #read in original polyphred.out
    
    my $id_contig = '0';
    my @contigs = MG::IO::Polyscan::getContig($polyscan,$id_contig);
    
    for my $contig (@contigs) {
	my $snps = MG::IO::Polyscan::Contig::getSNP_READSITES($contig);
	foreach my $pos (keys %{$snps}) {
	    $polyscan_positions->{$pos} = "polyscan_snp";
	    #print qq(polyscan_snp $pos\n);
	    
	}
	my $indels = MG::IO::Polyscan::Contig::getINDEL_READSITES($contig);
	foreach my $pos (keys %{$indels}) {
	    $polyscan_positions->{$pos} = "polyscan_indel";
	    #print qq(polyscan_indel $pos\n);
	}
    }

    return unless $polyscan_positions;
    return $polyscan_positions;

}



#~/svn/perl_modules/Genome/Model/Tools/Annotate 123> gmt annotate genotype-submission-report -root TEST2.BRCA1.BRCA2.TP53 -ace /gscmnt/238/medseq/human_misc/TCGA.OV.BRCA1.BRCA2.TP53.RT51195.091207/TCGA_Ovarian_Cancer_94_capture_pairs-0000675-Ensembl-56_37a/edit_dir/TCGA_Ovarian_Cancer_94_capture_pairs-0000675-Ensembl-56_37a.ace.1b.genomic.hs -gts /gscmnt/238/medseq/human_misc/TCGA.OV.BRCA1.BRCA2.TP53.RT51195.091207/TCGA_Ovarian_Cancer_94_capture_pairs-0000675-Ensembl-56_37a/edit_dir/TCGA_Ovarian_Cancer_94_capture_pairs-0000675-Ensembl-56_37a.genotype.submission.txt -indels /gscmnt/238/medseq/human_misc/TCGA.OV.BRCA1.BRCA2.TP53.RT51195.091207/TCGA_Ovarian_Cancer_94_capture_pairs-0000675-Ensembl-56_37a/edit_dir/rmeyer_100405.combo.indel.txt -polyphred /gscmnt/238/medseq/human_misc/TCGA.OV.BRCA1.BRCA2.TP53.RT51195.091207/TCGA_Ovarian_Cancer_94_capture_pairs-0000675-Ensembl-56_37a/edit_dir/rmeyer_100405.filtered.polyphred.out -polyscan /gscmnt/238/medseq/human_misc/TCGA.OV.BRCA1.BRCA2.TP53.RT51195.091207/TCGA_Ovarian_Cancer_94_capture_pairs-0000675-Ensembl-56_37a/edit_dir/rmeyer_100405.msa.polyscan -refseq /gscmnt/238/medseq/human_misc/TCGA.OV.BRCA1.BRCA2.TP53.RT51195.091207/TCGA_Ovarian_Cancer_94_capture_pairs-0000675-Ensembl-56_37a/edit_dir/TCGA_Ovarian_Cancer_94_capture_pairs-0000675-Ensembl-56_37a.c1.refseq.fasta 



1;
