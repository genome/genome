package Genome::Model::Tools::Annotate::GtsReport;

use strict;
use warnings;
use Genome;
use IPC::Run;

class Genome::Model::Tools::Annotate::GtsReport {
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
	output => {
	    type  =>  'String',
	    doc   =>  "provide the name for your output file. Default prints to screen",
	    is_optional  => 1,
	},

	], 
};

sub help_brief {
    return <<EOS
	gmt annotate gts-report -gts
EOS
}

sub help_synopsis {
    return <<EOS
	gmt annotate gts-report -gts
EOS
}

sub help_detail { 
    return <<EOS 
	This tool was design to read in a genotype submmision file, identify all the variants, annotate them, check for dbsnps and screen the nonsynonymous snps with sift and polyphen. However, it will work on any list of variants with this info chromosome,coordinate,sample,allele1,allele2. allele1 and allele2 may be substituted with ref_allele and variant allele but they should be in the positive orientation with respect to the reference. If the order of these variants is different from the default order based on the genotype submission file format then use the -order option also see the delimiter option if your file is not tab delimited. ref_id is not a manditory requirement. 
EOS
}

sub execute {


    my $self = shift;

    my $output = $self->output;
    if ($output) {
	open(OUTFILE,">$output") || die "\couldn't open the output file\n";
    }
    my $version = $self->version;
    my $organism = $self->organism;

    if ($output) {
	print OUTFILE qq(organism\tversion\tchromosome\tpos\tstop\tref\tvariant_allele\tsample\tvariant_type\tgenotype\tsource\tgene\ttranscript\tstrand\tTraans_stat\ttrv_type\tc_pos\taa\tpolyphen\tsift\tcons_score\tdomain\trs_id\tdbsnp_submittor\tdbsnp_alleles\tdbsnp_allele_match\n);
    } else {
	print qq(organism\tversion\tchromosome\tpos\tstop\tref\tvariant_allele\tsample\tvariant_type\tgenotype\tsource\tgene\ttranscript\tstrand\tTraans_stat\ttrv_type\tc_pos\taa\tpolyphen\tsift\tcons_score\tdomain\trs_id\tdbsnp_submittor\tdbsnp_alleles\tdbsnp_allele_match\n);
    }
    
    if ($organism eq "mouse" && $version eq "54_36p_v2") { $version = "54_37g_v2"; }
    
    my ($genotypes,$gt_counts) = &parse_gts($self);
    my ($annotation,$dbsnp,$SIFT_Prediction,$PPH_Prediction) = &annotate_variant_alleles($genotypes,$self);
    
    print qq(preparing your output\n);

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
			if ($variant_allele =~ /\+/ || $variant_allele =~ /\-/) {
			    my ($end,$ref,$var) = &indel_alleles($variant_allele,$pos,$chr,$organism);
			    $ref_base = $ref;
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
				if ($output) {
				    print OUTFILE qq($organism\t$version\t$chromosome\t$pos\t$stop\t$ref\t$variant_allele\t$sample\t$variant_type\t$genotype\t$source\t$gene\t$transcript\t$strand\t$Traans_stat\t$trv_type\t$c_pos\t$aa\t$polyphen\t$sift\t$cons_score\t$domain\t$rs_id\t$dbsnp_submittor\t$dbsnp_alleles\t$dbsnp_allele_match\n);
				} else {
				    print qq($organism\t$version\t$chromosome\t$pos\t$stop\t$ref\t$variant_allele\t$sample\t$variant_type\t$genotype\t$source\t$gene\t$transcript\t$strand\t$Traans_stat\t$trv_type\t$c_pos\t$aa\t$polyphen\t$sift\t$cons_score\t$domain\t$rs_id\t$dbsnp_submittor\t$dbsnp_alleles\t$dbsnp_allele_match\n);
				}
			    }
			}
		    }
		}
	    }
	}
    }
    if ($output) {
	print qq(your results have been printed in $output\n);
	close OUTFILE;
    }
    
}
sub write_genomic_prettybase {

    my ($genotypes) = @_; 

    print qq(writting the genomic_prettybase\n);
    open(GPB,">genomic_prettybase.txt");
    foreach my $chr (sort keys %{$genotypes}) {
	foreach my $pos (sort {$a<=>$b} keys %{$genotypes->{$chr}}) {
	    
	    #foreach my $variant_allele (sort keys %{$list->{$chr}->{$pos}}) {
	    
	    my $ref_base = $genotypes->{$chr}->{$pos}->{ref_base};
	    my $ref_sample = $genotypes->{$chr}->{$pos}->{ref_sample};

	    print GPB qq($pos\t$ref_sample\t$ref_base\t$ref_base\n);
	    foreach my $sample (sort keys %{$genotypes->{$chr}->{$pos}->{GT}}) {
		my $gt = $genotypes->{$chr}->{$pos}->{GT}->{$sample};
		my ($a1,$a2) = split(/\:/,$gt);
		print GPB qq($pos\t$sample\t$a1\t$a2\n);
		print qq(\nrunning report_msa_prettybase\n);
		my @command = ["report_msa_prettybase" , "genomic_prettybase.txt"];
		&ipc_run(@command);
		my $snp_counts_file = "snp_counts.txt";
		my $snp_summary_txt = "snp_summary.txt";
		my $snp_counts;


	    }
	}
    }
    close GPB;
}


sub annotate_variant_alleles {
    
    my ($genotypes,$self) = @_;
    
    my $version = $self->version;
    my $organism = $self->organism;
    
    if ($organism eq "mouse" && $version eq "54_36p_v2") { $version = "54_37g_v2"; }
    
    my $list = "GTS.annotation.list";
    open(LIST,">$list") || die "couldn't open annotation list";
    foreach my $chr (sort keys %{$genotypes}) {
	foreach my $pos (sort {$a<=>$b} keys %{$genotypes->{$chr}}) {

	    my $ref_base = $genotypes->{$chr}->{$pos}->{ref_base};
	    foreach my $variant_allele (sort keys %{$genotypes->{$chr}->{$pos}->{variant_allele}}) {
		
		if ($variant_allele =~ /\+/ || $variant_allele =~ /\-/) {
		    my ($stop,$ref,$var) = &indel_alleles($variant_allele,$pos,$chr,$organism);
		    
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
    print qq(running annotate transcript-variants\n);

    my @command = ["gmt" , "annotate" , "transcript-variants" , "--variant-file" , "$list" , "--output-file" , "$annotated_list" , "--reference-transcripts" , "NCBI-$organism.combined-annotation/$version"];

    &ipc_run(@command);

    my $dbsnp_out = "GTS.annotated.list.dbsnp";
    print qq(running annotate lookup-variants\n);
    @command = ["gmt" , "annotate" , "lookup-variants" , "-variant-file" , "$list" , "--output-file" , "$dbsnp_out", "--report-mode" , "full"];

    #@command = ["gmt" , "snp" , "get-dbsnps" , "--list" , "$list" , "--out" , "$dbsnp_out"];
    &ipc_run(@command);

    my ($annotation,$dbsnp,$sift_vars) = &parse_annotation($annotated_list,$dbsnp_out);
    my ($SIFT_Prediction,$PPH_Prediction);
    if ($sift_vars) {
	($SIFT_Prediction,$PPH_Prediction) = &run_sift_and_poylphen($sift_vars,$self);
    }
    return ($annotation,$dbsnp,$SIFT_Prediction,$PPH_Prediction);
}

sub run_sift_and_poylphen {

    my ($sift_vars,$self) = @_;

    my $version = $self->version;
    my $organism = $self->organism;
    if ($organism eq "mouse" && $version eq "54_36p_v2") { $version = "54_37g_v2"; }

    my $SIFT_Prediction;
    my $PPH_Prediction;
    
    foreach my $transcript (sort keys %{$sift_vars}) {
	open(OUT,">$transcript.nonsynonymous.snps.list");
	
	my $chr;
	foreach my $n (sort {$a<=>$b} keys %{$sift_vars->{$transcript}}) {
	    foreach my $var (sort keys %{$sift_vars->{$transcript}->{$n}}) {
		print OUT qq($var\n);
	    }
	}
	close OUT;
	
	my @command = ["rm" , "$transcript.protein.time.error"];
	&ipc_run(@command);

	@command = ["gmt" , "snp" , "screen-nonsynonymous-snps" , "--list" , "$transcript.nonsynonymous.snps.list" , "--transcript" , "$transcript" , "--organism" , "$organism" , "--version" , "$version"];
	
	print qq(will now screen nonsynonymous snps in $transcript\n);

	&ipc_run(@command);
	
	my $sift_resutls_file = "$transcript.protein.SIFTprediction";
	
	if (-f $sift_resutls_file) {
	    open(SIFT,$sift_resutls_file);
	    
	    #print qq(SIFT results are in $sift_resutls_file\n);
	    while (<SIFT>) {
		chomp;
		my $line = $_;
		my ($prot,$td,$n1,$n2,$n3,$n4) = split(/\s+/,$line);
		
		$SIFT_Prediction->{$transcript}->{$prot} = "$td $n1 $n2 $n3 $n4";
	    }
	    close SIFT;
	} else {
	    print qq(no SIFT results were found for $transcript\n);
	}
	
	
	my $polyphen_results_file = "$transcript.protein.PPHprediction";
	if (-f $polyphen_results_file) {
	    
	    #print qq(polypheb results are in $polyphen_results_file\n);
	    
	    open(PPH,$polyphen_results_file);
	    while (<PPH>) {
		chomp;
		my $line = $_;
		my ($snp_id,$status,$protein_id,$position,$aa1,$aa2,$prediction,$basis,$effect,$site,$region,$phat,$score_delta,$num_observ,$num_struct_init,$num_struct_filt,$pdb_id,$res_num,$chain_id,$ali_ide,$ali_len,$acc_normed,$sec_str,$map_region,$delta_volume,$delta_prop,$b_fact,$num_h_bonds,$het_cont_ave_num,$het_cont_min_dist,$inter_cont_ave_num,$inter_cont_min_dist,$sites_cont_ave_num,$sites_cont_min_dist) = split(/\t/,$line);
		my $prot = "$aa1$position$aa2";
		$PPH_Prediction->{$transcript}->{$prot} = "$prediction $score_delta";
	    }
	    close PPH;
	} else {
	    print qq(no polyphen results were found for $transcript\n);
	}
	
	#SIFT_Prediction "R132H DEoLETERIOUS 0.00 3.07 210 281"
	#PPH_Prediction "R132H probably damaging ? "
	
	@command = ["rm" , "$transcript.nonsynonymous.snps.list" , "$transcript.fasta" , "SIFT_$transcript.protein.fasta.out" , "$transcript.protein.PPHprediction" , "$transcript.protein.SIFTprediction" , "$transcript.protein.SIFTprediction.error" , "$transcript.protein.alignedfasta" , "$transcript.protein.clumped" , "$transcript.protein.clumped.consensuskey" , "$transcript.protein.clumped.error" , "$transcript.protein.fasta" , "$transcript.protein.fasta.out" , "$transcript.protein.fasta.query" , "$transcript.protein.fasta.query.globalX" , "$transcript.protein.fasta.query.globalX.error" , "$transcript.protein.fasta.query.out" , "$transcript.protein.fasta.query.selectedclumped" , "$transcript.protein.fasta.query.selectedclumped.error" , "$transcript.protein.fasta.query.selectedclumped.log" , "$transcript.protein.fasta.query.unfiltered" , "$transcript.protein.time.error"];
	
	&ipc_run(@command);
	
    }
    
    return($SIFT_Prediction,$PPH_Prediction);
}


sub parse_annotation {
    
    my ($annotation_file,$dbsnp_out) = @_;
    my $dbsnp;
    
    if ($dbsnp_out && -e $dbsnp_out) {
	#print qq(reading the dbsnp hit results\n);
	open(DBSNP,"$dbsnp_out") || die "couldn't open the dbsnp results";
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
	print qq(no dbsnps query result\n);
    }
    my $sift_vars;
    my $annotation;
    open(ANO,"$annotation_file") || die "couldn't open the annotated list\n";
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

sub ipc_run {
    
    my (@command) = @_;

    my ($in, $out, $err);
    my ($obj) = IPC::Run::run(@command, \$in, \$out, \$err);
    if ($err) {
	#print qq($err\n);
    }
    if ($out) {
	return ($out);
	#print qq($out\n);
    }
}

sub indel_alleles {

    my ($allele,$pos,$chr,$organism) = @_;

    my ($ref,$var,$end);
    if ($allele =~ /\+([\S]+)/) { ##INS
	$var = $1;
	$ref = "-";
	$end = $pos + 1;
    } elsif ($allele =~ /\-([\S]+)/) { ##DEL
	my $bases = $1;
	my $length = length($bases); 
	$end = ($pos - 1) + $length;
	$ref = &get_ref_base($pos,$end,$chr,$organism);
	$var = "-";
    }
    return ($end,$ref,$var);
}

sub parse_gts {

    my ($self) = @_;

    my $gts = $self->gts;
    my $delimiter = $self->delimiter;
    my $order = $self->order;
    my ($chr_n,$pos_n,$sample_n,$allele1_n,$allele2_n) = split(/\,/,$order);

    my $organism = $self->organism;
    my $ref_id = $self->ref_id;

    my ($genotypes,$gt_counts);
    open(GTS,$gts) || die "could't open the genotype submission file\n"; 
    print qq(\nExtracting variants from the genotype submission file\n);

    while (<GTS>) {
	chomp;
	my $line = $_;
	my ($chr,$pos,$sample,$allele1,$allele2) = (split(/$delimiter/,$line))[$chr_n - 1,$pos_n - 1,$sample_n - 1,$allele1_n - 1,$allele2_n - 1];

	$chr =~ s/C//;

	my $ref_base = &get_ref_base($pos,$pos,$chr,$organism);
	$genotypes->{$chr}->{$pos}->{ref_base}=$ref_base;
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

sub reverse_complement_allele {
    my ($allele_in) = @_;

    use Bio::Seq;

    unless ($allele_in =~ /[ACGT]/) { return ($allele_in); }
    my $seq_1 = new Bio::Seq(-seq => $allele_in);
    my $revseq_1 = $seq_1->revcom();
    my $rev1 = $revseq_1->seq;
    
    my $r_base = $rev1;
    
    return $r_base;
}

sub get_ref_base {

    my ($chr_start,$chr_stop,$chr_name,$organism) = @_;

    use Bio::DB::Fasta;

    my $RefDir;
    if ($organism eq "human"){
	$RefDir = "/gscmnt/sata180/info/medseq/biodb/shared/Hs_build36_mask1c/";
    } else {
	$RefDir = "/gscmnt/sata147/info/medseq/rmeyer/resources/MouseB37/";
    }

    my $refdb = Bio::DB::Fasta->new($RefDir);

    my $seq = $refdb->seq($chr_name, $chr_start => $chr_stop);
    $seq =~ s/([\S]+)/\U$1/;

    return $seq;
}

1;
