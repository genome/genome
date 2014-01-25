#!/usr/bin/env genome-perl
use strict;
use warnings;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use above 'Genome';

#Inputs
my $clinseq_model_id = '';
my $outdir = '';
my $bedtools_bin_dir = '';
my $deletion_regions_file = '';
my $control_regions_file = '';
my $filter = '';

GetOptions ('clinseq_model_id=i'=>\$clinseq_model_id, 'outdir=s'=>\$outdir, 'bedtools_bin_dir=s'=>\$bedtools_bin_dir, 
            'deletion_regions_file=s'=>\$deletion_regions_file, 'control_regions_file=s'=>\$control_regions_file,
            'filter=i'=>\$filter);

my $usage=<<INFO;
Example usage: 

  No filtering
  selectSnpSnvCaptureTargets.pl  --clinseq_model_id='2882726707'  --outdir=/gscmnt/sata132/techd/mgriffit/hgs/custom/all01/snp_array_design/result/pre_filtered/  --bedtools_bin_dir=/gscmnt/gc2142/techd/tools/bedtools/BEDTools-Version-2.14.3/bin/  --deletion_regions_file=/gscmnt/sata132/techd/mgriffit/hgs/custom/all01/snp_array_design/cnaseq.cnvhmm.losses.merged.final.txt  --control_regions_file=/gscmnt/sata132/techd/mgriffit/hgs/custom/all01/snp_array_design/control_regions.txt  --filter=0

  With filtering
  selectSnpSnvCaptureTargets.pl  --clinseq_model_id='2882726707'  --outdir=/gscmnt/sata132/techd/mgriffit/hgs/custom/all01/snp_array_design/result/post_filtered/  --bedtools_bin_dir=/gscmnt/gc2142/techd/tools/bedtools/BEDTools-Version-2.14.3/bin/  --deletion_regions_file=/gscmnt/sata132/techd/mgriffit/hgs/custom/all01/snp_array_design/cnaseq.cnvhmm.losses.merged.final.txt  --control_regions_file=/gscmnt/sata132/techd/mgriffit/hgs/custom/all01/snp_array_design/control_regions.txt  --filter=1

  Parameters
  --clinseq_model_id         ClinSeq model ID for the patient of interest.  SNVs will be selected from the underlying WGS and Exome somatic variation results
  --outdir                   Output directory for temp and results files
  --bedtools_bin_dir         Path to BEDtools binaries
  --deletion_regions_file    File containing deletion regions where heterozygous SNPs will be selected
  --control_regions_file     Control regions that have be pre-defined to exclude the deletion regions
  --filter                   Whether to apply coverage, VAF, and CNV diff cutoffs to the SNP positions accordingly (SNVs will not be filtered)

INFO

unless ($clinseq_model_id && $outdir && $bedtools_bin_dir && $deletion_regions_file && $control_regions_file && defined($filter)){
  print RED, "\n\nRequired parameter missing", RESET;
  print GREEN, "\n\n$usage", RESET;
  exit(1);
}

#Filtering cutoffs.  Apply only to SNPs in control and deletion regions.  SNVs will not be filtered
my $min_normal_read_coverage = 20;       #SNPs must have a read covearge level in normal greater than this
my $max_normal_read_coverage = 517;      #SNPs must have a read coverage level in normal less than this
my $het_vaf_diff_var = 10;               #SNPs must have a VAF within this distance from 50%.  (i.e. $vaf > 40 && $vaf < 60)
my $cnv_control_window_diff_var = 0.10;  #SNPs within control regions must have a CNV difference within this distance from 0 (ie. $diff > -0.1 && $diff < 0.1)
my $cnv_deletion_window_diff_var = 0.25; #SNPs within deletion region must have a CNV diff at least this large (i.e. $diff < -0.25)

my $target_probe_count = 5600;  #Target total number of positions
my $probes_per_deletion = 75;   #Target number of positions per deletion region defined


#Get somatic variation files from the ClinSeq model (both WGS and Exome)
my $clinseq_model = Genome::Model->get("id"=>$clinseq_model_id);
my $clinseq_build = $clinseq_model->last_succeeded_build;

my $wgs_somvar_build = $clinseq_build->wgs_build;
my $wgs_normal_refalign_build = $wgs_somvar_build->normal_build;
my $wgs_tumor_refalign_build = $wgs_somvar_build->tumor_build;
my $wgs_normal_bam = $wgs_somvar_build->normal_bam;
my $wgs_tumor_bam = $wgs_somvar_build->tumor_bam;

my $exome_somvar_build = $clinseq_build->exome_build;
my $exome_normal_refalign_build = $exome_somvar_build->normal_build;
my $exome_tumor_refalign_build = $exome_somvar_build->tumor_build;
my $exome_normal_bam = $exome_somvar_build->normal_bam;
my $exome_tumor_bam = $exome_somvar_build->tumor_bam;

#Obtain input files:
my $wgs_data_dir = $wgs_somvar_build->data_directory;
my $exome_data_dir = $exome_somvar_build->data_directory;

my $wgs_tier1_novel_snvs_file = $wgs_data_dir . "/effects/snvs.hq.novel.tier1.v2.bed";
my $wgs_tier2_novel_snvs_file = $wgs_data_dir . "/effects/snvs.hq.novel.tier2.v2.bed";
my $wgs_tier3_novel_snvs_file = $wgs_data_dir . "/effects/snvs.hq.novel.tier3.v2.bed";
my $wgs_tier1_known_snvs_file = $wgs_data_dir . "/effects/snvs.hq.previously_detected.tier1.v2.bed";
my $wgs_tier2_known_snvs_file = $wgs_data_dir . "/effects/snvs.hq.previously_detected.tier2.v2.bed";
my $wgs_tier3_known_snvs_file = $wgs_data_dir . "/effects/snvs.hq.previously_detected.tier3.v2.bed";
my $exome_tier1_novel_snvs_file = $exome_data_dir . "/effects/snvs.hq.novel.tier1.v2.bed";
my $exome_tier1_known_snvs_file = $exome_data_dir . "/effects/snvs.hq.previously_detected.tier1.v2.bed";
my $flt3_snvs_file = "/gscmnt/sata132/techd/mgriffit/hgs/custom/all01/snp_array_design/flt3_snps/all_flt3_het_snps.bed";

my $wgs_cnvs_file = $wgs_data_dir . "/variants/cnvs.hq";
my $wgs_normal_refalign_dir = $wgs_normal_refalign_build->data_directory;
my $wgs_tumor_refalign_dir = $wgs_tumor_refalign_build->data_directory;


#1.) Import SNVs from various files using a subroutine.  Define them as various classes.
print BLUE, "\n\nImporting SNV/InDel positions from various somatic variation results files:", RESET;
my %pos;
my $pos_ref = \%pos;
&importPositions('-pos_ref'=>$pos_ref, '-pos_file'=>$wgs_tier1_novel_snvs_file, '-class'=>"Tier1_SNV_Novel", '-region_name_prefix'=>"tier1_snv_novel_");
&importPositions('-pos_ref'=>$pos_ref, '-pos_file'=>$exome_tier1_novel_snvs_file, '-class'=>"Tier1_SNV_Novel", '-region_name_prefix'=>"tier1_snv_novel_");
&importPositions('-pos_ref'=>$pos_ref, '-pos_file'=>$wgs_tier1_known_snvs_file, '-class'=>"Tier1_SNV_Known", '-region_name_prefix'=>"tier1_snv_known_");
&importPositions('-pos_ref'=>$pos_ref, '-pos_file'=>$exome_tier1_known_snvs_file, '-class'=>"Tier1_SNV_Known", '-region_name_prefix'=>"tier1_snv_known_");
&importPositions('-pos_ref'=>$pos_ref, '-pos_file'=>$wgs_tier2_novel_snvs_file, '-class'=>"Tier2_SNV_Novel", '-region_name_prefix'=>"tier2_snv_novel_");
&importPositions('-pos_ref'=>$pos_ref, '-pos_file'=>$wgs_tier2_known_snvs_file, '-class'=>"Tier2_SNV_Known", '-region_name_prefix'=>"tier2_snv_known_");
&importPositions('-pos_ref'=>$pos_ref, '-pos_file'=>$wgs_tier3_novel_snvs_file, '-class'=>"Tier3_SNV_Novel", '-region_name_prefix'=>"tier3_snv_novel_");
&importPositions('-pos_ref'=>$pos_ref, '-pos_file'=>$wgs_tier3_known_snvs_file, '-class'=>"Tier3_SNV_Known", '-region_name_prefix'=>"tier3_snv_known_");
&importPositions('-pos_ref'=>$pos_ref, '-pos_file'=>$flt3_snvs_file, '-class'=>"FLT3_SNP", '-region_name_prefix'=>"flt3_snp");

my $snv_count = keys %pos;

#Note on annotations.  Some annotations may not be present in the refalign results (because the calls might actually be from a caller other than samtools)
#This means that we will have to look for annotations in more than one place. :(
#Go back to the original annotated SNP file and for the remaining SNPs, gather their ref_base and var_base values

#Now go back to the reference alignment results for the TUMOR and store extra annotations for each variant position loaded so far
#Get this extra annotation information from the post annotation file

#In particular we need to add: ref_base, var_base, gene_id, trans_id, var_type, aa_change
my $wgs_tumor_tier1_annot = "$wgs_data_dir/effects/snvs.hq.tier1.v1.annotated.top";
my $wgs_tumor_tier2_annot = "$wgs_data_dir/effects/snvs.hq.tier2.v1.annotated.top";
my $exome_tumor_tier1_annot = "$exome_data_dir/effects/snvs.hq.tier1.v1.annotated.top";
my $exome_tumor_tier2_annot = "$exome_data_dir/effects/snvs.hq.tier2.v1.annotated.top";
my $wgs_tumor_post_annot_file = "$wgs_tumor_refalign_dir/variants/filtered.variants.post_annotation";
my $wgs_normal_post_annot_file = "$wgs_normal_refalign_dir/variants/filtered.variants.post_annotation";
my @annot_files = ($wgs_tumor_tier1_annot, $wgs_tumor_tier2_annot, $exome_tumor_tier1_annot, $exome_tumor_tier2_annot, $wgs_tumor_post_annot_file, $wgs_normal_post_annot_file);

print BLUE, "\n\nAdding annotation values for each SNV from:", RESET;
foreach my $annot_file (@annot_files){
  print BLUE, "\n\t$annot_file", RESET;
  open (ANNO, "$annot_file") || die "\n\nCould not open annot file: $annot_file\n\n";
  while(<ANNO>){
    chomp($_);
    my @line = split("\t", $_);
    my $type = $line[5];
    my $chr = $line[0];
    my $start = $line[1];
    my $end = $line[2];
    my $coord = $chr . ":" . $start . "-" . $end;
    
    unless ($type eq "SNP"){
      next();
    }

    if ($pos{$coord}){
      #Only update annotation if it was not previously found
      if ($pos{$coord}{gene_id} eq "NA"){
        $pos{$coord}{ref_base} = $line[3];
        $pos{$coord}{var_base} = $line[4];
        $pos{$coord}{gene_id} = $line[6];
        $pos{$coord}{trans_id} = $line[7];
        $pos{$coord}{var_type} = $line[13];
        $pos{$coord}{aa_change} = $line[15];
      }
    }
  }
  close(ANNO);
}

#Store the distinct list of SNVs/InDels as a file for reference
my $distinct_snv_file = $outdir . "snvs.tsv";
open (SNV_OUT, ">$distinct_snv_file") || die "\n\nCould not open distinct tier1 SNV out file: $distinct_snv_file\n\n";
print SNV_OUT "chr\tstart\tend\tref_base\tvar_base\tgene_id\ttrans_id\tvar_type\taa_change\n";
foreach my $coord (sort keys %pos){
  print SNV_OUT "$pos{$coord}{chr}\t$pos{$coord}{start}\t$pos{$coord}{end}\t$pos{$coord}{ref_base}\t$pos{$coord}{var_base}\t$pos{$coord}{gene_id}\t$pos{$coord}{trans_id}\t$pos{$coord}{var_type}\t$pos{$coord}{aa_change}\n";
}
close(SNV_OUT);


#2.) Define the somatic deletion regions in the genome using hmm copy number segmentation
#    These regions were defined by starting with the hmmCopy number segments from the Clonality analysis of ALL1
#    Obi manually reviewed and merged these segments into a final list of target deletions
print BLUE, "\n\nImporting target regions (deletion regions and control regions)", RESET;
my %region_name_list;
my %dels;
open (DELS, "$deletion_regions_file") || die "\n\nCould not open deletion regions file\n\n";
my $header = 1;
while (<DELS>){
  if ($header){
    $header = 0;
    next();
  }
  chomp($_);
  my @line = split("\t", $_);
  my $chr = $line[0];
  my $start = $line[1];
  my $end = $line[2];

  #Make the target regions smaller to avoid SNPs being selected at their edges
  my $subtract_flank = 10000;
  my $size = $end - $start;
  if ($size > (($subtract_flank*2)+50000)){
    $start += $subtract_flank;
    $end -= $subtract_flank;
  }

  my $coord = $chr . ":" . $start . "-" . $end;
  $dels{$coord}{chr} = $chr;
  $dels{$coord}{start} = $start;
  $dels{$coord}{end} = $end;
  $dels{$coord}{size} = $line[3];
  $dels{$coord}{n_markers} = $line[4];
  $dels{$coord}{cn1} = $line[5];
  $dels{$coord}{adjusted_cn1} = $line[6];
  $dels{$coord}{cn2} = $line[7];
  $dels{$coord}{adjusted_cn2} = $line[8];
  $dels{$coord}{llr_somatic} = $line[9];
  $dels{$coord}{status} = $line[10];
  $dels{$coord}{merge_group} = $line[11];
  my $region_name = "deletion_region_"."$line[11]";
  $dels{$coord}{region_name} = $region_name;
  $region_name_list{$region_name} = $coord;
}
close (DELS);
my $target_deletion_count = keys %dels;
print BLUE, "\n\tImported $target_deletion_count deletion regions", RESET;


#3.) Define some 'control regions' that are not deleted
#    These regions were defined by considering several other types of regions
#    To identify the control regions, a series of problematic regions were excluded (Telomeres, centromeres, etc.).  These regions came from here:
#    /gscmnt/sata132/techd/mgriffit/reference_annotations/hg19/ideogram/hg19gaps.csv
#    All deletion regions from the previous step as well as all amplification regions were also excluded
#    The remaining segments of the genome will be imported as the control regions
my %controls;
open (CON, "$control_regions_file") || die "\n\nCould not open control regions file\n\n";
$header = 1;
my $control_region_count = 0;
while (<CON>){
  if ($header){
    $header = 0;
    next();
  }
  chomp($_);
  my @line = split("\t", $_);
  my $chr = $line[0];
  my $start = $line[1];
  my $end = $line[2];

  #Make the target regions smaller to avoid SNPs being selected at their edges.  Skip small control regions entirely
  my $subtract_flank = 10000;
  my $size = $end - $start;
  if ($size < (($subtract_flank*2)+50000)){
    next();
  }
  $control_region_count++;
  $start += $subtract_flank;
  $end -= $subtract_flank;

  my $coord = $chr . ":" . $start . "-" . $end;
  $controls{$coord}{chr} = $chr;
  $controls{$coord}{start} = $start;
  $controls{$coord}{end} = $end;
  my $region_name = "control_region_"."$control_region_count";
  $controls{$coord}{region_name} = $region_name;
  $region_name_list{$region_name} = $coord;
}
close (CON);
my $target_control_count = keys %controls;
print BLUE, "\n\tImported $target_control_count control regions", RESET;

my $region_name_count = keys %region_name_list;
print BLUE, "\n\tImported $region_name_count total distinct target regions (deletion regions + control regions)", RESET;

#Create a bed file for the deletion and control regions to allow a pre-filtering of SNP positions
my $target_regions_bed_file = $outdir . "target_regions.bed";
print BLUE, "\n\nCreating a bed file of all deletion and control regions:\n\t$target_regions_bed_file", RESET;
open (TBED, ">$target_regions_bed_file") || die "\n\nCould not open target regions bed for writing\n\n";
foreach my $coord (sort keys %dels){
  print TBED "$dels{$coord}{chr}\t$dels{$coord}{start}\t$dels{$coord}{end}\t+\n"; 
}
foreach my $coord (sort keys %controls){
  print TBED "$controls{$coord}{chr}\t$controls{$coord}{start}\t$controls{$coord}{end}\t+\n"; 
}
close(TBED);


#4.) Gather heterozygous SNPs detected in the normal genome that correspond to the regions above
#Identify SNPs in the 'filtered.variants.pre_annotation' by the reported change.  A SNP will not have a '-' as the ref or variant base
#Use the ambiguity code to limit the analysis only predicted heterozygous SNPs: Y, R, W, K, M, S
my $wgs_normal_snp_pre_annot_file = "$wgs_normal_refalign_dir/variants/filtered.variants.pre_annotation";
my $wgs_normal_snp_post_annot_file = "$wgs_normal_refalign_dir/variants/filtered.variants.post_annotation";
print BLUE, "\n\nSearching for heterozygous SNPs called in the normal tissue by the refalign pipeline here: $wgs_normal_snp_pre_annot_file", RESET;
my %normal_snps;

#Create a BED file of all possible SNPs called in the normal genome, using the 'filtered.variants.post_annotation'
print BLUE, "\n\nCreating a BED file of all SNPs called in the normal genome", RESET;
my $wgs_normal_snps_bed_file = $outdir . "all_normal_snps.bed";
open (SNP, "$wgs_normal_snp_pre_annot_file") || die "\n\nCould not open normal SNPs file: $wgs_normal_snp_pre_annot_file\n\n";
open (SNP_BED, ">$wgs_normal_snps_bed_file") || die "\n\nCould not open normal SNPs BED file for writing: $wgs_normal_snps_bed_file\n\n";
my $initial_het_snp_count = 0;
while(<SNP>){
  chomp($_);
  my @line = split("\t", $_);
  my $chr = $line[0];
  my $start = $line[1];
  my $end = $line[2];
  my $ref_base = $line[3];
  my $var_base = $line[4];

  #Skip Insertions and Deletions, leaving SNPs only
  if ($ref_base eq "-" || $var_base eq "-"){
    next();
  }
  #Skip the 'GL' chromosomes that will be difficult to tie into other data in the following steps:
  #Also skip variants on ChrY
  if ($chr =~ /^GL|^Y/){
    next();
  }
  #Skip all but the heterozygous ambiguity codes for the variant base
  unless ($var_base =~ /Y|R|W|K|M|S/){
    next();
  }
  my $coord = $chr . ":" . $start . "-" . $end;

  #Do not allow a SNP position to override a previous stored SNV
  if ($pos{$coord}){
    next();
  }
  
  $initial_het_snp_count++;
  print SNP_BED "$chr\t$start\t$end\t+\n";
}
close(SNP);
close(SNP_BED);
print BLUE, "\n\tIdentified $initial_het_snp_count heterozyous SNPs in the normal", RESET;

#Identify overlaps between the target regions and the SNP positions using BEDTools
print BLUE, "\n\nPre-filtering to only those normal SNPs within a target region (deletion or control region)", RESET;
my $bedtools_intersect_file = $outdir . "all_normal_snps.intersect.target_regions.bed";
my $bedtools_cmd = "$bedtools_bin_dir"."intersectBed -a $wgs_normal_snps_bed_file -b $target_regions_bed_file -f 1.0 -wa -wb > $bedtools_intersect_file";
print YELLOW, "\n\n\t$bedtools_cmd", RESET;
Genome::Sys->shellcmd(cmd => $bedtools_cmd, output_files=>["$bedtools_intersect_file"]);

#Parse the intersection file, import each normal SNP and organize by the target region it corresponds to
print BLUE, "\n\nParsing intersection and determining the target region each SNP lies within", RESET;
open (SNPBED, "$bedtools_intersect_file") || die "\n\nCould not open bedtools intersect file: $bedtools_intersect_file\n\n";
my $intersect_count = 0;
while(<SNPBED>){
  chomp($_);
  my @line = split("\t", $_);
  $intersect_count++;

  #SNP coords
  my $chr1 = $line[0];
  my $start1 = $line[1];
  my $end1 = $line[2];
  my $coord1 = $chr1 . ":" . "$start1" . "-" . "$end1";

  #Target region coords
  my $chr2 = $line[4];
  my $start2 = $line[5];
  my $end2 = $line[6];
  my $coord2 = $chr2 . ":" . "$start2" . "-" . "$end2";

  my $region_name;
  if ($dels{$coord2}){
    $region_name = $dels{$coord2}{region_name};
  }elsif($controls{$coord2}){
    $region_name = $controls{$coord2}{region_name};
  }else{
    print RED, "\n\nIntersection involving unknown target region coordinate ($coord2)\n\n", RESET;
    exit(1);
  }

  #Store candidate SNPs and organize by region name
  if ($normal_snps{$region_name}){
    my $snps_ref = $normal_snps{$region_name}{snps};
    $snps_ref->{$coord1}->{chr} = $chr1;
    $snps_ref->{$coord1}->{start} = $start1;
    $snps_ref->{$coord1}->{end} = $end1;
  }else{
    my %snps;
    $snps{$coord1}{chr} = $chr1;
    $snps{$coord1}{start} = $start1;
    $snps{$coord1}{end} = $end1;
    $normal_snps{$region_name}{snps} = \%snps;
  }
}
close(SNPBED);
print BLUE, "\n\t$intersect_count SNPs remain", RESET;


#Then pre-select up to N SNP positions for each target region (control or deletion region).  Select them evenly across each region.
#Thin the number of possible SNPs per region to a more tractable number by evenly selecting across each target region
my $snps_per_region = 250;

print BLUE, "\n\nThinning number of SNPs from very large regions to a more manageable number:", RESET;
my %selected_region_snps;
foreach my $region_name (sort keys %region_name_list){
  my $region_desc = $region_name_list{$region_name};
  my $region_snps = $normal_snps{$region_name}{snps};
  my $region_snps_count = keys %{$region_snps};
  print BLUE, "\n\tProcessing region: $region_name ($region_desc containing $region_snps_count possible SNPs)", RESET;

  #If the region has less than $snps_per_region SNPs, do nothing 
  if ($region_snps_count > $snps_per_region){
    my $x = $region_snps_count/$snps_per_region;
    #my $skip_value = sprintf("%.0f", $x);
    my $skip_value = int($x);
    print BLUE, "\n\t\tThinning to every $skip_value th other SNP", RESET;

    #Otherwise create a sorted array of the SNPs and thin them
    my @snps;
    foreach my $coord (sort {$region_snps->{$a}->{start} <=> $region_snps->{$b}->{start}} keys %{$region_snps}){
      push(@snps, $coord);
    }

    #Get the thinned SNPs list
    my %thinned_snps;
    for (my $i = 1; $i <= $region_snps_count; $i += $skip_value){
      my $coord = $snps[$i-1];
      $thinned_snps{$coord}=1;
    }

    #Now go through the hash of possible SNPs again and delete those that are not part of the thinned list
    foreach my $coord (sort {$region_snps->{$a}->{start} <=> $region_snps->{$b}->{start}} keys %{$region_snps}){
      unless($thinned_snps{$coord}){
        delete($region_snps->{$coord});
        next();
      }
    }
    my $new_region_snps_count = keys %{$region_snps};
    print BLUE, "\n\t\tRemaining SNPs in this region:  $new_region_snps_count", RESET;
  }else{
    print BLUE, "\n\t\tThinning not needed - skipping", RESET;
  }

  foreach my $coord (sort {$region_snps->{$a}->{start} <=> $region_snps->{$b}->{start}} keys %{$region_snps}){
    $selected_region_snps{$coord}{chr} = $region_snps->{$coord}->{chr};
    $selected_region_snps{$coord}{start} = $region_snps->{$coord}->{start};
    $selected_region_snps{$coord}{end} = $region_snps->{$coord}->{end};
    $selected_region_snps{$coord}{region_name} = $region_name;
  }
}
my $region_snps_count = keys %selected_region_snps;
print BLUE, "\n\nStored a new total of $region_snps_count SNPs for the deletion and control regions", RESET;

#Go back to the original annotated SNP file and for the remaining SNPs, gather their ref_base and var_base values
print BLUE, "\n\nAdding annotation values for the remaining normal SNPs in the candidate list from:\n\t$wgs_normal_snp_post_annot_file", RESET;
open (SNP, "$wgs_normal_snp_post_annot_file") || die "\n\nCould not open normal SNPs file: $wgs_normal_snp_post_annot_file\n\n";
while(<SNP>){
  chomp($_);
  my @line = split("\t", $_);
  my $type = $line[5];
  unless ($type eq "SNP"){
    next();
  }
  my $chr = $line[0];
  my $start = $line[1];
  my $end = $line[2];
  my $coord = $chr . ":" . $start . "-" . $end;
  if ($selected_region_snps{$coord}){
    $selected_region_snps{$coord}{ref_base} = $line[3];
    $selected_region_snps{$coord}{var_base} = $line[4];
    $selected_region_snps{$coord}{gene_id} = $line[6];
    $selected_region_snps{$coord}{trans_id} = $line[7];
    $selected_region_snps{$coord}{var_type} = $line[13];
    $selected_region_snps{$coord}{aa_change} = $line[15];
  }
}
close(SNP);


#5.) Compile a grand list of positions that includes both the heterozygous SNPs from the deletion and control regions as well as the SNV postions
print BLUE, "\n\nAdding SNPs from deletion and control regions to the list from SNV analysis", RESET;

foreach my $coord (keys %selected_region_snps){
  my $chr = $selected_region_snps{$coord}{chr};
  my $start = $selected_region_snps{$coord}{start};
  my $end = $selected_region_snps{$coord}{end};
  my $region_name = $selected_region_snps{$coord}{region_name};
  my $class;
  if ($region_name =~ /deletion/){
    $class = "Deletion_Het_SNP";
  }elsif($region_name =~ /control/){
    $class = "Control_Het_SNP";
  }else{
    print RED, "\n\nUnrecognized region name: $region_name\n\n", RESET;
    exit(1);
  }
  $pos{$coord}{chr} = $chr;
  $pos{$coord}{start} = $start;
  $pos{$coord}{end} = $end;
  $pos{$coord}{ref_base} = $selected_region_snps{$coord}{ref_base};
  $pos{$coord}{var_base} = $selected_region_snps{$coord}{var_base};
  $pos{$coord}{gene_id} = $selected_region_snps{$coord}{gene_id};
  $pos{$coord}{trans_id} = $selected_region_snps{$coord}{trans_id};
  $pos{$coord}{var_type} = $selected_region_snps{$coord}{var_type};
  $pos{$coord}{aa_change} = $selected_region_snps{$coord}{aa_change};
  $pos{$coord}{class} = $class;
  $pos{$coord}{region_name} = $region_name;
}

#Initialize some value that will be tracked for all positions
foreach my $coord (keys %pos){
  $pos{$coord}{normal_ref_rc} = 0;
  $pos{$coord}{normal_var_rc} = 0;
  $pos{$coord}{normal_cov} = 0;
  $pos{$coord}{normal_vaf} = 0;

  $pos{$coord}{tumor_ref_rc} = 0;
  $pos{$coord}{tumor_var_rc} = 0;
  $pos{$coord}{tumor_cov} = 0;
  $pos{$coord}{tumor_vaf} = 0;

  $pos{$coord}{cnv_window_coord} = "unknown";
  $pos{$coord}{cnv_window_read_coverage} = "unknown";
  $pos{$coord}{cnv_window_diff} = "unknown";
}
my $new_pos_count = keys %pos;
print BLUE, "\n\tCurrent list includes $new_pos_count positions", RESET;


#6.) Produce a variant file to be used for BAM read counting
#    File containing snvs in 1-based, 5-col format (chr, st, sp, var, ref). indels will be skipped
my $variant_file = $outdir . "all.variants.tsv";
my $variant_file_sorted = $outdir . "all.variants.sorted.tsv";
open (VAR, ">$variant_file") || die "\n\nCould not open variant file for writing\n\n";
foreach my $coord (sort keys %pos){
  print VAR "$pos{$coord}{chr}\t$pos{$coord}{start}\t$pos{$coord}{end}\t$pos{$coord}{ref_base}\t$pos{$coord}{var_base}\n";  
}
close(VAR);
print BLUE, "\n\nGet BAM read counts for the following variant file from the following BAMs:\n\tVAR FILE: $variant_file\n\tEXOME NORMAL BAM: $exome_normal_bam\n\tEXOME TUMOR BAM: $exome_tumor_bam\n\tWGS NORMAL BAM: $wgs_normal_bam\n\tWGS TUMOR BAM: $wgs_tumor_bam", RESET;

#Sort the output file on chr, start, end to improve efficiency of index lookups when accessing the BAM
my $sort_cmd = "sort -k 1,1 -k 2n,2n -k 3n,3n $variant_file > $variant_file_sorted";
print YELLOW, "\n\t$sort_cmd", RESET;
Genome::Sys->shellcmd(cmd => $sort_cmd, output_files=>["$variant_file_sorted"]);


#7.) Gather basic status about all of these positions from both the tumor and the normal samples.
#    - Get these stats from the WGS data not the Exome data:
#    - Variant and reference read counts from tumor and normal
#    - This will allow us to apply a minimum coverage cutoff to the positions we are going to select
#    - Variant allele frequency.  This will also be used to apply a cutoff (within a specified range expected for heterozygous)
#    - Get the CNV difference values for each position from the cnv windows file in the somatic variation results

#gmt analysis coverage bam-readcount
#Apply a minimum mapping quality score to the reads being counted using the option:  --min-quality-score (e.g. 20)

my $exome_normal_bamrc_outfile = &getBamReadCounts('-outdir'=>$outdir, '-outfile_name'=>'normal.exome.bam.readcounts.tsv', '-bam_file'=>$exome_normal_bam, '-variant_file'=>$variant_file_sorted, '-genome_build'=>'37lite', '-min_quality_score'=>'20');

my $exome_tumor_bamrc_outfile = &getBamReadCounts('-outdir'=>$outdir, '-outfile_name'=>'tumor.exome.bam.readcounts.tsv', '-bam_file'=>$exome_tumor_bam, '-variant_file'=>$variant_file_sorted, '-genome_build'=>'37lite', '-min_quality_score'=>'20');

my $wgs_normal_bamrc_outfile = &getBamReadCounts('-outdir'=>$outdir, '-outfile_name'=>'normal.wgs.bam.readcounts.tsv', '-bam_file'=>$wgs_normal_bam, '-variant_file'=>$variant_file_sorted, '-genome_build'=>'37lite', '-min_quality_score'=>'20');

my $wgs_tumor_bamrc_outfile = &getBamReadCounts('-outdir'=>$outdir, '-outfile_name'=>'tumor.wgs.bam.readcounts.tsv', '-bam_file'=>$wgs_tumor_bam, '-variant_file'=>$variant_file_sorted, '-genome_build'=>'37lite', '-min_quality_score'=>'20');

#Parse the bam read counts files to get reference and variant read counts.  Take the sum of counts from WGS and Exome, then calculate VAF
my %rc_files;
$rc_files{'1'}{file} = $wgs_normal_bamrc_outfile;
$rc_files{'1'}{tissue} = "normal";
$rc_files{'2'}{file} = $exome_normal_bamrc_outfile;
$rc_files{'2'}{tissue} = "normal";
$rc_files{'3'}{file} = $wgs_tumor_bamrc_outfile;
$rc_files{'3'}{tissue} = "tumor";
$rc_files{'4'}{file} = $exome_tumor_bamrc_outfile;
$rc_files{'4'}{tissue} = "tumor";

print BLUE, "\n\nParsing bam read count files:", RESET;
foreach my $c (sort {$a <=> $b} keys %rc_files){
  my $file = $rc_files{$c}{file};
  my $tissue = $rc_files{$c}{tissue};
  print BLUE, "\n\t$file", RESET;
  open (RC, "$file") || die "\n\nCould not open read count file: $file\n\n";
  while(<RC>){
    chomp($_);
    my @line = split("\t", $_);
    my $coord = $line[0] . ":" . $line[1] . "-" . "$line[1]";
    unless ($pos{$coord}){
      print RED, "\n\nFound an unrecognized coordinate ($coord) in the bam read counts results file:\n\t$file\n\n", RESET;
      exit(1);
    }
    if ($tissue eq "normal"){
      $pos{$coord}{normal_ref_rc} += $line[4];
      $pos{$coord}{normal_var_rc} += $line[5];
    }else{
      $pos{$coord}{tumor_ref_rc} += $line[4];
      $pos{$coord}{tumor_var_rc} += $line[5];
    }
  }
  close(RC);
}
#Calculate the WGS + Exome coverage and combined VAF
foreach my $coord (keys %pos){
  $pos{$coord}{normal_cov} = $pos{$coord}{normal_ref_rc} + $pos{$coord}{normal_var_rc};
  if($pos{$coord}{normal_cov}==0){
    $pos{$coord}{normal_vaf}=0;
  }else{
    $pos{$coord}{normal_vaf} = ($pos{$coord}{normal_var_rc} / $pos{$coord}{normal_cov})*100;
  }
  $pos{$coord}{tumor_cov} = $pos{$coord}{tumor_ref_rc} + $pos{$coord}{tumor_var_rc};
  if($pos{$coord}{tumor_cov}==0){
    $pos{$coord}{tumor_vaf}=0;
  }else{
    $pos{$coord}{tumor_vaf} = ($pos{$coord}{tumor_var_rc} / $pos{$coord}{tumor_cov})*100;
  }
}

#Determine the CNV difference value for each position.  Do this by intesecting the positions with CNV windows
#'intersectBed -a postions.bed -b cnv_windows.bed -f 1.0 -wa -wb'  
#The '-f 1.0' option should give the exons that are entirely overlapped by junctions
#The '-wa -wb' options, write the original coordinates for exons and junctions (as opposed to the merged coordinates).  
#Each overlaping pair will be reported as a seperate line

#Parse the CNV windows difference data and store
my %cnvs;
print BLUE, "\n\nImporting CNVs from file:\n\t$wgs_cnvs_file", RESET;
open (CNVS, "$wgs_cnvs_file") || die "\n\nCould not open CNVs file: $wgs_cnvs_file\n\n";
my $o = 0;
while(<CNVS>){
  if ($_ =~ /^\#|^CHR/){
    next();
  }
  $o++;
  chomp($_);
  my @line = split("\t", $_);
  my $chr = $line[0];
  my $start = $line[1];
  my $end = $start + (10000-1);
  my $window_coverage = $line[2] + $line[3];
  my $cnv_diff = $line[4];
  my $coord = $chr . ":" . "$start" . "-" . "$end";
  $cnvs{$coord}{order} = $o;
  $cnvs{$coord}{chr} = $chr;
  $cnvs{$coord}{start} = $start;
  $cnvs{$coord}{end} = $end;
  $cnvs{$coord}{window_read_coverage} = $window_coverage;
  $cnvs{$coord}{cnv_diff} = $cnv_diff;  
}
close(CNVS);

#Create a bed file of the CNV windows
print BLUE, "\n\nIntersecting CNV windows with target positions to get CNV difference values for each position", RESET;
my $windows_bed_file = $outdir . "cnv_windows.bed";
open (CNVBED, ">$windows_bed_file") || die "\n\nCould not open cnv windows bed for writing\n\n";
foreach my $coord (sort {$cnvs{$a}{order} <=> $cnvs{$b}{order}} keys %cnvs){
  print CNVBED "$cnvs{$coord}{chr}\t$cnvs{$coord}{start}\t$cnvs{$coord}{end}\t+\n";
}
close(CNVBED);

#Create a bed file of the positions
my $positions_bed_file = $outdir . "positions.bed";
open (POSBED, ">$positions_bed_file") || die "\n\nCould not open positions bed for writing\n\n";
foreach my $coord (sort keys %pos){
  print POSBED "$pos{$coord}{chr}\t$pos{$coord}{start}\t$pos{$coord}{end}\t+\n";  
}
close (POSBED);


#Identify overlaps between the windows and the positions using BEDTools
$bedtools_intersect_file = $outdir . "cnv_windows.intersect.positions.bed";
#$bedtools_cmd = "$bedtools_bin_dir"."intersectBed -a $positions_bed_file -b $windows_bed_file -f 1.0 -wa -wb > $bedtools_intersect_file";
$bedtools_cmd = "$bedtools_bin_dir"."intersectBed -a $positions_bed_file -b $windows_bed_file -wa -wb > $bedtools_intersect_file";
print YELLOW, "\n\n\t$bedtools_cmd", RESET;
Genome::Sys->shellcmd(cmd => $bedtools_cmd, output_files=>["$bedtools_intersect_file"]);

open (INT, "$bedtools_intersect_file") || die "\n\nCould not open bedtools intersect results file: $bedtools_intersect_file\n\n";
$intersect_count = 0;
while(<INT>){
  chomp($_);
  my @line = split("\t", $_);
  $intersect_count++;
  my $chr1 = $line[0];
  my $start1 = $line[1];
  my $end1 = $line[2];
  my $coord1 = $chr1 . ":" . "$start1" . "-" . "$end1";

  my $chr2 = $line[4];
  my $start2 = $line[5];
  my $end2 = $line[6];
  my $coord2 = $chr2 . ":" . "$start2" . "-" . "$end2";

  unless ($pos{$coord1} && $cnvs{$coord2}){
    print RED, "\n\nCould not match coord from intersect BED to original list of target positions or CNV windows\n\n", RESET;
    exit(1);
  }
  $pos{$coord1}{cnv_window_coord} = $coord2;
  $pos{$coord1}{cnv_window_read_coverage} = $cnvs{$coord2}{window_read_coverage};
  $pos{$coord1}{cnv_window_diff} = $cnvs{$coord2}{cnv_diff};
}
close(INT);

#Check counts to make sure an the window was found for every target position
my $positions_count = keys %pos;
unless ($positions_count == $intersect_count){
  print RED, "\n\nDid not identify a CNV window for every target position. Positions = $positions_count  CNV intersections = $intersect_count\n\n", RESET;
  exit(1);
}


#8.) Apply the SNP filters.  DO NOT filter any Tier1 Somatic SNVs
print BLUE, "\n\nApplying basic filters to SNP positions and removing those with low coverage, VAF that is inconsistent with heterozygous status, CNV diff, etc. (tier1 SNVs will not be filtered)", RESET;

#NOTE: To help determine what these cutoff values should be, I first ran the script without any filtering below and viewed distributions of the data in R
#In particular, I examined the distributions of 'normal_cov', 'normal_vaf', and 'cnv_window_diff'
if ($filter){
  print BLUE, "\n\tApplying the filtering step now ...", RESET;

  foreach my $coord (sort keys %pos){
    #Skip all SNVs as these can not be filtered
    unless ($pos{$coord}{class} eq "Control_Het_SNP" || $pos{$coord}{class} eq "Deletion_Het_SNP"){
      next();
    }
    #Apply a minimum read coverage cutoff (using the reads observed in WGS + Exome of Normal only) - only reads passing the minimum mapping quality cutoff are counted
    unless ($pos{$coord}{normal_cov} >= $min_normal_read_coverage){
      delete($pos{$coord});
      next();
    }
    #Apply a maximum read coverage cutoff.  Base this on the distribution of coverage observed for the Tier1 SNVs.  Use the median of this distribution + 1.5 times the interquartile range (i.e. standard definition of outlier)
    #For this data set this value is ~500 reads
    unless ($pos{$coord}{normal_cov} <= $max_normal_read_coverage){
      delete($pos{$coord});
      next();
    }

    #Apply a VAF min/max cutoff.  Base this on the distribution of VAFs observed for all normal heterozygous SNPs
    #The VAF value must be from the *normal* sample
    unless (($pos{$coord}{normal_vaf} >= (50-$het_vaf_diff_var)) && ($pos{$coord}{normal_vaf} <= (50+$het_vaf_diff_var))){
      delete($pos{$coord});
      next();
    }

    #For SNPs in control regions only.  Apply the same VAF min/max cutoff in the *tumor* sample
    #If these control region SNPs are really heterozygous and are really in a copy number neutral region, their VAF should always be close to 50%
    if ($pos{$coord}{class} eq "Control_Het_SNP"){
      unless (($pos{$coord}{tumor_vaf} >= (50-$het_vaf_diff_var)) && ($pos{$coord}{tumor_vaf} <= (50+$het_vaf_diff_var))){
        delete($pos{$coord});
        next();
      }
    }

    #Apply a CNV difference cutoff to those SNPs in the control regions so that they can not correspond to a CNV window that has a difference between tumor and normal
    if ($pos{$coord}{class} eq "Control_Het_SNP"){
      unless (($pos{$coord}{cnv_window_diff} >= (0-$cnv_control_window_diff_var)) && ($pos{$coord}{cnv_window_diff} <= (0+$cnv_control_window_diff_var))){
        delete($pos{$coord});
        next();
      }
    }

    #Apply a CNV difference cutoff to those SNPs in the deletion regions so that they must correspond to a CNV window that has a loss between tumor and normal
    if ($pos{$coord}{class} eq "Deletion_Het_SNP"){
      unless ($pos{$coord}{cnv_window_diff} <= (0-$cnv_deletion_window_diff_var)){
        delete($pos{$coord});
        next();
      }
    }
  }
}else{
  print YELLOW, "\n\tSkipping the filtering step", RESET;
}
$new_pos_count = keys %pos;
print BLUE, "\n\tCurrent list includes $new_pos_count positions", RESET;


#9.) Select heterozygous SNPs from within the deleted regions and control regions to fill the capacity of our capture reagent
#We have ~2500 SNV postions that must be included
#We have ~236 FLT3 SNP positions
#We have 27 deletion regions where we want ~75 SNP positions each
#Final design = 2500 SNV/InDel + 236 FLT3 specific + 2025 deletion regions + 500 control = 5261
my $snv_probe_count = $snv_count;

#For the deletions, attempt to select up to 12 positions per deletion
#Once again bin the positions by their regions
my %deletion_region_snps;
foreach my $coord (sort keys %pos){
  my $region_name = $pos{$coord}{region_name};
  my $class = $pos{$coord}{class};
  #Only thin the deletion SNPs
  unless ($class eq "Deletion_Het_SNP"){
    next();
  }
   if ($deletion_region_snps{$region_name}){
    my $snps_ref = $deletion_region_snps{$region_name}{snps};
    $snps_ref->{$coord}->{chr} = $pos{$coord}{chr};
    $snps_ref->{$coord}->{start} = $pos{$coord}{start};
    $snps_ref->{$coord}->{end} = $pos{$coord}{end};
  }else{
    my %snps;
    $snps{$coord}{chr} = $pos{$coord}{chr};
    $snps{$coord}{start} = $pos{$coord}{start};
    $snps{$coord}{end} = $pos{$coord}{end};
    $deletion_region_snps{$region_name}{snps} = \%snps;
  }
}

print BLUE, "\n\nThinning number of *deletion* SNPs to the final target number ($probes_per_deletion per deletion):", RESET;
my %selected_deletion_region_snps;
foreach my $region_name (sort keys %deletion_region_snps){
  my $region_desc = $region_name_list{$region_name};
  my $region_snps = $deletion_region_snps{$region_name}{snps};
  my $region_snps_count = keys %{$region_snps};
  print BLUE, "\n\tProcessing region: $region_name ($region_desc containing $region_snps_count possible SNPs)", RESET;

  #If the region has less than $snps_per_region SNPs, do nothing 
  if ($region_snps_count > $probes_per_deletion){
    my $x = $region_snps_count/$probes_per_deletion;
    my $skip_value = int($x);
    print BLUE, "\n\t\tThinning to every $skip_value th other SNP", RESET;

    #Otherwise create a sorted array of the SNPs and thin them
    my @snps;
    foreach my $coord (sort {$region_snps->{$a}->{start} <=> $region_snps->{$b}->{start}} keys %{$region_snps}){
      push(@snps, $coord);
    }

    #Get the thinned SNPs list
    my %thinned_snps;
    for (my $i = 1; $i <= $region_snps_count; $i += $skip_value){
      my $coord = $snps[$i-1];
      $thinned_snps{$coord}=1;
    }

    #Now go through the hash of possible SNPs again and delete those that are not part of the thinned list
    foreach my $coord (sort {$region_snps->{$a}->{start} <=> $region_snps->{$b}->{start}} keys %{$region_snps}){
      unless($thinned_snps{$coord}){
        delete($region_snps->{$coord});
        next();
      }
    }
    my $new_region_snps_count = keys %{$region_snps};
    print BLUE, "\n\t\tRemaining SNPs in this region:  $new_region_snps_count", RESET;
  }else{
    print BLUE, "\n\t\tThinning not needed - skipping", RESET;
  }

  foreach my $coord (sort {$region_snps->{$a}->{start} <=> $region_snps->{$b}->{start}} keys %{$region_snps}){
    $selected_deletion_region_snps{$coord} = 1;
  }
}
my $final_deletion_probe_count = keys %selected_deletion_region_snps;
my $target_control_probe_count = $target_probe_count - ($snv_probe_count + $final_deletion_probe_count);


#Now go through the original positions list and remove deletion SNPs that were not selected
foreach my $coord (sort keys %pos){
  my $region_name = $pos{$coord}{region_name};
  my $class = $pos{$coord}{class};
  #Only thin the deletion SNPs
  unless ($class eq "Deletion_Het_SNP"){
    next();
  }
  unless ($selected_deletion_region_snps{$coord}){
    delete($pos{$coord});
  }
}
$new_pos_count = keys %pos;
print BLUE, "\n\n\tSelected $final_deletion_probe_count positions in deletion regions", RESET;
print BLUE, "\n\tCurrent list includes $new_pos_count positions", RESET;


#For the controls, select evenly across all control regions
print BLUE, "\n\nThinning number of *control* SNPs to the final target number ($target_control_probe_count total across genome):", RESET;
my %control_region_snps;
foreach my $coord (sort keys %pos){
  my $region_name = $pos{$coord}{region_name};
  my $class = $pos{$coord}{class};
  #Only thin the control SNPs
  unless ($class eq "Control_Het_SNP"){
    next();
  }
  $control_region_snps{$coord} = 1;
}
my $control_region_snps_count = keys %control_region_snps;
print BLUE, "\n\tStarting with a total of $control_region_snps_count control region SNPs", RESET;
if ($control_region_snps_count > $target_control_probe_count){

  my $x = $control_region_snps_count/$target_control_probe_count;
  my $skip_value = int($x);
  print BLUE, "\n\t\tThinning to every $skip_value th other SNP", RESET;

  #Otherwise create a sorted array of the SNPs and thin them
  my @snps;
  foreach my $coord (sort keys %control_region_snps){
    push(@snps, $coord);
  }

  #Get the thinned SNPs list
  my %thinned_snps;
  for (my $i = 1; $i <= $control_region_snps_count; $i += $skip_value){
    my $coord = $snps[$i-1];
    $thinned_snps{$coord}=1;
  }

  #Now go through the hash of possible SNPs again and delete those that are not part of the thinned list
  foreach my $coord (sort keys %control_region_snps){
    unless($thinned_snps{$coord}){
      delete($control_region_snps{$coord});
      next();
    }
  }
  my $new_region_snps_count = keys %control_region_snps;
  print BLUE, "\n\t\tRemaining SNPs in this region:  $new_region_snps_count", RESET;

  foreach my $coord (sort keys %pos){
    my $class = $pos{$coord}{class};
    #Only thin the control SNPs
    unless ($class eq "Control_Het_SNP"){
      next();
    }
    unless ($control_region_snps{$coord}){
      delete($pos{$coord});
    }
  }
}
$control_region_snps_count = keys %control_region_snps;
$new_pos_count = keys %pos;
print BLUE, "\n\n\tSelected $control_region_snps_count positions in control regions", RESET;
print BLUE, "\n\tCurrent list includes $new_pos_count positions", RESET;


#10.) Produce a final spreadsheet summarizing all positions selected
my $final_target_positions_file = $outdir . "final_target_positions.tsv";
print BLUE, "\n\nPrinting final list of candidate positions to: $final_target_positions_file", RESET;
open (TSV, ">$final_target_positions_file") || die "\n\nCould not open final target positions outfile for writing: $final_target_positions_file\n\n";
print TSV "chr\tstart\tend\tref_base\tvar_base\tgene_id\ttrans_id\tvar_type\taa_change\tclass\tregion_name\tnormal_ref_rc\tnormal_var_rc\tnormal_cov\tnormal_vaf\ttumor_ref_rc\ttumor_var_rc\ttumor_cov\ttumor_vaf\tcnv_window_coord\tcnv_window_read_coverage\tcnv_window_diff\n";
foreach my $coord (sort keys %pos){
  print TSV "$pos{$coord}{chr}\t$pos{$coord}{start}\t$pos{$coord}{end}\t$pos{$coord}{ref_base}\t$pos{$coord}{var_base}\t$pos{$coord}{gene_id}\t$pos{$coord}{trans_id}\t$pos{$coord}{var_type}\t$pos{$coord}{aa_change}\t$pos{$coord}{class}\t$pos{$coord}{region_name}\t$pos{$coord}{normal_ref_rc}\t$pos{$coord}{normal_var_rc}\t$pos{$coord}{normal_cov}\t$pos{$coord}{normal_vaf}\t$pos{$coord}{tumor_ref_rc}\t$pos{$coord}{tumor_var_rc}\t$pos{$coord}{tumor_cov}\t$pos{$coord}{tumor_vaf}\t$pos{$coord}{cnv_window_coord}\t$pos{$coord}{cnv_window_read_coverage}\t$pos{$coord}{cnv_window_diff}\n";
}
close(TSV);

print "\n\n";
exit();


########################################################################################################################################
#Import SNVs/Indels from a positions bed file                                                                                          #
########################################################################################################################################
sub importPositions{
  my %args = @_;
  my $pos_ref = $args{'-pos_ref'};
  my $pos_file = $args{'-pos_file'};
  my $class = $args{'-class'};
  my $region_name_prefix = $args{'-region_name_prefix'};

  print BLUE, "\n\tGetting positions from $pos_file", RESET;
  my $pos_count = 0;

  open (SNVS, "$pos_file") || die "\n\nCould not open SNV file\n\n";
  while(<SNVS>){
    chomp($_);
    my @line = split("\t", $_);
    my $chr = $line[0];
    #Skip the 'GL' chromosomes that will be difficult to tie into other data in the following steps:
    #Also skip variants on ChrY
    if ($chr =~ /^GL|^Y/){
      next();
    }
    $pos_count++;
    my $start = $line[1]+1;
    my $end = $line[2];
    my $ref_var = $line[3];
    my $ref_base;
    my $var_base;
    if ($ref_var =~ /(\S+)\/(\S+)/){
      $ref_base = $1;
      $var_base = $2;
    }else{
      print RED, "\n\nCould not resolve ref/var base from string: $ref_var\n\n", RESET;
      exit(1);
    }

    #my $ref_base = $line[3];
    my $coord = "$chr" . ":". $start . "-" . "$end";
    my $region_name = "$region_name_prefix"."$pos_count";
    $pos{$coord}{chr} = $chr;
    $pos{$coord}{start} = $start;
    $pos{$coord}{end} = $end;
    $pos{$coord}{ref_base} = $ref_base;
    $pos{$coord}{var_base} = $var_base;
    $pos{$coord}{gene_id} = "NA";
    $pos{$coord}{trans_id} = "NA";
    $pos{$coord}{var_type} = "NA";
    $pos{$coord}{aa_change} = "NA";
    $pos{$coord}{class} = $class;
    $pos{$coord}{region_name} = $region_name;
  }
  close (SNVS);
  $pos_count = keys %pos;
  print BLUE, "\n\t\tFound $pos_count unique positions so far", RESET;

  return();
}


########################################################################################################################################
#Get BAM read counts for all positions                                                                                                 #
########################################################################################################################################
sub getBamReadCounts{
  my %args = @_;
  my $outdir = $args{'-outdir'};
  my $outfile_name = $args{'-outfile_name'};
  my $bam_file = $args{'-bam_file'};
  my $variant_file = $args{'-variant_file'};
  my $genome_build = $args{'-genome_build'};
  my $min_quality_score = $args{'-min_quality_score'};

  print BLUE, "\n\n\tChecking BAM results for $outfile_name", RESET;
  #Check to see if the output bam read counts file is already present.
  #If it is not present, run the BAM read counts command on all positions the input variant file
  my $outfile = "$outdir"."$outfile_name";
  my $temp_bam_outfile = $outfile . ".tmp";
  my $stdout_file = "$outdir"."$outfile_name".".stdout";
  my $stderr_file = "$outdir"."$outfile_name".".stderr";

  if (-e $outfile){

    #If the BAM read counts file is present:
    #1.) Grab the results for all positions that are in the input variant file.
    #2.) Store these in memory.
    my %variants;
    open (VAR, "$variant_file") || die "\n\nCould not open variants file: $variant_file\n\n";
    my $o = 0;
    while(<VAR>){
      $o++;
      chomp($_);
      my @line = split("\t", $_);
      my $chr = $line[0];
      my $start = $line[1];
      my $end = $line[2];
      my $coord = $chr . ":" . $start . "-" . $end;
      $variants{$coord}{chr} = $chr;
      $variants{$coord}{start} = $start;
      $variants{$coord}{end} = $end;
      $variants{$coord}{ref_base} = $line[3];
      $variants{$coord}{var_base} = $line[4];
      $variants{$coord}{line} = "$_\n";
      $variants{$coord}{order} = $o;
    }
    close(VAR);

    my %bam_results;
    open (BAM, "$outfile") || die "\n\nCould not open BAM results file: $outfile\n\n";
    while(<BAM>){
      chomp($_);
      my @line = split("\t", $_);
      my $chr = $line[0];
      my $start = $line[1];
      my $coord = $chr . ":" . $start . "-" . $start;

      #Only store BAM results for coordinates in the target variant list
      if ($variants{$coord}){
        $bam_results{$coord}{line} = "$_\n";
      }
    }
    close (BAM);

    #3.) Check for missing variants by comparing the complete variant list to the current bam results
    my %missing_variants;
    foreach my $coord (keys %variants){
      unless ($bam_results{$coord}){
        $missing_variants{$coord} = 1;
      }
    }
    my $missing_variant_count = keys %missing_variants;

    if ($missing_variant_count > 0){
      print YELLOW, "\n\tFound $missing_variant_count missing variants comparing the variant list to the BAM read counts results - updating\n", RESET;

      #4.) Make a new temp variants file with the positions that do not have BAM read counts
      my $temp_var_file = $variant_file . ".tmp";
      my $temp_var_file_sorted = $variant_file . ".sorted.tmp";
      open (VAR, ">$temp_var_file") || die "\n\nCould not open missing variants temp file: $temp_var_file\n\n";
      foreach my $coord (sort keys %missing_variants){
        print VAR "$variants{$coord}{line}";
      }
      close(VAR);

      #Sort the result temp variants file
      my $sort_cmd = "sort -k 1,1 -k 2n,2n -k 3n,3n $temp_var_file > $temp_var_file_sorted";
      print YELLOW, "\n\t$sort_cmd", RESET;
      Genome::Sys->shellcmd(cmd => $sort_cmd, output_files=>["$temp_var_file_sorted"]);

      #5.) Run the BAM read counts command on only these positions
      my $bamrc_cmd = "gmt analysis coverage bam-readcount  --bam-file=$bam_file  --output-file=$temp_bam_outfile  --variant-file=$temp_var_file_sorted  --genome-build=$genome_build  --min-quality-score=$min_quality_score  1>$stdout_file  2>$stderr_file";
      print YELLOW, "\n\t$bamrc_cmd", RESET;
      Genome::Sys->shellcmd(cmd => $bamrc_cmd, output_files=>["$temp_bam_outfile"]);

      #6.) Parse the result and merge into the result already stored in memory
      open (BAM, "$temp_bam_outfile") || die "\n\nCould not open temp BAM results file: $temp_bam_outfile\n\n";
      while(<BAM>){
        chomp($_);
        my @line = split("\t", $_);
        my $chr = $line[0];
        my $start = $line[1];
        my $coord = $chr . ":" . $start . "-" . $start;
        $bam_results{$coord}{line} = "$_\n";
      }
      close (BAM);

      #7.) Make sure a BAM result was found for every variant in the original file
      my $variant_count = keys %variants;
      my $bam_results_count = keys %bam_results;
      unless ($variant_count == $bam_results_count){
        print RED, "\n\nStill do not have a BAM result for all variants\n\n", RESET;
        exit();
      }

      #8.) Merge the old BAM read counts file with the new one and write over the old BAM results file with the updated results
      print BLUE, "\n\tUpdating old BAM file to remove entries for variants not in the variants file and add entries for variants that were missing in the BAM results file", RESET;
      open (BAM, ">$outfile") || die "\n\nCould not open BAM file for updating: $outfile\n\n";
      foreach my $coord (sort keys %bam_results){
        print BAM "$bam_results{$coord}{line}";
      }
      close (BAM);

      #9.) Clean up the temp files
      my $rm_cmd = "rm -f $temp_var_file $temp_var_file_sorted $temp_bam_outfile";
      Genome::Sys->shellcmd(cmd => $rm_cmd);
    }
  }else{
    my $bamrc_cmd = "gmt analysis coverage bam-readcount  --bam-file=$bam_file  --output-file=$outfile  --variant-file=$variant_file  --genome-build=$genome_build  --min-quality-score=$min_quality_score  1>$stdout_file  2>$stderr_file";
    print YELLOW, "\n\t$bamrc_cmd", RESET;
    Genome::Sys->shellcmd(cmd => $bamrc_cmd, output_files=>["$outfile"]);
  }
  return($outfile);
}

