#!/usr/bin/env genome-perl
use strict;
use warnings;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use above 'Genome';
use Genome::Model::ClinSeq::Util qw(:all);

my $working_dir = "/gscmnt/sata132/techd/mgriffit/luc/rnaseq_vs_snv/";
my $original_files_dir = "/gscmnt/sata132/techd/mgriffit/luc/rnaseq_vs_snv/original/";
my $liftover_files_dir = "/gscmnt/sata132/techd/mgriffit/luc/rnaseq_vs_snv/liftover/";
my $distinct_coords_file_hg18 = "/gscmnt/sata132/techd/mgriffit/luc/rnaseq_vs_snv/distinct.coords.hg18.bed";
my $distinct_coords_file_hg19 = "/gscmnt/sata132/techd/mgriffit/luc/rnaseq_vs_snv/distinct.coords.hg19.bed";
my $distinct_coords_file_unmapped = "/gscmnt/sata132/techd/mgriffit/luc/rnaseq_vs_snv/distinct.coords.unmapped.bed";
#http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz
my $chain_file = "/gscmnt/sata132/techd/mgriffit/luc/rnaseq_vs_snv/hg18ToHg19.over.chain";
my $snv_summary_files_dir_hg18 = "/gscmnt/sata132/techd/mgriffit/luc/rnaseq_vs_snv/summary_snv_hg18/";
my $snv_summary_files_dir_hg19 = "/gscmnt/sata132/techd/mgriffit/luc/rnaseq_vs_snv/summary_snv_hg19/";
my $readcounts_files_dir_hg18 = "/gscmnt/sata132/techd/mgriffit/luc/rnaseq_vs_snv/readcounts_hg18/"; #For WGS readcounts only
my $readcounts_files_dir_hg19 = "/gscmnt/sata132/techd/mgriffit/luc/rnaseq_vs_snv/readcounts_hg19/"; #For RNAseq readcounts only
my $readcounts_files_merge_dir_hg19 = "/gscmnt/sata132/techd/mgriffit/luc/rnaseq_vs_snv/readcounts_merge_hg19/";
my $bam_locations_dir = "/gscmnt/sata132/techd/mgriffit/luc/rnaseq_vs_snv/bam_locations/";
my $bam_locations_dir_individual = "/gscmnt/sata132/techd/mgriffit/luc/rnaseq_vs_snv/bam_locations/individual/";
my $all_bams_file = $bam_locations_dir."all_bams.txt";
my $entrez_ensembl_data = &loadEntrezEnsemblData();
my $rnaseq_data_dir = "/gscmnt/sata132/techd/mgriffit/luc/rnaseq_vs_snv/rnaseq_data/";

#Make a new copy of the original file
unless (-e $original_files_dir && -d $original_files_dir){
  mkdir($original_files_dir);
}
#Make a new copy of the original file
unless (-e $liftover_files_dir && -d $liftover_files_dir){
  mkdir($liftover_files_dir);
}
#Make a new copy of the original file
unless (-e $snv_summary_files_dir_hg18 && -d $snv_summary_files_dir_hg18){
  mkdir($snv_summary_files_dir_hg18);
}
unless (-e $snv_summary_files_dir_hg19 && -d $snv_summary_files_dir_hg19){
  mkdir($snv_summary_files_dir_hg19);
}
#Make a new copy of the original file
unless (-e $readcounts_files_dir_hg18 && -d $readcounts_files_dir_hg18){
  mkdir($readcounts_files_dir_hg18);
}
unless (-e $readcounts_files_dir_hg19 && -d $readcounts_files_dir_hg19){
  mkdir($readcounts_files_dir_hg19);
}
unless (-e $readcounts_files_merge_dir_hg19 && -d $readcounts_files_merge_dir_hg19){
  mkdir($readcounts_files_merge_dir_hg19);
}
unless (-e $rnaseq_data_dir && -d $rnaseq_data_dir){
  mkdir($rnaseq_data_dir);
}

#Load a list of SNV files
print BLUE, "\n\nLoad data:", RESET;

my %infiles;
my $infile_list = "/gscmnt/sata132/techd/mgriffit/luc/rnaseq_vs_snv/input_files.txt";
open (IN, "$infile_list") || die "\n\nCould not open infile: $infile_list\n\n";
my $c = 0;
while(<IN>){
  chomp($_);
  my $path = $_;
  my $result = `wc -l $path`;
  my $line_count;
  my $original_file_name;
  my $patient;

  if ($result =~ /^(\d+)/){
    $line_count = $1;
  }else{
    print RED, "\n\nCould not get line count\n\n", RESET;
    exit();
  }
  if ($path =~ /\/(hc.*)$/){
    $original_file_name = $1;
  }else{
    print RED, "\n\nCould not determine file name from path: $path", RESET;
    exit();
  }
  if ($path =~ /(LUC\d+)/){
    $patient = $1;
  }else{
    print RED, "\n\nCould not determine patient name from path: $path", RESET;
    exit();
  }

  $c++;
  $infiles{$c}{patient} = $patient;
  $infiles{$c}{original_path} = $path;
  $infiles{$c}{line_count} = $line_count;
  $infiles{$c}{original_file_name} = $original_file_name;
  my $new_file_name = $patient . ".tier1.somatic.annot";
  $infiles{$c}{new_file_name} = $new_file_name;
  $infiles{$c}{new_path} = $original_files_dir . $new_file_name;
  $infiles{$c}{liftover_path} = $liftover_files_dir . $new_file_name;
  $infiles{$c}{summary_path_hg18} = $snv_summary_files_dir_hg18 . $new_file_name;
  $infiles{$c}{summary_path_hg19} = $snv_summary_files_dir_hg19 . $new_file_name;
  $infiles{$c}{readcounts_path_hg18} = $readcounts_files_dir_hg18 . $new_file_name;
  $infiles{$c}{readcounts_path_hg19} = $readcounts_files_dir_hg19 . $new_file_name;
  my $cp_cmd = "cp $infiles{$c}{original_path} $infiles{$c}{new_path}";
  print YELLOW, "\n\t$cp_cmd", RESET;
  Genome::Sys->shellcmd(cmd => $cp_cmd);

}
close(IN);

#Build a list of all unique positions encounted in all SNV files 
#Chromosome names should be formated as per UCSC: chr1
print BLUE, "\n\nImport unique hg18 coords", RESET;
my %coords_hg18;
foreach my $c (keys %infiles){
  my $path = $infiles{$c}{new_path};
  open (IN, "$path") || die "\n\nCould not open input file: $path\n\n";
  while(<IN>){
    chomp($_);
    my @line = split("\t", $_);
    my $chr = $line[0];
    my $start = $line[1];
    my $end = $line[2];
    if ($_ =~ /\w+\s+\d+\s+\d+/){
      my $cs = "$chr".":"."$start-$end";
      $coords_hg18{$cs}{chromosome} = "chr$chr";
      $coords_hg18{$cs}{start} = $start;
      $coords_hg18{$cs}{end} = $end;
    }
  }
  close(IN);
}
my $coord_count_hg18 = keys %coords_hg18;
print BLUE, "\n\tFound $coord_count_hg18 distinct coord positions", RESET;

#Write all distinct hg18 coords to a file in BED format for use with lift-over
print BLUE, "\n\nWriting a BED file of distinct hg18 coords for liftOver purposes", RESET;
open (OUT, ">$distinct_coords_file_hg18") || die "\n\nCould not open output file: $distinct_coords_file_hg18\n\n";
foreach my $cs (sort keys %coords_hg18){
  #adjust hg18 coords to be 0 based for liftOver step
  my $start = $coords_hg18{$cs}{start} - 1;
  my $end = $coords_hg18{$cs}{end};
  print OUT "$coords_hg18{$cs}{chromosome}\t$start\t$end\t$cs\n";
}
close(OUT);

#print Dumper %infiles;
#print Dumper %coords_hg18;

#Use liftOver to convert coordinates
#liftOver test.hg18.bed hg18ToHg19.over.chain test.hg19.bed test.unmapped.bed
print BLUE, "\n\nRunning liftOver", RESET;
my $liftover_cmd = "liftOver $distinct_coords_file_hg18 $chain_file $distinct_coords_file_hg19 $distinct_coords_file_unmapped";
print YELLOW, "\n\t$liftover_cmd\n", RESET;
Genome::Sys->shellcmd(cmd => $liftover_cmd);

#Parse the liftover BED file to get the converted hg19 coords
my %coords_hg19;
open(IN, "$distinct_coords_file_hg19") || die "\n\nCould not open input file: $distinct_coords_file_hg19\n\n";
while(<IN>){
  chomp($_);
  my @line = split("\t", $_);
  my $cs = $line[3];
  my $chr = $line[0];
  my $start = $line[1];
  my $end = $line[2];

  #In the following hash, the 'cs' key is still the same as it was in hg19, but the coords are as calculated by liftOver
  $coords_hg19{$cs}{chromosome} = "$chr";
  $coords_hg19{$cs}{start} = $start;
  $coords_hg19{$cs}{end} = $end;
}
close(IN);
my $coord_count_hg19 = keys %coords_hg19;
print BLUE, "\n\tFound $coord_count_hg19 distinct coord positions that could be mapped to hg19", RESET;


#Now Go through the original files and make new versions, swapping out the hg18 coords with hg19
print BLUE, "\n\nCreate new SNV files with coordinates converted to hg19", RESET;
foreach my $c (keys %infiles){
  my $path = $infiles{$c}{new_path};
  my $liftover_path = $infiles{$c}{liftover_path};
  open (IN, "$path") || die "\n\nCould not open input file: $path\n\n";
  open (OUT, ">$liftover_path") || die "\n\nCould not open input file: $liftover_path\n\n";
  my $header = 1;
  while(<IN>){
    chomp($_);
    my @line = split("\t", $_);
    if ($header){
      print OUT "$_\n";
      $header = 0;
      next();
    }
    my $chr = $line[0];
    my $start = $line[1];
    my $end = $line[2];
    my $cs = "$chr".":"."$start-$end";
    if ($coords_hg19{$cs}){
      #Convert coords
      my $new_chr = $coords_hg19{$cs}{chromosome};
      my $new_start = $coords_hg19{$cs}{start};
      my $new_end = $coords_hg19{$cs}{end};
      if ($new_chr =~ /^chr(.*)$/){
        $new_chr = $1;
      }else{
        print RED, "\n\nChr format not understood: $new_chr\n\n", RESET;
        exit();
      }

      #Undo the 0-base conversion need during the liftover step
      $line[0] = $new_chr;
      $line[1] = $new_end;
      $line[2] = $new_end;

      #Print out the modified line
      local $" = "\t";
      print OUT "@line\n";

    }else{
      print YELLOW, "\n\nUnmapped entry:\n$_", RESET;
    }
  }
  close(IN);
  close(OUT);
}

#Thin down the SNV files to remove undeeded columns - repeat this process for both the original hg18 files and the new hg19 files
print BLUE, "\n\nCreate new summarized SNV files with mapped gene names - hg18", RESET;
#Get a complete list of mutated genes, and the number of times each gene was mutated.   Use mapped gene names to allow them to be paired up with the RNA-seq analysis
my %mutant_genes;
my %mutant_patients;
foreach my $c (keys %infiles){
  my $path = $infiles{$c}{new_path};
  my $summary_path_hg18 = $infiles{$c}{summary_path_hg18};
  my $patient = $infiles{$c}{patient};
  #coord	        gene_name	mapped_gene_name	aa_changes	ref_base	var_base
  #7:36698835-36698835	AOAH	        AOAH	                p.H109L,p.H77L	T               A	
  open (IN, "$path") || die "\n\nCould not open input file: $path\n\n";
  open (OUT, ">$summary_path_hg18") || die "\n\nCould not open output file: $summary_path_hg18\n\n";
  my $header = 1;
  while(<IN>){
    chomp($_);
    my @line = split("\t", $_);
    if ($header){
      $header = 0;
      print OUT "coord\tgene_name\tmapped_gene_name\taa_changes\tref_base\tvar_base\n";
      next();
    }
    my $coord = "$line[0]:$line[1]-$line[1]";
    my $gene_name = $line[6];
    my $aa_changes = $line[15];
    my $ref_base = $line[3];
    my $var_base = $line[4];
    #Fix gene name
    my $mapped_gene_name = &fixGeneName('-gene'=>$gene_name, '-entrez_ensembl_data'=>$entrez_ensembl_data, '-verbose'=>0);
    print OUT "$coord\t$gene_name\t$mapped_gene_name\t$aa_changes\t$ref_base\t$var_base\n";

    my $sample = "$patient"."_T";
    $mutant_patients{$sample}{$mapped_gene_name}{status}=1;
    $mutant_patients{$sample}{$mapped_gene_name}{count}++;

    #Store the mutations in a unified mutation hash
    if ($mutant_genes{$mapped_gene_name}){
      my $patients = $mutant_genes{$mapped_gene_name}{patients};
      $patients->{$patient}=1;
      $mutant_genes{$mapped_gene_name}{patient_count}++;
    }else{
      $mutant_genes{$mapped_gene_name}{gene_name} = $gene_name;
      $mutant_genes{$mapped_gene_name}{patient_count} = 1;
      my %patients;
      $patients{$patient}=1;
      $mutant_genes{$mapped_gene_name}{patients} = \%patients;
    }
  
  }
  close (IN);
  close (OUT);
}

print BLUE, "\n\nCreate new summarized SNV files with mapped gene names - hg19", RESET;
foreach my $c (keys %infiles){
  my $path = $infiles{$c}{liftover_path};
  my $summary_path_hg19 = $infiles{$c}{summary_path_hg19};
  #coord	        gene_name	mapped_gene_name	aa_changes	ref_base	var_base
  #7:36698835-36698835	AOAH	        AOAH	                p.H109L,p.H77L	T               A	
  open (IN, "$path") || die "\n\nCould not open input file: $path\n\n";
  open (OUT, ">$summary_path_hg19") || die "\n\nCould not open output file: $summary_path_hg19\n\n";
  my $header = 1;
  while(<IN>){
    chomp($_);
    my @line = split("\t", $_);
    if ($header){
      $header = 0;
      print OUT "coord\tgene_name\tmapped_gene_name\taa_changes\tref_base\tvar_base\n";
      next();
    }
    my $coord = "$line[0]:$line[1]-$line[1]";
    my $gene_name = $line[6];
    my $aa_changes = $line[15];
    my $ref_base = $line[3];
    my $var_base = $line[4];
    #Fix gene name
    my $mapped_gene_name = &fixGeneName('-gene'=>$gene_name, '-entrez_ensembl_data'=>$entrez_ensembl_data, '-verbose'=>0);
    print OUT "$coord\t$gene_name\t$mapped_gene_name\t$aa_changes\t$ref_base\t$var_base\n";
  }
  close (IN);
  close (OUT);
}


#Set up the data paths files
my %bam_data;
my $bam_data_header;
my $header = 1;
open (BAM, "$all_bams_file") || die "\n\nCould not open bams data file: $all_bams_file\n\n";
while(<BAM>){
  chomp($_);
  my $line_record = $_;
  my @line = split("\t", $_);
  if ($header){
    $bam_data_header = $_;
    $header = 0;
    next();
  }
  my $patient = $line[0];
  my $sample_type = $line[1];
  my $data_type = $line[2];
  my $sample_data_type = "$sample_type"."_"."$data_type";
  if ($bam_data{$patient}){
    my $bams_ref = $bam_data{$patient}{bams};
    $bams_ref->{$sample_data_type}->{sample_type} = $sample_type;
    $bams_ref->{$sample_data_type}->{data_type} = $data_type;
    $bams_ref->{$sample_data_type}->{line_record} = $line_record;
  }else{
    my %tmp;
    $tmp{$sample_data_type}{sample_type} = $sample_type;
    $tmp{$sample_data_type}{data_type} = $data_type;
    $tmp{$sample_data_type}{line_record} = $line_record;
    $bam_data{$patient}{bams} = \%tmp;
  }
}
close(BAM);

#print Dumper %bam_data;

#Now write out individual BAM files, one for WGS data (hg18) and one for RNAseq data (hg19)
print BLUE, "\n\nCreate individual bam locations files for each patient", RESET;
foreach my $c (sort {$infiles{$a}->{patient} cmp $infiles{$b}->{patient}} keys %infiles){
  my $patient = $infiles{$c}{patient};
  my $bams = $bam_data{$patient}{bams};
  #print Dumper $bams;

  print BLUE, "\n\t$patient", RESET;
  if ($bams->{'Tumor_WGS'}){
    print BLUE, "\n\t\tFound WGS tumor BAM for $patient", RESET;
    my $wgs_bams_file = $bam_locations_dir_individual . "$patient"."_WGS_hg18.txt";
    $infiles{$c}{wgs_bam_locations_file} = $wgs_bams_file;
    open (OUT, ">$wgs_bams_file") || die "\n\nCould not open individual bams file for output: $wgs_bams_file\n\n";
    print OUT "$bam_data_header\n";
    foreach my $b (keys %{$bams}){
      if ($bams->{$b}->{data_type} eq 'WGS'){
        print OUT "$bams->{$b}->{line_record}\n";
      }
    }
    close(OUT);
  }else{
    print YELLOW, "\n\t\tNo Tumor_WGS bam defined for $patient", RESET;
  }
  if ($bams->{'Tumor_RNAseq'}){
    print BLUE, "\n\t\tFound RNAseq tumor BAM for $patient", RESET;
    my $rnaseq_bams_file = $bam_locations_dir_individual . "$patient"."_RNAseq_hg19.txt";
    $infiles{$c}{rnaseq_bam_locations_file} = $rnaseq_bams_file;
    open (OUT, ">$rnaseq_bams_file") || die "\n\nCould not open individual bams file for output: $rnaseq_bams_file\n\n";
    print OUT "$bam_data_header\n";
    foreach my $b (keys %{$bams}){
      if ($bams->{$b}->{data_type} eq 'RNAseq'){
        print OUT "$bams->{$b}->{line_record}\n";
      }
    }
    close(OUT);
  }else{
    print YELLOW, "\n\t\tNo Tumor_RNAseq bam defined for $patient", RESET;
  }
}


#Get BAM read counts for all positions in each file
#Unfortunately, the same hg18 vs. hg19 problem comes up again for the BAM read count step.  We need to use the original hg18 coordinates for the WGS data and then the hg19 coordinates for the RNAseq data
print BLUE, "\n\nGet BAM read counts for all BAMs all positions", RESET;
foreach my $c (sort {$infiles{$a}->{patient} cmp $infiles{$b}->{patient}} keys %infiles){
  my $patient = $infiles{$c}{patient};
  print BLUE, "\n\t$patient", RESET;

  #WGS BAMS first
  if ($infiles{$c}{wgs_bam_locations_file}){
    my $hg18_path = $infiles{$c}{summary_path_hg18};
    my $wgs_data_paths_file = $infiles{$c}{wgs_bam_locations_file};
    my $readcounts_path_hg18 = $infiles{$c}{readcounts_path_hg18};
    #my $wgs_readcounts_cmd = "/gscmnt/sata206/techd/git/genome/lib/perl/Genome/Model/ClinSeq/original-scripts/snv/getBamReadCounts.pl  --positions_file=$hg18_path  --data_paths_file=$wgs_data_paths_file  --output_file=$readcounts_path_hg18  --no_fasta_check=1";
    
    #Call BamReadCounts tools from genome model clin-seq - since it uses Bio::Db:Sam we need a hack to deal with genome versions...
    my $wgs_readcounts_cmd = "genome5.10.1 model clin-seq get-bam-read-counts  --positions-file=$hg18_path  --data-paths-file=$wgs_data_paths_file  --output-file=$readcounts_path_hg18  --no-fasta-check=1";

    print BLUE, "\n\tWGS", RESET;
    #print YELLOW, "\n\t$wgs_readcounts_cmd", RESET;
    if (-e $readcounts_path_hg18){
      print YELLOW, "\n\t\tFile already exists - skipping", RESET;
    }else{
      #System call to ClinSeq GetBamReadCounts tool
      #print YELLOW, "\n\t$wgs_readcounts_cmd", RESET;
      #Genome::Sys->shellcmd(cmd => $wgs_readcounts_cmd);

      #Or call directly using a Genome API call...
      Genome::Model::ClinSeq::Command::GetBamReadCounts->execute('positions_file'=>$hg18_path, 'data_paths_file'=>$wgs_data_paths_file, 'output_file'=>$readcounts_path_hg18, 'no_fasta_check'=>1);

    }

  }

  #Now RNAseq BAMS
  if ($infiles{$c}{rnaseq_bam_locations_file}){
    my $hg19_path = $infiles{$c}{summary_path_hg19};
    my $rnaseq_data_paths_file = $infiles{$c}{rnaseq_bam_locations_file};
    my $readcounts_path_hg19 = $infiles{$c}{readcounts_path_hg19};
    #my $rnaseq_readcounts_cmd = "/gscmnt/sata206/techd/git/genome/lib/perl/Genome/Model/ClinSeq/original-scripts/snv/getBamReadCounts.pl  --positions_file=$hg19_path  --data_paths_file=$rnaseq_data_paths_file  --output_file=$readcounts_path_hg19  --no_fasta_check=1";
    my $rnaseq_readcounts_cmd = "genome5.10.1 model clin-seq get-bam-read-counts  --positions-file=$hg19_path  --data-paths-file=$rnaseq_data_paths_file  --output-file=$readcounts_path_hg19  --no-fasta-check=1";

    print BLUE, "\n\tRNAseq", RESET;
    if (-e $readcounts_path_hg19){
      print YELLOW, "\n\t\tFile already exists - skipping", RESET;
    }else{
      #print YELLOW, "\n\t$rnaseq_readcounts_cmd", RESET;
      #Genome::Sys->shellcmd(cmd => $rnaseq_readcounts_cmd);

      #Or call directly using a Genome API call...
      Genome::Model::ClinSeq::Command::GetBamReadCounts->execute('positions_file'=>$hg19_path, 'data_paths_file'=>$rnaseq_data_paths_file, 'output_file'=>$readcounts_path_hg19, 'no_fasta_check'=>1);

    }
  }
}

#Get the WGS readcount results, convert coordinates from hg18 to hg19 and merge with the RNAseq readcount results
print BLUE, "\n\nConvert coordinates for WGS readcounts and merging with RNAseq readcounts results", RESET;
foreach my $c (sort {$infiles{$a}->{patient} cmp $infiles{$b}->{patient}} keys %infiles){
  my $patient = $infiles{$c}{patient};
  print BLUE, "\n\t$patient", RESET;

  #If both RNAseq and WGS files are ready, do the convert and merge
  if (-e $infiles{$c}{readcounts_path_hg18} && -e $infiles{$c}{readcounts_path_hg19}){
    print BLUE, "\n\t\tProceeding with merge", RESET;

    #Make a directory for this patient
    my $working_dir = "$readcounts_files_merge_dir_hg19"."$patient"."/";
    mkdir($working_dir);

    my $merged_file = $working_dir ."$patient".".readcounts.hg19.tsv";
    $infiles{$c}{merged_file} = $merged_file;
    $infiles{$c}{merged_file_dir} = $working_dir;

    my %merged;

    #Go through the hg19 RNA file - get the coords of interest and line records
    my $header_line1_string;
    my $header1 = 1;
    open (RNA, "$infiles{$c}{readcounts_path_hg19}") || die "\n\nCould not open readcounts file: $infiles{$c}{readcounts_path_hg19}\n\n";
    while(<RNA>){
      chomp($_);
      my @line = split("\t", $_);
      my $line_record = $_;
      if ($header1){
        $header1 = 0;
        $header_line1_string = $_;
        next();
      }
      my $coord = $line[0];
      $merged{$coord}{rnaseq_record} = $line_record;
    }
    close(RNA);

    #Now go through hg18 WGS file - convert coords to hg19, match to previously stored RNA record, and add on WGS extra columns
    my $header_line2_string;
    my @header_line2;
    my $header2 = 1;
    open (WGS, "$infiles{$c}{readcounts_path_hg18}") || die "\n\nCould not open readcounts file: $infiles{$c}{readcounts_path_hg18}\n\n";
    while(<WGS>){
      chomp($_);
      my @line = split("\t", $_);
      my $line_record = $_;
      if ($header2){
        $header2 = 0;
        @header_line2 = @line;
        $header_line2_string = $_;
        next();
      }
      my $coord = $line[0];
      if ($coords_hg19{$coord}){
        my $new_chr = $coords_hg19{$coord}{chromosome};
        my $new_start = $coords_hg19{$coord}{start};
        my $new_end = $coords_hg19{$coord}{end};
        if ($new_chr =~ /^chr(.*)$/){
          $new_chr = $1;
        }else{
          print RED, "\n\nChr format not understood: $new_chr\n\n", RESET;
          exit();
        }
        my $new_coord = "$new_chr".":"."$new_end-$new_end";
        
        if ($merged{$new_coord}){
          $merged{$new_coord}{wgs_record} = $line_record;
        }else{
          print RED, "\n\nCould not match coords: $coord vs $new_coord", RESET;
        }
      }else{
        print RED, "\n\nCould not find new coords for: $coord", RESET;
      }
    }
    close(WGS);

    #Print out the joined header and data lines
    #print Dumper %merged;
    open(OUT, ">$merged_file") || die "\n\nCould not open merge file for output: $merged_file\n\n";
    my $col_count = scalar(@header_line2);
    my @merge_columns = @header_line2[6 .. $col_count-1];
    my $merge_columns_string = join("\t", @merge_columns);
    print OUT "$header_line1_string\t$merge_columns_string\n";
    foreach my $coord (sort keys %merged){
      my $wgs_line_record = $merged{$coord}{wgs_record};
      my @wgs_line_record = split("\t", $wgs_line_record);
      my $col_count = scalar(@wgs_line_record);
      my @merge_columns = @wgs_line_record[6 .. $col_count-1];
      my $merge_columns_string = join("\t", @merge_columns);
      print OUT "$merged{$coord}{rnaseq_record}\t$merge_columns_string\n";
    }
    close(OUT);

  }else{
    print YELLOW, "\n\t\tOne of the two files is missing - skipping", RESET;
  }
}

#Get the RNAseq expression data
#/gscmnt/sata206/techd/git/genome/lib/perl/Genome/Model/ClinSeq/original-scripts/rnaseq/outlierGenesAbsolute.pl  --cufflinks_dir=/gscmnt/gc8001/info/model_data/2881616601/build117345038/expression/  --working_dir=/gscmnt/sata132/techd/mgriffit/luc/rnaseq_vs_snv/rnaseq_data/luc7/
print BLUE, "\n\nGetting RNAseq expression data", RESET;
foreach my $c (sort {$infiles{$a}->{patient} cmp $infiles{$b}->{patient}} keys %infiles){
  my $patient = $infiles{$c}{patient};
  $infiles{$c}{rnaseq_data_name} = $patient."_T";
  print BLUE, "\n\t$patient", RESET;
  if ($infiles{$c}{rnaseq_bam_locations_file}){
    my $rnaseq_locations_file = $infiles{$c}{rnaseq_bam_locations_file};
    open(IN, "$rnaseq_locations_file") || die "\n\nCould not open rnaseq_locations_file: $rnaseq_locations_file";
    my $header = 1;
    my $expression_dir;
    my $alignments_dir;
    while(<IN>){
      chomp($_);
      my @line = split("\t", $_);
      if ($header){
        $header = 0;
        next();
      }
      $expression_dir = "$line[4]"."expression/";
      $alignments_dir = "$line[4]"."alignments/";
      $infiles{$c}{rnaseq_expression_dir} = $expression_dir;
      $infiles{$c}{rnaseq_alignments_dir} = $alignments_dir;
    }
    close(IN);

    my $rnaseq_results_dir = "$rnaseq_data_dir"."$patient"."_T/";
    $infiles{$c}{rnaseq_data_dir} = $rnaseq_results_dir;
    my $rnaseq_data_cmd = "/gscmnt/sata206/techd/git/genome/lib/perl/Genome/Model/ClinSeq/original-scripts/rnaseq/outlierGenesAbsolute.pl  --cufflinks_dir=$expression_dir  --working_dir=$rnaseq_results_dir";

    if (-e $rnaseq_results_dir && -d $rnaseq_results_dir){
      print YELLOW, "\n\t\tRNAseq results dir already exists - skipping this step", RESET;
      $infiles{$c}{rnaseq_expression_matrix} = "$rnaseq_results_dir"."isoforms_merged/isoforms.merged.fpkm.expsort.tsv";
    }else{
      print BLUE, "\n\t\tGetting RNAseq expression matrix files", RESET;
      #print YELLOW, "\n\t$rnaseq_data_cmd", RESET;
      mkdir($rnaseq_results_dir);
      Genome::Sys->shellcmd(cmd => $rnaseq_data_cmd);
      $infiles{$c}{rnaseq_expression_matrix} = "$rnaseq_results_dir"."isoforms_merged/isoforms.merged.fpkm.expsort.tsv";
    }
  }
}


#Generate figure and statistics for the RNAseq to WGS variant allele frequency comparisons
#/gscmnt/sata206/techd/git/genome/lib/perl/Genome/Model/ClinSeq/original-scripts/snv/WGS_vs_Exome_vs_RNAseq_VAF_and_FPKM.R  /gscmnt/sata132/techd/mgriffit/hgs/test/ /gscmnt/sata132/techd/mgriffit/hgs/all1/snv/wgs_exome/snvs.hq.tier1.v1.annotated.compact.readcounts.tsv /gscmnt/sata132/techd/mgriffit/hgs/all1/rnaseq/tumor/absolute/isoforms_merged/isoforms.merged.fpkm.expsort.tsv
print BLUE, "\n\nSummarize SNV WGS and RNAseq variant allele frequencies", RESET;
foreach my $c (sort {$infiles{$a}->{patient} cmp $infiles{$b}->{patient}} keys %infiles){
  my $patient = $infiles{$c}{patient};
  print BLUE, "\n\t$patient", RESET;
  if ($infiles{$c}{rnaseq_expression_matrix} && $infiles{$c}{merged_file_dir} && $infiles{$c}{merged_file}){
    my $summary_dir = "$infiles{$c}{merged_file_dir}"."summary/";
    my $summary_cmd = "/gscmnt/sata206/techd/git/genome/lib/perl/Genome/Model/ClinSeq/original-scripts/snv/WGS_vs_Exome_vs_RNAseq_VAF_and_FPKM.R  $summary_dir $infiles{$c}{merged_file} $infiles{$c}{rnaseq_expression_matrix}";
    if (-e $summary_dir && -d $summary_dir){
      print YELLOW, "\n\t\tSummary stats dir already exists - skipping this step", RESET;
    }else{
      print BLUE, "\n\t\tGenerating figures and statistics for the RNAseq to WGS variant allele frequency comparisons", RESET;
      print YELLOW, "\n\t$summary_cmd", RESET;
      mkdir($summary_dir);
      Genome::Sys->shellcmd(cmd => $summary_cmd);
      $infiles{$c}{vaf_summary_dir} = $summary_dir;
    }
  }
}




#Inject the LUC9 Normal sample into the hash of samples
my $c_count = keys %infiles;
$c = $c_count+1;
$infiles{$c}{rnaseq_expression_matrix} = "/gscmnt/sata132/techd/mgriffit/luc/rnaseq_vs_snv/rnaseq_data/LUC9_N/isoforms_merged/isoforms.merged.fpkm.expsort.tsv";
$infiles{$c}{rnaseq_data_name} = "LUC9_N";
$infiles{$c}{rnaseq_data_dir} = "/gscmnt/sata132/techd/mgriffit/luc/rnaseq_vs_snv/rnaseq_data/LUC9_N/";
$infiles{$c}{patient} = "LUC9";
$infiles{$c}{rnaseq_alignments_dir} = "/gscmnt/gc7001/info/model_data/2880885516/build116087250/alignments/";


#Build an expression matrix for all RNA-seq libraries that are currently available...
print BLUE, "\n\nBuild expression matrix for all RNAseq libraries", RESET;
my $first_file = 1;
my %exp;
my %gene_map;
my %name_list;
foreach my $c (sort {$infiles{$a}->{rnaseq_data_name} cmp $infiles{$b}->{rnaseq_data_name}} keys %infiles){
  my $rnaseq_data_name = $infiles{$c}{rnaseq_data_name};
  my $patient = $infiles{$c}{patient};
  print BLUE, "\n\t$rnaseq_data_name", RESET;
  if ($infiles{$c}{rnaseq_expression_matrix}){
    $name_list{$rnaseq_data_name}=$patient;
    my $fpkm_file = $infiles{$c}{rnaseq_expression_matrix};

    my $header = 1;
    my %columns;
    open (FPKM, "$fpkm_file") || die "\n\nCould not open FPKM file: $fpkm_file\n\n";
    while(<FPKM>){
      chomp($_);
      my @line = split("\t", $_);
      if ($header){
        $header = 0;
        my $p = 0;
        foreach my $col (@line){
          $columns{$col}{position} = $p;
          $p++;
        }
        next();
      }
      if ($first_file){
        my $tracking_id = $line[$columns{'tracking_id'}{position}];
        my $mapped_gene_name = $line[$columns{'mapped_gene_name'}{position}];
        my $fpkm = $line[$columns{'FPKM'}{position}];
        my $rank = $line[$columns{'rank'}{position}];
        $gene_map{$tracking_id}{mapped_gene_name} = $mapped_gene_name;
        $exp{$tracking_id}{$rnaseq_data_name}{fpkm} = $fpkm;
        $exp{$tracking_id}{$rnaseq_data_name}{rank} = $rank;
      }else{
        my $tracking_id = $line[$columns{'tracking_id'}{position}];
        my $fpkm = $line[$columns{'FPKM'}{position}];
        my $rank = $line[$columns{'rank'}{position}];
        $exp{$tracking_id}{$rnaseq_data_name}{fpkm} = $fpkm;
        $exp{$tracking_id}{$rnaseq_data_name}{rank} = $rank;
      }
    }
    close(FPKM);
    $first_file = 0;
  }
}
my $gene_count = keys %exp;
my $name_count = keys %name_list;
print BLUE, "\nFound expression values for $gene_count genes for $name_count libraries", RESET;

my @name_list = sort keys %name_list;
my $name_list_string = join("\t", @name_list);

my $fpkm_matrix_file = "$rnaseq_data_dir"."fpkm_matrix.tsv";
my $rank_matrix_file = "$rnaseq_data_dir"."rank_matrix.tsv";
my $categorical_file = "$rnaseq_data_dir"."expression_categorical.tsv";

open(FPKM_MATRIX, ">$fpkm_matrix_file") || die "\n\nCould not open: $fpkm_matrix_file\n\n";
open(RANK_MATRIX, ">$rank_matrix_file") || die "\n\nCould not open: $rank_matrix_file\n\n";
open(CAT, ">$categorical_file") || die "\n\nCould not open: $categorical_file\n\n";

my $matrix_header = "tracking_id\tmapped_gene_name\tmutation_status\tmutation_count\tmutated_patient_list\t$name_list_string";
print FPKM_MATRIX "$matrix_header\n";
print RANK_MATRIX "$matrix_header\n";
my $categorical_header = "tracking_id\tmapped_gene_name\tpatient\tsample\tmutation_status\tmutation_count\tFPKM\tRANK";
print CAT "$categorical_header\n";

foreach my $g (sort keys %exp){
  my $mapped_gene_name = $gene_map{$g}{mapped_gene_name};
  my $mutation_count = 0;
  my $mutation_status = 0;
  my $patient_list = "na";
  if ($mutant_genes{$mapped_gene_name}){
    $mutation_status = 1;
    $mutation_count = $mutant_genes{$mapped_gene_name}{patient_count};
    my $pt = $mutant_genes{$mapped_gene_name}{patients};
    my @pt = keys %{$pt};
    $patient_list = join(",", @pt);
  }
  my @fpkms;
  my @ranks;
  foreach my $p (sort keys %name_list){
    push (@fpkms, $exp{$g}{$p}{fpkm});
    push (@ranks, $exp{$g}{$p}{rank});

    #Print out categorical entry to file
    my $mutation_status = 0;
    my $mutation_count = 0;
    if ($mutant_patients{$p}{$mapped_gene_name}){
      $mutation_status = $mutant_patients{$p}{$mapped_gene_name}{status};
      $mutation_count = $mutant_patients{$p}{$mapped_gene_name}{count};
    }
    print CAT "$g\t$mapped_gene_name\t$name_list{$p}\t$p\t$mutation_status\t$mutation_count\t$exp{$g}{$p}{fpkm}\t$exp{$g}{$p}{rank}\n";

  }
  my $fpkms_string = join("\t", @fpkms);
  my $ranks_string = join("\t", @ranks);
  
  my $exp_fpkm_string = "$g\t$mapped_gene_name\t$mutation_status\t$mutation_count\t$patient_list\t$fpkms_string";
  my $exp_rank_string = "$g\t$mapped_gene_name\t$mutation_status\t$mutation_count\t$patient_list\t$ranks_string";

  print FPKM_MATRIX "$exp_fpkm_string\n";
  print RANK_MATRIX "$exp_rank_string\n";

}
close(FPKM_MATRIX);
close(RANK_MATRIX);
close(CAT);


print BLUE, "\n\nSummarizing Tophat alignment dirs", RESET;
foreach my $c (sort {$infiles{$a}->{patient} cmp $infiles{$b}->{patient}} keys %infiles){
  my $patient = $infiles{$c}{patient};
  print BLUE, "\n\t$patient", RESET;

  my $rnaseq_results_dir = $infiles{$c}{rnaseq_data_dir};
  my $rnaseq_alignments_dir = $infiles{$c}{rnaseq_alignments_dir};
  my $tophat_summary_dir = $rnaseq_results_dir ."tophat_summary/";

  if (-e $tophat_summary_dir && -d $tophat_summary_dir){
    print YELLOW, "\n\tTophat summary dir already exists - skipping this step", RESET;    
  }else{
    mkdir($tophat_summary_dir);
    my $summary_cmd = "/gscmnt/sata206/techd/git/genome/lib/perl/Genome/Model/ClinSeq/original-scripts/qc/tophatAlignmentSummary.pl  --reference_fasta_file='/gscmnt/sata420/info/model_data/2857786885/build102671028/all_sequences.fa'  --reference_annotations_dir='/gscmnt/sata132/techd/mgriffit/reference_annotations/hg19/'  --working_dir=$tophat_summary_dir  --tophat_alignment_dir='$rnaseq_alignments_dir'  --verbose=1";
    print YELLOW, "\n\t$summary_cmd", RESET;
    Genome::Sys->shellcmd(cmd => $summary_cmd);
  }
}


print "\n\n";





exit();


