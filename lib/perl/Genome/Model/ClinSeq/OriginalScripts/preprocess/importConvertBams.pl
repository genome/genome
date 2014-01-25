#!/usr/bin/env genome-perl

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Benchmark;
use above 'Genome';

my $usage=<<INFO;

./importConvertBams.pl 

e.g.

./importConvertBams.pl  --working_dir=/gscmnt/gc2142/techd/alexa_seq/input_data/  --library_name='H_KH-540079-41-11865-cDNA-2-lib1' [trim="16,95"]

Use extra option --bam_file to specify a single bam file directly instead of looking up using the library name (need for external data)

Use extra option --trim if you want to trim the reads.  e.g. --trim="16,95" will trim the first 15 bases and the last 5 bases from both R1 and R2

Use --p to specify multiple CPUs for bzip2 compression

INFO

my $working_dir = '';
my $library_name = '';
my $bam_file = '';    #Start with a BAM file instead if the data was generated external to WASHU
my $p = '';
my $trim = '';

GetOptions ('working_dir=s'=>\$working_dir, 'library_name=s'=>\$library_name, 'bam_file=s'=>\$bam_file, 'p=i'=>\$p, 'trim=f'=>\$trim);

unless ($working_dir && $library_name){
  print RED, "\n\nParameters missing\n\n", RESET;
  print GREEN, "$usage", RESET;
  exit();
}

my $f=0;
my $l=0;
if ($trim){
  if ($trim =~ /(\d+)\,(\d+)/){
     $f = $1;
     $l = $2;
  }else{
    print RED, "\n\nFormat of trim specification is not correct: $trim (should be [firstbase,lastbase])\n\n", RESET;
    exit();
  }
}

my $bzip_bin = "/gscmnt/gc2142/techd/tools/pbzip2-1.1.4/pbzip2";

print BLUE, "\n\nProcessing library: $library_name", RESET;

#Get bam file paths for this library name
#genome instrument-data list solexa --show id,flow_cell_id,lane,index_sequence,sample_name,library_name,clusters,read_length,bam_path --filter library_name='H_KH-540079-41-11865-cDNA-2-lib1'
my $bam_search = "genome instrument-data list solexa --show id,flow_cell_id,lane,index_sequence,sample_name,library_name,clusters,read_length,bam_path --filter library_name='$library_name'";
my @result = `$bam_search`;
my %bams;
my $lane_list = '';

if ($bam_file){
  my $lane = 1;
  my $flowcell = "unknown";
  my $index;
  $bams{$library_name}{path} = $bam_file;
  $bams{$library_name}{lane} = $lane;
  if ($bam_file =~ /^\/\w+\/(\w+)/){
    $flowcell = $1;
  }
  $bams{$library_name}{flowcell} = $flowcell;

  $lane_list .= "\n\t$flowcell (lane $lane)";

  unless (-e $bam_file){
    print RED, "\n\nBAM file: $bam_file does not appear to exist - aborting", RESET;
    exit();
  }
}else{
  foreach my $line (@result){
    chomp($line);
    my ($path, $base, $flowcell, $lane, $index);
    if ($line =~ /(\S+\.bam)/){
      $path = $1;
    }else{
      next();
    }
    if ($line =~ /(\w+)\.bam/){
      $base = $1;
    }
    if ($base =~ /gerald\_(\w+)\_(\d+)$/){
      $flowcell = $1;
      $lane = $2;
    }elsif($base =~ /gerald\_(\w+)\_(\d+)\_(\w+)$/){
      $flowcell = $1;
      $lane = $2;
      $index = $3;
    }else{
      print RED, "\n\nCould not determine flowcell and lane from base name: $base\n\n", RESET;
    }	
    unless ($path && $base && $flowcell && $lane){
      print RED, "\n\nCould not determine base name, path, flowcell and lane from bam path: $line\n\n", RESET;
      exit();
    }
    $bams{$base}{path} = $path;
    $bams{$base}{flowcell} = $flowcell;
    $bams{$base}{lane} = $lane;
    $lane_list .= "\n\t$flowcell (lane $lane)";
  }
}
my $lane_count = keys %bams;
print BLUE, "\n\tFound $lane_count lane(s) of data: $lane_list", RESET;

#Get start time
my $t0 = new Benchmark;

#Write out a bash script that:
#1.) Creates a library analysis dir within the working dir
my $lib_dir = "$working_dir"."$library_name";
unless (-e $lib_dir && -d $lib_dir){
  print BLUE, "\n\nMaking directories", RESET;
  mkdir ($lib_dir);
}

#3.) Converts the BAM files to fastq files and writes these to the data dir, then compress them to save space
foreach my $name (sort keys %bams){
  my $bampath = $bams{$name}{path};
  my $flowcell = $bams{$name}{flowcell};
  my $lane = $bams{$name}{lane};
  my $data_dir = "$lib_dir/$flowcell"."_"."$lane/";

  print BLUE, "\n\nConverting BAM to fastq and the compressing: $flowcell $lane $bampath", RESET;
  unless (-e $data_dir && -d $data_dir){
    mkdir($data_dir);
  }
  my $r1_fastq = "$data_dir"."s_"."$lane"."_1_sequence.txt";
  my $r1_fastq_temp = "$r1_fastq".".tmp";
  my $r1_fastq_bz = "$r1_fastq".".bz2";
  my $r2_fastq = "$data_dir"."s_"."$lane"."_2_sequence.txt";
  my $r2_fastq_temp = "$r2_fastq".".tmp";
  my $r2_fastq_bz = "$r2_fastq".".bz2";
  $bams{$name}{r1_fastq} = $r1_fastq;
  $bams{$name}{r2_fastq} = $r2_fastq;
  $bams{$name}{r1_fastq_bz} = $r1_fastq_bz;
  $bams{$name}{r2_fastq_bz} = $r2_fastq_bz;
  $bams{$name}{data_dir} = $data_dir;
  $bams{$name}{lane} = $lane;
  $bams{$name}{flowcell} = $flowcell;
  $bams{$name}{line_count} = 0;

  if (-e $r1_fastq_bz && -e $r2_fastq_bz){
    print YELLOW, "\n\n\tCompressed fastq files already present, existing input files for this lane will be used\n\n", RESET;
  }else{
    print BLUE, "\n\nConverting BAM to FASTQ using picard:", RESET;
    #gmt picard sam-to-fastq --input=$bampath --fastq=$r1_fastq --fastq2=$r2_fastq --fragment-fastq=/dev/null --maximum-memory=12  --maximum-permgen-memory=256
    #OR
    #java -Xmx12g -XX:MaxPermSize=256m  -cp $ENV{GENOME_SW_LEGACY_JAVA}/samtools/picard-tools-1.22/sam-1.22.jar:$ENV{GENOME_SW_LEGACY_JAVA}/samtools/picard-tools-1.22/picard-1.22.jar:/gsc/scripts/opt/genome/snapshots/stable/genome-1035/lib/perl/Genome/Model/Tools/Picard/SamToFastq/GCSamToFastq.jar edu.wustl.genome.samtools.GCSamToFastq  INPUT='$bampath' FASTQ='$r1_fastq' SECOND_END_FASTQ='$r2_fastq' FRAGMENT_FASTQ='/dev/null'
    my $rm_cmd = "rm -f $r1_fastq $r2_fastq $r1_fastq_bz $r2_fastq_bz $r1_fastq_temp $r2_fastq_temp";
    my $cmd_picard = "gmt picard sam-to-fastq --input=$bampath --fastq=$r1_fastq --fastq2=$r2_fastq --fragment-fastq=/dev/null --maximum-memory=12  --maximum-permgen-memory=256";
    print BLUE, "\n$cmd_picard\n\n", RESET;
    my $exit_code_picard = system($cmd_picard);
    print " (Exit code = $exit_code_picard)";
    if ($exit_code_picard){
      print RED, "\n\nERROR: Picard process reported exit code: $exit_code_picard - deleting fastq files", RESET;
      system($rm_cmd);
      next();
    }

    #Is specified by the user, trim the reads as desired
    if ($trim && $f && $l){
      my $r1_cmd_fastx = "cat $r1_fastq | fastx_trimmer -f $f -l $l -Q33 > $r1_fastq_temp";
      my $r1_cmd_mv = "mv -f $r1_fastq_temp $r1_fastq";
      my $r2_cmd_fastx = "cat $r2_fastq | fastx_trimmer -f $f -l $l -Q33 > $r2_fastq_temp";
      my $r2_cmd_mv = "mv -f $r2_fastq_temp $r2_fastq";
      print BLUE, "\n$r1_cmd_fastx\n\n", RESET;
      my $exit_code_fastx = system($r1_cmd_fastx);
      if ($exit_code_fastx){
        print RED, "\n\nERROR: Fastx process reported exit code: $exit_code_fastx - deleting fastq files", RESET;
        system($rm_cmd);
        next();
      }
      print BLUE, "\n$r1_cmd_mv\n\n", RESET;
      system($r1_cmd_mv);

      print BLUE, "\n$r2_cmd_fastx\n\n", RESET;
      $exit_code_fastx = system($r2_cmd_fastx);
      if ($exit_code_fastx){
        print RED, "\n\nERROR: Fastx process reported exit code: $exit_code_fastx - deleting fastq files", RESET;
        system($rm_cmd);
        next();
      }
      print BLUE, "\n$r2_cmd_mv\n\n", RESET;
      system($r2_cmd_mv);
    }

    #Before compressing, get and store the number of lines in this pair of files
    my $line_count = `wc -l $r1_fastq`; chomp($line_count);
    if ($line_count =~ /^(\d+)/){$line_count=$1;}else{print RED, "\n\nLine count on file failed: $r1_fastq\n\n", RESET; exit();}
    $bams{$name}{line_count} = $line_count;

    print BLUE, "\n\nCompressing FASTQ files", RESET;
    my $cmd_bzip = "$bzip_bin -f $r1_fastq $r2_fastq";
    if ($p){
      $cmd_bzip = "$bzip_bin -p$p -f $r1_fastq $r2_fastq";
    }
    print BLUE, "\n\t$cmd_bzip\n", RESET;
    my $exit_code_bzip = system($cmd_bzip);
    print " (Exit code = $exit_code_bzip)";
    if ($exit_code_bzip){
      print RED, "\n\nERROR: Bzip process reported exit code: $exit_code_bzip - deleting fastq_files", RESET;
      system($rm_cmd);
      next();
    }
  }
}


#Divide very large single lanes into pieces to improve performance downstream
#HiSeq lanes may have up to ~100-150 million reads in a single lane - divide into smaller pieces according to a target size
foreach my $name (sort keys %bams){
  my $flowcell = $bams{$name}{flowcell};
  my $lane = $bams{$name}{lane};
  my $data_dir = $bams{$name}{data_dir};
  my $r1_fastq_bz = $bams{$name}{r1_fastq_bz};
  my $r2_fastq_bz = $bams{$name}{r2_fastq_bz};

  print BLUE, "\n\nDividing large fastq lanes into pieces: $flowcell $lane\n\t$r1_fastq_bz\n\t$r2_fastq_bz", RESET;

  #Check if split files appear to already be present for this fastq lane
  my $test_path = "$data_dir"."s_"."$lane"."0_1_sequence.txt.bz2";
  if (-e $test_path){
    print YELLOW, "\n\n\tCompressed split fastq files already present, existing input files for this lane will be used\n\n", RESET;
    next();
  }

  #Get the line count of this fastq pair
  print BLUE, "\n\n\tCounting lines in input files", RESET;
 
  my $line_count1;
  if ($bams{$name}{line_count}){
    $line_count1 = $bams{$name}{line_count};
  }else{
    $line_count1 = `bzcat $r1_fastq_bz | wc -l`; chomp($line_count1);
    if ($line_count1 =~ /^(\d+)/){$line_count1=$1;}else{print RED, "\n\nLine count on file 1 failed\n\n", RESET; exit();}
  }
  print BLUE, "\n\t\tFound $line_count1 lines in file1 - assuming same number in file2", RESET;

  #Only allow a lane to be divided into a maximum of $max_bins pieces
  #If the read count is > $target_size*$max_bins, adjust the target size (by adding 1% of target size)
  my $target_bin_size = 15000000;
  #my $target_bin_size = 15000;
  my $max_bins = 10;
  my $read_count = ($line_count1/4); #Fastq files have one read per 4 lines

  #If there are a small number of reads, dont split into piece
  if ($read_count <= $target_bin_size){
    print YELLOW, "\n\n\tFound $read_count total reads - since this is less than the target bin size ($target_bin_size), skipping...", RESET;
    next();
  }else{
    print BLUE, "\n\n\tFound $read_count total reads - going to divide file into up to $max_bins of $target_bin_size reads", RESET;
  }
  my $increment = $target_bin_size * 0.01;
  if ($read_count > ($target_bin_size*$max_bins)){
    while ($read_count > ($target_bin_size*$max_bins)){
      $target_bin_size += $increment;
    }
  }

  #Determine the actual bin size that would evenly divide the total reads into bins without leaving a small amount left over

  #What number of bins gets closest to the target bin size
  #print YELLOW, "\n\nDetermine number of bins", RESET;
  #print YELLOW, "\n\ti\ttarget_bin_size\tbin_read_count\tabs_diff\tdiff\tremainder", RESET;
  my %diffs;
  for (my $i = 1; $i <= $max_bins; $i++){
    my $bin_read_count = sprintf("%d", ($read_count/$i));
    my $diff = sprintf("%d", ($bin_read_count - $target_bin_size)); #If -ve, the current bin is smaller than the target bin.  If +ve, the current bin is larger than the target bin.
    my $abs_diff = abs($diff);
    my $remainder = $read_count - ($i * $bin_read_count);
    #print YELLOW, "\n\t$i\t$target_bin_size\t$bin_read_count\t$diff\t$abs_diff\t$remainder", RESET;
    $diffs{$i}{diff} = $diff;
    $diffs{$i}{abs_diff} = $abs_diff;
    $diffs{$i}{bin_read_count} = $bin_read_count;
    $diffs{$i}{remainder} = $remainder;
  }
  my $o = 0;
  my $starting_bin_count;
  my $starting_bin_size;
  my $starting_bin_diff;
  foreach my $i (sort {$diffs{$a}->{abs_diff} <=> $diffs{$b}->{abs_diff}} keys %diffs){
    $o++;
    if ($o == 1){
      $starting_bin_count = $i;
      $starting_bin_size = $diffs{$i}{bin_read_count};
      $starting_bin_diff = $diffs{$i}{diff};
    }
  }

  #print YELLOW, "\n\nSelect: starting_bin_count = $starting_bin_count\tstarting_bin_size = $starting_bin_size\t$starting_bin_diff = $starting_bin_diff", RESET;

  #Using the starting bin size above, figure out the final bin size that will not leave a small remainder
  #Do this by incrementing $starting_bin_size until: ($read_count/$starting_bin_size) is less than $starting_bin_count
  while(($read_count/$starting_bin_size) > $starting_bin_count){
    my $x = $read_count/$starting_bin_size;
    #print YELLOW, "\n\t$read_count/$starting_bin_size = $x", RESET;
    $starting_bin_size++;
  }
  #Add one more for good measure
  $starting_bin_size += 1;
  #print YELLOW, "\n\nFinal select: starting_bin_size = $starting_bin_size", RESET;
  
  #Using this bin size, how many reads would actually be contained in each bin?
  #print YELLOW, "\n\nFinal bins:", RESET;
  my $read_pool = $read_count;
  my $i = 0;
  while($read_pool > $starting_bin_size){
    $i++;
    #print YELLOW "\n\t$i\t$starting_bin_size", RESET;
    $read_pool -= $starting_bin_size;
  }
  $i++;
  #print YELLOW "\n\t$i\t$read_pool", RESET;

  #Since fastq files have four lines per read, the actual divisor will be the $starting_bin_size * 4
  my $n_lines = $starting_bin_size * 4;
  #print YELLOW, "\n\nMAX LINES PER FILE: $n_lines", RESET;
  print BLUE, "\n\t\tDecided on $starting_bin_count bins of $starting_bin_size reads each", RESET;

  #Execute split command
  #my $split_cmd_r1 = "bzcat $r1_fastq_bz | head -n 400000 | split -l $n_lines -d -a 1 - $data_dir"."R1_";
  my $split_cmd_r1 = "bzcat $r1_fastq_bz | split -l $n_lines -d -a 1 - $data_dir"."R1_";
  print BLUE, "\n\n\tExecuting: $split_cmd_r1", RESET;
  Genome::Sys->shellcmd(cmd => $split_cmd_r1);

  #my $split_cmd_r2 = "bzcat $r2_fastq_bz | head -n 400000 | split -l $n_lines -d -a 1 - $data_dir"."R2_";
  my $split_cmd_r2 = "bzcat $r2_fastq_bz | split -l $n_lines -d -a 1 - $data_dir"."R2_";
  print BLUE, "\n\n\tExecuting: $split_cmd_r2", RESET;
  Genome::Sys->shellcmd(cmd => $split_cmd_r2);

  #Get the list of files created ... rename and compress them
  my $ls_cmd = "ls $data_dir"."R*";
  my @split_files = `$ls_cmd`;
  chomp(@split_files);
  foreach my $split_file (@split_files){
    if ($split_file =~ /(\S+\/)R(\d+)\_(\d+)$/){
      my $base_path = $1;
      my $read = $2;
      my $sublane = $3;
      my $old_file = "$split_file";
      my $new_file = "$base_path"."s_"."$lane$sublane"."_"."$read"."_sequence.txt";
      my $mv_cmd = "mv -f $old_file $new_file";
      print BLUE, "\n\n\t$mv_cmd", RESET;
      Genome::Sys->shellcmd(cmd => $mv_cmd);

      my $cmd_bzip = "$bzip_bin -f $new_file";
      if ($p){
        $cmd_bzip = "$bzip_bin -p$p -f $new_file";
      }
      print BLUE, "\n\t$cmd_bzip", RESET;
      system($cmd_bzip);

    }else{
      print RED, "\n\nCould not understand file path of split file: $split_file\n\n", RESET;
      exit();
    }
  }
}

#6.) Reports the total run time for this process
my $t1 = new Benchmark;
my $td = timediff($t1, $t0);
my $string = timestr($td);
my $seconds = '';
if ($string =~ /(\d+)\s+wallclock/){
  $seconds = $1;
}
print YELLOW, "\n\nImport, convert and compress run for $library_name took $seconds seconds", RESET;
print "\n\n";

exit();



