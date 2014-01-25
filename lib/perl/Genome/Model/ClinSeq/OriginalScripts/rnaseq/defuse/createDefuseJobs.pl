#!/usr/bin/env genome-perl

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Benchmark;

my $usage=<<INFO;

./createDefuseJobs.pl 

e.g.

./createDefuseJobs.pl  --defuse_bin=/gscmnt/gc2142/techd/tools/defuse/defuse-0.4.2/scripts/defuse.pl  --defuse_config=/gscmnt/gc2142/techd/tools/defuse/defuse-0.4.2/scripts/config.txt  --working_dir=/gscmnt/gc2142/techd/analysis/defuse/results/Colon/  --library_name='H_KH-540079-41-11865-cDNA-2-lib1'

Use extra option --bam_file to specify a single bam file directly instead of looking up using the library name (need for external data)

If there are multiple libraries for a single sample, use:  --sample_name

e.g.

./createDefuseJobs.pl  --defuse_bin=/gscmnt/gc2142/techd/tools/defuse/defuse-0.4.2/scripts/defuse.pl  --defuse_config=/gscmnt/gc2142/techd/tools/defuse/defuse-0.4.2/scripts/config.txt  --working_dir=/gscmnt/gc2142/techd/analysis/defuse/results/Colon/  --sample_name='H_GP-13-0890-01A-01R-0451-09'

INFO

my $working_dir = '';
my $sample_name = '';
my $library_name = '';
my $bam_file = '';    #Start with a BAM file instead if the data was generated external to WASHU
my $defuse_bin = '';
my $defuse_config = '';

GetOptions ('working_dir=s'=>\$working_dir, 'sample_name=s'=>\$sample_name, 'library_name=s'=>\$library_name, 'bam_file=s'=>\$bam_file, 'defuse_bin=s'=>\$defuse_bin, 'defuse_config=s'=>\$defuse_config);

unless ($working_dir && ($library_name || $sample_name) && $defuse_bin && $defuse_config){
  print RED, "\n\nParameters missing\n\n", RESET;
  print GREEN, "$usage", RESET;
  exit();
}
if ($sample_name && $library_name){
  print RED, "\n\nSpecificy library_name OR sample_name, not both", RESET;
  exit();
}
unless ($working_dir =~ /\/$/){
  $working_dir .= "/";
}
unless (-e $working_dir && -d $working_dir){
  print RED, "\n\nWorking dir not valid: $working_dir\n\n", RESET;
  exit();
}
unless (-e $defuse_bin){
  print RED, "\n\nDefuse script not found: $defuse_bin\n\n", RESET;
  exit();
}
unless (-e $defuse_config){
  print RED, "\n\nDefuse config file not found: $defuse_config\n\n", RESET;
  exit();
}
if ($library_name){
  print "\n\nProcessing library: $library_name";
}elsif($sample_name){
  print "\n\nProcessing sample: $sample_name";
}



#Get bam file paths for this library name
#genome instrument-data list solexa --show id,flow_cell_id,lane,index_sequence,sample_name,library_name,clusters,read_length,bam_path,archive_path --filter library_name='H_KH-540079-41-11865-cDNA-2-lib1'
#genome instrument-data list solexa --show id,flow_cell_id,lane,index_sequence,sample_name,library_name,clusters,read_length,bam_path,archive_path --filter sample_name='H_GP-13-0890-01A-01R-0451-09'
my $bam_search;
if ($library_name){
  $bam_search = "genome instrument-data list solexa --show id,flow_cell_id,lane,index_sequence,sample_name,library_name,clusters,read_length,bam_path,archive_path --filter library_name='$library_name'";
}elsif($sample_name){
  $bam_search = "genome instrument-data list solexa --show id,flow_cell_id,lane,index_sequence,sample_name,library_name,clusters,read_length,bam_path,archive_path --filter sample_name='$sample_name'";
}
my @result = `$bam_search`;
my %bams;
my $lane_list = '';
my $group_name;
if ($library_name){
  $group_name = $library_name;
}elsif($sample_name){
  $group_name = $sample_name;
}

if ($bam_file){
  my $lane = 1;
  my $flowcell = "unknown";
  $bams{$group_name}{path} = $bam_file;
  $bams{$group_name}{file_type} = "bam";
  $bams{$group_name}{lane} = $lane;
  if ($bam_file =~ /^\/\w+\/(\w+)/){
    $flowcell = $1;
  }
  $bams{$group_name}{flowcell} = $flowcell;

  $lane_list .= "\n\t$flowcell (lane $lane)";

  unless (-e $bam_file){
    print RED, "\n\nBAM file: $bam_file does not appear to exist - aborting", RESET;
    exit();
  }
}else{
  foreach my $line (@result){
    chomp($line);
    my ($path, $base, $flowcell, $lane, $file_type);
    if ($line =~ /(\S+\.bam)/){
      $path = $1;
      $file_type = "bam";
      if ($line =~ /(\w+)\.bam/){
        $base = $1;
      }
    }elsif($line =~ /(\S+\.tar\.gz)/){
      $path = $1;
      $file_type = "fastq";
      if ($line =~ /(\w+)\.tar\.gz/){
        $base = $1;
      }
    }else{
      next();
    }
    if ($base =~ /gerald\_(\w+)\_(\d+)$/){
      $flowcell = $1;
      $lane = $2;
    }elsif ($base =~ /sequence\_(\w+)\_(\d+)$/){
      $flowcell = $1;
      $lane = $2;
    }

    unless ($path && $base && $flowcell && $lane){
      print RED, "\n\nCould not determine base name, path, flowcell and lane from bam path: $line\n\n", RESET;
      exit();
    }
    $bams{$base}{path} = $path;
    $bams{$base}{file_type} = $file_type;
    $bams{$base}{flowcell} = $flowcell;
    $bams{$base}{lane} = $lane;
    $lane_list .= "\n\t$flowcell (lane $lane)";
  }
}
my $lane_count = keys %bams;
print "\n\tFound $lane_count lane(s) of data: $lane_list";

#Get start time
my $t0 = new Benchmark;

#Write out a bash script that:
#1.) Creates a library analysis dir within the working dir
my $group_dir = "$working_dir"."$group_name";
unless (-e $group_dir && -d $group_dir){
  print "\n\nMaking directories";
  mkdir ($group_dir);
}

#2.) Creates a data dir within the library analysis dir
my $data_dir = "$group_dir/data/";
unless (-e $data_dir && -d $data_dir){
  mkdir ($data_dir);
}

#3.) Converts the BAM files to fastq files and writes these to the data dir, then compress them to save space
foreach my $name (sort keys %bams){
  my $filepath = $bams{$name}{path};
  my $flowcell = $bams{$name}{flowcell};
  my $lane = $bams{$name}{lane};
  my $file_type = $bams{$name}{file_type};
  my $r1_fastq_archive = "$data_dir"."s_"."$lane"."_1_sequence.txt";
  my $r2_fastq_archive = "$data_dir"."s_"."$lane"."_2_sequence.txt";
  my $r1_fastq = "$data_dir"."$flowcell". "_"."lane$lane"."_1.fastq";
  my $r2_fastq = "$data_dir"."$flowcell". "_"."lane$lane"."_2.fastq";
  my $r1_fastq_gz = "$r1_fastq".".gz";
  my $r2_fastq_gz = "$r2_fastq".".gz";
  $bams{$name}{r1_fastq} = $r1_fastq;
  $bams{$name}{r2_fastq} = $r2_fastq;
  $bams{$name}{r1_fastq_gz} = $r1_fastq_gz;
  $bams{$name}{r2_fastq_gz} = $r2_fastq_gz;
  my $rm_cmd = "rm -f $r1_fastq $r2_fastq $r1_fastq_gz $r2_fastq_gz";
  if (-e $r1_fastq_gz && -e $r2_fastq_gz){
    print "\n\n\tCompressed fastq files already present, existing input files for this lane will be used\n\n";
  }else{

    if ($file_type eq "bam"){
      print "\n\nConverting BAM to FASTQ using picard\n\n";
      #gmt picard sam-to-fastq --input=$filepath --fastq=$r1_fastq --fastq2=$r2_fastq --fragment-fastq=/dev/null --maximum-memory=12  --maximum-permgen-memory=256
      #OR
      #java -Xmx12g -XX:MaxPermSize=256m  -cp $ENV{GENOME_SW_LEGACY_JAVA}/samtools/picard-tools-1.22/sam-1.22.jar:$ENV{GENOME_SW_LEGACY_JAVA}/samtools/picard-tools-1.22/picard-1.22.jar:/gsc/scripts/opt/genome/snapshots/stable/genome-1035/lib/perl/Genome/Model/Tools/Picard/SamToFastq/GCSamToFastq.jar edu.wustl.genome.samtools.GCSamToFastq  INPUT='$filepath' FASTQ='$r1_fastq' SECOND_END_FASTQ='$r2_fastq' FRAGMENT_FASTQ='/dev/null'
      my $cmd_picard = "gmt picard sam-to-fastq --input=$filepath --fastq=$r1_fastq --fastq2=$r2_fastq --fragment-fastq=/dev/null --maximum-memory=12  --maximum-permgen-memory=256";
      my $exit_code_picard = system($cmd_picard);
      print " (Exit code = $exit_code_picard)";
      if ($exit_code_picard){
        print "\n\nERROR: Picard process reported exit code: $exit_code_picard - deleting fastq files";
        system($rm_cmd);
        next();
      }
      print "\n\nCompressing FASTQ files";
      my $cmd_gzip = "gzip -f --fast $r1_fastq $r2_fastq";
      my $exit_code_gzip = system($cmd_gzip);
      print " (Exit code = $exit_code_gzip)";
      if ($exit_code_gzip){
        print "\n\nERROR: Gzip process reported exit code: $exit_code_gzip - deleting fastq_files";
        system($rm_cmd);
        next();
      }
    }elsif($file_type eq "fastq"){
      print "\n\nCopying fastq archive from archive dir, unpacking, renaming and compressing files\n\n";      

      #Copy the tar archive to the data dir
      my $cp_cmd = "cp $filepath $data_dir";
      print "\n\n$cp_cmd";
      system($cp_cmd);

      #Unpack the tar containing R1 and R2 fastq files
      my $tar_cmd = "tar -zxvf $data_dir"."$name".".tar.gz  --directory=$data_dir";
      print "\n\n$tar_cmd\n";
      system($tar_cmd);

      #Rename files
      my $mv_cmd1 = "mv $r1_fastq_archive $r1_fastq";
      my $mv_cmd2 = "mv $r2_fastq_archive $r2_fastq";
      print "\n\n$mv_cmd1";
      system($mv_cmd1);
      print "\n\n$mv_cmd2";
      system($mv_cmd2);

      #Compress
      my $cmd_gzip = "gzip -f --fast $r1_fastq $r2_fastq";
      print "\n\n$cmd_gzip";
      my $exit_code_gzip = system($cmd_gzip);
      print " (Exit code = $exit_code_gzip)";
      if ($exit_code_gzip){
        print "\n\nERROR: Gzip process reported exit code: $exit_code_gzip - deleting fastq_files";
        system($rm_cmd);
        next();
      }

      #Delete the tar archive
      my $cmd_rm_tar = "rm -f $data_dir"."$name".".tar.gz";
      print "\n\n$cmd_rm_tar";
      system($cmd_rm_tar);

    }else{
      print RED, "\n\nFile type not determined\n\n", RESET;
      system($rm_cmd);
      exit();
    }
  }
}





#4.) Runs the defuse command
my $defuse_cmd = "$defuse_bin -c $defuse_config -d $data_dir -o $group_dir -n $group_name -s lsf -p 100";
print "\n\nRunning defuse command:\n$defuse_cmd\n\n";
my $exit_code_defuse = system($defuse_cmd);
print " (Exit code = $exit_code_defuse)";

#5.) If Defuse was successful - delete the fastq files, otherwise save them for a rerun
if ($exit_code_defuse){
  print "\n\nERROR: Defuse process reported exit code: $exit_code_defuse";
}else{
  foreach my $name (sort keys %bams){
    print "\n\nDeleting fastq files:";
    my $rm_cmd = "rm -f $bams{$name}{r1_fastq} $bams{$name}{r2_fastq} $bams{$name}{r1_fastq_gz} $bams{$name}{r2_fastq_gz}";
    print "\n\n$rm_cmd"; 
    system($rm_cmd);
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
print "\n\nDefuse run for $group_name took $seconds seconds";
print "\n\n";



exit();



