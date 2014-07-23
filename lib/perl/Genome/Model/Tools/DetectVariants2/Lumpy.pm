#lumpy.pm 
#this program is ment to run the bam files specified.  it will run both a paired end and a single read file before depositing the data into a histo file in the lumpy-sv_results dir

#/gscuser/mfulton/lumpy-sv/bin/lumpy -mw 4 -tt 0.0 -pe bam_file:/gscuser/mfulton/lumpy-sv/data/pe.pos_sorted.bam,histo_file:/gscuser/mfulton/lumpy-sv_results/pe.pos_sorted.new.histo,mean:499.8528,stdev:49.5842750493,read_length:150,min_non_overlap:150,discordant_z:4,back_distance:20,weight:1,id:1,min_mapping_threshold:20 -sr bam_file:/gscuser/mfulton/lumpy-sv/data/sr.pos_sorted.bam,back_distance:20,weight:1,id:2,min_mapping_threshold:20 > /gscuser/mfulton/lumpy-sv_results/results_pesr.bedpe


package Genome::Model::Tools::DetectVariants2::Lumpy;

use warnings;
use strict;
 
use Genome;
use File::Basename;
 
 my @FULL_CHR_LIST = (1..22, 'X', 'Y', 'MT');
 
 class Genome::Model::Tools::DetectVariants2::Lumpy{
     is => 'Genome::Model::Tools::DetectVariants2::Detector',
};


sub bam_split{
    my $orig_bam = shift;
    my $split_loc = Genome::Sys->create_temp_directory(
);
    my $command = "bamtools split -in $orig_bam -stub $split_loc/test1 -tag RG";
    print "I'm going to run $command ";
    my $return = Genome::Sys->shellcmd(
       cmd => $command,
); 
    print "\n\n the split bam files are in: $split_loc \n\n";   
    my @file_list = glob("$split_loc/*");
    print Data::Dumper::Dumper(\@file_list); 
 return @file_list;
}


sub pe_alignment{
   my $orig_bam = shift;   
   my $pe_bam = Genome::Sys->create_temp_file_path(
   ); 
   my $pe_alignments ="samtools view -u -F 0x0002 $orig_bam |  samtools view -u -F 0x0100 - | samtools view -u -F 0x0004 - | samtools view -u -F 0x0008 - | samtools view -b -F 0x0400 - > $pe_bam"; 
   print "200 - I'm going to run $pe_alignments ";
   my $pe_split  = Genome::Sys->shellcmd(
      cmd => $pe_alignments,
      allow_zero_size_output_files => 1,
      );
   return $pe_bam;
}

sub sr_alignment{
    my $orig_bam = shift;
    my $sr_bam = Genome::Sys->create_temp_file_path(
    );
   my $sr_alignment = "samtools view -h $orig_bam | /gscuser/mfulton/lumpy-sv/scripts/extractSplitReads_BwaMem -i stdin | samtools view -Sb - > $sr_bam";
   print "300 - I'm going to run $sr_alignment";
    my $sr_split = Genome::Sys->shellcmd(
       cmd => $sr_alignment,
       allow_zero_size_output_files => 1,
       );
    return $sr_bam;
}

sub pe_cmd_arrangement{
#-pe bam_file:disc-lib1.bam,histo_file:lib1.histo,mean:563.83,stdev:17.839,read_length:100,min_non_overlap:100,discordant_z:4,back_distance:20,weight:1,id:10,min_mapping_threshold:10  
    
    my $new_bam = shift;
    my $all_perm = shift;
    my $self = shift;

    my @perm = split('//',$all_perm);

    my $pe_histo = Genome::Sys->create_temp_file_path();  
   
    print "\n\nThe values in @ perm are: @perm ,\n\n The values in new_bam are: $new_bam, \n\n the histo file is located in: $pe_histo\n\n";
 
    print "I'm going to calculate the mean and stdv.\n";
    my @st_mn = &mean_stdv_reader($new_bam,$pe_histo);
    
    print "\nI just recieved @st_mn\n";
    

    my $mean = $st_mn[0];
    my $std = $st_mn[1];    
    my $pe_loc = &pe_alignment($new_bam);
    my $pe_text = $perm[1];
    my $pe_cmd = "-pe bam_file:$pe_loc,histo_file:$pe_histo,mean:$mean,stdev:$std,read_length:150,$pe_text";   
    
    print "I'm about to return pe_cmd whic contains $pe_cmd \n";

   return $pe_cmd;
}


sub sr_cmd_arrangement{
#This command should return a string that looks something like this:
#-sr bam_file:split-lib2.bam,back_distance:20,min_mapping_threshold:10,weight:1,id:40,min_clip:20
    
    my $new_bam = shift;
    my $all_perm = shift;
    my $self = shift;
    my $sr_loc = &sr_alignment($new_bam);
    my @perm = split('//',$all_perm);

    my $sr_text = $perm[2];

    my $sr_cmd = " -sr bam_file:$sr_loc,$sr_text";    

}
sub mean_stdv_reader{
  my $new_bam = shift;  
  my $pe_histo = shift;
  my $export_loc = "/gscuser/mfulton/Practice/mean_stdv.txt";
  my @mn_stdv = qq(samtools view $new_bam | tail -n+100000 | /gscuser/mfulton/lumpy-sv/scripts/pairend_distro.py -r1 100 -X 4 -N 10000 -o $pe_histo);
  
  print "400 - the mean and stdv command reads: @mn_stdv \n";
  
  require IPC::System::Simple;
  my $ms_output = IPC::System::Simple::capture(@mn_stdv);
  
  print "\n\n\n500 - $ms_output";
  
  foreach ($ms_output){  
  my @values = split;
  my $stdv1 = pop(@values);
  my $mean1 = pop(@values);
  my @mean2 = split(':',$mean1);  
  my $mean = pop(@mean2);
  my @stdv2 = split(':',$stdv1);
  my $stdv = pop(@stdv2);
  my @stdv_mean = ($mean,$stdv);

  return @stdv_mean;
 }
}


sub _open_params {

    my $self = shift;

    my $all_perm = $self->params;
    my @perm = split('//',$all_perm);

    my $lump_text = $perm[0];
    print "$lump_text";
    my @perm = split(',',$lump_text);
    $lump_text = join(' ',@perm);
    
    my @per = split(':',$lump_text);
    $lump_text = join(' ',@per);    

    print "$lump_text";
    my $executable_path = "~/lumpy-sv/bin/lumpy";
    my $output_files = $self->_sv_staging_output;

    my @sur_cmd = ("$executable_path $lump_text "," > $output_files");
    
    print "\nI'm returning the command for the beginning and end of the command line.  I'm returning: @sur_cmd \n";

    return @sur_cmd;
}

sub _detect_variants{
    my $self = shift;
 
    print "\nI'm going to detect the parameters.\n";  

    my $all_perm = $self->params;
    
    print "\n\n$all_perm\n\n";

    my @perm = split('//',$all_perm);
 
    my @cmd = $self->_open_params();
   
    foreach my $cms (@cmd){
    print "\n$cms";
    }

    my $orig_bam = $self->aligned_reads_input;
    
    print "\nI'm going to split the original bam file.\n";
    
    my @new_bam = &bam_split($orig_bam);

    print "the new bam file locations are in: @new_bam \n\n";
 
    my @pe_cmds;
    my @sr_cmds;
    
    foreach my $cur_bam (@new_bam){
        print "\nI'm going to detect the paired end commands for each file.\n";

        my $pe_cmd = &pe_cmd_arrangement($cur_bam,$all_perm);

        print "\nI'm going to detect the split read commands for each file.\n";

        my $sr_cmd = &sr_cmd_arrangement($cur_bam,$all_perm);

        print "\n\n@pe_cmds\n\n"; 
        print "\n\n@pe_cmds\n\n"; 
        
        push (@pe_cmds,"$pe_cmd"); 
        push (@sr_cmds,"$sr_cmd");
     }

    print "\nThe PE array\n"; 
    print "\n\n@pe_cmds\n\n"; 

    print "\nThe SR array\n"; 
    print "\n\n@sr_cmds\n\n"; 
    
    my $pe_cmd = join(",@pe_cmds");   
    my $sr_cmd = join(",@sr_cmds");   
 
    splice @cmd, 1, 0, "@pe_cmds", "@sr_cmds";  

    my $cmmd = "@cmd";

    print Data::Dumper::Dumper(\@cmd);   
 

    foreach my $cms (@cmd){
    print "$cms";
    }

    my $run = Genome::Sys->shellcmd(
        cmd => $cmmd,
        output_files => [$self->_sv_staging_output],
        allow_zero_size_output_files => 1,
       );

}      
