#lumpy.pesr.pm 
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



sub pe_retrieve{

 #-----------------------------------------------------------------
 #     retrieve the pe bam info and save it to a string
 #-----------------------------------------------------------------
  
   print "\n I'm going to get the info for the pe bam file. \n";

   my $pe_bam =  Genome::Sys->create_temp_file_path(
        
        ); 
        
        return $pe_bam;
}

sub sr_retrieve{

#-----------------------------------------------------------------
#     retrieve the sr bam info and save it to a string
#-----------------------------------------------------------------  
   
    print "\n Now I'm going to get the info for the sr bam file. \n";

   my $sr_bam =  Genome::Sys->create_temp_file_path(

        );
      
        print "done\n\n";
        
        return $sr_bam;
} 

sub pe_alignment{
 #-----------------------------------------------------------------
 #                       pe_alignment work
 #-----------------------------------------------------------------

   my $orig_bam = shift;
   my $pe_bam = shift;

   print  "I'm going to write the pe samtools string. \n";

   my $pe_alignments ="samtools view -u -F 0x0002 $orig_bam |  samtools view -u -F 0x0100 - | samtools view -u -F 0x0004 - | samtools view -u -F 0x0008 - | samtools view -b -F 0x0400 - > $pe_bam"; 


       print "\n I wrote the pe samtools string. \n";

    my $pe_split  = Genome::Sys->shellcmd(
         cmd => $pe_alignments,
         allow_zero_size_output_files => 1,
        );
}

sub sr_alignment{
#-----------------------------------------------------------------
#                       sr_alignment
#-----------------------------------------------------------------

    my $orig_bam = shift;
    my $sr_bam = shift;

    print "\n I'm getting ready to write the sr_alignment samtools command! \n";
 
    my $sr_alignment = "samtools view -h $orig_bam | /gscuser/mfulton/lumpy-sv/scripts/extractSplitReads_BwaMem -i stdin | samtools view -Sb - > $sr_bam";
   
    print "\n I wrote it so now I'm going to run it using shellcmd. \n";

    my $sr_split = Genome::Sys->shellcmd(
       cmd => $sr_alignment,
       allow_zero_size_output_files => 1,
       );
}


sub mean_stdv_reader{
#-----------------------------------------------------------------
#    This will calculate the mean and standard 
#    deviation and then saves it to a text file.
#-----------------------------------------------------------------
  my $orig_bam = shift;  
  my $pe_histo = shift;

  print "\n I'm going to look up the location for the mean and stdv.\n";

  my $export_loc = "/gscuser/mfulton/Practice/mean_stdv.txt";
  
  my @mn_stdv = qq(samtools view $orig_bam | tail -n+100000 | /gscuser/mfulton/lumpy-sv/scripts/pairend_distro.py -r1 100 -X 4 -N 10000 -o $pe_histo);

  print "\n I'm going to run the command now.\n";

  require IPC::System::Simple;
  
  my $ms_output = IPC::System::Simple::capture(@mn_stdv);

  print "$ms_output";
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

#------------------------------------------------------
#              paramater reader
#------------------------------------------------------

sub read_params{
  
  my $self = shift;   

  my $all_params = $self->params;
   
  my @params_list = split('//',$all_params);
  
  my $lumpy_params = join(' >> ',@params_list);
  
  print "$lumpy_params";
  
  print "\n\n\n";
 
  return @params_list;
   
}


#------------------------------------------------------
#            lumpy  paramater reader
#------------------------------------------------------

sub lumpy_params_{
  
    my @params1 = shift;

    my $params_l = shift(@params1);

    my @params2 = split(',',$params_l);
 
    my $lumpy_params = join(' ',@params2);

    my @params3 = split(':',$lumpy_params);
 
    my $lumpy_params2 = join(' ',@params3);

    return $lumpy_params2;

}


#------------------------------------------------------
#                pe paramater reader
#-----------------------------------------------------

sub pe_param_reader{

     my @params1 = shift;

     my $params_l = shift(@params1);
 
     return $params_l;

}

#-------------------------------------------------------
#               sr parameter reader
#------------------------------------------------------

sub sr_param_reader{

     my @params1 = shift;
     
     my $params_l = shift(@params1);
      
     return $params_l;
}


sub _detect_variants {

    my $self = shift;

     print "\nI'm getting ready to recieve the file! \n"; 
        
   my $orig_bam = $self->aligned_reads_input;
   #my $orig_bam = "~/lumpy-sv/data/pe_sr.pos_sorted.bam";    
        print "\nI recieved the file! \n";

             
    
    my $pe_histo =  Genome::Sys->create_temp_file_path(
       );
 
    my $pe_loc = &pe_retrieve;

    my $sr_loc = &sr_retrieve;

    my @st_mn = &mean_stdv_reader($orig_bam,$pe_histo);
    
    print "@st_mn";
 
    my $mean = shift(@st_mn);

    my $std = pop(@st_mn);

    &pe_alignment($orig_bam,$pe_loc);
    
    &sr_alignment($orig_bam,$sr_loc);

    my @perm =$self->read_params;
   
    print "the array perm is @perm";
    print "\n\n\n";

    my $lump_text = &lumpy_params_(@perm);

    print "the lumpy paramaters are $lump_text \n\n";

    shift(@perm);

    my $pe_text = &pe_param_reader(@perm);
   
    print "the paired end paramaters are $pe_text \n\n";

    shift(@perm);
    
    my $sr_text = &sr_param_reader(@perm);
  
    print "the split read paramaters are $sr_text \n\n";    


#-----------------------------------------------------------------
#                      running lumpy
#-----------------------------------------------------------------

#path to lumpy program

    my $executable_path = "/gscuser/mfulton/lumpy-sv/bin/lumpy";

# out_put to designated spot

    my $output_files = $self->_sv_staging_output;


    my $pe_info = "bam_file:$pe_loc,histo_file:$pe_histo,mean:$mean,stdev:$std,read_length:150,$pe_text";

#"bam_file:/gscuser/mfulton/lumpy-sv/data/sr.pos_sorted.bam

    my $sr_info = "bam_file:$sr_loc,$sr_text";


    my $cmd = "$executable_path $lump_text -pe $pe_info -sr $sr_info > $output_files";

 print "\n I'm going to execute the shell command. \n";

    my $return = Genome::Sys->shellcmd(
        cmd => $cmd,
        output_files => [$self->_sv_staging_output],
        allow_zero_size_output_files => 1,
       );
    }

