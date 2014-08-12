#lumpy.pm 
#this program is ment to run the bam files specified.  it will run both a paired end and a single read file before depositing the data into a histo file in the lumpy-sv_results dir

package Genome::Model::Tools::DetectVariants2::Lumpy;
use warnings;
use strict;
 
use Genome;
use File::Basename;
use IPC::System::Simple;
 
 my @FULL_CHR_LIST = (1..22, 'X', 'Y', 'MT');
 
class Genome::Model::Tools::DetectVariants2::Lumpy{
     is => 'Genome::Model::Tools::DetectVariants2::Detector',
};

sub _detect_variants{
    my $self = shift;
    my $orig_bam = $self->aligned_reads_input;
    my @new_bam = bam_split($orig_bam);
    my @pe_cmds;
    my @sr_cmds;    
    foreach my $cur_bam (@new_bam){
        if ($self->pe_param){
            my $pe_cmd = $self->pe_cmd_arrangement($cur_bam);
            push (@pe_cmds,"$pe_cmd"); 
        }
        if ($self->sr_param){
            my $sr_cmd = $self->sr_cmd_arrangement($cur_bam);
            push (@sr_cmds,"$sr_cmd");
        }
     }
    if (defined $self->control_aligned_reads_input){
       my $bam2 = $self->control_aligned_reads_input;
       my @new_bam2 = &bam_split($bam2);
       foreach my $cur_bam (@new_bam2){
            if ($self->pe_param){
                my $pe_cmd = $self->pe_cmd_arrangement($cur_bam);
                push (@pe_cmds,"$pe_cmd"); 
            }
            if ($self->sr_param){
                my $sr_cmd = $self->sr_cmd_arrangement($cur_bam);
                push (@sr_cmds,"$sr_cmd");
            }
       }
    }
    my $pe_cmd = join(",@pe_cmds");   
    my $sr_cmd = join(",@sr_cmds");   
    
    my @cmd = $self->_open_params();
    splice @cmd, 1, 0, "@pe_cmds", "@sr_cmds";  
    my $cmmd = "@cmd";
    
    my $run = Genome::Sys->shellcmd(
        cmd => $cmmd,
        output_files => [$self->_sv_staging_output],
        allow_zero_size_output_files => 1,
       );
}

sub bam_split{
    my $orig_bam = shift;
    my $split_loc = Genome::Sys->create_temp_directory();
    my $command = "bamtools split -in $orig_bam -stub $split_loc/test1 -tag RG";
    my $return = Genome::Sys->shellcmd(
       cmd => $command,
    ); 
    my @file_list = glob("$split_loc/*");
    return @file_list;
}

sub pe_alignment{
    my $orig_bam = shift;   
    my $pe_bam = Genome::Sys->create_temp_file_path(); 
    my $pe_alignments ="samtools view -b -F 1294 $orig_bam -o $pe_bam"; 
    my $pe_split  = Genome::Sys->shellcmd(
      cmd => $pe_alignments,
      allow_zero_size_output_files => 1,
      );
    return $pe_bam;
}

sub sr_alignment{
    my $self = shift;
    my $orig_bam = shift;
    my $sr_bam = Genome::Sys->create_temp_file_path();
    my $extraction_script = $self->lumpy_script_for_extract_split_reads_bwamem();
    my $sr_alignment = "samtools view -h $orig_bam | $extraction_script -i stdin | java -Xmx8g -XX:MaxPermSize=256m -cp /gsc/scripts/lib/java/samtools/picard-tools-1.82/SamFormatConverter.jar net.sf.picard.sam.SamFormatConverter I=/dev/stdin O=$sr_bam";
    my $sr_split = Genome::Sys->shellcmd(
       cmd => $sr_alignment,
       allow_zero_size_output_files => 1,
       );
    return $sr_bam;
}

sub pe_cmd_arrangement{
  
    my $self = shift;
    my $current_split_bam = shift;
    
    my $pe_loc = &pe_alignment($current_split_bam);
    my %st_mn = $self->mean_stdv_reader($current_split_bam);
    my $mean = $st_mn{mean};
    my $std = $st_mn{stdv};    
    my $pe_histo = $st_mn{histo};  
    my $pe_text = $self->pe_param;
    
    my $pe_cmd = " -pe bam_file:$pe_loc,histo_file:$pe_histo,mean:$mean,stdev:$std,read_length:150,$pe_text";   
    return $pe_cmd;
}

sub sr_cmd_arrangement{
    my $self = shift;  
    my $new_bam = shift;
    my $sr_loc = $self->sr_alignment($new_bam);
    my $sr_cmd =$self->sr_arrange($sr_loc); 
    return $sr_cmd;
}

sub sr_arrange{
    my $self = shift;
    my $sr_loc = shift;
    my $sr_text = $self->sr_param;
    my $sr_cmd = " -sr bam_file:$sr_loc,$sr_text"; 
    return $sr_cmd;
}

sub mean_stdv_reader{
    my $self = shift;
    my $new_bam = shift;  
    my $pe_histo = Genome::Sys->create_temp_file_path();  
    my $temp_dir = Genome::Sys->create_temp_file_path();
    
    my $export_loc = "$temp_dir/mean_stdv.txt";
    my $pe_extraction_script = $self->lumpy_script_for_pairend_distro();
    my @mn_stdv = qq(samtools view $new_bam | tail -n+100 | $pe_extraction_script -r1 100 -X 4 -N 10000 -o $pe_histo);
    my $ms_output = IPC::System::Simple::capture(@mn_stdv);
    
    if ($ms_output =~ m/mean:([-+]?[0-9]*\.?[0-9]*)\s+stdev:([-+]?[0-9]*\.?[0-9]*)/){
        my $mean = $1;
        my $stdv = $2;
        my %stdv_mean = (
            mean => $mean,
            stdv => $stdv,
            histo=> $pe_histo,
        );
        return %stdv_mean;
    }
    else {
        die"ERROR couldn't find mean and stdev";
    }
 }
sub sr_param{
    my $self = shift;
    my %params =$self->params_hash();
    return $params{'sr'};
}

sub pe_param{
    my $self = shift;
    my %params =$self->params_hash();
    return $params{'pe'};
}

sub lumpy_param{
    my $self = shift;
    my %params =$self->params_hash();
    return $params{'lp'};
}

sub params_hash{
    my $self = shift;
    my $unparsed_params = $self->params;
    my @params = split('//',$unparsed_params);
    my %parameters;
    foreach my $place (@params){
        if ($place =~ m/^\-([a-z]{2}),(.*)$/){
            $parameters{$1} = $2;
        }
        else {
            die sprintf("You specified the parameters incorrectly. Unparsed parametere were: (%s)  The malformed parameters were: (%s)",
                $unparsed_params, $place);
        }
    }
    return %parameters;
}

sub _open_params {
    my $self = shift;
    my $lump_text =$self->lumpy_param;
    $lump_text =~ s/,/ /g;
    $lump_text =~ s/:/ /g;    
    print "$lump_text";
    my $executable_path = "~/lumpy-sv/bin/lumpy";
    my $output_files = $self->_sv_staging_output;
    my @sur_cmd = ("$executable_path $lump_text "," > $output_files");
    return @sur_cmd;
}

sub lumpy_directory{
    my $self = shift;
    my $version = $self->version();

    return _lumpy_directory($version); 
}

sub _lumpy_directory {
    my $version = shift;  
    return File::Spec->catdir(File::Spec->rootdir,"usr","lib","lumpy"."$version"); 
}

sub lumpy_command{
    my $self = shift;
    return File::Spec->catfile($self->lumpy_directory(),"bin","lumpy");
}

sub lumpy_scripts_directory{
    my $self = shift;
    return File::Spec->catfile($self->lumpy_directory(),'scripts');
}

sub lumpy_script_for {
    my $self = shift;
    my $script_name = shift;

    die "no script name given" if not $script_name;
    my $script_location = File::Spec->catfile($self->lumpy_scripts_directory(), "$script_name");

    die "script does not exist $script_location" if not -e $script_location; 
    return $script_location;
}

sub lumpy_script_for_extract_split_reads_bwamem{
    my $self = shift;
    return $self->lumpy_script_for("extractSplitReads_BwaMem");
}

sub lumpy_script_for_pairend_distro{
    my $self = shift;
    return $self->lumpy_script_for("pairend_distro.py");
}

sub has_version {
    my $self = shift;
    my $version = shift; 
    if (-d _lumpy_directory($version)){
        return 1;
    }
    else {
        return 0;
    }
}







