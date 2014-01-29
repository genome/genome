package Genome::Model::Tools::HmpShotgun::AlignMetagenomes;

use strict;
use warnings;

use Genome;
use Workflow;
use IO::File;

use Data::Dumper;

class Genome::Model::Tools::HmpShotgun::AlignMetagenomes {
    is  => ['Command'],
    has => [
        working_directory => {
            is  => 'String',
            is_input => '1',
            doc => 'The working directory.',
        },
        reference_name => {
        	is  => 'String',
            is_input => '1',
            doc => 'The reference sequence.',
        },
        instrument_data_id => {
        	is  => 'String',
            is_input => '1',
            doc => 'The reads to align.',
        },
        aligned_file => {
        	is  => 'String',
            is_output => '1',
            is_optional => '1',
            doc => 'The resulting alignment.',
        },
        
	],
    has_param => [
           lsf_resource => {
           default_value => 'select[model!=Opteron250 && type==LINUX64] rusage[mem=4000]',
           },
    ],
};

sub help_brief {
    'Align reads against a given metagenomic reference.';
}

sub help_detail {
    return <<EOS
    Align reads against a given metagenomic reference.
EOS
}

sub execute {
    my $self = shift;

    $self->dump_status_messages(1);
    $self->debug_message(">>>Running AlignMetagenomes at ".UR::Context->current->now);
    #my $model_id = $self->model_id;
    $self->debug_message("Ref seq: ".$self->reference_name);
    #$self->debug_message("Reads: ".$self->reads_file);
    
  #  my $align_basename = File::Basename::basename($self->reads_file);
    
#    my $working_directory = $self->working_directory."/alignments/".$align_basename."/";
#    unless (-e $working_directory) {
#    	Genome::Sys->create_directory("$working_directory");
#    }
    
    $self->debug_message("Working directory: ".$self->working_directory);
    
    #my $alignment_file = $working_directory."/alignment_file.bam";
    #$self->aligned_file($alignment_file);
    #$self->debug_message("<<<Completed AlignMetagenomes for testing at at ".UR::Context->current->now);
    #return 1;
    
    #expected output files
    #Move these to resolver methods in a build object or something similar
#    my $aligner_output_file = $working_directory."/aligner_output.txt";
#    my $unaligned_reads_file = $working_directory."/unaligned.txt";
#    my $alignment_file = $working_directory."/alignment_file.bam";
# 
#    $self->aligned_file($alignment_file);
#    
#    #check to see if those files exist
#    my @expected_output_files = ( $aligner_output_file, $unaligned_reads_file, $alignment_file );
#    my $rv_check = Genome::Sys->are_files_ok(input_files=>\@expected_output_files);
#    
#    if (defined($rv_check)) {
#	    if ($rv_check == 1) {
#	    	#shortcut this step, all the required files exist.  Quit.
#	    	$self->debug_message("Skipping this step.  If you would like to regenerate these files, remove them and rerun.");
#	   	    $self->debug_message("<<<Completed alignment at ".UR::Context->current->now);
#	   	    return 1;
#	    }
#	}
     
     
     
     
     my %alignment_params = (
            instrument_data_id => $self->instrument_data_id,
            reference_name     => $self->reference_name,
            aligner_name       => "bwa",
            aligner_version    => "0.5.4",
            aligner_params     => "-t4",
            force_fragment     => 0,
            samtools_version   => "r453",
            picard_version     =>  "r104",
        );
     
      my $alignment = Genome::InstrumentData::Alignment->create(%alignment_params);
   
       if (!$alignment || $@) {
            $self->error_message($@);
            $self->error_message('Failed to create an alignment object with params: '. Data::Dumper::Dumper(%alignment_params) );
       }
           
    
       $alignment->find_or_generate_alignment_data;
         	
       $self->debug_message("\n************ Alignment object: ".Dumper($alignment) );
         	
       $self->debug_message($alignment->output_dir."/all_sequences.bam");  	
       $self->aligned_file($alignment->output_dir."/all_sequences.bam");
     
     
#    my $aligner = Genome::Model::Tools::Bwa::AlignReads->create(dna_type=>'dna', 
#    															align_options=>' -t 4 ', 
#    															ref_seq_file=>$self->reference_sequence_file,
#    															files_to_align_path=>$self->reads_file,
#    															aligner_output_file=>$aligner_output_file,
#    															unaligned_reads_file=>$unaligned_reads_file,
#    															alignment_file=>$alignment_file,
#    															);
#    															
#    $self->debug_message("Aligning at ".UR::Context->current->now);
#    my $rv_aligner = $aligner->execute;
#    
#   
#    if ($rv_aligner != 1) {
#    	$self->error_message("Aligner failed.  Return value: $rv_aligner");
#    	return;
#    }
#           	
#    #sort the alignment file
#    my $sorted_file = $working_directory."/alignment_file.sorted.bam";
#    my $sorter = Genome::Model::Tools::Sam::SortBam->create(file_name=>$alignment_file,
#    														name_sort=>0,
#    														output_file=>$sorted_file);
#    														
#   	my $rv_sort = $sorter->execute;
#   	if ($rv_sort != 1) {
#   		$self->error_message("Sort failed.  Return value: $rv_sort");
#    	return;
#   	}
#   	
#   	my $mv_cmd = "mv $sorted_file $alignment_file";
#   	my $rv_mv = Genome::Sys->shellcmd(cmd=>$mv_cmd);
#   	if ($rv_mv != 1) {
#   		$self->error_message("Move of sorted file failed.  Return value: $rv_sort");
#    	return;	
#   	}
#    	
#    Genome::Sys->mark_files_ok(input_files=>\@expected_output_files);
#    
    $self->debug_message("<<<Completed alignment at ".UR::Context->current->now);
    
    return 1;
}
1;
