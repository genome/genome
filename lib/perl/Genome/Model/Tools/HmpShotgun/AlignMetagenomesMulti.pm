package Genome::Model::Tools::HmpShotgun::AlignMetagenomesMulti;

use strict;
use warnings;

use Genome;
use Workflow;
use IO::File;

use Data::Dumper;

class Genome::Model::Tools::HmpShotgun::AlignMetagenomesMulti {
    is  => ['Command'],
    has => [
        working_directory => {
            is  => 'String',
            is_input => '1',
            is_output => '1',
            doc => 'The working directory.',
        },
        reads_and_references => {
            is  => 'String',
            is_input => '1',
            doc => 'pipe delimited list of reads and ref seq files to align against',
        },
        generate_concise => {
            is  => 'String',
            is_input => '1',
            is_optional => '1',
            doc => 'flag to indicate whether or not to generate concise output in this run',
            default_value => '0',
        },
        aligned_file => {
            is  => 'String',
            is_output => '1',
            is_optional => '1',
            doc => 'The resulting alignment.',
        },
        unaligned_file => {
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
    $self->status_message(">>>Running AlignMetagenomes at ".UR::Time->now);
    #my $model_id = $self->model_id;
    $self->status_message("Reads and reference: ".$self->reads_and_references);

	#split between reads and ref on @ sign...
    my @list = split(/\@/,$self->reads_and_references);

    my $reads_file = $list[0];
    my $reference_sequence = $list[1];
    
    #split between | on reads if paired end
    my $reads_basename;
    my @reads_list = split(/\|/,$reads_file);
    if ( scalar(@reads_list) == 2 ) {
    	my $name1 = File::Basename::basename($reads_list[0]);
    	my $name2 = File::Basename::basename($reads_list[1]);
    	$reads_basename = $name1."_paired_".$name2;
    } elsif ( scalar(@reads_list) == 1) {
    	$reads_basename = File::Basename::basename($reads_list[0]);
    } else {
    	$self->status_message("Found an invalid number of input read files.  Need 1 or 2. Quitting.");
    	return;
    }
    
    $self->status_message("Reads: ".$reads_file);
    $self->status_message("Reference: ".$reference_sequence);
    
    #my $refseq_basename = File::Basename::basename($reference_sequence);
    my @refseq_path_dirs = split(/\//,$reference_sequence);
    my $refseq_basename = $refseq_path_dirs[-2]; 
    
    #my $refseq_dirname = File::Basename::dirname($reference_sequence);
    #$self->status_message("Refseq directory: ".$refseq_dirname);
 
 	#switch here on all whether or not to generate a concise alignment only
 	my $subdirectory = "alignments_top_hit";
 	my $alignment_options = " -t4 ";
 	my $alignment_file_name = "aligned.bam";
 	my $top_hits = 1;
 	if ( $self->generate_concise ) {
 		$subdirectory = "alignments_multiple_hits";
 		$alignment_file_name = "concise.txt";
 		$alignment_options = " -t4 -n5 ";
 		$top_hits = "10000";
 	}
 
    my $parent_directory = $self->working_directory."/$subdirectory/";
    my $working_directory = $self->working_directory."/$subdirectory/".$reads_basename."_aligned_against_".$refseq_basename;
    unless (-e $working_directory) {
    	Genome::Sys->create_directory("$working_directory");
    }
    
    #$self->aligned_file($alignment_file);
    #$self->status_message("<<<Completed AlignMetagenomes for testing at at ".UR::Time->now);
    #return 1;
    
    #expected output files
    #Move these to resolver methods in a build object or something similar
    my $alignment_file = $working_directory."/".$alignment_file_name;
    my $aligner_output_file = $working_directory."/aligner_output.txt";
    my $unaligned_reads_file = $working_directory."/unaligned.txt";

    $self->unaligned_file($unaligned_reads_file);
    
    #check to see if those files exist
    my @expected_output_files = ( $aligner_output_file, $alignment_file );
    my $rv_check = Genome::Sys->are_files_ok(input_files=>\@expected_output_files);
    
    if (defined($rv_check)) {
	    if ($rv_check == 1) {
	    	#shortcut this step, all the required files exist.  Quit.
	    	$self->status_message("Skipping this step.  Alignments exist for reads file $reads_basename against reference sequence $refseq_basename. If you would like to regenerate these files, remove them and rerun.");
                $self->aligned_file($alignment_file);
    	        $self->working_directory($parent_directory);
    	        $self->status_message("<<<Completed alignment at ".UR::Time->now);
                return 1;
	    } 
    } 
 

	my $aligner;

    if ($self->generate_concise) { 
    	$self->status_message("Aligning with Concise option at ".UR::Time->now);
    	$self->status_message("Reads file: ".$reads_file);
    	$aligner = Genome::Model::Tools::Bwa::AlignReadsMulti->create(dna_type=>'dna', 
    									align_options=>$alignment_options, 
    									ref_seq_file=>$reference_sequence,
    									files_to_align_path=>$reads_file,
    									aligner_output_file=>$aligner_output_file,
    	         						        unaligned_reads_file=>$unaligned_reads_file,
    									concise_file=>$alignment_file,
                                                                        temp_directory=>$working_directory,
                                                                        top_hits=>$top_hits,
            							);
    } else {				
    	$self->status_message("Aligning with standard options at ".UR::Time->now);
    	$aligner = Genome::Model::Tools::Bwa::AlignReads->create(dna_type=>'dna', 
    									align_options=>$alignment_options, 
    									ref_seq_file=>$reference_sequence,
    									files_to_align_path=>$reads_file,
    									aligner_output_file=>$aligner_output_file,
    	         						unaligned_reads_file=>$unaligned_reads_file,
    									alignment_file=>$alignment_file,
                                        temp_directory=>$working_directory,
                                        top_hits=>$top_hits,
            							);
    }
    															
    $self->status_message("Aligning at ".UR::Time->now);
    my $rv_aligner = $aligner->execute;

    if ($rv_aligner != 1) {
              $self->error_message("Aligner failed.  Return value: $rv_aligner");
              return;
    }
   
    if ($self->generate_concise) {
   		#set outputs for next step in workflow
    	$self->aligned_file($alignment_file);
    	$self->working_directory($parent_directory);
    	Genome::Sys->mark_files_ok(input_files=>\@expected_output_files);
    	$self->status_message("<<<Completed alignment at ".UR::Time->now);
    	return 1;
    }
    
	#if not concise, keep going...
    ##sort the alignment file

    my $sorted_file = $working_directory."/aligned_coord_sorted.bam";
    my $sorter = Genome::Model::Tools::Sam::SortBam->create(file_name=>$alignment_file,
							    name_sort=>0,
    							output_file=>$sorted_file);
   	
   my $rv_sort = $sorter->execute;
   if ($rv_sort != 1) {
   		$self->error_message("Sort failed.  Return value: $rv_sort");
    	return;
   }

	$self->status_message("Moving $sorted_file into $alignment_file.");
   	my $mv_cmd = "mv $sorted_file $alignment_file";
   	my $rv_mv = Genome::Sys->shellcmd(cmd=>$mv_cmd);
   	if ($rv_mv != 1) {
   		$self->error_message("Move of sorted file failed.  Return value: $rv_sort");
   	return;	
   	}
    	
    $self->aligned_file($alignment_file);
    $self->working_directory($parent_directory);
    Genome::Sys->mark_files_ok(input_files=>\@expected_output_files);
    
    $self->status_message("<<<Completed alignment at ".UR::Time->now);
    
    return 1;
}
1;
