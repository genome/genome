package Genome::Model::Tools::HmpShotgun::FilterResults;

use strict;
use warnings;

use Genome;
use Workflow;
use IO::File;

class Genome::Model::Tools::HmpShotgun::FilterResults {
    is  => ['Command'],
    has => [
        working_directory => {
            is  => 'String',
            is_input => '1',
            doc => 'The working directory.',
        },
        reference1_top_hit_alignment_file => {
        	is  => 'String',
                is_input => '1',
                doc => 'Reads aligned to the first reference.',
        },
        reference2_top_hit_alignment_file => {
        	is  => 'String',
                is_input => '1',
                doc => 'Reads aligned to the second reference.',
        },
        paired_end1_concise_file => {
        	    is  => 'String',
                is_input => '1',
                doc => 'The conscise file of multiple hits.',
        },
        paired_end2_concise_file => {
        	    is  => 'String',
                is_input => '1',
                doc => 'The conscise file of multiple hits.',
        },
        taxonomy_file => {
        		is  => 'String',
                is_input => '1',
                doc => 'Taxonomy file.',
        },
        sam_header => {
        		is  => 'String',
                is_input => '1',
                doc => 'Ref seq sam header.',
        },

        filtered_alignment_file => {
        	    is  => 'String',
                is_output => '1',
                is_optional => '1',
                doc => 'The resulting filtered alignment.',
        },
        read_count_file => {
        	    is  => 'String',
                is_output => '1',
                is_optional => '1',
                doc => 'a report file.',
        },
        other_hits_file => {
        	    is  => 'String',
                is_output => '1',
                is_optional => '1',
                doc => 'a report file.',
        },
        lsf_resource => {
                is_param => 1,
                value => "-R 'select[mem>64000 && model!=Opteron250 && type==LINUX64] span[hosts=1] rusage[mem=64000]' -M 64000000",
                #value => "-R 'select[mem>30000 && model!=Opteron250 && type==LINUX64] span[hosts=1] rusage[mem=30000]' -M 30000000",
        },
        lsf_queue => {
                is_param => 1,
                value => 'bigmem'
                #value => 'long',
                }	
        ],
        has_param => [
                       #lsf_resource => {
           #     default_value => 'select[model!=Opteron250 && type==LINUX64] rusage[mem=12000] -M 24000000',
           #},
        ],
};

sub help_brief {
    'Run the reference coverage report.';
}

sub help_detail {
    return <<EOS
    Runs the reference coverage report.
EOS
}

sub execute {
    my $self = shift;

    $self->dump_status_messages(1);
    $self->status_message(">>>Running FilterResults at ".UR::Time->now);
    #my $model_id = $self->model_id;
    $self->status_message("Aligned Bam File for refseq1: ".$self->reference1_top_hit_alignment_file);
    $self->status_message("Aligned Bam File for refseq2: ".$self->reference2_top_hit_alignment_file);
    
    $self->status_message("Paired end 1 concise file: ".$self->paired_end1_concise_file);
    $self->status_message("Paired end 2 concise file: ".$self->paired_end2_concise_file);
    
    $self->status_message("Header: ".$self->sam_header);
    
    $self->status_message("Taxonomy file: ".$self->taxonomy_file);
    
    my $working_directory = $self->working_directory."/alignments_filtered/";
    my $report_directory = $self->working_directory."/reports/";
    
    unless (-e $report_directory) {
    	Genome::Sys->create_directory($report_directory);
    }

    my $other_hits_file = $report_directory."/other_hits.txt";
    my $read_count_file = $report_directory."/reads_per_contig.txt";
    my $filtered_alignment_file_no_header = $working_directory."/combined_no_header.sam";
    my $filtered_alignment_file_unsorted_bam = $working_directory."/combined_unsorted.bam";
    my $filtered_alignment_file = $working_directory."/combined.sam";
    my $filtered_alignment_file_bam = $working_directory."/combined.bam";

    my $genus_file = $report_directory."/genus.txt";
    my $phyla_file = $report_directory."/phyla.txt";
    
    my @expected_output_files = ($read_count_file,$genus_file,$phyla_file,$filtered_alignment_file_bam);
    
    my $rv_check = Genome::Sys->are_files_ok(input_files=>\@expected_output_files);
    if ($rv_check) {
        $self->filtered_alignment_file($filtered_alignment_file_bam);
        $self->read_count_file($read_count_file);
        $self->other_hits_file("some_other_hits_file_tbd");
 
    	$self->status_message("Expected output files exist.  Skipping processing.");
    	$self->status_message("<<<Completed FilterResults at ".UR::Time->now);
    	return 1;
    }

#/gscmnt/sata409/research/mmitreva/sabubuck/HMP_CLINICAL_SAMPLES_JAN_2010/SCRIPTS/merge_unique_top_hit_sam.no_phyla_counts.pl 
   
    #my $cmd = "/gsc/var/tmp/perl-5.10.0/bin/perl64 /gscmnt/sata409/research/mmitreva/sabubuck/HMP_CLINICAL_SAMPLES_JAN_2010/SCRIPTS/merge_unique_top_hit_sam_opt.pl -sam1 ".$self->reference1_top_hit_alignment_file." -sam2 ". $self->reference2_top_hit_alignment_file."  -taxonomy ".$self->taxonomy_file." -sam_concise1 ".$self->paired_end1_concise_file." -sam_concise2 ".$self->paired_end2_concise_file." -sam_combined_out $filtered_alignment_file_no_header  -read_count_output $read_count_file ";
    
    my $cmd = "/gsc/var/tmp/perl-5.10.0/bin/perl64 /gscmnt/sata409/research/mmitreva/sabubuck/HMP_CLINICAL_SAMPLES_JAN_2010/SCRIPTS/merge_unique_top_hit_sam.no_phyla_counts.pl -sam1 ".$self->reference1_top_hit_alignment_file." -sam2 ". $self->reference2_top_hit_alignment_file."  -taxonomy ".$self->taxonomy_file." -sam_concise1 ".$self->paired_end1_concise_file." -sam_concise2 ".$self->paired_end2_concise_file." -sam_combined_out $filtered_alignment_file_no_header  -read_count_output $read_count_file  -phyla_output $phyla_file  -genus_output $genus_file";
    
     														
    $self->status_message("FilterResults cmd: $cmd");

    $self->status_message("Running filter at ".UR::Time->now);
    my $rv_filter = Genome::Sys->shellcmd(cmd=>$cmd);
  
    if ( $rv_filter != 1) {
        $self->error_message("<<<Failed FilterResults on filter script.  Return value: $rv_filter");
    } 
    $self->status_message("Completed filter at ".UR::Time->now);
    
    my @input_files = ( $self->sam_header, $filtered_alignment_file_no_header );
    my $rv_cat = Genome::Sys->cat(input_files=>\@input_files,output_file=>$filtered_alignment_file); 
 
    if ( $rv_cat != 1) {
        $self->error_message("<<<Failed FilterResults on header cat.  Return value: $rv_cat");
    } 
    $self->status_message("Completed cat.");
 
    $self->status_message("Converting from sam to bam file: $filtered_alignment_file to $filtered_alignment_file_unsorted_bam");
    my $picard_path = "/gsc/scripts/lib/java/samtools/picard-tools-1.07/";
    my $cmd_convert = "java -Xmx2g -cp $picard_path/SamFormatConverter.jar net.sf.picard.sam.SamFormatConverter VALIDATION_STRINGENCY=SILENT I=$filtered_alignment_file O=$filtered_alignment_file_unsorted_bam";  
    #my $cmd_convert = "samtools view -bS $ > $merged_alignment_files_per_refseq_sam";
    my $rv_convert = Genome::Sys->shellcmd(cmd=>$cmd_convert);											 
            
    if ($rv_convert != 1) {
        $self->error_message("<<<Failed FilterResults on sam to bam conversion.  Return value: $rv_convert");
        return;
    }

    $self->status_message("Starting bam sort.");
    my $cmd_sort = Genome::Model::Tools::Sam::SortBam->create(file_name=>$filtered_alignment_file_unsorted_bam, output_file=>$filtered_alignment_file_bam);
    my $rv_sort = $cmd_sort->execute;
 
    if ($rv_sort != 1) {
        $self->error_message("<<<Failed FilterResults on sam to bam conversion.  Return value: $rv_convert");
        return;
    }

    $self->filtered_alignment_file($filtered_alignment_file_bam);
    $self->read_count_file($read_count_file);
    $self->other_hits_file("boo");
    
    Genome::Sys->mark_files_ok(input_files=>\@expected_output_files);
    
    
    #clean up
    $self->status_message("Removing intermediate files.");
    unlink($filtered_alignment_file);
    unlink($filtered_alignment_file_unsorted_bam);
    unlink($filtered_alignment_file_no_header);
    
    $self->status_message("<<<Completed FilterResults at ".UR::Time->now);
    
    return 1;
}

sub find_top_hit {
	my $self = shift;
	return 1;
}

1;
