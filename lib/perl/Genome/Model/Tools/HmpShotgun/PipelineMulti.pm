package Genome::Model::Tools::HmpShotgun::PipelineMulti;

use strict;
use warnings;

use Genome;
use Workflow;

#run with
#gmt hmp-shotgun pipeline-multi --reads-files='/gscmnt/sata409/research/mmitreva/jpeck/testdata/t1.fna|/gscmnt/sata409/research/mmitreva/jpeck/testdata/t2.fna' --reference-sequences=/gscmnt/sata409/research/mmitreva/jpeck/refseq_metagenome1/all_sequences.fa,/gscmnt/sata409/research/mmitreva/jpeck/refseq_metagenome2/all_sequences.fa --working-directory=/gscmnt/sata409/research/mmitreva/jpeck/multi_build0 --workflow-log-directory=/gscmnt/sata409/research/mmitreva/jpeck/multi_build0/workflow_logs --regions-file=/gscmnt/sata409/research/mmitreva/jpeck/refseq_metagenome1/combined_ref_cov_regions.txt --reads-and-references=1 --generate-concise=1

class Genome::Model::Tools::HmpShotgun::PipelineMulti {
    is => ['Workflow::Operation::Command'],
    workflow => sub { Workflow::Operation->create_from_xml(\*DATA); },
    has => [
        workflow_log_directory => {
                    is => 'String',
                    doc => 'The directory where the workflow logs (LSF output) should be dumped.' ,
        },
        cleanup => { 
                    is => 'Boolean',
                    is_optional => '1',
                    default_value => '0',
                    doc => 'A clean up flag.  Will remove intermediate files if set. Default = 0, no cleanup.',
        },
    ]
};

sub help_synopsis{
    my $self = shift;
    return "TBD";
}

sub pre_execute {
    my $self = shift;

    #make required directories if they don't exist
    my $working_dir = $self->working_directory;

    $self->_operation->log_dir($self->workflow_log_directory);

    $self->status_message("Launching HMP Metagenomic Pipeline");
    $self->status_message("Using working directory:".$working_dir);
    $self->status_message("Workflow log directory:".$self->workflow_log_directory);
    $self->status_message("Delete intermediate files on completion: ".$self->cleanup);
    
    $self->status_message("Reads file string: ".$self->reads_files);
    my @reads_files = split(/,/ , $self->reads_files);
    my $list_string = join("\n",@reads_files);
    $self->status_message("Reads files: \n".$list_string); 
	
    $self->status_message("Reference sequence files string: ".$self->reference_sequences);
    my @reference_sequences = split(/,/ , $self->reference_sequences);
    my $list_refseqs = join("\n",@reference_sequences);
    $self->status_message("Ref seq files: \n".$list_refseqs); 
    
    #for paired end, top hit alignments
    my @reads_and_references;
    #for fragments, mulit hit alignments
    my @reads_and_references_frags;
    for my $read_item (@reads_files) {
        for my $refseq_item (@reference_sequences) {
            push (@reads_and_references,$read_item."@".$refseq_item);
            #break out the frags for multi alignment
            my @frag_reads = split(/\|/,$read_item);
            for my $frag_read_item (@frag_reads) {
            	push (@reads_and_references_frags,$frag_read_item."@".$refseq_item);
            }
        }
    } 
    
   

    $self->status_message("Paired end Reads and References");
    $self->status_message(join("\n",@reads_and_references) );
    
    $self->status_message("Fragment Reads and References");
    $self->status_message(join("\n",@reads_and_references_frags) );
    
    
    
    #$self->reads_file(\@reads_and_references);
    $self->reads_and_references(\@reads_and_references);
    $self->fragment_reads_and_references(\@reads_and_references_frags);
	
    $self->status_message("Creating required directories.");
    Genome::Sys->create_directory("$working_dir");
    Genome::Sys->create_directory("$working_dir/alignments_filtered");
    Genome::Sys->create_directory("$working_dir/alignments_top_hit");
    Genome::Sys->create_directory("$working_dir/alignments_multiple_hits");
    Genome::Sys->create_directory("$working_dir/logs");
    Genome::Sys->create_directory("$working_dir/workflow_logs");
    Genome::Sys->create_directory("$working_dir/tmp");
    Genome::Sys->create_directory("$working_dir/reports");

    $self->status_message("Pre-execute of Pipeline complete.");

    return 1;
}

sub post_execute {
    my $self = shift;
    my $working_dir = $self->working_directory;

    my $cleanup = $self->cleanup;

    if ($cleanup) {
        $self->status_message("Cleaning up intermediate files.");
     
    } else {
        $self->status_message("Leaving intermediate files behind.");
    }

	$self->status_message("Done.");
    return 1;
}

1;
__DATA__
<?xml version='1.0' standalone='yes'?>

<workflow name="HMP Metagenomic Pipeline">

  <link fromOperation="input connector" fromProperty="working_directory"            toOperation="AlignTopHit" toProperty="working_directory" /> 
  <link fromOperation="input connector" fromProperty="reads_and_references"	        toOperation="AlignTopHit" toProperty="reads_and_references" />
  <link fromOperation="AlignTopHit"     fromProperty="aligned_file"                 toOperation="MergeAlignments" toProperty="alignment_files" />
  <link fromOperation="AlignTopHit"     fromProperty="unaligned_file"                 toOperation="MergeAlignments" toProperty="unaligned_files" />
  <link fromOperation="AlignTopHit"     fromProperty="working_directory"            toOperation="MergeAlignments" toProperty="working_directory" />

  <link fromOperation="input connector" fromProperty="working_directory"            toOperation="AlignMulti" toProperty="working_directory" /> 
  <link fromOperation="input connector" fromProperty="fragment_reads_and_references"	toOperation="AlignMulti" toProperty="reads_and_references" />
  <link fromOperation="input connector" fromProperty="generate_concise"	        	toOperation="AlignMulti" toProperty="generate_concise" />
  <link fromOperation="AlignMulti"      fromProperty="aligned_file"                 toOperation="MergeAlignmentsMulti" toProperty="concise_files" />
  <link fromOperation="AlignMulti"      fromProperty="working_directory"            toOperation="MergeAlignmentsMulti" toProperty="working_directory" />
   
  <link fromOperation="MergeAlignments"	fromProperty="reference1_aligned_file"    	toOperation="FilterResults" toProperty="reference1_top_hit_alignment_file" />
  <link fromOperation="MergeAlignments"	fromProperty="reference2_aligned_file"    	toOperation="FilterResults" toProperty="reference2_top_hit_alignment_file" />
 
  <link fromOperation="MergeAlignmentsMulti"	fromProperty="paired_end1_concise_file"      	toOperation="FilterResults" toProperty="paired_end1_concise_file" />
  <link fromOperation="MergeAlignmentsMulti"	fromProperty="paired_end2_concise_file"      	toOperation="FilterResults" toProperty="paired_end2_concise_file" />
   
  <link fromOperation="input connector" fromProperty="sam_header"	        toOperation="FilterResults" toProperty="sam_header" />
  <link fromOperation="input connector" fromProperty="working_directory"            toOperation="FilterResults" toProperty="working_directory" /> 
  <link fromOperation="input connector" fromProperty="taxonomy_file"                toOperation="FilterResults" toProperty="taxonomy_file" /> 
  <link fromOperation="FilterResults"   fromProperty="filtered_alignment_file"      toOperation="RefCov" toProperty="aligned_bam_file" /> 
  <link fromOperation="FilterResults"   fromProperty="read_count_file"      		toOperation="RefCov" toProperty="read_count_file" /> 
  <link fromOperation="FilterResults"   fromProperty="other_hits_file"      		toOperation="RefCov" toProperty="other_hits_file" /> 

  <link fromOperation="input connector" fromProperty="working_directory"            toOperation="RefCov" toProperty="working_directory" /> 
  <link fromOperation="input connector" fromProperty="regions_file"                 toOperation="RefCov" toProperty="regions_file" /> 
  <link fromOperation="RefCov"          fromProperty="combined_file"                toOperation="Report" toProperty="align_final_file" /> 
  
  <link fromOperation="input connector" fromProperty="working_directory"            toOperation="Report" toProperty="working_directory" />
  <link fromOperation="Report"       	fromProperty="final_file"                   toOperation="output connector" toProperty="final_file" />

<operation name="AlignTopHit" parallelBy="reads_and_references">
    <operationtype commandClass="Genome::Model::Tools::HmpShotgun::AlignMetagenomesMulti" typeClass="Workflow::OperationType::Command">
    </operationtype>
</operation>

<operation name="AlignMulti" parallelBy="reads_and_references">
    <operationtype commandClass="Genome::Model::Tools::HmpShotgun::AlignMetagenomesMulti" typeClass="Workflow::OperationType::Command">
    </operationtype>
</operation>

<operation name="MergeAlignments">
    <operationtype commandClass="Genome::Model::Tools::HmpShotgun::MergeAlignments" typeClass="Workflow::OperationType::Command" />
</operation>

<operation name="MergeAlignmentsMulti">
    <operationtype commandClass="Genome::Model::Tools::HmpShotgun::MergeAlignmentsMulti" typeClass="Workflow::OperationType::Command" />
</operation>

<operation name="FilterResults">
    <operationtype commandClass="Genome::Model::Tools::HmpShotgun::FilterResults" typeClass="Workflow::OperationType::Command" />
</operation>

<operation name="RefCov">
    <operationtype commandClass="Genome::Model::Tools::HmpShotgun::RefCov" typeClass="Workflow::OperationType::Command" />
</operation>

<operation name="Report">
    <operationtype commandClass="Genome::Model::Tools::HmpShotgun::Report" typeClass="Workflow::OperationType::Command" />
</operation>

<operationtype typeClass="Workflow::OperationType::Model">
    <inputproperty>working_directory</inputproperty>
    <inputproperty>generate_concise</inputproperty>
    <inputproperty>reads_and_references</inputproperty>
    <inputproperty>fragment_reads_and_references</inputproperty>
    <inputproperty>reference_sequences</inputproperty>
    <inputproperty>reads_files</inputproperty>
    <inputproperty>taxonomy_file</inputproperty>
    <inputproperty>regions_file</inputproperty>
    <inputproperty>sam_header</inputproperty>
    <outputproperty>final_file</outputproperty>
</operationtype>

</workflow>
