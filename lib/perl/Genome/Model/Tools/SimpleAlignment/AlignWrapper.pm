package Genome::Model::Tools::SimpleAlignment::AlignWrapper;

use strict;
use warnings;

use Genome;
use Workflow;
use IO::File;

use Data::Dumper;

class Genome::Model::Tools::SimpleAlignment::AlignWrapper {
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
    $self->debug_message(">>>Running AlignWrapper at ".UR::Context->current->now);
    $self->debug_message("Ref seq: ".$self->reference_name);
    $self->debug_message("Working directory: ".$self->working_directory);
   
    #get these from pp  
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
    
    #resolve this from alignment object    
    my $alignment_file = $alignment->output_dir."/all_sequences.bam";
    $self->debug_message("Setting aligned file to: $alignment_file");  	
    $self->aligned_file($alignment_file);
     
    $self->debug_message("<<<Completed alignment at ".UR::Context->current->now);
    
    return 1;
}
1;
