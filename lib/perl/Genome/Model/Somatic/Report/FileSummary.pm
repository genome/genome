package Genome::Model::Somatic::Report::FileSummary;

use strict;
use warnings;

use Genome;
use Path::Class::Dir;

class Genome::Model::Somatic::Report::FileSummary {
    is => 'Genome::Model::Report',
};

sub description {
    my $self = shift();

    return 'A summary of the somatic pipeline run, including linecounts of every output file in the pipeline';
}

sub _add_to_report_xml {
    my $self = shift();

    my $doc = $self->_xml;
    my $files_node = $doc->createElement('files');
    
    my $line_counts = $self->get_file_counts;
    for my $file_name (keys %$line_counts) {
        my $file_node = $files_node->addChild($doc->createElement('file') );
        $file_node->addChild( $doc->createAttribute('file-name', $file_name) );
        $file_node->addChild( $doc->createAttribute('count', $line_counts->{$file_name}) );
    }
    
    $self->_main_node->addChild($files_node);
    
    return 1;
}

sub get_file_counts {
    my $self = shift;
    my $line_counts;
    
    my $can_get_files = $self->build->can('somatic_workflow_input');
    if($can_get_files) {
        $self->build->newest_workflow_instance->ordered_child_instances;
    }

    for my $property_name ($self->files_to_report) {
        
        my $file_name;
        if($can_get_files) {
            $file_name = $self->build->somatic_workflow_input($property_name);
        }

        if ($file_name) {
            if (-e $file_name) {
                $line_counts->{$property_name} = `wc -l < $file_name`;
                chomp $line_counts->{$property_name};
            } else {
                $line_counts->{$property_name} = "File not found";
            }
        } else {
            $line_counts->{$property_name} = '-';
        }
    }
    
    return $line_counts;
}

sub files_to_report {
qw( 
sniper_snp_output
sniper_indel_output
breakdancer_output_file
snp_filter_output
filter_ceu_yri_output
dbsnp_output
loh_output_file
loh_fail_output_file
annotate_output_snp
ucsc_output_snp
ucsc_unannotated_output_snp
indel_lib_filter_single_output
indel_lib_filter_multi_output
annotate_output_indel
ucsc_output_indel
ucsc_unannotated_output_indel
tier_1_snp_file
tier_2_snp_file
tier_3_snp_file
tier_4_snp_file
tier_1_indel_file
tier_2_indel_file
tier_3_indel_file
tier_4_indel_file
tier_1_snp_high_confidence_file
tier_2_snp_high_confidence_file
tier_3_snp_high_confidence_file
tier_4_snp_high_confidence_file
);
}


1;

