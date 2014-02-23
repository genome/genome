package Genome::Model::Tools::Sam::MarkDuplicates;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;
use File::Basename;

class Genome::Model::Tools::Sam::MarkDuplicates {
    is  => 'Genome::Model::Tools::Sam',
    has => [
        file_to_mark => {
            is  => 'String',
            doc => 'The bam file to mark ',
        },
        marked_file => {
            is  => 'String',
            doc => 'The resulting marked file',
        },
        metrics_file => {
            is  => 'String',
            doc => 'The resulting metrics file with deduplication statistics',
        },
        remove_duplicates => {
            is  => 'Integer',
            doc => 'Denoting whether the output file should have duplicates removed.  Default is 1, duplicates will be removed.',
            default_value => 1,
	        is_optional => 1
        },
        max_jvm_heap_size => {
            is  => 'Integer',
            doc => 'The size in gigabytes of the Java Virtual Machine maximum memory allocation.',
            default_value => 2,
	        is_optional => 1,
        },
        max_permgen_size => {
            is => 'Integer',
            doc => 'the maximum memory (Mbytes) to use for the "permanent generation" of the Java heap (e.g., for interned Strings)',
            is_optional => 1,
            default_value => 256,
        },
        assume_sorted => {
            is  => 'Integer',
            doc => 'Assume the input file is coordinate order sorted.  Default is 1, true.',
	        default_value => 1,
	        is_optional => 1,
        },
        validation_stringency => {
            is => 'String',
            doc => 'Controls how strictly to validate a SAM file being read.',
            is_optional => 1,
            default_value => 'SILENT',
            valid_values => Genome::Model::Tools::Picard->__meta__->property('validation_stringency')->valid_values,
        },
        log_file => {
            is  => 'String',
            doc => 'The stdout of the mark duplicates tool',
	        is_optional => 1
        },
        tmp_dir => {
            is  => 'String',
            doc => 'The temporary working directory.  Provide this if you are marking duplicates on a whole genome bam file.',
	        is_optional => 1
        },
        dedup_version => {
            is => 'String',
            doc => 'The version of picard to use for deduplication',
            is_optional => 1,
        },
        dedup_params => {
            is => 'String',
            doc => 'All parameters of picard to use for deduplication',
            is_optional => 1,
        },
        include_comment => {
            is => 'Text',
            doc => 'comment to include as a @CO in the BAM header',
            is_optional => 1,
        },
    ],
    doc => 'This module is now just a wrapper for `gmt picard mark-duplicates` but preserves the interface and outputs of this tool, which provided duplicate functionality',
};

sub help_brief {
    'Tool to mark or remove duplicates from BAM or SAM files.';
}

sub help_detail {
    return <<EOS
    Tool to mark or remove duplicates from BAM or SAM files used by the pipeline.  This is just another interface to `gmt picard mark-duplicates`.
EOS
}

sub execute {
    my $self = shift;

    my $input_file = $self->file_to_mark;
    my $result = $self->marked_file; 
    
    $self->debug_message("Attempting to mark duplicates." );
   
    unless (-e $input_file)  {
       $self->error_message("Source file $input_file not found!");
       return;
    }
    
    if (-e $result )  {
       $self->error_message("The target file already exists at: $result . Please remove this file and rerun to generate a new merged file.");
       return;
    }
    
    #merge those Bam files...BAM!!!
    my $now = UR::Context->current->now;
    $self->debug_message(">>> Beginning mark duplicates at $now");

    my %params = (
        input_file             => $self->file_to_mark,
        output_file            => $self->marked_file,
        metrics_file           => $self->metrics_file,
        maximum_memory         => $self->max_jvm_heap_size,
        maximum_permgen_memory => $self->max_permgen_size,
        log_file               => $self->log_file,
        temp_directory         => $self->tmp_dir,
    );

    $params{include_comment} = $self->include_comment if $self->include_comment;
    $params{use_version} = $self->dedup_version if $self->dedup_version;

    my %markdup_params = (
        remove_duplicates => $self->remove_duplicates,
        assume_sorted     => $self->assume_sorted,
    );

    my @valid_markdup_params = qw(max_sequences_for_disk_read_ends_map); #list here each time adding new parameters

    my $dedup_params = $self->dedup_params;
    if ($dedup_params) {
        $dedup_params =~ s/^\s*//;
        my %given_params = split /\s+|\=/, $dedup_params;
        for my $given_param (keys %given_params) {
            my $preset_param = lc($given_param);
            if ($self->can($preset_param)) {
                if ($self->$preset_param) {
                    $self->warning_message("Overwrite $given_param to ". $given_params{$given_param}); 
                    $markdup_params{$preset_param} = $given_params{$given_param};
                }
                else {
                    $self->error_message("No preset value for picard markdup parameter $given_param");
                    return;
                }
            }
            else {
                unless (grep{lc($given_param) eq $_}@valid_markdup_params) {
                    $self->error_message("parameter $given_param is not a valid picard markdup parameter or implement it into Genome::Model::Tools::Picard::MarkDuplicates");
                    return;
                }
                unless ($given_params{$given_param}) {
                    $self->error_message("No value provided for picard markdup parameter $given_param");
                    return;
                }
                $markdup_params{lc($given_param)} = $given_params{$given_param};
            }
        }
    }
    %params = (%params, %markdup_params);

    my $picard_cmd = Genome::Model::Tools::Picard::MarkDuplicates->create(%params);
    my $md_rv = $picard_cmd->execute();
    
    $self->debug_message("Mark duplicates return value: $md_rv");
    if ($md_rv != 1) {
        $self->error_message("Mark duplicates error!  Return value: $md_rv");
    } 
    else {
        $self->debug_message("Success.  Duplicates marked in file: $result");
    }

    $now = UR::Context->current->now;
    $self->debug_message("<<< Completing mark duplicates at $now.");
    
    return 1;
}


1;
