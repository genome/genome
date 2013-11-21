package Genome::Model::SomaticValidation::Command::PrepareIndelValidationInput;

use strict;
use warnings;

use Genome;
use File::Spec;

class Genome::Model::SomaticValidation::Command::PrepareIndelValidationInput {
    is => 'Genome::Command::Base',
    has => [
        build => {
            is => 'Genome::Model::Build::SomaticValidation',
            id_by => 'build_id',
        },
        build_id => {
            is => 'Text',
            is_input => 1,
            doc => 'build id of SomaticValidation model',
        },
    ],
    has_transient_calculated => [
        _large_indel_directory => {
            calculate_from => 'build',
            calculate => q| return File::Spec->join($build->data_directory, 'validation/large_indel') |,
        },
        _small_indel_directory => {
            calculate_from => 'build',
            calculate => q| return File::Spec->join($build->data_directory, 'validation/small_indel') |,
        },
        _large_indel_output => {
            calculate_from => '_large_indel_directory',
            calculate => q| return File::Spec->join($_large_indel_directory, 'indels_to_validate.bed') |,
            is_output => 1,
        },
        _small_indel_output => {
            calculate_from => '_small_indel_directory',
            calculate => q| return File::Spec->join($_small_indel_directory, 'indels_to_validate.bed') |,
            is_output => 1,
        },
        _small_indel_annotation_output => {
            calculate_from => '_small_indel_directory',
            calculate => q| return File::Spec->join($_small_indel_directory, 'indels_to_validate.annotation') |,
            is_output => 1,
        },
    ],
    has_param => [
        lsf_queue => {
            default => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER_ALT},
        },
    ],
};

sub execute {
    my $self = shift;

    return 1 if $self->build and not $self->build->normal_sample;

    $self->_create_output_directories();

    $self->_create_input_beds;

    return 1;
}

# Consolidate indels from both DV2 detection in this build and from indel_variant_list, splitting the list into small and large indels
sub _create_input_beds {
    my $self = shift;
    my $build = $self->build;

    my ($large_fh, $large_file) = Genome::Sys->create_temp_file;
    my ($small_fh, $small_file) = Genome::Sys->create_temp_file;
    my ($small_annotation_fh, $small_annotation_file) = Genome::Sys->create_temp_file;
    my @fh = ($large_fh, $small_fh, $small_annotation_fh);
    my @input_files = ($large_file, $small_file, $small_annotation_file);
    my @output_files = ($self->_large_indel_output, $self->_small_indel_output, $self->_small_indel_annotation_output);

    # Consolidate the indels into large and small files
    my $detected_indels = $self->_detected_indels;
    if ($detected_indels) {
        my $ifh = Genome::Sys->open_file_for_reading($detected_indels);
        while (my $line = $ifh->getline) {
           $self->_process_one_indel($line, $large_fh, $small_fh, $small_annotation_fh); 
        }
    }

    if ($build->indel_variant_list) {
        my $indel_list = $build->indel_variant_list->output_dir . "/indels.hq.bed";
        unless (-s $indel_list) {
            die $self->error_message("Build has an indel_variant_list set but there is no indel file at the expected location $indel_list");
        }

        my $ifh = Genome::Sys->open_file_for_reading($indel_list);
        while (my $line = $ifh->getline) {
           $self->_process_one_indel($line, $large_fh, $small_fh, $small_annotation_fh); 
        }
    }

    # Sort!
    while (@fh) {
        my $input_file = shift @input_files;
        my $output_file = shift @output_files;
        my $fh = shift @fh;
        $fh->close;
        my $sort_cmd = Genome::Model::Tools::Joinx::Sort->create(
            input_files => [$input_file],
            output_file => $output_file,
        );
        unless ($sort_cmd->execute) {
            die $self->error_message("Failed to sort indel output $input_file into $output_file");
        }
    }

    return 1;
}

# Parse an indel and put it in the correct bin by size
sub _process_one_indel {
    my ($self, $line, $large_fh, $small_fh, $small_annotation_fh) = @_;

    chomp($line);
    my ($chr, $start, $stop, $ref_var, @everything_else) = split(/\t/, $line);
    my ($size, $annostart, $annostop, $type);
    my ($ref,$var) = split(/\//,$ref_var);
    if ($ref eq '-' || $ref eq '0' || $ref eq '*') { #ins
        #count number of bases inserted
        $ref = "-";
        $size = length($var);
        $annostart = ($start);
        $annostop = ($stop + 1);
        $type = 'INS';
    }
    elsif ($var eq '-' || $var eq '0' || $var eq '*') { #del
        $var = "-";
        $size = length($ref);
        $annostart = ($start + 1);
        $annostop = ($stop);
        $type = 'DEL';
    }
    else {
        $self->error_message("Line $line has wrong insertion or deletion nomenclature. Either ref or var should be 0 or -\n");
        next;
    }

    if ( $size > 0 && $size <= 2) {
        #Add 1 bp padding to bed because we just want to look at regions
        $start-=1;
        $stop+=1;
        $small_fh->print(join("\t", ($chr,$start,$stop,"$ref/$var") ) . "\n");
        $small_annotation_fh->print(join("\t", ($chr,$annostart,$annostop,$ref,$var) ) . "\n");
    }
    elsif ($size > 2) {
        $large_fh->print(join("\t", ($chr, $start, $stop, "$ref/$var", $type) ) . "\n");
    } else {
        die $self->error_message("Weird size found ($size) for indel found on line $line");
    }

    return 1;
}

sub _create_output_directories {
    my $self = shift;
    for my $dir ($self->_small_indel_directory, $self->_large_indel_directory) {
        Genome::Sys->create_directory($dir);
    }

    return 1;
}

# file representing indels detected in the DV2 portion of this build
sub _detected_indels {
    my $self = shift;
    my $file = $self->build->data_directory . "/variants/indels.hq.bed";
    unless (-e $file) {
        return;
    }
    return $file;
}

1;
