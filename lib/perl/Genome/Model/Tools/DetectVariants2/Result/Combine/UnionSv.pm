package Genome::Model::Tools::DetectVariants2::Result::Combine::UnionSv;

use warnings;
use strict;

use Genome;
use Tie::File;
use File::Basename;

class Genome::Model::Tools::DetectVariants2::Result::Combine::UnionSv{
    is  => 'Genome::Model::Tools::DetectVariants2::Result::Combine',
    doc => 'Union svs into one file',
};

sub _needs_symlinks_followed_when_syncing { 0 };
sub _working_dir_prefix { 'union-sv' };
sub resolve_allocation_disk_group_name { $ENV{GENOME_DISK_GROUP_MODELS} };
sub allocation_subdir_prefix { 'union_sv' };
sub _variant_type { 'svs' };

sub _combine_variants {
    my $self = shift;
    my $base_name = 'svs.hq';

    my ($dir_a, $dir_b) = ($self->input_directory_a, $self->input_directory_b);
    my @files = map{$_ .'/'. $base_name}($dir_a, $dir_b);
    my $output_file = $self->temp_staging_directory . '/' . $base_name;

    if (-z $files[0] and -z $files[1]) {
        $self->warning_message("0 size of $base_name from both input dir. Probably for testing of small bams");
        `touch $output_file`;
    }
    else {
        my $input_files = join ',', @files;
        my $union_command = Genome::Model::Tools::Breakdancer::MergeFiles->create(
            input_files => $input_files,
            output_file => $output_file,
        );
    
        unless ($union_command->execute) {
            $self->error_message("Error executing union command");
            die $self->error_message;
        }
    }

    $self->debug_message("Now make copy to sv merge outputs");
    
    for my $dir ($dir_a, $dir_b) {
        $self->debug_message("dir is $dir");
        my $subclass_name = $self->_get_sr_name($dir, 'subclass_name');
        if ($subclass_name =~ /Combine\:\:UnionSv/) {
            $self->debug_message("This is a union directory: $dir");
            $self->_copy_file($dir);
            next;
        }

        #If not combine result, it must have a detector_name
        my $dir_type;
        my $detector_name = $self->_get_sr_name($dir, 'detector_name');
        if ($detector_name =~ /Squaredancer/) {
            $dir_type = 'squaredancer.';
        }
        elsif ($detector_name =~ /Breakdancer/) { 
            my $param_list = $self->_get_sr_name($dir, 'detector_params');
            next unless $param_list;
                        
            if ($param_list =~ /\-t/) {
                $dir_type = 'Inter.';
            }
            elsif ($param_list =~ /\-o/) {
                $dir_type = 'Intra.';
            }
            else {
                $self->warning_message("Failed to determine dir_type for params ($param_list).");
                next;
            }
        }
        else {
            $self->error_message("Unknown detector: $detector_name");
            die;
        }
        $self->_copy_file($dir, $dir_type);
    }        

    #MG needs a final merge.annot.somatic file to only list somatic events.
    map{$self->_get_combine_somatic_list($_)}qw(annot file);
    
    return 1;
}

sub _get_combine_somatic_list {
    my ($self, $type) = @_;
    my @merge_files;
    my $base_name;

    if ($type eq 'annot') {
        $base_name   = 'svs.hq.merge.annot.somatic'; 
        @merge_files = glob($self->temp_staging_directory . "/*.file.annot");
    }
    elsif ($type eq 'file') {
        $base_name   = 'svs.merge.file.somatic';
        @merge_files = glob($self->temp_staging_directory . "/*.merge.file");
    }
    else {
        die $self->error_message('Wrong type');
    }

    if (@merge_files) {
        my $output_file  = $self->temp_staging_directory . '/' . $base_name;
        my $out_fh    = Genome::Sys->open_file_for_writing($output_file) 
            or die "Failed to open $output_file for writing\n";
        my ($temp_fh, $temp_file) = Genome::Sys->create_temp_file;

        my $cat_cmd = "cat @merge_files";
        my $cat_fh = IO::File->new($cat_cmd . "|");
        binmode $cat_fh;
        my $header;

        while (my $line = $cat_fh->getline) {
            if ($line =~ /^#/) {#no header
                $header = $line unless $header;
                next;
            }
            my @columns = split /\s+/, $line;
            unless ($columns[11] =~ /normal/) {  #no normal
                shift @columns;   #remove the useless first column
                $line = join "\t", @columns;
                $temp_fh->print($line."\n");
            }
        }
        map{$_->close}($temp_fh, $cat_fh);

        my $sort = Genome::Model::Tools::Snp::Sort->create(
            snp_file    => $temp_file,
            output_file => $output_file,
        );
        $self->warning_message('Failed to execute sort command') unless $sort->execute;

        my @header_names = split /\s+/, $header;
        shift @header_names;
        $header = '#' . join "\t", @header_names;
    
        tie my @array, 'Tie::File', $output_file or die "Failed to tie $output_file";
        unshift @array, $header;
        untie @array;
    }
    else {
        $self->warning_message("There are no $type file to concat");
    }

    return 1;
}


sub _get_sr_name {
    my ($self, $dir, $property) = @_;
    my $sr = $self->_get_sr($dir);
    my $name = $sr->$property;
    unless ($name) {
        $self->error_message("Failed to get $property for sr: ".$sr->id);
        die;
    }
    return $name;
}


sub _get_sr {
    my ($self, $dir) = @_;
    my $dir_target = -l $dir ? readlink($dir) : $dir;

    unless ($dir_target) {
        $self->warning_message("Failed to read the target from symlink ($dir).");
        return;
    }

    $dir_target =~ s/\/$//;
    my $result_id = (split('-', $dir_target))[-1];
    my $result    = Genome::SoftwareResult->get($result_id);

    unless ($result) {
        $self->error_message("Failed to get result for result ID ($result_id)");
        die;
    }

    return $result;
}


sub _copy_file {
    my ($self, $dir, $dir_type) = @_;
    my @file_names = map{$self->_variant_type.'.merge.'.$_}qw(fasta file file.annot out);

    for my $i (0..$#file_names) {
        my $file_type = $file_names[$i];
        if ($dir_type) {
            my $dest_name = $dir_type . $file_type;
            my $target = $dir .'/'. $file_type;
            my $dest   = $self->temp_staging_directory ."/$dest_name";
            unless (-e $target) {
                $self->warning_message("Target file: $target not existing");
                next;
            }
            unless (Genome::Sys->copy_file($target, $dest)) {
                $self->warning_message("Failed to copy $target to $dest");
            }
        }
        else { # for the case in sub union directory
            my @targets = glob($dir .'/*'. $file_type);
            unless (@targets) {
                $self->warning_message("$file_type can not be found in $dir");
                next;
            }
            unless (@targets == 2) {
                $self->warning_message("Expect 2 $file_type not ". scalar @targets);
            }
            for my $target (@targets) {
                my $dest_name = basename($target);
                my $dest      = $self->temp_staging_directory ."/$dest_name";
                unless (Genome::Sys->copy_file($target, $dest)) {
                    $self->warning_message("Failed to copy $target to $dest");
                }
            }
        }
    }
    return 1;
}


sub _validate_output {
    my $self = shift;
    my $variant_type = $self->_variant_type;
    my $out_file     = $self->temp_staging_directory.'/'.$variant_type.'.hq';

    for my $file ($out_file) {
        unless (-e $out_file) {
            die $self->error_message("Fail to find valid output file: $out_file");
        }
    }
    return 1;
}

1;
