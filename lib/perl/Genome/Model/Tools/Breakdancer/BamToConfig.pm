package Genome::Model::Tools::Breakdancer::BamToConfig;

use strict;
use Genome;
use File::Basename;
use File::Copy;
use File::chdir;

class Genome::Model::Tools::Breakdancer::BamToConfig {
    is  => 'Genome::Model::Tools::Breakdancer',
    has => [
        normal_bam => {
            is  => 'String',
            doc => 'input normal bam file',
            is_input    => 1,
            is_optional => 1,
        },
        tumor_bam  => {
            is  => 'String',
            doc => 'input tumor bam file',
            is_input => 1,
        },
        params     => {
            is  => 'String',
            doc => 'bam2cfg parameters',
            default_value => '-g -h',
            is_input      => 1,
            is_optional   => 1,
        },
        output_file => {
            is  => 'String',
            doc => 'The output breakdancer config file',
            is_output => 1,
            is_input  => 1,
        },
    ],
};

sub help_brief {
    'Tool to create breakdancer config file';
}

sub help_detail {
    return <<EOS
    Tool to create breakdancer config file
EOS
}

sub execute {
    my $self = shift;

    my $out_file = $self->output_file;
    my $out_dir  = dirname $out_file;

    if (-s $out_file) {
        $self->status_message("breakdancer config file $out_file existing. Skip this step");
        return 1;
    }

    unless (-d $out_dir) {
        $self->warning_message("output dir $out_dir not existing. Now try to make it");
        File::Path::make_path($out_dir);
        die "Failed to make out_dir $out_dir\n" unless -d $out_dir;
    }

    unless (Genome::Sys->validate_directory_for_write_access($out_dir)) {
        $self->error_message("$out_dir can not be written to");
        die;
    }

    unless (Genome::Sys->validate_file_for_writing($out_file)) {
        die "output file $out_file can not be written\n";
    }

    my @bam_list;
    for my $type ('tumor', 'normal') {
        my $property = $type . '_bam';
        push @bam_list, $self->$property if $self->$property;
    }
    
    my $bam_string = join ' ', @bam_list;

    my $cfg_cmd = $self->breakdancer_config_command; 
    $cfg_cmd .= ' ' . $self->params . ' ' . $bam_string . ' > '. $out_file;
    $self->debug_message("Breakdancer command: $cfg_cmd");

    {
        local $CWD = $out_dir;  #change current work dir to out_dir so *.insert.histogram can be written there
        my $rv = Genome::Sys->shellcmd(
            cmd => $cfg_cmd,
            input_files  => \@bam_list,
            output_files => [$self->output_file],
            allow_zero_size_output_files => 1,
        );
        unless ($rv) {
            $self->error_message("Running breakdancer config failed using command: $cfg_cmd");
            die;
        }
        $self->status_message("bam2cfg finished ok. The insert.histogram files are created in $out_dir");
    }
    return 1;
}


1;

