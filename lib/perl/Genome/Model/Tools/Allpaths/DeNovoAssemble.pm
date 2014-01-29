package Genome::Model::Tools::Allpaths::DeNovoAssemble;

use strict;
use warnings;

use Genome;
use File::Path qw(make_path);

my $PICARD_TOOLS_DIR="$ENV{GENOME_SW_LEGACY_JAVA}/samtools/picard-tools-1.52";

class Genome::Model::Tools::Allpaths::DeNovoAssemble {
    is => 'Genome::Model::Tools::Allpaths::Base',
    has => [
        pre => {
            is => 'Text',
            doc => 'The output directory prefix.' 
        },
        ploidy => {
            is => 'Text',
            doc => 'Ploidy',
        },
        in_group_file => {
            is => 'Text',
            doc => 'in_group_file',
        },
        in_libs_file => {
            is => 'Text',
            doc => 'in_libs_file',
        },
        run => {
            is => 'Text',
            doc => 'name of the run',
            default_value => 'run',
        },
        sub_dir => {
            is => 'Text',
            doc => 'name of subdirectory',
            default_value => 'test',
        },
        overwrite => {
            is => 'Boolean',
            doc => 'should existing results be overwritten',
            default_value => 1,
        },
        reference_name => {
            is => 'String',
            doc => 'name of the reference',
            default_value => 'sample',
        },
        haploidify => {
            is => 'Boolean',
            doc => 'Run Allpaths with haploidify option',
            default_value => 0,
        },
        max_memory_gb => {
            is => 'Number',
            doc => 'Maximum memory that allpaths should use',
            default_value => 0,
        },
        min_contig => {
            is => 'Number',
            doc => 'Minimum length of contigs to report',
            is_optional => 1,
        },
        ca_min_unique => {
            is => 'Number',
            is_optional => 1,
        },
    ],
};

sub help_brief {
    'run ALLPATHS de novo assembler';
}

sub help_detail {
    return;
}

sub execute {
    my $self = shift;

    if (not $self->pre or not -d $self->pre) {
        $self->error_message("Output directory prefix does not exist");
        return;
    }

    my $output_dir = $self->pre."/".$self->reference_name;
    if (! -d $output_dir ) {
        make_path($output_dir);
    }

    # Prepare
    # -need group file (generate based on inputs)
    # -need library file (generate based on inputs)
    # -separate inputs by library
    # -must be at least 2 paired libs, one short, one long
    # -may have add'l long frag lib
    # -may have add'l long jumping paired lib

    if (! -d $output_dir."/data") {
        make_path($output_dir."/data");
    }

    if (Genome::DataSource::GMSchema->has_default_handle) {
        $self->debug_message("Disconnecting GMSchema default handle.");
        Genome::DataSource::GMSchema->disconnect_default_handle();
    }

    my $prepare_cmd = 'ulimit -s 100000 && '. $self->executable_for_version("PrepareAllPathsInputs.pl").' PICARD_TOOLS_DIR='.$PICARD_TOOLS_DIR.' DATA_DIR='.$output_dir.'/data PLOIDY='.$self->ploidy.' IN_GROUPS_CSV='.$self->in_group_file.' IN_LIBS_CSV='.$self->in_libs_file;

    $self->debug_message("Run PrepareAllPathsInput");
    Genome::Sys->shellcmd(cmd => $prepare_cmd); 
    if ($? != 0) {
        $self->error_message("Failed to run PrepareAllPathsInput: $@");
        return;
    }

    my $overwrite;
    if ($self->overwrite) {
        $overwrite="True";
    }
    else {
        $overwrite = "False";
    }
    my $haploidify;
    if ($self->haploidify) {
        $haploidify="True";
    }
    else {
        $haploidify="False";
    }
    my $cmd = 'ulimit -s 100000 && '.$self->executable_for_version("RunAllPathsLG").' PRE='.$self->pre.' REFERENCE_NAME='.$self->reference_name.' DATA_SUBDIR=data RUN='.$self->run.' SUBDIR='.$self->sub_dir.' TARGETS=standard OVERWRITE='.$overwrite.' HAPLOIDIFY='.$haploidify.' MAX_MEMORY_GB='.$self->max_memory_gb;

    if ($self->min_contig) {
        $cmd = "$cmd MIN_CONTIG=".$self->min_contig;
    }

    if ($self->ca_min_unique) {
        $cmd = "$cmd CA_MIN_UNIQUE=".$self->ca_min_unique;
    }

    $self->debug_message("Run ALLPATHS de novo");
    Genome::Sys->shellcmd(cmd => $cmd);
    if ( $? != 0) {
        $self->error_message("Failed to run ALLPATHS de novo: $@");
        return;
    }
    $self->debug_message("Run ALLPATHS de novo...OK");

    return 1;
}

1;

