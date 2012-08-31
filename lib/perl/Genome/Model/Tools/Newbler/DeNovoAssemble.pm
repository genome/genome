package Genome::Model::Tools::Newbler::DeNovoAssemble;

use strict;
use warnings;

use Genome;
use Data::Dumper 'Dumper';

class Genome::Model::Tools::Newbler::DeNovoAssemble {
    is => 'Genome::Model::Tools::Newbler',
    has => [
        version => {
            is => 'Text',
            doc => 'Version of newbler to use',
            valid_values => ['mapasm454_source_03152011', 'mapasm454_source_04282011'],
            #version that allows assembling of fastq files only
        },
        output_directory => {
            is => 'Text',
            doc => 'Directory to output assembly to',
        },
        input_files => {
            is => 'Text',
            is_many => 1,
            doc => 'Input fastq files for assembly',
        },
    ],
    #addl' params: -ace -large -rip -scaffold -finish -a 1 -pairt -cdna -stopjoin -vt -vs -cpu -het
    has_optional => [
        consed => {
            is => 'Boolean',
            doc => 'Creates ace file in edit_dir for consed use',
        },
        rip => {
            is => 'Boolean',
            doc => 'Allow splitting of reads',
        },
        ace => {
            is => 'Boolean',
            doc => 'Creates ace file',
        },
        vt => {
            is => 'Text',
            doc => 'To trim primers,adapters, polyA tails',
        },
        vs => {
            is => 'Text',
            doc => 'To remove reads that match cloning vectors',
        },
        fe => {
            is => 'Text',
            doc => 'Text file of contaminant read names',
        },
    ],
};

sub help_brief {
    'Tool to run newbler assembler and other newbler tools';
}

sub help_detail {
    return <<"EOS"
gmt newbler de-novo-assemble --version mapasm454_source_03152011 --fastq-files reads.fastq,reads2.fastq
gmt newbler de-novo-assemble --version mapasm454_source_03152011 --fastq-files reads.fastq,reads2.fastq --rip
gmt newbler de-novo-assemble --version mapasm454_source_03152011 --fastq-files reads.fastq,reads2.fastq --rip --consed
EOS
}

sub execute {
    my $self = shift;

    my $command;
    if ( not $command = $self->_build_assemble_command ) {
        $self->error_message( "Failed to build assemble command" );
        return;
    }

    if ( system( "$command" ) ) {
        $self->error_message( "Failed to run de-novo-assemble command: $command" );
        return;
    }

    return 1;
}

sub _build_assemble_command {
    my $self = shift;

    my $assembler = $self->path_to_version_run_assembly;
    unless( $assembler ) {
        $self->error_message( "Failed to get path to newbler assembler for version: ".$self->version );
        return;
    }

    if ( $self->ace and $self->consed ) {
        $self->error_message( "--ace and --consed are mutually execlusive in number, please specify one of the two" );
        return;
    }

    if ( $self->vt ) {
        unless ( -s $self->vt ) {
            $self->error_message( "Failed to find trim file specified in --vt option: ".$self->vt );
            return;
        }
    }
    if( $self->vs ) {
        unless( -s $self->vs ) {
            $self->error_message( "Failed to find filter reads file specified in --vs option: ".$self->vs );
            return;
        }
    }
    if ( $self->fe ) {
        unless( -s $self->fe ) {
            $self->error_message( "Failed to find contaminant reads file specified in --fe option: ".$self->fe );
            return;
        }
    }

    my $cmd = $assembler.' -o '.$self->output_directory.' -force'; #force is needed for newbler to put files is existing directories .. it does not remove any existing files is that dir
    $cmd .= ' -consed' if $self->consed;
    $cmd .= ' -rip' if $self->rip;
    $cmd .= ' -ace' if $self->ace;
    $cmd .= ' -vt '.$self->vt if $self->vt;
    $cmd .= ' -vs '.$self->vs if $self->vs;
    $cmd .= ' -fe '.$self->fe if $self->fe;

    my $input_fastqs;
    for my $input_fastq ( $self->input_files ) {
        $input_fastqs .= " $input_fastq";
    }

    $cmd .= $input_fastqs;

    return $cmd;
}

1
