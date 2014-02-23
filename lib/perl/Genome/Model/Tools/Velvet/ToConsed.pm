package Genome::Model::Tools::Velvet::ToConsed;

use strict;
use warnings;

use Genome;
use Cwd 'abs_path';
use File::Basename;
use Data::Dumper;

#TODO - make an assembly-directory param and everything else optional 
class Genome::Model::Tools::Velvet::ToConsed {
    is           => 'Command::V2',
    has          => [
        assembly_directory => {
            is      => 'String',
            doc     => 'Main assembly directory above edit_dir and chromat_dir',
        },
    ],
    has_optional => [
        fastq_file  => {
            is      => 'String',
            doc     => 'Input fastq file to start the velvet assembly',
            is_mutable => 1,
        },
        afg_file    => {
            is      => 'String', 
            doc     => 'input velvet_asm.afg file path(s)',
            is_mutable => 1,
        },
        velvet_seq_file => {
            is      => 'String',
            doc     => 'velvet generated Sequences file',
            is_mutable => 1,
        },
        out_acefile => {
            is      => 'String', 
            doc     => 'name for output acefile, default is ./velvet_asm.ace',
            is_mutable => 1,
        },
        chunk_size  => {
            is      => 'Integer',
            doc     => 'number of fastq sequences each chunk',
            default => 10000,
        },
	make_traces      => {
	    is      => 'Boolean',
	    doc     => 'Make a trace file for each read .. this is NOT recommanded for assemblies with 10,000 or more reads',
	    default => 0,
	},
        no_ace      => {
            is      => 'Boolean',
            doc     => 'Do not make ace file',
        },
        no_phdball  => {
            is      => 'Boolean',
            doc     => 'Do not make phdball file',
        },
        time_stamp  => {
            is      => 'String',
            doc     => 'Time stamp for reads in ace and phd ball files',
        },
    ],
};

sub help_brief {
    'This tool converts velvet assembly to acefile, then convert input read fastq file into phdball and scf files',
}

sub help_detail {
    return <<EOS
gmt velvet to-consed --assembly-directory /gscmnt/gc9999/assembly/e_coli --no-phdball
gmt velvet to-consed --assembly-directory /gscmnt/gc9999/assembly/e_coli --no-ace
gmt velvet to-consed --assembly-directory /gscmnt/gc9999/assembly/e_coli --fastq-file /gscmnt/sata001/assembly/e_coli_2/input.fastq

EOS
}

sub execute {
    my $self = shift;
    my $time = $self->time_stamp ? $self->time_stamp : localtime;
    
    return unless $self->validate_inputs;

    #make ace file
    unless ( $self->no_ace ) {
        my $to_ace  = Genome::Model::Tools::Velvet::ToAce->create(
            assembly_directory => $self->assembly_directory,
            afg_file    => $self->afg_file,
            out_acefile => $self->out_acefile,
            seq_file    => $self->velvet_seq_file,
            time        => $time,
            );
        if ( not $to_ace->execute ) {
            $self->error_message( "Failed to create ace file" );
            return;
        }
        $self->debug_message( "Successfully made ace file" );
    }

    #make phdball file
    if ( not $self->no_phdball ) {
        my %to_phdscf_params = (
            fastq_file => $self->fastq_file,
            ball_file  => $self->assembly_directory.'/edit_dir/phd.ball',
            base_fix   => 1,
            time       => $time,
            chunk_size => $self->chunk_size,
            );

        #conditon for creating scf files is if scf_dir (chromat_dir) is specified in gmt fastq to-phdball-scf
        if ( $self->make_traces ) {
            $to_phdscf_params{scf_dir} = $self->assembly_directory.'/chromat_dir';
        }
        
        my $to_phdscf = Genome::Model::Tools::Fastq::ToPhdballScf::Chunk->create( %to_phdscf_params );
        if ( not $to_phdscf->execute ) {
            $self->error_message( "Failed to create phdball file" );
            return;
        }
        $self->debug_message( "Successfully made phdball file" );
    }

    return 1;
}

sub validate_inputs {
    my $self = shift;

    #validate assembly dir .. make edit_dir
    unless ( -d $self->assembly_directory ) {
        $self->error_message( "Invalid assembly directory: ".$self->assembly_directory.' does not exist' );
        return;
    }

    unless ( -d $self->assembly_directory.'/edit_dir' ) {
        Genome::Sys->create_directory( $self->assembly_directory.'/edit_dir' );
        unless ( -d $self->assembly_directory.'/edit_dir' ) {
            $self->error_message( "Failed to create edit_dir in assembly dir, check directory permission ".$self->assembly_directory );
            return;
        }
    }
    
    #validate input fasta, afg and seq files
    my %param_file_names = $self->param_and_default_file_names;
    for my $param_name ( keys %param_file_names ) {
        if ( $self->$param_name ) {
            if ( not -s $self->$param_name ) {
                $self->error_message( "Failed to find user supplied file: ".$self->$param_name );
                return;
            }
        }
        else {
            my $file = $self->assembly_directory.'/'.$param_file_names{$param_name};
            $self->$param_name( $file );
            if ( not -s $file ) {
                $self->error_message( "Failed to find file: ".$file );
                return;
            }
        }
    }

    #validate out ace file
    if ( not $self->out_acefile ) {
        $self->out_acefile( $self->assembly_directory.'/edit_dir/velvet_asm.ace' );
    }

    return 1;
}

sub param_and_default_file_names {
    return (
        'fastq_file'      => 'input.fastq',
        'afg_file'        => 'velvet_asm.afg',
        'velvet_seq_file' => 'Sequences',
    );
}

1;

