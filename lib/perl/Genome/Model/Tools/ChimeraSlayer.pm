package Genome::Model::Tools::ChimeraSlayer;

use strict;
use warnings;

use Genome;
use Data::Dumper 'Dumper';

class Genome::Model::Tools::ChimeraSlayer {
    is => 'Command',
    has => [
        exec_dir => {
            is => 'Text',
            doc => 'Directory to execute program in',
            is_optional => 1,
        },
        query_NAST => {
            is => 'Text',
            doc => 'Multi-fasta file containing query sequences in alignment format',
            is_optional => 1,
        },
        db_NAST => {
            is => 'Text',
            is_optional => 1,
            doc => 'Database in NAST format (default: file installed 2010-07-07)',
        },
        db_FASTA => {
            is => 'Text',
            is_optional => 1,
            doc => 'Database in fasta format (default: file installed 2010-07-07)',
        },
        n => {
            is => 'Number',
            is_optional => 1,
            doc => 'number of top matching database sequences to compare to (default 15)',
        },
        R => {
            is => 'Number',
            is_optional => 1,
            doc => 'min divergence ratio (default 1.007)',
        },
        P => {
            is => 'Number',
            is_optional => 1,
            doc => 'min percent identity among matching sequences (default 90)',
        },
        M => {
            is => 'Number',
            is_optional => 1,
            doc => 'match score (default 5)',
        },
        N => {
            is => 'Number',
            is_optional => 1,
            doc => 'mismatch penalty (default -4)',
        },
        Q => {
            is => 'Number',
            is_optional => 1,
            doc => 'min query coverage by matching database sequence (default 70)',
        },
        T => {
            is => 'Number',
            is_optional => 1,
            doc => 'maximum traverses of the multiple alignment (default 1)',
        },
        windowSize => {
            is => 'Number',
            is_optional => 1,
            doc => 'window size (default 50)',
        },
        windowStep => {
            is => 'Number',
            is_optional => 1,
            doc => 'winow step (default 5)',
        },
        minBS => {
            is => 'Number',
            is_optional => 1,
            doc => 'minimum bootstrap support for calling chimera (default 90)',
        },
        S => {
            is => 'Number',
            is_optional => 1,
            doc => 'percent of SNPs to sample on each side of breakpoint for computing bootstrap support (default 10)',
        },
        num_parents_test => {
            is => 'Number',
            is_optional => 1,
            doc => 'number of potential parents to test for chimeras (default 3)',
        },
        MAX_CHIMERA_PARENT_PER_ID => {
            is => 'Number',
            is_optional => 1,
            doc => 'Chimera/Parent alignments with perID above this are considered non-chimeras (default 100)',
        },
        printFinalAlignments => {
            is => 'Boolean',
            is_optional => 1,
            doc => 'shows alignment between query sequence and pair of candidate chimera parents (default: off)',
        },
        printCSalignments => {
            is => 'Boolean',
            is_optional => 1,
            doc => 'print ChimeraSlayer alignments in ChimeraSlayer output (default: off)',
        }
    ],
};

sub help_brief {
    'Tool to run chimera detector: chimera_slayer',
}

sub help_detail {
    return <<"EOS"
gmt chimera-slayer --query-NAST chims.NAST --exec-dir /gscmnt/100/mc16s --printCSalignments
gmt chimera-slayer --query-NAST chims.NAST --exec-dir /gscmnt/103/mc16s --S 10 --windowSize 50
EOS
}

sub execute {
    my $self = shift;

    unless( $self->exec_dir and -d $self->exec_dir ) {
        $self->error_message("Failed to find exec_dir dir: ".$self->exec_dir);
        return;
    }

    unless( $self->query_NAST and -s $self->query_NAST ) {
        $self->error_message("Failed to find query_NAST file or file is zero size: ".$self->query_NAST);
        return;
    }
    my $cmd = $self->command_string;
    if ( not $cmd ) {
        $self->error_message("Failed to build command string to run chimera slayer");
        return;
    }

    $self->debug_message("Running chimera slayer with command:\n $cmd");
    #print $cmd."\n"; return 1;

    my $rv = eval { Genome::Sys->shellcmd( cmd => $cmd ) };
    if( not $rv ) {
        $self->error_message("Failed to run chimera slayer with command: $cmd");
        return;
    }

    return 1;
}

sub command_string {
    my $self = shift;

    my $cmd = $self->path_to_chimera_slayer;

    my $db_nast_file = ( $self->db_NAST ) ? $self->db_NAST : $self->version_db_nast_file;
    if ( not -s $db_nast_file ) {
        $self->error_message("Failed to find file or file is zero size: $db_nast_file");
        return;
    }
    my $db_fasta_file = ( $self->db_FASTA ) ? $self->db_FASTA : $self->version_db_fasta_file;
    if ( not -s $db_fasta_file ) {
        $self->error_message("Failed to find file or file is zero size: $db_fasta_file");
        return;
    }
    
    #required by gmt
    $cmd .= ' --query_NAST '.$self->query_NAST;
    $cmd .= ' --db_NAST '.$db_nast_file;
    $cmd .= ' --db_FASTA '.$db_fasta_file;
    $cmd .= ' --exec_dir '.$self->exec_dir;

    #optional with defaults values already hardcoded in
    for my $param_name ( $self->optional_params_with_defaults ) {
        if ( $self->$param_name ) {
            $cmd .= " --$param_name ".$self->$param_name;
        }
    }

    #boolean
    $cmd .= ' --printCSalignments' if $self->printCSalignments;
    $cmd .= ' --printFinalAlignments' if $self->printFinalAlignments;

    return $cmd;
}

sub optional_params_with_defaults {
    return qw/
        n
        R
        P
        M
        N
        Q
        T
        windowSize
        windowStep
        minBS
        S
        num_parents_test
        MAX_CHIMERA_PARENT_PER_ID
    /;
}

sub version {
    return '2010-07-07';
}

sub path_to_db_files {
    return '/gscmnt/gc4096/info/reference_sequences/chimera-detector-16SrRNA/';
}

sub version_db_nast_file {
    return $_[0]->path_to_db_files.'/'.$_[0]->version.'/rRNA16S.gold.NAST_ALIGNED.fasta';
}

sub version_db_fasta_file {
    return $_[0]->path_to_db_files.'/'.$_[0]->version.'/rRNA16S.gold.fasta';
}

sub path_to_chimera_slayer {
    my $self = shift;

    my $script = $ENV{GENOME_SW} . '/broad/ChimeraSlayer/ChimeraSlayer.pl';
    unless( -x $script ) {
        $self->error_message("Failed to find script or script is not executable: $script");
        return;
    }

    return $script;
}

1;
