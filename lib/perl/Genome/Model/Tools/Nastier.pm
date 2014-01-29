package Genome::Model::Tools::Nastier;

use strict;
use warnings;

use Genome;
use Data::Dumper 'Dumper';

class Genome::Model::Tools::Nastier {
    is => 'Command',
    has => [
        query_FASTA => {
            is => 'Text',
            doc => 'query database in FASTA format (seqs to NAST-align)',
            is_optional => 1,
        },
        db_NAST => {
            is => 'Text',
            is_optional => 1,
            doc => 'reference database in NAST format (default: file installed 2010-07-07)',
        },
        db_FASTA => {
            is => 'Text',
            is_optional => 1,
            doc => 'reference database in FASTA format (default: file installed 2010-07-07)',
        },
        num_top_hits => {
            is => 'Number',
            is_optional => 1,
            doc => 'number of top hits to use for profile-alignment (default: 10)',
        },
        Evalue => {
            is => 'Number',
            is_optional => 1,
            doc => 'Evalue cutoff for top hits (default: 1e-50)',
        },
        output_file => {
            is => 'Text',
            is_optional => 1,
            doc => 'File to output results to, otherwise prints to screen',
        },
    ],
};

sub help_brief {
    'Tool to run chimera detector: nastier';
}

sub help_detail {
    return <<"EOS"
gmt nastier --query-FASTA chims.NAST.fasta --Evalue 1e-50 --output-file chims.NAST.align.cs
EOS
}

sub execute {
    my $self = shift;

    unless ( $self->query_FASTA and -s $self->query_FASTA ) {
        $self->error_message("Failed to find query FASTA file or file is zero size: ".$self->query_FASTA);
        return;
    }

    my $cmd = $self->build_command_string;
    if ( not $cmd ) {
        $self->error_message("Failed to build command string to run Nastier");
        return;
    }

    $self->debug_message("Running nastier with command:\n$cmd");

    my $rv = eval { Genome::Sys->shellcmd( cmd => $cmd ) };
    if( not $rv ) {
        $self->error_message("Failed to run nastier with command: $cmd");
        return;
    }

    return 1;
}

sub build_command_string {
    my $self = shift;

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

    my $cmd = $self->path_to_nastier;
    #mandidatory
    $cmd .= ' --query_FASTA '.$self->query_FASTA;
    $cmd .= ' --db_NAST '.$db_nast_file;
    $cmd .= ' --db_FASTA '.$db_fasta_file;
    #optional with default already set in orig script
    $cmd .= ' --num_top_hits '.$self->num_top_hits if $self->num_top_hits;
    $cmd .= ' --Evalue '.$self->Evalue if $self->Evalue;
    #optional
    $cmd .= ' > '.$self->output_file if $self->output_file; #otherwise prints to screen

    return $cmd;
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

sub path_to_nastier {
    my $self = shift;

    my $script = $ENV{GENOME_SW} . '/broad/NAST-iEr/run_NAST-iEr.pl';
    if ( not -x $script ) {
        $self->error_message("Failed to find script of script is not executable: $script");
        return;
    }
    return $script;
}

1;
