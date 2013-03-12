package Genome::Model::Tools::Allpaths::EfastaToFasta;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Allpaths::EfastaToFasta {
    #is => 'Genome::Model::Tools::Allpaths::Base',
    is => 'Command::V2',
    has => [
        assembly_directory => {
            is => 'Text',
            doc => 'Path to allpaths assembly.',
        },
        efasta_in => {
            is => 'Text',
            is_optional => 1,
            doc => 'Allpaths created assembly.efasta file',
        },
        clean_fasta_out => {
            is => 'Text',
            is_optional => 1,
            doc => 'Output fasta with ambiguous bases removed',
        },
    ],
};

sub help_brief {
    return 'Produce metrics for allpaths assemblies'
}

sub help_detail {
    return;
}

sub __errors__ {
    my $self = shift;

    my @errors = $self->SUPER::__errors__(@_);
    return @errors if @errors;

    if ( $self->assembly_directory ) {
        if ( not -d $self->assembly_directory ) {
            push @errors, UR::Object::Tag->create(
                type       => 'invalid',
                properties => [qw/ assembly_directory /],
                desc       => 'The assembly_directory is not a directory!', 
            );
            return @errors;
        }
    }
    
    if ( $self->efasta_in ) {
        if ( not -s $self->efasta_in ) {
            push @errors, UR::Object::Tag->create(
                type       => 'invalid',
                properties => [qw/ efasta_in /],
                desc       => 'The efasta file is missing or empty',
            );
            return @errors;
        }
    } else {
        my $scaffolds_efasta = $self->resolve_scaffolds_efasta;
        if ( not $scaffolds_efasta ) {
            push @errors, UR::Object::Tag->create(
                type => 'invalid',
                properties => [qw/ assembly_directory /],
                desc => $self->error_message,
            );
        }
    }

    return @errors;
}

sub execute {
    my $self = shift;

    my $scaffolds_efasta = ( $self->efasta_in )
        ? $self->efasta_in
        : $self->resolve_scaffolds_efasta ;
    
    my $reader = Genome::Model::Tools::Sx::PhredEnhancedSeqReader->create(
        file => $scaffolds_efasta
    );    
    
    my $clean_fasta_out = ( $self->clean_fasta_out )
        ? $self->clean_fasta_out
        : $self->resolve_clean_fasta_out ;

    unlink $clean_fasta_out;

    my $writer = Genome::Model::Tools::Sx::PhredWriter->create(
        file => $clean_fasta_out
    );

    while( my $seq = $reader->read ) {
        $writer->write( $seq );
    }

    return 1;
}

sub resolve_clean_fasta_out  {
    my $self = shift;

    return $self->assembly_directory.'/clean.final.assembly.fasta';
}

sub resolve_scaffolds_efasta {
    my $self = shift;

    my @files = glob($self->assembly_directory.'/*/data/*/ASSEMBLIES/*/final.assembly.efasta');
    if ( not @files ) {
        $self->error_message('No scaffold efasta file found in '.$self->assembly_directory);
        return;
    }
    elsif ( @files > 1 ) {
        $self->error_message("More than one scaffold assembly file found!\n".join("\n", @files));
        return;
    }

    return $files[0];
}

1;
