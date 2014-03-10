package Genome::Model::Tools::Soap::CreateContigsBasesFile;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Soap::CreateContigsBasesFile {
    is => 'Genome::Model::Tools::Soap::Base',
    has => [
	    assembly_directory => {
            is => 'Text',
            is_optional => 1,
            doc => 'Soap assembly directory',
        },
        input_file => {
            is => 'Text',
            is_optional => 1,
            doc => 'The input scaffold sequences file from SOAP. Deafult is the gap filled scaffold file in the edit_dir ("gapfill"), then the scaffold sequences file ("*.scafSeq") in the assembly directory.',
        },
	    output_file => {
            is => 'Text',
            is_optional => 1,
            doc => 'Output contigs file. Default is "contigs.bases" in the assembly directory',
        },
        min_contig_length => {
            is => 'Integer',
            doc => 'Minimum contig length',
        },
    ],
};

sub help_brief {
    'Create contigs fasta file from soap scaffold sequences file.';
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
                type => 'invalid',
                properties => [qw/ assembly_directory /],
                desc => 'The assembly_directory is not a directory!',
            );
            return @errors;
        }
        my $create_edit_dir = $self->create_edit_dir;
        return if not $create_edit_dir;
        my $input_file;
        for my $method (qw/ _resolve_gapfill_file _resolve_scaffold_sequence_file /) {
            $input_file = $self->$method;
            last if -s $input_file;
        }
        $self->input_file($input_file);
        $self->output_file( $self->_resolve_contigs_bases_file ) if not defined $self->output_file;
    }
    elsif ( not $self->output_file ) { 
        push @errors, UR::Object::Tag->create(
            type => 'invalid',
            properties => [qw/ output_file /],
            desc => 'No output file given and no assembly_directory given to determine the output file!',
        );
    }

    my $file = $self->input_file;
    if ( not $file ) {
        push @errors, UR::Object::Tag->create(
            type => 'invalid',
            properties => [ 'input_file' ],
            desc => "A gapfill or scaffold sequences input file is required!",
        );
        return @errors;
    }
    if ( not -s $file ) {
        push @errors, UR::Object::Tag->create(
            type => 'invalid',
            properties => [ 'input_file' ],
            desc => "File ($file) does not have any size!",
        );
        return @errors;
    }

    return @errors;
}

sub execute {
    my $self = shift;
    $self->debug_message('Create contigs bases file...');

    $self->debug_message('Input sequences file: '.$self->input_file);
    my $in = Genome::Model::Tools::Sx::PhredReader->create(file => $self->input_file);
    return if not $in;
    $self->debug_message('Output contigs bases file: '.$self->output_file);
    my $out = Genome::Model::Tools::Sx::PhredWriter->create(file => $self->output_file);
    return if not $out;

    my $supercontig_number = 0;
    while ( my $seq = $in->read ) {
        #filter out scaffolds less than min_contig_length
        next unless length $seq->{seq} >= $self->min_contig_length;
        my $contig_number = 0;
        foreach my $bases ( split (/N+/, $seq->{seq}) ) {
            #filter out contigs less than min contig length
            next unless length $bases >= $self->min_contig_length;
            my $id = 'Contig'.$supercontig_number.'.'.++$contig_number;
            $out->write({ id => $id, seq => $bases });
        }
        $supercontig_number++;
    }
    $out->close;

    $self->debug_message('Create contigs bases file...OK');
    return 1;
}


1;
