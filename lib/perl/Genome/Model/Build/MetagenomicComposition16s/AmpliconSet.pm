package Genome::Model::Build::MetagenomicComposition16s::AmpliconSet;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::MetagenomicComposition16s::AmpliconSet {
    is => 'UR::Object',
    has => [
        name => { is => 'Text', },
        primers => { is => 'Text', is_many => 1, is_optional => 1, },
        file_base_name => { is => 'Text', },
        directory => { is => 'Text', },
        # fasta
        fasta_dir => { calculate => q| return $self->directory.'/fasta'; |, },
        processed_fasta_file => { calculate => q| return $self->_file_for('processed_fasta'); |, },
        processed_qual_file => { calculate => q| return $_[0]->_file_for('processed_qual'); |, },
        chimera_free_fasta_file => { calculate => q| return $_[0]->_file_for('chimera_free_fasta'); |, },
        chimera_free_qual_file => { calculate => q| return $_[0]->_file_for('chimera_free_qual'); |, },
        oriented_fasta_file => { calculate => q| return $_[0]->_file_for('oriented_fasta'); |, },
        oriented_qual_file => { calculate => q| return $_[0]->_file_for('oriented_qual'); |, },
        # classifier
        classifier => { is => 'Text', },
        classification_dir => { calculate => q| return $self->directory.'/classification'; |, },
        classification_file => { calculate => q| return $self->_file_for('classification'); |, },
        # chimera
        chimera_dir => {  calculate => q| return $self->directory.'/chimera'; |, }, 
        chimera_file => {  calculate => q| return $self->_file_for('chimera'); |, }, 
    ],
    has_optional => [
        oriented_qual_file => { is => 'Text', },
        _amplicon_iterator => { is => 'Code', },
    ],
};

#< Amplicons >#
sub has_amplicons {
    my $self = shift;

    # Does the iterator exist?
    return 1 if $self->{_amplicon_iterator};

    # Do the fasta/qual files exist?
    my $fasta_file = $self->processed_fasta_file;
    my $qual_file = $self->processed_qual_file;
    return 1 if -e $fasta_file and -e $qual_file;

    return;
}

sub next_amplicon {
    my $self = shift;
    my $amplicon_iterator = $self->amplicon_iterator;
    return if not $amplicon_iterator;
    return $amplicon_iterator->();
}

sub amplicon_iterator {
    my $self = shift;

    return $self->{_amplicon_iterator} if $self->{_amplicon_iterator};

    my %input = $self->amplicon_iterator_input_fasta_and_qual;
    return if not %input;

    my $reader =  Genome::Model::Tools::Sx::PhredReader->create(%input);
    if ( not  $reader ) {
        $self->error_message('Failed create phred reader');
        return;
    }

    my $classification_file = $self->classification_file;
    my ($classification_io, $classification_line);
    if ( -s $classification_file ) {
        $classification_io = eval{ Genome::Sys->open_file_for_reading($classification_file); };
        if ( not $classification_io ) {
            $self->error_message('Failed to open classification file: '.$classification_file);
            return;
        }
        $classification_line = $classification_io->getline;
        chomp $classification_line;
    }

    my $amplicon_iterator = sub{
        my $seq = $reader->read;
        return unless $seq;  #<-- HERER

        my %amplicon = (
            name => $seq->{id},
            reads => [ $seq->{id} ],
            reads_processed => 1,
            seq => $seq,
        );

        return \%amplicon if not $classification_line;

        my @classification = split(';', $classification_line); # 0 => id | 1 => ori
        if ( not defined $classification[0] ) {
            Carp::confess('Malformed classification line: '.$classification_line);
        }
        if ( $seq->{id} ne $classification[0] ) {
            return \%amplicon;
        }

        $classification_line = $classification_io->getline;
        chomp $classification_line if $classification_line;

        $amplicon{classification} = \@classification;
        return \%amplicon;
    };

    return $self->{_amplicon_iterator} = $amplicon_iterator;
}
#<>#

#< FILES/READERS/WRITERS >#
sub _file_for {
    my ($self, $type) = @_;

    Carp::confess("No type given to get fasta (qual) file") unless defined $type;
    my %types_and_props = (
        processed_fasta => [qw/ fasta_dir processed.fasta /],
        processed_qual => [qw/ fasta_dir processed.fasta.qual /],
        chimera_free_fasta => [qw/ fasta_dir chimera_free.fasta /],
        chimera_free_qual => [qw/ fasta_dir chimera_free.fasta.qual /],
        oriented_fasta => [qw/ fasta_dir oriented.fasta /],
        oriented_qual => [qw/ fasta_dir oriented.fasta.qual /],
        chimera => [qw/ chimera_dir chimera /],
        classification => [ 'classification_dir', $self->classifier ],
    );
    Carp::confess("Invalid type ($type) given to get file") unless $types_and_props{$type};
    my ($method, $ext) = @{$types_and_props{$type}};

    return sprintf(
        '%s/%s%s.%s',
        $self->$method,
        $self->file_base_name,
        ( $self->name eq '' ? '' : '.'.$self->name ),
        $ext,
    );
}

sub seq_reader_for {
    my ($self, $type) = @_;
    
    Carp::confess("No type given to get seq reader") unless defined $type;

    my $fasta_file = $self->_file_for($type.'_fasta');
    my $qual_file = $self->_file_for($type.'_qual');
    return unless -e $fasta_file and -e $qual_file; # ok

    my %params = (
        file => $fasta_file,
        qual_file => $qual_file,
    );
    my $reader =  Genome::Model::Tools::Sx::PhredReader->create(%params);
    if ( not  $reader ) {
        $self->error_message("Failed to create $type seq reader for amplicon set name (".$self->name.')');
        return;
    }

    return $reader;
}

sub seq_writer_for {
    my ($self, $type) = @_;

    Carp::confess("No type given to get fasta and qual writer") unless defined $type;

    my $fasta_file = $self->_file_for($type.'_fasta');
    my $qual_file = $self->_file_for($type.'_qual');
    unlink $fasta_file, $qual_file;

    my $writer =  Genome::Model::Tools::Sx::PhredWriter->create(
        file => $fasta_file,
        qual_file => $qual_file,
    );
    unless ( $writer ) {
        $self->error_message("Failed to create $type seq writer for amplicon set name (".$self->name.')');
        return;
    }

    return $writer;
}

sub amplicon_iterator_input_fasta_and_qual {
    my $self = shift;

    my $fasta_file = $self->chimera_free_fasta_file;
    my $qual_file = $self->chimera_free_qual_file;
    if ( not -e $fasta_file or not -e $qual_file ) {
        $fasta_file = $self->processed_fasta_file;
        $qual_file = $self->processed_qual_file;
        if ( not -e $fasta_file or not -e $qual_file ) {
            return; # no amplicons to iterate
        }
    }

    return ( file => $fasta_file, qual_file => $qual_file, );
}
#<>#

1;

