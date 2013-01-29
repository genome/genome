package Genome::Model::Build::MetagenomicComposition16s::AmpliconSet;

use strict;
use warnings;

use Genome;

require File::Basename;

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
        processed_fastq_file => { calculate => q| return $self->_file_for('processed_fastq'); |, },
        chimera_free_fasta_file => { calculate => q| return $_[0]->_file_for('chimera_free_fasta'); |, },
        chimera_free_qual_file => { calculate => q| return $_[0]->_file_for('chimera_free_qual'); |, },
        oriented_fasta_file => { calculate => q| return $_[0]->_file_for('oriented_fasta'); |, },
        oriented_qual_file => { calculate => q| return $_[0]->_file_for('oriented_qual'); |, },
        # classifier
        classifier => { is => 'Text', },
        classification_dir => { calculate => q| return $self->directory.'/classification'; |, },
        classification_file => { calculate => q| return $self->_file_for('classification'); |, },
        chimera_free_classification_file => { calculate => q| return $self->_file_for('chimera_free_classification'); |, },
        # chimera
        chimera_detector => { is => 'Text', is_optional => 1, },
        chimera_dir => {  calculate => q| return $self->directory.'/chimera'; |, }, 
        chimera_file => {  calculate => q| return $self->_file_for('chimera'); |, }, 
    ],
    has_optional => [
        oriented_qual_file => { is => 'Text', },
        _amplicon_iterator => { is => 'Code', },
    ],
};

#< Amplicons >#
sub next_amplicon {
    my $self = shift;
    if ( not $self->_amplicon_iterator ) {
        $self->_amplicon_iterator( $self->amplicon_iterator );
    }
    return $self->_amplicon_iterator->();
}

sub amplicon_iterator {
    # This iterator will ALWAYS include chimeric amplicons!
    my $self = shift;

    # Sequence
    my $fasta_file = $self->oriented_fasta_file;
    my $qual_file = $self->oriented_qual_file;
    if ( not -e $fasta_file or not -e $qual_file ) {
        $fasta_file = $self->processed_fasta_file;
        $qual_file = $self->processed_qual_file;
        if ( not -e $fasta_file or not -e $qual_file ) {
            return;
        }
    }

    my $reader =  Genome::Model::Tools::Sx::PhredReader->create( 
        file => $fasta_file,
        qual_file => $qual_file, 
    );
    if ( not  $reader ) {
        $self->error_message('Failed create phred reader');
        return;
    }

    # Classifcation
    #  Try chimera free classification first, then classification file
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

    # Chimera
    my $chimera_file = $self->chimera_file;
    my ($chimera_reader, $chimera);
    if ( -s $chimera_file ) {
        $chimera_reader = $self->chimera_detector_reader_class->create(input => $chimera_file);
        if ( not $chimera_reader ) {
            $self->error_message('Failed to get chimera reader!');
            return;
        }
        $chimera = $chimera_reader->read;
    }

    # Iterator
    my $amplicon_iterator = sub{
        my $seq = $reader->read;
        return unless $seq;

        my %amplicon = (
            id => $seq->{id},
            name => $seq->{id},
            seq => $seq,
        );

        if ( $classification_line ) {
            my @classification = split(';', $classification_line); # 0 => id | 1 => ori
            if ( not defined $classification[0] ) {
                Carp::confess('Malformed classification line: '.$classification_line);
            }
            if ( $seq->{id} eq $classification[0] ) {
                $amplicon{classification_line} = $classification_line;
                $amplicon{classification} = \@classification;

                $classification_line = $classification_io->getline;
                chomp $classification_line if $classification_line;
            }
        }

        if ( $chimera and $amplicon{id} eq $chimera->{id} ) {
            $amplicon{chimera_result} = $chimera;
            $chimera = $chimera_reader->read;
        }

        return \%amplicon;
    };

    return $amplicon_iterator;
}
#<>#

#< FILES/READERS/WRITERS >#
sub base_name_for {
    my ($self, $type) = @_;

    my $file = $self->_file_for($type);
    my $base_name = File::Basename::basename($file);
    Carp::confess("Failed to get base name for file! $file") if not $base_name;

    return $base_name;
}

sub _file_for {
    my ($self, $type) = @_;

    Carp::confess("No type given to get fasta (qual) file") unless defined $type;
    my %types_and_props = (
        processed_fasta => [qw/ fasta_dir processed.fasta /],
        processed_qual => [qw/ fasta_dir processed.fasta.qual /],
        processed_fastq => [qw/ fasta_dir processed.fastq /],
        chimera_free_fasta => [qw/ fasta_dir oriented.chimera_free.fasta /],
        chimera_free_qual => [qw/ fasta_dir oriented.chimera_free.fasta.qual /],
        oriented_fasta => [qw/ fasta_dir oriented.fasta /],
        oriented_qual => [qw/ fasta_dir oriented.fasta.qual /],
        chimera => [qw/ chimera_dir chimera /],
        classification => [ 'classification_dir', $self->classifier ],
        chimera_free_classification => [ 'classification_dir', 'chimera_free.'.$self->classifier ],
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

sub chimera_detector_reader_class {
    my $self = shift;
    return 'Genome::Model::Tools::'.Genome::Utility::Text::string_to_camel_case(join(' ', split('-', $self->chimera_detector))).'::Reader';

}
#<>#

1;

