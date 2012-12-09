package Genome::Model::Tools::AmpliconAssembly::Collate;

use strict;
use warnings;

use Genome;

use Bio::SeqIO;
use Data::Dumper 'Dumper';

class Genome::Model::Tools::AmpliconAssembly::Collate {
    is => 'Genome::Model::Tools::AmpliconAssembly',
};
#< Helps >#
sub help_detail {
    return 'This command will get the fasta (and quality) for each amplicon, and combine them into one file.  The types of fastas (and qualities) retrieved are: '.join(', ', Genome::Model::Tools::AmpliconAssembly::Set->amplicon_fasta_types);
}

sub help_synopsis {
}

#< Command >#
sub sub_command_sort_position { 50; }

sub execute {
    my $self = shift;

    my $amplicons = $self->get_amplicons
        or return;
    
    my @amplicon_fasta_types = Genome::Model::Tools::AmpliconAssembly::Set->amplicon_fasta_types;

    $self->_open_fasta_and_qual_writers(@amplicon_fasta_types)
        or return;

    for my $amplicon ( @$amplicons ) {
        for my $type ( @amplicon_fasta_types ) {
            $self->_collate_amplicon_fasta_and_qual($amplicon, $type)
                or return;
        }
    }

    return 1;
}

#< FHs >#
sub _key_for_fasta_writer {
    return '_'.$_[0].'_fasta_writer';
}

sub _key_for_qual_writer {
    return '_'.$_[0].'_qual_writer';
}

sub _open_fasta_and_qual_writers {
    my ($self, @types) = @_;

    for my $type ( @types ) {
        my $fasta_file = $self->amplicon_assembly->fasta_file_for_type($type);
        unlink $fasta_file if -e $fasta_file;
        $self->{ _key_for_fasta_writer($type) } = Bio::SeqIO->new(
            '-file' => '>'.$fasta_file,
            '-format' => 'fasta',
        )
            or return; # bad

        my $qual_file = $self->amplicon_assembly->qual_file_for_type($type);
        unlink $qual_file if -e $qual_file;
        $self->{ _key_for_qual_writer($type) } = Bio::SeqIO->new(
            '-file' => '>'.$qual_file,
            '-format' => 'qual',
        )
            or return; # bad
    }

    return 1;
}

#< Collating the Amplicon Fastas >#
sub _collate_amplicon_fasta_and_qual {
    my ($self, $amplicon, $type) = @_;

    my $method = $self->amplicon_assembly->amplicon_bioseq_method_for_type($type);
    unless ( $method ) { # bad
        $self->error_message("Can't determine method for getting bioseqs for type ($type)");
        return;
    }
    
    my @bioseqs = $amplicon->$method
        or return 1; # ok

    for my $bioseq ( @bioseqs ) {
        $self->{ _key_for_fasta_writer($type) }->write_seq($bioseq);
        $self->{ _key_for_qual_writer($type) }->write_seq($bioseq);
    }

    return 1;
}

1;

#$HeadURL$
#$Id$
