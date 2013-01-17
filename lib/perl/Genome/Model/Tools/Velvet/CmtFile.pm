package Genome::Model::Tools::Velvet::CmtFile;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Velvet::CmtFile {
    is => 'Genome::Model::Tools::Velvet::Base',
    has => [
        sequencing_technologies => {
            is => 'Text',
            is_many => 1,
            doc => 'Technology used to sequence data',
        },
        current_finishing_status => {
            is => 'Text',
            doc => 'Finishing status of this assembly project',
            is_optional => 1,
            default => 'High Quality Draft',
        },
        finishing_goal => {
            is => 'Text',
            doc => 'Assembly finishing goal',
            is_optional => 1,
            default => 'High Quality Draft',
        },
        assembly_directory => {
            is => 'Text',
            doc => 'Main assembly directory .. above edit_dir',
        },
        min_contig_length => {
            is => 'Number',
            doc => 'Minimum contig length to process',
            is_optional => 1,
        }
    ],
};

sub help_brief {
    'Tool to create contigs.cmt file used by submissions group';
}

sub help_detail {
    return <<EOS
    gmt velvet cmt-file --version 1.1.04 --assembly-directory /gscmnt/101/e_coli --min-contig-length 200
EOS
}

sub execute {
    my $self = shift;

    my $sequencing_technology;
    if( not $sequencing_technology = $self->_resolve_sequencing_technology ) {
        $self->error_message('Failed to resolve sequencing technology');
        return;
    }

    my $coverage;
    if( not $coverage = $self->_calculate_genome_coverage ) {
        $self->error_message('Failed to calculate coverage');
        return;
    }

    my $assembler_version = ( $self->version )
        ? $self->version
        : 'Unknown' ;

    unlink $self->contigs_cmt_file;
    my $out_fh = Genome::Sys->open_file_for_writing( $self->contigs_cmt_file );

    my $output = "StructuredCommentPrefix\t##Genome-Assembly-Data-START##\n".
                 "Finishing Goal\t".$self->finishing_goal."\n".
                 "Current Finishing Status\t".$self->current_finishing_status."\n".
                 "Assembly Method\tVelvet v. ".$assembler_version."\n".
                 "Genome Coverage\t".$coverage."x\n".
                 "Sequencing Technology\t". $sequencing_technology."\n".
                 "StructuredCommentSuffix\t##Genome-Assembly-Data-END##\n";
    $out_fh->print( $output );
    $out_fh->close;

    return 1;
}

sub _resolve_sequencing_technology {
    my $self = shift;

    my @valid_technologies = $self->valid_sequencing_technologies;
    my @used_technologies;
    for my $technology ( $self->sequencing_technologies ) {
        $technology = 'unknown' if not grep { /$technology/i } @valid_technologies;
        $technology = 'illumina' if lc $technology eq 'solexa';
        push @used_technologies, $technology if not grep { /$technology/i} @used_technologies;
    }
    
    return join(',', map {ucfirst $_} sort @used_technologies );
}

sub valid_sequencing_technologies {
    return qw/ sanger solexa illumina 454 /;
}

sub _calculate_genome_coverage {
    my $self = shift;

    my $genome_size = $self->_assembly_genome_size;
    if ( not $genome_size ) {
        $self->error_message('Could not calculate assembly genome size');
        return;
    }

    my $input_bases = $self->_input_bases;
    if ( not $input_bases ) {
        $self->error_message('Could not calculate input bases count');
        return;
    }

    return int( $input_bases / $genome_size );
}

sub _assembly_genome_size {
    my $self = shift;

    return if not -s $self->contigs_bases_file;

    my $genome_size;
    my $reader = Genome::Model::Tools::Sx::PhredReader->create(
        file => $self->contigs_bases_file,
    );
    while( my $read = $reader->read ) {
        if ( $self->min_contig_length ) {
            next if length $read->{seq} < $self->min_contig_length;
        }
        $genome_size += length $read->{seq};
    }

    return $genome_size;
}

sub _input_bases {
    my $self = shift;

    return if not -s $self->input_collated_fastq_file;

    my $input_bases;
    my $reader = Genome::Model::Tools::Sx::FastqReader->create(
        file => $self->input_collated_fastq_file,
    );
    while( my $read = $reader->read ) {
        $input_bases += length $read->{seq};
    }
    
    return $input_bases;
}

1;
