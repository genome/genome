package Genome::Model::Tools::Fasta::TrimQuality;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Fasta::TrimQuality {
    is  => 'Genome::Model::Tools::Fasta',
    has => [	 
    min_trim_quality => {
        type => 'Integer',
        doc => 'Minimum quality value cutoff (10)',
        default => 10,
    },
    min_trim_length => {
        type => 'Integer',
        doc => 'Minimum clipped read length (100)',
        default => 100,
    },		 
    ],
};

use Data::Dumper;
require File::Copy;

sub help_brief {
    return 'Trims by FASTA and Quality files by minimum quality and length';
}

sub help_detail { 
    return <<EOS 
EOS
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_)
        or return;

    unless ( $self->have_qual_file ) {
        $self->error_message("A qual file (fasta file name w/ .qual) is required");
        return;
    }

    return $self;
}

sub execute {
    my $self = shift;

    $self->chdir_fasta_directory 
        or return;

    my $command = sprintf(
        'trim3 %s -m %s -q %s',# -x 10',
        $self->fasta_base,
        $self->min_trim_length,
        $self->min_trim_quality,
    );

    if ( system $command ) {
        ($self->error_message("trim3 failed.") and return);
    }

    # Makes a <fasta>.clip and <fasta>.clip.qua, move to file names
    # FASTA
    my $fasta_bak = sprintf('%s.preclip', $self->fasta_base);
    File::Copy::copy($self->fasta_base, $fasta_bak)
        or ($self->error_message(sprintf('Can\'t copy %s to %s: %s', $self->fasta_base, $fasta_bak, $!)) 
	    and return);
    unlink $self->fasta_base;
    my $fasta_clip = sprintf('%s.clip', $self->fasta_base);
    File::Copy::copy($fasta_clip, $self->fasta_base)
        or ($self->error_message( sprintf('Can\'t copy output file (%s) to %s: %s', $fasta_clip, $self->fasta_base, $!) )
	    and return);
    unlink $fasta_clip;

    # QUAL
    my $qual_bak = sprintf('%s.preclip', $self->qual_base);
    File::Copy::copy($self->qual_base, $qual_bak)
        or ($self->error_message( sprintf('Can\'t copy %s to %s: %s', $self->qual_base, $qual_bak, $!) ) and return);
    unlink $self->qual_base;
    my $qual_clip = sprintf('%s.clip.qual', $self->fasta_base);
    File::Copy::copy($qual_clip, $self->qual_base)
        or ($self->error_message( sprintf('Can\'t copy output qual file (%s) to %s: %s', $qual_clip, $self->qual_base, $!) ) and return);
    unlink $qual_clip;

    $self->debug_message("No sequences made the quality and length cut.") unless -s $self->fasta_base;

    $self->chdir_cwd 
        or return;

    return 1;
}

1;

#$HeadURL$
#$Id$
