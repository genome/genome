package Genome::Model::Tools::Sam::StatsForPos;

use strict;
use warnings;

use Genome;
use Command;

class Genome::Model::Tools::Sam::StatsForPos {
    is => 'Genome::Model::Tools::Sam',
    has => [
    bam_file => {
        doc => 'The input bam file WHICH MUST BE INDEXED',
        is => 'String',
    },
    position1 => {
        doc => 'The position of interest',
        is => 'String',
    },
    position2 => {
        doc => 'The position of interest',
        is => 'String',
        is_optional=>1,
    },

    chromosome => {
        doc => 'The position of interest',
        is => 'String',
    },
    ],
};

sub help_brief {
    'fill this out later';
}

sub help_detail {
    return "fill this out later";
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);



    return $self;
}

sub execute {
    my $self = shift;

    unless (-s $self->bam_file) {
        $self->error_message('Map file '. $self->map_file .' not found or has zero size.');
        return;
    }
    my $bam_file = $self->bam_file;
    my $bai_file = $bam_file .'.bai';
    if (-e $bai_file) {
        my $bam_mtime = (stat($bam_file))[9];
        my $bai_mtime = (stat($bai_file))[9];
        if ($bam_mtime > $bai_mtime) {
            unless (unlink $bai_file) {
                die('Failed to remove old bai file'. $bai_file);
            }
        }
    }
    unless (-e $bai_file) {
        my $index = Genome::Model::Tools::Sam::IndexBam->execute(bam_file => $bam_file);
        unless ($index->result) {
            print "Failed to index $bam_file\n" and return;
        }
    }

    unless ($self->position1) {
        $self->error_message('SUPPLY A POSITION!');
        return;
    }
    unless ($self->chromosome) {
        $self->error_message('SUPPLY A CHROMOSOME!');
        return;
    }
    unless($self->position2) {
        $self->position2($self->position1);
    }
    my $cmd = "samtools view " . $self->bam_file . " " . $self->chromosome . ":" . $self->position1 . "-" . $self->position2;
    print $cmd;
    return 1;
}
1;
