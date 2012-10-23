package Genome::Model::Tools::Fastq::ToPhdballScf;

use strict;
use warnings;

use Genome;
use Bio::SeqIO;
use Bio::Seq::Quality;


class Genome::Model::Tools::Fastq::ToPhdballScf {
    is           => 'Genome::Model::Tools::Fastq',
    has_optional => [
        time    => {
            is  => 'String',
            doc => 'time stamp inside phd file, often need sync with timestamp in acefile',
        },
        scf_dir => {
            is  => 'String',
            doc => 'The path to place scf trace files. This is also a flag to make scf traces',
        },
        ball_file => {
            is  => 'String',
            doc => 'phd ball file name with path',
            default => './phd.ball',
        },
        id_range => {
            is  => 'String',
            doc => 'useful for split huge fastq file(velvet) into parts, must be something like 200-3999',
        },
        base_fix => {
            is  => 'Boolean',
            doc => 'fix fastq bases, change non-ATCG bases to A (solexa way)',
        },
    ],
};


sub help_brief {
    "convert from fastq files to phd ball and scf files"
}


sub help_detail {                           
    return <<EOS 

EOS
}


sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);

    if ($self->id_range) {
        if ($self->id_range =~ /^(\d+)\-(\d+)$/) {
            unless ($1 < $2) {
                $self->error_message("Small number first in id_range: ".$self->id_range);
                return;
            }
        }
        else {
            $self->error_message("id range format must be like: \"number-number\" ".$self->id_range);
            return;
        }
    }       
    
    my $ball_file = $self->ball_file;
    if (-s $ball_file) {
        $self->warning_message("ball file $ball_file already exists, Now remove");
        unlink $ball_file;
        if (-s $ball_file) {
            $self->error_message("Failed to delete existing $ball_file");
            return;
        }
    }

    if ($self->scf_dir) {
        my $scf_dir = $self->scf_dir;

        unless (-d $scf_dir) {
            mkdir $scf_dir, 0777;
            unless (-d $scf_dir) {
                $self->error_message("Fail to make dir: $scf_dir");
                return;
            }
        }
    }
    return $self;
}


sub execute {
    my $self  = shift;
    
    my $id_range = $self->id_range;
    my $fq_io    = $self->get_fastq_reader($self->fastq_file);       
    my $ball_io  = $self->get_ball_appender($self->ball_file, 'phd');

    my ($new_id, $max_id) = split /-/, $id_range if $id_range;

    while (my $fq = $fq_io->next_seq) {
        my $seq = $fq->seq;
        $seq =~ s/[^AaTtCcGg]/A/g if $self->base_fix;
        
        my $id   = $id_range ? $new_id : $fq->id;
        my $qual = $fq->qual;
        map{$_ -= 31}@$qual if $self->solexa_fastq; #solexa fastq qual = Qphred + 64. Bio::SeqIO fastq quality = sanger fastq quality = Qphred + 33
        
        my %params = (
            -seq  => $seq,
            -qual => $qual,
            -id   => $id,
            -force_flush => 1,
            -trace       => [map{$_*10}(0..$fq->length-1)], 
        );
        
        my $phd_swq = Bio::Seq::Quality->new(%params);  
        $phd_swq->chromat_file($id);
        $phd_swq->time($self->time) if $self->time;
        $ball_io->write_seq($phd_swq);  
        
        if ($self->scf_dir) {
            my $scf_file = $self->scf_dir.'/'.$id;
            my $gz_file  = $scf_file.'.gz';
        
            if (-s $gz_file) {
                $self->warning_message("$gz_file already exist, now delete it");
                unlink $gz_file;
            }
                
            my $scf_swq = Bio::Seq::Quality->new(%params);
            my $scf_io  = $self->get_scf_writer($scf_file, 'scf');
        
            $scf_io->write_seq(-target => Bio::Seq::SequenceTrace->new(-swq => $scf_swq));
            `gzip $scf_file 2>&1`;
        }
        $new_id++ if $id_range;
    }
    
    if ($id_range) {
        $self->warning_message("max_id $max_id in id_range not equal to last assigned id $new_id")
            unless $new_id - 1 == $max_id;
    }

    return 1;
}


sub get_ball_appender {
    return shift->_get_bioseq(@_, '>>');
}


sub get_scf_writer {
    return shift->_get_bioseq(@_, '>');
}

1;
