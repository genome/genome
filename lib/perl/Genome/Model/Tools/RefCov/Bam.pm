package Genome::Model::Tools::RefCov::Bam;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::RefCov::Bam {
    has => [
        bam_file => {
            is => 'String',
            doc => 'The BAM file to load',
        },
        bai_file => {
            is_calculated => 1,
            calculate_from => ['bam_file'],
            calculate => sub {
                my $bam_file = shift;
                return $bam_file .'.bai';
            },
        }
    ],
    has_optional => {
        bio_db_bam => { },
        bio_db_index => { },
        header => { },
        _chr_to_tid_hash_ref => {},
        _chr_to_length_hash_ref => {},
    },
};


sub create {
    my $class = shift;
    my %params = @_;
    unless ($] > 5.010) {
        die "Bio::DB::Sam requires perl 5.10!";
    }
    require Bio::DB::Sam;
    my $self = $class->SUPER::create(@_);
    $self->_load;
    return $self;
}

sub _load {
    my $self = shift;
    my $bai_file = $self->bai_file;
    my $bai_file_sz = -s $bai_file;
    if ( defined $bai_file_sz and $bai_file_sz == 0 ) {
        unlink $bai_file;
    }
    my $bam  = Bio::DB::Bam->open( $self->bam_file );
    unless ($bam) {
        die('Failed to open BAM file '. $self->bam_file);
    }
    $self->bio_db_bam($bam);
    $self->header($bam->header);
    
    my $bam_file = $self->bam_file;
    my $bam_mtime = (stat($bam_file))[9];
    
    my $bai_mtime = (stat($bai_file))[9];
    if ((defined($bam_mtime) && defined($bai_mtime)) && $bam_mtime > $bai_mtime) {
        die('BAM file '. $bam_file .' appears to be newer than BAM index '. $bai_file);
    }
    # I'm not sure this symlink business ever worked... deprecated June 14th, 2012 by Jason Walker
    if (0) { 
        my @bam_symlinks;
        while (-l $bam_file) {
            push @bam_symlinks, $bam_file;
            $bam_file = readlink($bam_file);
        }
        my $bam_mtime = (stat($bam_file))[9];
        my @symlinks;
        if (-e $bai_file) {
            while (-l $bai_file) {
                push @symlinks, $bai_file;
                $bai_file = readlink($bai_file);
            }
            my $bai_mtime = (stat($bai_file))[9];
            if ((defined($bam_mtime) && defined($bai_mtime)) && $bam_mtime > $bai_mtime) {
                unless (unlink $bai_file) {
                    die('Failed to remove old bai file'. $bai_file);
                }
            }
        }
        unless (-e $bai_file) {
            Bio::DB::Bam->index_build($self->bam_file);
            if (@symlinks) {
                @symlinks = reverse(@symlinks);
                my $to_file = $bai_file;
                for my $from_file (@symlinks) {
                    unless (symlink($to_file,$from_file)) {
                        die('Failed to create symlink '. $from_file .' => '. $to_file);
                    }
                    $to_file = $from_file;
                }
            }
        }
    }
    my $index  = Bio::DB::Bam->index($self->bam_file,1);
    unless ($index) {
        die('Failed to find, load, or generate index for BAM file '. $self->bam_file);
    }
    $self->bio_db_index($index);
    return 1;
}

sub tid_for_chr {
    my $self = shift;
    my $chr = shift;
    my $target_name_index = $self->chr_to_tid_hash_ref;
    my $tid = $target_name_index->{$chr};
    unless (defined $tid) { die('Failed to get tid for chromosome '. $chr); }
    return $tid;
}

sub chr_to_tid_hash_ref {
    my $self = shift;

    unless  ($self->_chr_to_tid_hash_ref) {
        my $header = $self->header();
        my $targets = $header->n_targets();
        my $target_names = $header->target_name();
        my %target_name_index;
        my $i = 0;
        for my $target_name (@{ $target_names }) {
            $target_name_index{$target_name} = $i++;
        }
        # Make sure our index is not off
        unless ($targets == $i) {
            die 'Expected '. $targets .' targets but counted '. $i .' indices!';
        }
        $self->_chr_to_tid_hash_ref(\%target_name_index);
    }
    return $self->_chr_to_tid_hash_ref;
}

sub chr_to_length_hash_ref {
    my $self = shift;

    unless  ($self->_chr_to_length_hash_ref) {
        my $chr_to_tid = $self->chr_to_tid_hash_ref;
        my $header = $self->header;
        my $target_len = $header->target_len;
        my %chr_to_length;
        for my $chr (keys %{$chr_to_tid}) {
            my $tid = $chr_to_tid->{$chr};
            my $length = $target_len->[$tid];
            $chr_to_length{$chr} = $length;
            $self->_chr_to_length_hash_ref(\%chr_to_length);
        }
    }
    return $self->_chr_to_length_hash_ref;
}


1;
