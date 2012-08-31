package Genome::Model::Tools::RefCov::ROI::Bam;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::RefCov::ROI::Bam {
    is => ['Genome::Model::Tools::RefCov::ROI::FileI'],
    has => [
        _refcov_bam => { is_optional => 1, },
        _target_names => { is_optional => 1, },
        _target_len => { is_optional => 1, },
        _tid => { is_optional => 1, },
    ],
};

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    unless ($self) { return; }
    my $refcov_bam = Genome::Model::Tools::RefCov::Bam->create(
        bam_file => $self->file,
    );
    $self->_refcov_bam($refcov_bam);

    my $header = $refcov_bam->header();
    my $target_names = $header->target_name;
    my $target_len = $header->target_len;
    $self->_target_names($target_names);
    $self->_target_len($target_len);
    $self->_tid(0);
    if ($self->load_all || $self->region_index_substring) {
        $self->_read_file;
    }
    return $self;
}

sub next_region {
    my $self = shift;
    if ($self->load_all) {
        unless ($self->_all_regions) {
            my $regions = $self->all_regions;
            $self->_all_regions($regions);
        }
        return shift(@{$self->_all_regions});
    } else {
        my $tid = $self->_tid;
        $self->_tid($tid+1);
        my $chr = $self->_target_names->[$tid];
        my $length = $self->_target_len->[$tid];
        unless (defined($chr) && defined($length)) {
            return;
        }
        my $region = $self->_parse_line($chr,$length);
        unless ($region) { return; }
        if ($self->make_objects) {
            $region = Genome::Model::Tools::RefCov::ROI::Region->create(%{$region});
            return $region;
        } else {
            $region->{length} = (($region->{end} - $region->{start}) + 1);
            $region->{id} = $region->{chrom} .':'. $region->{start} .'-'. $region->{end};
            return $region;
        }
    }
}


sub _read_file {
    my $self = shift;
    my $chr_to_length_hash_ref = $self->_refcov_bam->chr_to_length_hash_ref;
    for my $chr (keys %{$chr_to_length_hash_ref}) {
        my $length = $chr_to_length_hash_ref->{$chr};
        my $region = $self->_parse_line($chr,$length);
        $self->_add_region($region);
    }
    return 1;
}

sub _load_chromosomes {
    my $self = shift;
    my $chr_to_tid_hash_ref = $self->_refcov_bam->chr_to_tid_hash_ref;
    my @chromosomes = keys %{$chr_to_tid_hash_ref};
    $self->_chromosomes(\@chromosomes);
}

sub _parse_line {
    my $self = shift;
    my $chr = shift;
    my $end = shift;
    my $strand = '+';
    my $start = 1;

    unless (defined($chr) && defined($start) && defined($end)) {
        return;
    }
    my $wingspan = $self->wingspan;
    if ($wingspan) {
        $start -= $wingspan;
        $end += $wingspan;
    }
    my %region = (
        name => $chr .':'. $start .'-'. $end,
        chrom => $chr,
        start => $start,
        end => $end,
        strand => '+',
    );
    return \%region;
}

1;
