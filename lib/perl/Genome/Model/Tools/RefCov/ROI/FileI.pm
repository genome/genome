package Genome::Model::Tools::RefCov::ROI::FileI;

use strict;
use warnings;

use Genome;

my $DEFAULT_REGION_INDEX_SUBSTRING = 0;
my $DEFAULT_MAKE_OBJECTS = 0;
my $DEFAULT_LOAD_ALL = 0;

class Genome::Model::Tools::RefCov::ROI::FileI {
    has => [
        file => {
            is => 'String',
            doc => 'The file path of the defined regions/intervals',
        },
    ],
    has_optional => {
        make_objects => {
            is => 'Boolean',
            default_value => $DEFAULT_MAKE_OBJECTS,
        },
        load_all => {
            is => 'Boolean',
            default_value => $DEFAULT_LOAD_ALL,
        },
        region_index_substring => {
            is => 'Integer',
            default_value => $DEFAULT_REGION_INDEX_SUBSTRING,
        },
        wingspan => {
            is => 'Integer',
            doc => 'An integer distance to add to each end of a region.',
        },
        _fh => { },
        _chromosomes => { },
        _all_regions => { },
    }
};

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    unless ($self) { return; }
    if ($class eq 'Genome::Model::Tools::RefCov::ROI::Bam') { return $self; }
    my $fh = IO::File->new($self->file,'r');
    unless ($fh) { die('Failed to load file: '. $self->file); }
    $self->_fh($fh);
    if ($self->load_all || $self->region_index_substring) {
        $self->_read_file;
    }
    return $self;
}

sub next_region {
    my $self = shift;
    if ($self->load_all || $self->region_index_substring ) {
        unless ($self->_all_regions) {
            my $regions = $self->all_regions;
            $self->_all_regions($regions);
        }
        return shift(@{$self->_all_regions});
    } else {
        my $line = $self->_fh->getline;
        unless ($line) { $self->_fh->close; return; }
        my $region = $self->_parse_line($line);
        unless ($region) { next; }
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

sub chromosomes {
    my $self = shift;
    unless ($self->_chromosomes) {
        $self->_load_chromosomes;
    }
    return @{$self->_chromosomes};
}

sub overlaps_regions {
    my $self = shift;
    my ($chr,$start,$stop) = @_;
    my $start_substr = substr($start, 0, $self->region_index_substring) || 0;
    my $stop_substr = substr($stop, 0, $self->region_index_substring) || 0;

    my %overlapping_regions;
    for (my $position_key = $start_substr; $position_key <= $stop_substr; $position_key++) {
        my $region_key = $chr .':'. $position_key;
        if ($self->{indexed_regions}->{$region_key}) {
            my @region_list = split(/\n/, $self->{indexed_regions}->{$region_key});
            foreach my $region (@region_list) {
                (my $region_start, my $region_stop) = split(/\t/, $region);
                # This determines overlap, and it identifies reads that contain a region
                if ( ($start >= $region_start && $start <= $region_stop) ||
                         ($stop >= $region_start && $stop <= $region_stop) ||
                             # A read that spans the region
                             ($start < $region_start && $stop > $region_stop)
                         ) {
                    $overlapping_regions{$region} = 1;
                }
            }
        }
    }
    my @overlapping_regions = keys %overlapping_regions;
    if (@overlapping_regions) {
        return \@overlapping_regions;
    }
    return 0;
}

sub chromosome_regions {
    my $self = shift;
    my $chrom = shift;
    unless ($self->{chrom_regions}->{$chrom}) {
        return;
    }
    return @{$self->{chrom_regions}->{$chrom}};
}

sub all_regions {
    my $self = shift;
    my @chromosomes = $self->chromosomes;
    my @regions;
    for my $chrom (@chromosomes) {
        push @regions, $self->chromosome_regions($chrom);
    }
    return \@regions;
}

sub _add_region {
    my $self = shift;
    my $region = shift;
    unless ($region && ref($region) eq 'HASH') {
        die ('Must supply a HASH reference to method _add_region');
    }
    unless (defined($region->{start}) && defined($region->{end}) && defined($region->{chrom})) {
        die('HASH ref must contain start, end, and chrom attributes!');
    }
    my $chrom = $region->{chrom};
    my $start = $region->{start};
    my $stop = $region->{end};
    my $start_substr = substr($start, 0, $self->region_index_substring) || 0;
    my $stop_substr = substr($stop, 0, $self->region_index_substring) || 0;
    for (my $position_key = $start_substr; $position_key <= $stop_substr; $position_key++) {
        my $region_key = $chrom . ':' . $position_key;
        $self->{indexed_regions}->{$region_key} .=  $start."\t". $stop ."\n";
    }
    if ($self->make_objects) {
        # This is only necessary to perform arithmetic operations
        $region = Genome::Model::Tools::RefCov::ROI::Region->create(%{$region});
    }
    push @{$self->{chrom_regions}->{$chrom}}, $region;
    return 1;
}

sub _read_file {
    my $self = shift;
    my $fh = $self->_fh;
    while (my $line = $fh->getline) {
        my $region = $self->_parse_line($line);
        unless ($region) { next; }
        $self->_add_region($region);
    }
    $fh->close;
    return 1;
}

sub _parse_line {
    die ('_parse_line is an abstract method.  Please implement in '. __PACKAGE__);
}
sub _load_chromosomes {
die ('_load_chromosomes is an abstract method.  Please implement in '. __PACKAGE__);
}

1;
