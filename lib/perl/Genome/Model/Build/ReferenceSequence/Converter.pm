package Genome::Model::Build::ReferenceSequence::Converter;

use Genome;
use warnings;
use strict;
use Sys::Hostname;
use List::Util qw( first );
use List::MoreUtils qw( each_array );
use Data::Dumper;

class Genome::Model::Build::ReferenceSequence::Converter {
    is => ['Genome::SoftwareResult'],
    has => [
        source_reference_build => {
            is => 'Genome::Model::Build::ReferenceSequence',
            id_by => 'source_reference_build_id',
        },
        destination_reference_build => {
            is => 'Genome::Model::Build::ReferenceSequence',
            id_by => 'destination_reference_build_id',
        },
    ],
    has_input => [
        source_reference_build_id => {
            is => 'Number',
            doc => 'the reference to use by id',
        },
        destination_reference_build_id => {
            is => 'Number',
            doc => 'the reference to use by id',
        },
    ],
    has_metric => [
        algorithm => { # README - If adding an algorithm, please add to the list of valid algorithms
            is => 'Text',
            valid_values => [qw/ convert_chrXX_contigs_to_GL chop_chr prepend_chr lift_over no_op drop_extra_contigs chr_md5/],
            doc => 'method to use to convert from the source to the destination',
        },
        resource => {
            is => 'Text',
            doc => 'additional resource to facilitate conversion if the algorithm requires (e.g. lift_over chain file)',
            is_optional => 1,
        },
    ],
    has_transient_optional => [
        _chr_map => {
            is => 'HASH',
            doc => 'A hash ref with the source chr as the key and the destination chr as the value',
        },
    ],
};

sub exists_for_references {
    my $class = shift;
    my $source_reference = shift;
    my $destination_reference = shift;

    my $ref = $class->_faster_get(source_reference_build => $source_reference, destination_reference_build => $destination_reference);
    return defined($ref);
}

sub convert_bed {
    my $class = shift;
    my ($source_bed, $source_reference, $destination_bed, $destination_reference) = @_;

    my $self = $class->_faster_get(source_reference_build => $source_reference, destination_reference_build => $destination_reference);
    unless($self) {
        $class->error_message('Could not find converter from ' . $source_reference->__display_name__ . ' to ' . $destination_reference->__display_name__  . '. (See `genome model reference-sequence converter list` for available conversions.)');
        return;
    }

    if($self->__errors__) {
        $class->error_message('Loaded converter could not be used due to errors: ' . join(' ',map($_->__display_name__, $self->__errors__)));
        return;
    }

    if($self->is_per_position_algorithm($self->algorithm)) {
        $self->parse_and_write_bed($source_bed, $destination_bed, $self->algorithm, undef);
    } else {
        #operate on the whole BED file at once
        my $algorithm = $self->algorithm;
        $self->$algorithm($source_bed, $source_reference, $destination_bed, $destination_reference);
    }

   return $destination_bed;
}

sub is_per_position_algorithm {
    my $self = shift;
    my ($algorithm) = @_;

    return !grep($_ eq $algorithm, 'lift_over', 'no_op', 'drop_extra_contigs');
}

sub convert_position {
    my $self = shift;
    my ($algorithm, $chrom, $start, $stop) = @_;

    unless(_coords_are_valid($chrom, $start, $stop)) {
        $self->error_message('Missing one or more of chrom, start, stop. Got: (' . ($chrom || '') . ', ' . ($start || '') . ', ' . ($stop || '') . ').');
        return;
    }

    my ($new_chrom, $new_start, $new_stop) =  $self->$algorithm($chrom, $start, $stop);
    unless(_coords_are_valid($new_chrom, $new_start, $new_stop)) {
        $self->error_message('Could not convert one or more of chrom, start, stop. Got: (' . ($new_chrom || '') . ', ' . ($new_start || '') . ', ' . ($new_stop || '') . ').');
        return;
    }

    return ($new_chrom, $new_start, $new_stop);
}

sub convert_chrXX_contigs_to_GL {
    my $self = shift;
    my ($chrom, $start, $stop) = $self->chop_chr(@_);

    if($chrom =~ /\d+_(GL\d+)R/) {
        $chrom = $1 . '.1';
    } elsif ($chrom =~ /Un_gl(\d+)/i) {
        $chrom = 'GL' . $1 . '.1';
    }elsif ( $chrom =~ /\d+_GL(\d+)_random$/i) {
        $chrom = "GL" .  $1 . '.1';
    }

    return ($chrom, $start, $stop);
}

sub chop_chr {
    my $self = shift;
    my ($chrom, $start, $stop) = @_;

    $chrom =~ s/^chr//;

    return ($chrom, $start, $stop);
}

sub prepend_chr {
    my $self = shift;
    my ($chrom, $start, $stop) = @_;

    unless($chrom =~ /^chr/) {
        $chrom = 'chr' . $chrom;
    }

    return ($chrom, $start, $stop);
}

sub chr_md5 {
    my $self = shift;
    my ($chrom, $start, $stop) = @_;

    my $chr_map = $self->get_or_create_chr_map_hash_ref();
    my $new_chrom = $chr_map->{$chrom};
    unless ($new_chrom) {
        $self->fatal_message('Failed to convert '. $chrom .':'. $start .'-'. $stop .' with chromosome map: '. Data::Dumper::Dumper($chr_map));
    }
    return ($new_chrom, $start, $stop);
}

sub lift_over {
    my $self = shift;
    my ($source_bed, $source_reference, $destination_bed, $destination_reference) = @_;

    my $chain_file = $self->resource;
    unless($chain_file and -s $chain_file) {
        $self->error_message('No chain file resource found for LiftOver.');
        die $self->error_message;
    }
    my $lift_over_cmd = Genome::Model::Tools::LiftOver->create(
        source_file => $source_bed,
        chain_file => $chain_file,
        destination_file => $destination_bed,
    );
    unless($lift_over_cmd->execute()) {
        die $self->error_message('LiftOver failed: ' . $lift_over_cmd->error_message);
    }

    return $destination_bed;
}

sub no_op {
    my $self = shift;
    my ($source_bed, $source_reference, $destination_bed, $destination_reference) = @_;

    #possibly useful for recording that two reference sequences are completely equivalent except in name
    Genome::Sys->copy_file($source_bed, $destination_bed);
    return $destination_bed;
}

sub parse_and_write_bed {
    my $self = shift;
    my ($source_bed, $destination_bed, $position_algorithm, $site_test) = @_;

    my $source_fh = Genome::Sys->open_file_for_reading($source_bed);
    my $destination_fh = Genome::Sys->open_file_for_writing($destination_bed);
    
    while(my $line = <$source_fh>) {
        chomp $line;
        my ($chrom, $start, $stop, @extra) = split("\t", $line);
        unless (_coords_are_valid($chrom, $start, $stop)) {
            $self->debug_message('Not converting non-entry line %s', $line);
            $destination_fh->say($line);
            next;
        }
        my ($new_chrom, $new_start, $new_stop) = ($chrom, $start, $stop);
        if ($position_algorithm) {
            ($new_chrom, $new_start, $new_stop) = $self->convert_position($position_algorithm, $chrom, $start, $stop);
        }
        my $new_line = join("\t", $new_chrom, $new_start, $new_stop, @extra) . "\n";
        print $destination_fh $new_line if !$site_test || $site_test->($new_chrom, $new_start, $new_stop);
    }
    $source_fh->close;
    $destination_fh->close;
}

sub drop_extra_contigs {
    my $self = shift;
    my ($source_bed, $source_reference, $destination_bed, $destination_reference) = @_;

    my $chromosome_set = Set::Scalar->new(@{$self->destination_reference_build->chromosome_array_ref});
    my $algorithm = $self->resource;
    if($algorithm) {
        unless($self->validate_algorithm_name($algorithm) && $self->is_per_position_algorithm($algorithm)) {
            die "Per-position algorithm from resource property is not a valid converter\n";
        }
    }

    my $site_test = sub {
        my ($chrom) = @_;
        return $chromosome_set->has($chrom);
    };
    $self->parse_and_write_bed($source_bed, $destination_bed, $algorithm, $site_test);
}

sub validate_algorithm_name {
    my $self = shift;
    my ($algorithm_name) = @_;

    return first { $_ eq $algorithm_name } @{$self->__meta__->property_meta_for_name('algorithm')->valid_values};
}

sub _coords_are_valid {
    my ($chrom, $start, $stop) = @_;
    return (
        defined $chrom and $chrom ne q{}
            and defined $start and $start ne q{}
            and defined $stop and $stop ne q{}
            and $start <= $stop
    );
}

sub get_or_create_chr_map_hash_ref {
    my $self = shift;

    return $self->_chr_map if $self->_chr_map;
    my %chr_map;

    my $source_seqdict = $self->source_reference_build->get_or_create_seqdict_hash_ref;
    my $destination_seqdict = $self->destination_reference_build->get_or_create_seqdict_hash_ref;

    my @sorted_source_chrs = sort { $source_seqdict->{$a}{md5} cmp $source_seqdict->{$b}{md5}} keys %{$source_seqdict};
    my @sorted_destination_chrs = sort { $destination_seqdict->{$a}{md5} cmp $destination_seqdict->{$b}{md5}} keys %{$destination_seqdict};
    my $each_array = each_array(@sorted_source_chrs,@sorted_destination_chrs);
    while ( my ($source_chr,$destination_chr) = $each_array->()) {
        my $source_md5 = $source_seqdict->{$source_chr}{md5};
        my $destination_md5 = $destination_seqdict->{$destination_chr}{md5};
        if ($source_md5 eq $destination_md5) {
            $chr_map{$source_chr} = $destination_chr;
        } else {
            $self->fatal_message('Unable to convert by chr_md5 due to unequal chromosome md5: '. $source_md5 .' is not equal to '. $destination_md5 .' for chromosomes '. $source_chr .' and '. $destination_chr);
        }
    }
    $self->_chr_map(\%chr_map);
    return $self->_chr_map;
}

1;
