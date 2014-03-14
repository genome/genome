package Genome::InstrumentData::Command::Import::WorkFlow::SplitBamByReadGroup;

use strict;
use warnings;

use Genome;

require IO::File;
require List::MoreUtils;

class Genome::InstrumentData::Command::Import::WorkFlow::SplitBamByReadGroup { 
    is => 'Command::V2',
    has_input => [
        bam_path => {
            is => 'Text',
            doc => 'The path of the unsorted bam to sort.',
        }
    ],
    has_output => [ 
        read_group_bam_paths => {
            is => 'Text',
            is_many => 1,
            doc => 'The paths of the read group bams.',
        },
    ],
    has_optional_transient => [
        headers => { is => 'Array', },
        read_groups_and_headers => { is => 'Hash', default => {} },
        removed_reads_cnt => { is => 'Number', },
        _read_group_fhs => { is => 'HASH', default => {} },
    ],
    has_optional_calculated => [
        read_group_ids => {
            calculate_from => [qw/ read_groups_and_headers /],
            calculate => q( return keys %$read_groups_and_headers ), 
        },
    ],
};

sub execute {
    my $self = shift;
    $self->debug_message('Spilt bam by read group...');

    my $set_headers_and_read_groups = $self->_set_headers_and_read_groups;
    return if not $set_headers_and_read_groups;

    my @read_group_ids = $self->read_group_ids;

    my $write_reads_ok = $self->_write_reads();
    return if not $write_reads_ok;

    my $verify_read_count_ok = $self->_verify_read_count;
    return if not $verify_read_count_ok;

    $self->debug_message('Spilt bam by read group...done');
    return 1;
}

sub _set_headers_and_read_groups {
    my $self = shift;

    my $helpers = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get;
    my $headers = $helpers->load_headers_from_bam($self->bam_path);
    return if not $headers;
    
    my $read_group_headers = delete $headers->{'@RG'};
    my (@read_group_ids, $read_groups_and_headers);
    if ( $read_group_headers ) {
        $read_groups_and_headers = $helpers->read_groups_from_headers($read_group_headers);
        return if not $read_groups_and_headers;
    }

    $self->headers($headers);
    $self->read_groups_and_headers($read_groups_and_headers);

    return 1;
}


sub _write_reads {
    my ($self) = @_;
    $self->debug_message('Write reads...');

    my $bam_path = $self->bam_path;
    $self->debug_message("Bam path: $bam_path");
    my $bam_fh = IO::File->new("samtools view $bam_path |");
    if ( not $bam_fh ) {
        $self->error_message('Failed to open file handle to samtools command!');
        return;
    }

    my $removed_reads_cnt = 0;
    my $previous_line;
    my $previous_id;
    my $previous_read_group_id;
    while ( my $line = $bam_fh->getline ) {
        my @tokens = split(/\t/, $line);
        my $id = $tokens[0];

        $line =~ m/\sRG:Z:(.*?)\s/;
        my $read_group_id = $1;
        $read_group_id //= 'unknown';

        unless($previous_line) {
            $previous_line = $line;
            $previous_id = $id;
            $previous_read_group_id = $read_group_id;
            next;
        }

        if($id eq $previous_id and $previous_read_group_id eq $read_group_id) {
            my $fh = $self->_fh_for_read_group_and_pairedness($read_group_id, 'paired');
            $fh->print($previous_line, $line);
            undef $previous_line;
            undef $previous_id;
            undef $previous_read_group_id;
        } else {
            my $fh = $self->_fh_for_read_group_and_pairedness($previous_read_group_id, 'singleton');
            $fh->print($previous_line);

            $previous_line = $line;
            $previous_id = $id;
            $previous_read_group_id = $read_group_id;
        }
    }

    if($previous_line) {
        my $fh = $self->_fh_for_read_group_and_pairedness($previous_read_group_id, 'singleton');
        $fh->print($previous_line);
    }

    for my $fh ( $bam_fh, values %{$self->_read_group_fhs} ) {
        $fh->close;
    }

    $self->debug_message('Removed reads count: '.$removed_reads_cnt);
    $self->removed_reads_cnt($removed_reads_cnt);

    $self->debug_message('Write reads...done');
    return 1;
}

sub _verify_read_count {
    my $self = shift;
    $self->debug_message('Verify read count...');

    my $helpers = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get;

    my @read_group_bam_paths = $self->read_group_bam_paths;
    my @validated_read_group_bam_paths;
    my $read_count = 0;
    for my $read_group_bam_path ( @read_group_bam_paths ) {
        my $flagstat = $helpers->run_flagstat($read_group_bam_path);
        return if not $flagstat;
        next if $flagstat->{total_reads} == 0; # skip
        $read_count += $flagstat->{total_reads};
        push @validated_read_group_bam_paths, $read_group_bam_path;
    }

    if ( not @validated_read_group_bam_paths ) {
        $self->error_message('No read group bams passed validation!');
        return;
    }

    $self->read_group_bam_paths(\@validated_read_group_bam_paths);

    my $original_flagstat = $helpers->load_or_run_flagstat($self->bam_path);
    return if not $original_flagstat;

    $self->debug_message('Original bam read count: '.$original_flagstat->{total_reads});
    $self->debug_message('Removed reads count: '.$self->removed_reads_cnt);
    $self->debug_message('Read group bams read count: '.$read_count);

    if ( $original_flagstat->{total_reads} - $self->removed_reads_cnt != $read_count ) {
        $self->error_message('Original and read group bam read counts do not match!');
        return;
    }

    $self->debug_message('Verify read count...done');
    return 1;
}

sub _fh_for_read_group_and_pairedness {
    my ($self, $read_group, $pairedness) = @_;

    my $fhs = $self->_read_group_fhs;
    my $key = join('*', $read_group, $pairedness);
    unless(exists $fhs->{$key}) {
        $fhs->{$key} = $self->_open_fh_for_read_group_and_pairedness($read_group, $pairedness);
        $self->_read_group_fhs($fhs);
    }

    return $fhs->{$key};
}

sub _open_fh_for_read_group_and_pairedness {
    my ($self, $read_group_id, $pairedness) = @_;

    my $read_group_bam_path = $self->bam_path;
    $read_group_bam_path =~ s/\.bam$//;
    $read_group_bam_path = join('.', $read_group_bam_path, $read_group_id, $pairedness, 'bam');
    my $fh = IO::File->new("| samtools view -S -b -o $read_group_bam_path -");
    if ( not $fh ) {
        $self->error_message('Failed to open file handle to samtools command!');
        return;
    }

    $self->_write_headers_for_read_group($fh, $read_group_id);

    my @read_group_bam_paths = $self->read_group_bam_paths;
    push @read_group_bam_paths, $read_group_bam_path;
    $self->read_group_bam_paths(\@read_group_bam_paths);

    return $fh;
}

sub _write_headers_for_read_group {
    my ($self, $fh, $read_group_id) = @_;

    my $headers = $self->headers;
    Carp::confess('No headers to write to read group bams!') if not $headers;

    my $helpers = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get;
    my $headers_as_string = $helpers->headers_to_string($headers);
    return if not $headers_as_string;

    $fh->print( $headers_as_string );

    my $read_groups_and_headers = $self->read_groups_and_headers;
    unless(exists $read_groups_and_headers->{$read_group_id}) {
        $read_groups_and_headers->{$read_group_id} = join("\t", 'CN:NA',);
    }

    $fh->print(
        join("\t", '@RG', 'ID:'.$read_group_id, $read_groups_and_headers->{$read_group_id})."\n"
    );

    return 1;
}


1;

