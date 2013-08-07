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
        read_groups_and_headers => { is => 'Hash', },
        removed_reads_cnt => { is => 'Number', },
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
    $self->status_message('Spilt bam by read group...');

    my $set_headers_and_read_groups = $self->_set_headers_and_read_groups;
    return if not $set_headers_and_read_groups;

    my @read_group_ids = $self->read_group_ids;
    if ( not @read_group_ids or @read_group_ids == 1 ) {
        $self->status_message('Spilting bam by read group is NOT necessary. There is only one read group [or none] in headers.');
        $self->read_group_bam_paths([ $self->bam_path ]);
        return 1;
    }

    my $read_group_fhs = $self->_open_file_handles_for_each_read_group_bam(@read_group_ids);
    return if not $read_group_fhs;

    my $write_headers_ok = $self->_write_headers_to_read_group_bams($read_group_fhs);
    return if not $write_headers_ok;

    my $write_reads_ok = $self->_write_reads($read_group_fhs);
    return if not $write_reads_ok;

    my $verify_read_count_ok = $self->_verify_read_count;
    return if not $verify_read_count_ok;

    $self->status_message('Spilt bam by read group...done');
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

sub _open_file_handles_for_each_read_group_bam {
    my ($self, @read_group_ids) = @_;
    $self->status_message('Open file handle for each read group bam...');

    Carp::confess('No read group ids to open bams!') if not @read_group_ids;

    $self->status_message('Read group count: '.@read_group_ids);
    my (%read_group_fhs, @read_group_bam_paths);
    for my $read_group_id ( @read_group_ids ){
        my $read_group_bam_path = $self->bam_path;
        $read_group_bam_path =~ s/\.bam$//;
        $read_group_bam_path .= '.'.$read_group_id.'.bam';
        push @read_group_bam_paths, $read_group_bam_path;
        my $fh = IO::File->new("| samtools view -S -b -o $read_group_bam_path -");
        if ( not $fh ) {
            $self->error_message('Failed to open file handle to samtools command!');
            return;
        }
        $read_group_fhs{$read_group_id} = $fh;
    }
    $self->read_group_bam_paths(\@read_group_bam_paths);

    $self->status_message('Open file handle for each read group bam...done');
    return \%read_group_fhs;
}

sub _write_headers_to_read_group_bams {
    my ($self, $read_group_fhs) = @_;

    Carp::confess('No read group fhs to write headers to read group bams!') if not $read_group_fhs;

    my $headers = $self->headers;
    Carp::confess('No headers to write to read group bams!') if not $headers;

    my $read_groups_and_headers = $self->read_groups_and_headers;
    Carp::confess('No read groups and headers to write headers to read group bams!') if not $read_groups_and_headers;

    my $helpers = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get;
    my $headers_as_string = $helpers->headers_to_string($headers);
    return if not $headers_as_string;

    for my $read_group_id ( keys %$read_group_fhs ) {
        $read_group_fhs->{$read_group_id}->print( $headers_as_string );
        $read_group_fhs->{$read_group_id}->print(
            join("\t", '@RG', 'ID:'.$read_group_id, $read_groups_and_headers->{$read_group_id})."\n"
        );
    }

    return 1;
}

sub _write_reads {
    my ($self, $read_group_fhs) = @_;
    $self->status_message('Write reads...');

    my $bam_path = $self->bam_path;
    $self->status_message("Bam path: $bam_path");
    my $bam_fh = IO::File->new("samtools view $bam_path |");
    if ( not $bam_fh ) {
        $self->error_message('Failed to open file handle to samtools command!');
        return;
    }

    my $removed_reads_cnt = 0;
    while ( my $line = $bam_fh->getline ) {
        my @tokens = split(/\t/, $line);
        print "$tokens[0] $tokens[1]\n";
        if ( $tokens[1] & 0x100 ) { # secondary alignment
            $removed_reads_cnt++;
            next;
        }
        $line =~ m/\sRG:Z:(.*?)\s/;
        my $read_group_id = $1;
        $read_group_id //= 'unknown';
        $read_group_fhs->{$read_group_id}->print($line);
    }

    for my $fh ( $bam_fh, values %$read_group_fhs ) {
        $fh->close;
    }

    $self->status_message('Removed reads count: '.$removed_reads_cnt);
    $self->removed_reads_cnt($removed_reads_cnt);

    $self->status_message('Write reads...done');
    return 1;
}

sub _verify_read_count {
    my $self = shift;
    $self->status_message('Verify read count...');

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

    $self->status_message('Original bam read count: '.$original_flagstat->{total_reads});
    $self->status_message('Removed reads count: '.$self->removed_reads_cnt);
    $self->status_message('Read group bams read count: '.$read_count);

    if ( $original_flagstat->{total_reads} - $self->removed_reads_cnt != $read_count ) {
        $self->error_message('Original and read group bam read counts do not match!');
        return;
    }

    $self->status_message('Verify read count...done');
    return 1;
}


1;

