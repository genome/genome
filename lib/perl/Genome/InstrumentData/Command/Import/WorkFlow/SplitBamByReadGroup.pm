package Genome::InstrumentData::Command::Import::WorkFlow::SplitBamByReadGroup;

use strict;
use warnings;

use Genome;

require IO::File;
require List::MoreUtils;
use Params::Validate ':types';

class Genome::InstrumentData::Command::Import::WorkFlow::SplitBamByReadGroup { 
    is => 'Command::V2',
    roles => [qw/
        Genome::InstrumentData::Command::Import::WorkFlow::Role::WithWorkingDirectory
        Genome::InstrumentData::Command::Import::WorkFlow::Role::RemovesInputFiles
    /],
    has_input => [
        bam_path => {
            is => 'FilePath',
            doc => 'The path of the unsorted bam to sort.',
        }
    ],
    has_output => [ 
        output_bam_paths => {
            is => 'FilePath',
            is_many => 1,
            doc => 'The paths of the read group bams.',
        },
    ],
    has_optional_transient => [
        headers => { is => 'Array', },
        read_groups_and_tags => { is => 'HASH', default => {}, },
        old_and_new_read_group_ids => { is => 'HASH', default => {}, },
        _read_group_fhs => { is => 'HASH', default => {} },
    ],
};

sub execute {
    my $self = shift;
    $self->debug_message('Spilt bam by read group...');

    my $set_headers_and_read_groups = $self->_set_headers_and_read_groups;

    my $write_reads_ok = $self->_process_reads;
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
    die $self->error_message('Failed to get headers!') if not $headers;
    $self->headers($headers);
    
    my $read_group_headers = delete $headers->{'@RG'} || [];
    my $read_groups_and_tags = $helpers->read_groups_and_tags_from_headers($read_group_headers);
    die $self->error_message('Failed to get read groups and tags from headers!') if not $read_groups_and_tags;
    $self->read_groups_and_tags($read_groups_and_tags);

    return 1;
}

sub _process_reads {
    my ($self) = @_;
    $self->debug_message('Processing reads...');

    my $bam_path = $self->bam_path;
    $self->debug_message("Bam path: $bam_path");
    my $bam_fh = IO::File->new("samtools view $bam_path |");
    if ( not $bam_fh ) {
        die $self->error_message('Failed to open file handle to samtools command!');
    }

    my $previous_tokens;
    my $previous_read_group_id;
    while ( my $line = $bam_fh->getline ) {
        chomp $line;
        my @tokens = split(/\t/, $line);
        my $read_group_id = $self->_get_rg_id_from_sam_tokens(\@tokens);

        unless($previous_tokens) {
            $previous_tokens = \@tokens;
            $previous_read_group_id = $read_group_id;
            next;
        }

        if($tokens[0] eq $previous_tokens->[0] and $previous_read_group_id eq $read_group_id) {
            $self->_write_reads_based_on_read_group_and_type(
                rg_id => $read_group_id,
                type => 'paired',
                reads => [ $previous_tokens, \@tokens ],
            );
            undef $previous_tokens;
            undef $previous_read_group_id;
        } else {
            $self->_write_reads_based_on_read_group_and_type(
                rg_id => $previous_read_group_id,
                type => _read1_or_read2($previous_tokens->[1]),
                reads => [ $previous_tokens, ],
            );
            $previous_tokens = \@tokens;
            $previous_read_group_id = $read_group_id;
        }
    }

    if ($previous_tokens) {
        $self->_write_reads_based_on_read_group_and_type(
            rg_id => $previous_read_group_id,
            type => _read1_or_read2($previous_tokens->[1]),
            reads => [ $previous_tokens, ],
        );
    }

    for my $fh ( $bam_fh, values %{$self->_read_group_fhs} ) {
        $fh->close;
    }

    $self->debug_message('Processing reads...done');
    return 1;
}

sub _get_rg_id_from_sam_tokens {
    my ($self, $tokens) = @_;

    my $rg_tag_idx = List::MoreUtils::firstidx { $_ =~ m/^RG:/ } @$tokens;
    my $rg_tag;
    if ( defined $rg_tag_idx ) {
        return (split(':', $tokens->[$rg_tag_idx]))[2];
    }
    else {
        return 'unknown';
    }
}

sub _read1_or_read2 {
    $_[0] & 0x80 ? 'read2' : 'read1'; # unlabeled reads will be read1
}

sub _verify_read_count {
    my $self = shift;
    $self->debug_message('Verify read count...');

    my $helpers = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get;

    my @output_bam_paths = $self->output_bam_paths;
    my @validated_output_bam_paths;
    my $read_count = 0;
    for my $read_group_bam_path ( @output_bam_paths ) {
        my $flagstat = $helpers->run_flagstat($read_group_bam_path);
        return if not $flagstat;
        next if $flagstat->{total_reads} == 0; # skip
        $read_count += $flagstat->{total_reads};
        push @validated_output_bam_paths, $read_group_bam_path;
    }

    if ( not @validated_output_bam_paths ) {
        die $self->error_message('No read group bams passed validation!');
    }

    $self->output_bam_paths(\@validated_output_bam_paths);

    my $original_flagstat = $helpers->load_or_run_flagstat($self->bam_path);
    return if not $original_flagstat;

    $self->debug_message('Original bam read count: '.$original_flagstat->{total_reads});
    $self->debug_message('Read group bams read count: '.$read_count);

    if ( $original_flagstat->{total_reads} != $read_count ) {
        die $self->error_message('Original and split by read group bam read counts do not match!');
    }

    $self->debug_message('Verify read count...done');
    return 1;
}

sub _write_reads_based_on_read_group_and_type {
    my ($self, %params) = @_;

    my $fhs = $self->_read_group_fhs;
    my $key = join('*', $params{rg_id}, $params{type});
    unless(exists $fhs->{$key}) {
        $fhs->{$key} = $self->_open_fh_for_read_group_and_type($params{rg_id}, $params{type});
        $self->_read_group_fhs($fhs);
    }

    for my $read_tokens ( @{$params{reads}} ) {
        $self->_update_read_group_for_sam_tokens_based_on_type($read_tokens, $params{type});
        $fhs->{$key}->print( join("\t", @$read_tokens)."\n" );
    }

    return 1;
}

sub _update_read_group_for_sam_tokens_based_on_type {
    my ($self, $tokens, $type) = @_;

    my $rg_tag_idx = List::MoreUtils::firstidx { $_ =~ m/^RG:/ } @$tokens;
    my $rg_tag;
    if ( defined $rg_tag_idx ) {
        $rg_tag = splice(@$tokens, $rg_tag_idx, 1);
    }
    else {
        $rg_tag = 'RG:Z:unknown';
        $rg_tag_idx = $#$tokens;
    }

    my $rg_id = (split(':', $rg_tag))[2];
    my $new_rg_id = $self->old_and_new_read_group_ids->{$rg_id}->{$type};
    $rg_tag = join(':', 'RG',  'Z', $new_rg_id);
    splice(@$tokens, $rg_tag_idx, 0, $rg_tag); # Add RG tag back

    return 1;
}

sub _open_fh_for_read_group_and_type {
    my ($self, $read_group_id, $type) = @_;

    my $read_group_bam_path = $self->get_working_bam_path_with_new_extension($self->bam_path, $read_group_id, $type);
    my $samtools_cmd = "| samtools view -S -b -o $read_group_bam_path -";
    $self->debug_message("Opening fh for $read_group_bam_path $type with:\n$samtools_cmd");
    my $fh = IO::File->new($samtools_cmd);
    if ( not $fh ) {
        die $self->error_message('Failed to open file handle to samtools command!');
    }

    $self->_write_headers_for_read_group($fh, $read_group_id, $type);

    my @output_bam_paths = $self->output_bam_paths;
    push @output_bam_paths, $read_group_bam_path;
    $self->output_bam_paths(\@output_bam_paths);

    return $fh;
}

sub _write_headers_for_read_group {
    my ($self, $fh, $rg_id, $type) = @_;

    # Add mapping for old RG id to new RG UUID
    my $old_and_new_read_group_ids = $self->old_and_new_read_group_ids;
    if ( not exists $old_and_new_read_group_ids->{$rg_id} ) {
        $old_and_new_read_group_ids->{$rg_id} = {
            paired => UR::Object::Type->autogenerate_new_object_id_uuid,
            read1 => UR::Object::Type->autogenerate_new_object_id_uuid,
            read2 => UR::Object::Type->autogenerate_new_object_id_uuid,
        };
    }

    my $helpers = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get;
    my $headers_as_string = $helpers->headers_to_string( $self->headers );
    return if not $headers_as_string;
    $fh->print( $headers_as_string );
    $fh->print( $self->_header_for_read_group_and_type($rg_id, $type) );

    return 1;
}

sub _header_for_read_group_and_type {
    my ($self, $rg_id, $type) = @_;

    my $read_groups_and_tags = $self->read_groups_and_tags;
    if ( not exists $read_groups_and_tags->{$rg_id} ) {
        # Add RG tags for groups that are not in the header. This includes the 'unknown' group
        $read_groups_and_tags->{$rg_id} = { CN => 'NA', };
    }
    delete $read_groups_and_tags->{$rg_id}->{ID} if exists $read_groups_and_tags->{$rg_id}->{ID};
    my %rg_tags = %{$read_groups_and_tags->{$rg_id}};
    my @tag_names = sort keys %rg_tags;

    return join(
        "\t", '@RG',
        'ID:'.$self->old_and_new_read_group_ids->{$rg_id}->{$type},
        map { join(':', $_, $rg_tags{$_}) } @tag_names
    )."\n";
}

1;

