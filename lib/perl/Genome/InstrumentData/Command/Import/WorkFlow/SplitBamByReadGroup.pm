package Genome::InstrumentData::Command::Import::WorkFlow::SplitBamByReadGroup;

use strict;
use warnings;

use Genome;

require IO::File;
require List::MoreUtils;
use Params::Validate ':types';

class Genome::InstrumentData::Command::Import::WorkFlow::SplitBamByReadGroup { 
    is => 'Command::V2',
    has_input => [
        bam_path => {
            is => 'Text',
            doc => 'The path of the unsorted bam to sort.',
        }
    ],
    has_output => [ 
        output_bam_paths => {
            is => 'Text',
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
    has_optional_calculated => [
        read_group_ids => {
            calculate_from => [qw/ read_groups_and_tags /],
            calculate => q( return keys %$read_groups_and_tags ), 
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
    my $read_groups_and_tags = $helpers->read_groups_and_tags_from_headers($read_group_headers);
    return if not $read_groups_and_tags;

    # Add unknown rg
    $read_groups_and_tags->{unknown} = { ID => 'unknown', CN => 'NA', };
    
    # Add instdata uuid for each rg
    my %old_and_new_read_group_ids;
    for my $rg_tags ( values %$read_groups_and_tags ) {
        my $rg_id = delete $rg_tags->{ID};
        $old_and_new_read_group_ids{$rg_id} = {
            paired => UR::Object::Type->autogenerate_new_object_id_uuid,
            #singleton => UR::Object::Type->autogenerate_new_object_id_uuid,
        };
        $old_and_new_read_group_ids{$rg_id}->{singleton} = $old_and_new_read_group_ids{$rg_id}->{paired};
    }
    $self->old_and_new_read_group_ids(\%old_and_new_read_group_ids);

    $self->headers($headers);
    $self->read_groups_and_tags($read_groups_and_tags);

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

    my $previous_tokens;
    my $previous_read_group_id;
    while ( my $line = $bam_fh->getline ) {
        my @tokens = split(/\t/, $line);
        my $read_group_id = $self->_get_rg_id_from_sam_tokens(\@tokens);

        unless($previous_tokens) {
            $previous_tokens = \@tokens;
            $previous_read_group_id = $read_group_id;
            next;
        }

        if($tokens[0] eq $previous_tokens->[0] and $previous_read_group_id eq $read_group_id) {
            my $fh = $self->_write_reads_based_on_read_group_and_pairedness(
                rg_id => $read_group_id,
                pairedness => 'paired',
                reads => [ $previous_tokens, \@tokens ],
            );
            undef $previous_tokens;
            undef $previous_read_group_id;
        } else {
            my $fh = $self->_write_reads_based_on_read_group_and_pairedness(
                rg_id => $previous_read_group_id,
                pairedness => 'singleton',
                reads => [ $previous_tokens, ],
            );
            $previous_tokens = \@tokens;
            $previous_read_group_id = $read_group_id;
        }
    }

    if($previous_tokens) {
        my $fh = $self->_write_reads_based_on_read_group_and_pairedness(
            rg_id => $previous_read_group_id,
            pairedness => 'singleton',
            reads => [ $previous_tokens, ],
        );
    }

    for my $fh ( $bam_fh, values %{$self->_read_group_fhs} ) {
        $fh->close;
    }

    $self->debug_message('Write reads...done');
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
        $self->error_message('No read group bams passed validation!');
        return;
    }

    $self->output_bam_paths(\@validated_output_bam_paths);

    my $original_flagstat = $helpers->load_or_run_flagstat($self->bam_path);
    return if not $original_flagstat;

    $self->debug_message('Original bam read count: '.$original_flagstat->{total_reads});
    $self->debug_message('Read group bams read count: '.$read_count);

    if ( $original_flagstat->{total_reads} != $read_count ) {
        $self->error_message('Original and split by read group bam read counts do not match!');
        return;
    }

    $self->debug_message('Verify read count...done');
    return 1;
}

sub _write_reads_based_on_read_group_and_pairedness {
    my ($self, %params) = @_;

    my $fhs = $self->_read_group_fhs;
    my $key = join('*', $params{rg_id}, $params{pairedness});
    unless(exists $fhs->{$key}) {
        $fhs->{$key} = $self->_open_fh_for_read_group_and_pairedness($params{rg_id}, $params{pairedness});
        $self->_read_group_fhs($fhs);
    }

    for my $read_tokens ( @{$params{reads}} ) {
        $self->_update_read_group_for_sam_tokens_based_on_paired_endness($read_tokens, $params{pairedness});
        $fhs->{$key}->print( join("\t", @$read_tokens) );
    }

    return 1;
}

sub _update_read_group_for_sam_tokens_based_on_paired_endness {
    my ($self, $tokens, $pairedness) = @_;

    my $rg_tag_idx = List::MoreUtils::firstidx { $_ =~ m/^RG:/ } @$tokens;
    my $rg_tag;
    if ( defined $rg_tag_idx ) {
        $rg_tag = splice(@$tokens, $rg_tag_idx, 1);
    }
    else {
        $rg_tag = 'RG:Z:unknown';
        $rg_tag_idx = $#$tokens;
    }

    my @rg_tokens = split(':', $rg_tag);
    my $new_rg_id = $self->old_and_new_read_group_ids->{$rg_tokens[2]}->{$pairedness};
    $rg_tag = 'RG:Z:'.$new_rg_id;
    splice(@$tokens, $rg_tag_idx, 0, $rg_tag); # Add RG tag back

    return 1;
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

    $self->_write_headers_for_read_group($fh, $read_group_id, $pairedness);

    my @output_bam_paths = $self->output_bam_paths;
    push @output_bam_paths, $read_group_bam_path;
    $self->output_bam_paths(\@output_bam_paths);

    return $fh;
}

sub _write_headers_for_read_group {
    my ($self, $fh, $read_group_id, $pairedness) = @_;

    my $headers = $self->headers;
    Carp::confess('No headers to write to read group bams!') if not $headers;

    my $helpers = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get;
    my $headers_as_string = $helpers->headers_to_string($headers);
    return if not $headers_as_string;

    $fh->print( $headers_as_string );

    my $read_groups_and_tags = $self->read_groups_and_tags;
    my %rg_tags = %{$read_groups_and_tags->{$read_group_id}};
    my @tag_names = sort keys %rg_tags;
    $fh->print( 
        join(
            "\t", '@RG',
            'ID:'.$self->old_and_new_read_group_ids->{$read_group_id}->{$pairedness},
            map { join(':', $_, $rg_tags{$_}) } @tag_names
        )."\n"
    );

    return 1;
}


1;

