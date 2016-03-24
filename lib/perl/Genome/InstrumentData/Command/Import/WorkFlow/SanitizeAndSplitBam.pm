package Genome::InstrumentData::Command::Import::WorkFlow::SanitizeAndSplitBam;

use strict;
use warnings;

use Genome;

require IO::File;
require List::MoreUtils;
use Params::Validate ':types';

class Genome::InstrumentData::Command::Import::WorkFlow::SanitizeAndSplitBam {
    is => 'Command::V2',
    roles => [qw/
        Genome::InstrumentData::Command::Import::WorkFlow::Role::WithWorkingDirectory
        Genome::InstrumentData::Command::Import::WorkFlow::Role::RemovesInputFiles
    /],
    has_input => [
        bam_path => {
            is => 'FilePath',
            doc => 'The path of the unsorted bam to sort.',
        },
        library => {
            is => 'Genome::Library',
            doc => 'Library to use in the bam headers.',
        },
    ],
    has_output => [ 
        output_bam_paths => {
            is => 'FilePath',
            is_many => 1,
            doc => 'The paths of the read group bams.',
        },
    ],
    has_optional_transient => [
        old_and_new_read_group_ids => { is => 'HASH', default => {}, },
        _read_group_fhs => { is => 'HASH', default => {} },
    ],
};

sub execute {
    my $self = shift;
    $self->debug_message('Spilt bam by read group...');

    my $write_reads_ok = $self->_process_reads;
    return if not $write_reads_ok;

    my $verify_read_count_ok = $self->_verify_read_count;
    return if not $verify_read_count_ok;

    $self->debug_message('Spilt bam by read group...done');
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

    my ($template_name, @reads);
    while ( my $line = $bam_fh->getline ) {
        chomp $line;
        my @tokens = split(/\t/, $line);

        if ( not $template_name ) {
            $template_name = $tokens[0];
        }

        if ( $tokens[0] eq $template_name ) {
            # group reads by template name
            push @reads, \@tokens;
            next;
        }

        $self->_write_reads(@reads);
        $template_name = $tokens[0];
        @reads = ( \@tokens );
    }

    if ( @reads ) {
        $self->_write_reads(@reads);
    }

    for my $fh ( $bam_fh, values %{$self->_read_group_fhs} ) {
        $fh->close;
    }

    $self->debug_message('Processing reads...done');
    return 1;
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

sub _write_reads {
    my ($self, $reads) = @_;

    my $rg_id = $self->_read_group_id_for_reads($reads);
    my ($read1, $read2) =_separate_reads($reads);
    my $type = _determine_type($read1, $read2);

    my $fhs = $self->_read_group_fhs;
    my $key = join('*', $rg_id, $type);
    unless(exists $fhs->{$key}) {
        $fhs->{$key} = $self->_open_fh_for_read_group_and_type($rg_id, $type);
        $self->_read_group_fhs($fhs);
    }

    for my $read_tokens ( grep { defined } ( $read1, $read2 ) ) {
        # Revcomp!
        $read_tokens->[1] = _correct_sequence_and_qualities($read_tokens);
        # Remove alignment information
        splice @$read_tokens, 2, 7, (qw/ * 0 0 * * 0 0 /);
        # Remove ALL tags, add new RG tag
        splice @$read_tokens, 11;
        push @$read_tokens, 'RG:Z:'.$self->old_and_new_read_group_ids->{$rg_id}->{$type};
        #print( join( "\t", @$read_tokens)."\n");
        $fhs->{$key}->print( join( "\t", @$read_tokens)."\n");
    }

    return 1;
}

sub _read_group_id_for_reads {
    my ($self, @reads) = @_;

    my @rg_ids;
    for my $read ( @reads ) {
        my $rg_tag_idx = List::MoreUtils::firstidx { $_ =~ m/^RG:/ } @$read;
        next if defined $rg_tag_idx;
        push @rg_ids, (split(':', $read->[$rg_tag_idx]))[2];
    }

    return 'unknown' if not @rg_ids;
    return (List::MoreUtils::uniq(@rg_ids))[0];
}

sub _separate_reads {
    # Separate reads into 1/2 removing supplementary but including the secondary alignments
    my (@read1s, @read2s);
    for my $read ( @_ ) {
        # Remove Supplementary 2048 0x800
        next if $read->[1] & 0x800;
        if ( $read->[1] & 0x80 ) {
            # read 2
            push @read2s, $read; # flagged read2
        }
        else {
            push @read1s, $read; # unflagged reads will be read1
        }
    }

    return ( $read1s[0], $read2s[0] );
}

sub _determine_type_and_set_flags {
    # Determine the type paired/singleton and set the flags accordingly
    my $type;
    # paired will get:
    # 64/128 first/second in pair
    # 1      read paired
    # 4      read unmapped
    # 8      mate unmapped
    if ( defined $_[0] and defined $_[1] ) {
        $type = 'paired';
        $_[0]->[1] = 77;
        $_[1]->[1] = 141;
    }
    # sigleton will get:
    # 4  first in pair
    # 64 read unmapped
    elsif ( defined $_[0] ) {
        $type = 'singleton';
        $_[0]->[1] = 68;
    }
    elsif ( defined $_[1] ) {
        $type = 'singleton';
        $_[1]->[1] = 68;
    }
    # Gotta send in reads!
    else {
        die 'No reads given to _determine_type_and_set_flags!';
    }
    $type;
}

sub _correct_sequence_and_qualities {
    $_[0]->[9] = uc $_[0]->[9];
    if ( $_[0]->[1] & 0x10 ) {
        $_[0]->[9] = reverse $_[0]->[9];
        $_[0]->[9] =~ tr/ATCG/TAGC/;
        $_[0]->[10] = reverse $_[0]->[10];
    }
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

    $self->_create_new_read_group_ids($read_group_id);
    $self->_write_header($fh, $self->old_and_new_read_group_ids->{$read_group_id}->{$type});

    my @output_bam_paths = $self->output_bam_paths;
    push @output_bam_paths, $read_group_bam_path;
    $self->output_bam_paths(\@output_bam_paths);

    return $fh;
}

sub _create_new_read_group_ids {
    my ($self, $rg_id) = @_;

    return 1 if exists $self->old_and_new_read_group_ids->{$rg_id};
    $self->old_and_new_read_group_ids->{$rg_id} = {
        paired => UR::Object::Type->autogenerate_new_object_id_uuid,
        singleton => UR::Object::Type->autogenerate_new_object_id_uuid,
    };
}

sub _write_header {
    my ($self, $fh, $rg_id) = @_;

    # Header Line
    $fh->print( join("\t", '@HD', 'VN:1.4', 'SO:queryname')."\n" );

    # Read Group
    $fh->print(
        join(
            "\t", '@RG',
            'ID:'.$rg_id,
            'LB:'.$self->library->name,
            'SM:'.$self->library->sample->name,
        )."\n"
    );
}

1;

