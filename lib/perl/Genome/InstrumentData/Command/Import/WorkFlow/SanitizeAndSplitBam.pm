package Genome::InstrumentData::Command::Import::WorkFlow::SanitizeAndSplitBam;

use strict;
use warnings;

use Genome;

require IO::File;
require List::MoreUtils;
require List::Util;
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
    $self->debug_message('Sanitize and split bam...');

    my $metrics = $self->_process_reads;
    return if not $metrics;

    my $verify_read_count_ok = $self->_verify_read_count($metrics);
    return if not $verify_read_count_ok;

    $self->debug_message('Sanitize and split bam...done');
    return 1;
}

sub _process_reads {
    my ($self) = @_;
    $self->debug_message('Processing reads...');

    my $bam_path = $self->bam_path;
    $self->debug_message("Bam path: $bam_path");
    my $bam_fh = IO::File->new("samtools view $bam_path |");
    if ( not $bam_fh ) {
        $self->fatal_message('Failed to open file handle to samtools command!');
    }

    my ($template_name, @reads);
    my %metrics = ( input => 0, output => 0 );
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

        $metrics{input} += @reads;
        $metrics{output} += $self->_write_reads(@reads);

        $template_name = $tokens[0];
        @reads = ( \@tokens );
    }

    if ( @reads ) {
        $metrics{input} += @reads;
        $metrics{output} += $self->_write_reads(@reads);
    }

    for my $fh ( $bam_fh, values %{$self->_read_group_fhs} ) {
        $fh->close;
    }

    $self->debug_message('Processing reads...done');
    return \%metrics;
}

sub _verify_read_count {
    my ($self, $metrics) = @_;
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

    $self->debug_message('Original bam reads:   '.$original_flagstat->{total_reads});
    $self->debug_message('Reads processed:      '.$metrics->{input});
    $self->debug_message('Reads Excluded:       '.($metrics->{input} - $metrics->{output}));
    $self->debug_message('Reads Written:        '.$metrics->{output});
    $self->debug_message('Read group bam reads: '.$read_count);

    if ( $original_flagstat->{total_reads} != $metrics->{input} ) {
        $self->fatal_message(
            'Conflicting read counts from flagstat and input metrics! %s <=> %s',
            $original_flagstat->{total_reads}, $self->_metrics->{input},
        );
    }

    if ( $read_count != $metrics->{output}  ) {
        $self->fatal_message(
            'Conflicting read counts from output flagstats and output metrics! %s <=> %s',
            $read_count, $metrics->{output},
        );
    }

    $self->debug_message('Verify read count...done');
    return 1;
}

sub _write_reads {
    my $self = shift;

    my ($read1, $read2) =_separate_reads(@_);
    my @separated_reads = grep { defined } ( $read1, $read2 );
    return 0 if not @separated_reads;

    my $type = _determine_type_and_set_flags($read1, $read2);
    my $rg_id = _read_group_id_for_reads(@separated_reads);
    my $fh = $self->_get_fh_for_read_group_and_type($rg_id, $type);
    for my $read_tokens ( @separated_reads ) {
        # Sanitize!
        _sanitize_read($read_tokens);
        # Add RG tag
        push @$read_tokens, 'RG:Z:'.$self->old_and_new_read_group_ids->{$rg_id}->{$type};
        $fh->print( join( "\t", @$read_tokens)."\n");
    }

    return @separated_reads;
}

sub _read_group_id_for_reads {
    my @rg_ids;
    for ( @_ ) {
        my $rg_tag = List::Util::first { $_ =~ m/^RG:/ } @$_;
        next if not defined $rg_tag;
        push @rg_ids, (split(':', $rg_tag))[2];
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

sub _sanitize_read {
    # Upcase sequence
    $_[0]->[9] = uc $_[0]->[9];
    # If listed as revcomp, revcomp seq and qual
    if ( $_[0]->[1] & 0x10 ) {
        $_[0]->[9] = reverse $_[0]->[9];
        $_[0]->[9] =~ tr/ATCG/TAGC/;
        $_[0]->[10] = reverse $_[0]->[10];
    }
    # Remove alignment information
    splice @{$_[0]}, 2, 7, (qw/ * 0 0 * * 0 0 /);
    # Remove ALL tags
    splice @{$_[0]}, 11;
}

sub _get_fh_for_read_group_and_type {
    my ($self, $read_group_id, $type) = @_;

    my $key = join('*', $read_group_id, $type);
    unless ( exists $self->_read_group_fhs->{$key} ) {
        my $read_group_bam_path = $self->get_working_bam_path_with_new_extension($self->bam_path, $read_group_id, $type);
        my $samtools_cmd = "| samtools view -S -b -o $read_group_bam_path -";
        $self->debug_message("Opening fh for $read_group_bam_path $type with:\n$samtools_cmd");
        my $fh = IO::File->new($samtools_cmd);
        if ( not $fh ) {
            $self->fatal_message('Failed to open file handle to samtools command!');
        }
        $self->_read_group_fhs->{$key} = $fh;

        $self->_create_new_read_group_ids($read_group_id);
        $self->_write_header($fh, $self->old_and_new_read_group_ids->{$read_group_id}->{$type});

        my @output_bam_paths = $self->output_bam_paths;
        push @output_bam_paths, $read_group_bam_path;
        $self->output_bam_paths(\@output_bam_paths);
    }

    return $self->_read_group_fhs->{$key};
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

