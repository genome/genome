package Genome::InstrumentData::IntermediateAlignmentResult::Bwa;

use Genome;
use File::Copy qw/mv/;
use File::Basename;
use File::Slurp qw/read_file/;
use Carp qw/confess/;

use warnings;
use strict;

class Genome::InstrumentData::IntermediateAlignmentResult::Bwa {
    is=>['Genome::InstrumentData::IntermediateAlignmentResult'],
};

sub output_file_prefix {
    my $self = shift;
    my $prefix = $self->input_file;
    my $input_pass = $self->input_pass || '';
    $prefix .= ".$input_pass" if $input_pass;
    return $prefix;
}

sub sai_file {
    my $self = shift;
    return $self->output_dir.'/'.$self->output_file_prefix.".sai";
}

sub log_file {
    my $self = shift;
    return $self->output_dir.'/'.$self->output_file_prefix.".log";
}

sub md5sum {
    my $self = shift;
    my $path = $self->output_dir.'/'.$self->output_file_prefix . '.sai.md5';
    return read_file($path);
}

sub _run_aligner {
    my $self = shift;

    my $fasta_file = $self->aligner_index->full_consensus_path('fa');
    my $tmp_dir = $self->temp_scratch_directory;
    my $input_file = "$tmp_dir/" . $self->input_file;

    my $bam_flag = "";
    if ($input_file =~ /\.bam/) {
        $bam_flag = "-b" . ($self->input_pass||'');
    }

    my $output_file_prefix = $self->output_file_prefix;
    my $sai_file = "$tmp_dir/$output_file_prefix.sai";
    my $log_file = "$tmp_dir/$output_file_prefix.log";

    my $bwa = Genome::Model::Tools::Bwa->path_for_bwa_version($self->aligner_version);

    my @args = (
        "aln",
        $self->aligner_params || '',
        $bam_flag,
        $fasta_file,
        $input_file,
        "1> $sai_file 2>> $log_file"
    );

    my $cmd = join(" ", $bwa, @args);

    # disconnect the db handle before this long-running event
    if (Genome::DataSource::GMSchema->has_default_handle) {
        $self->debug_message("Disconnecting GMSchema default handle.");
        Genome::DataSource::GMSchema->disconnect_default_dbh();
    }

    Genome::Sys->shellcmd(
        cmd          => $cmd,
        input_files  => [ $fasta_file, $input_file ],
        output_files => [ $sai_file, $log_file ],
        skip_if_output_is_present => 0,
    );

    unless ($self->_verify_bwa_aln_did_happen(sai_file => $sai_file, log_file => $log_file)) {
        confess $self->error_message("bwa aln did not seem to successfully take place for " . $fasta_file);
    }

    my $sai_bytes = -s $sai_file;
    $self->debug_message("Created $sai_file ($sai_bytes bytes)");

    my $md5_file = $self->_create_md5($sai_file);
    unless($md5_file) {
        confess $self->error_message("Failed to create md5 file for sai file $sai_file");
    }

    $self->_promote_to_staging($sai_file, $log_file, $md5_file);
}

# TODO: put this somewhere else, maybe SR::Stageable
sub _create_md5 {
    my ($self, $file) = @_;

    my $md5_file = "$file.md5";
    my $cmd      = "md5sum $file > $md5_file";

    my $rv  = Genome::Sys->shellcmd(
        cmd                        => $cmd,
        input_files                => [$file],
        output_files               => [$md5_file],
        skip_if_output_is_present  => 0,
    );
    $self->error_message("Fail to run: $cmd") and return unless $rv == 1;
    return $md5_file;
}

sub _promote_to_staging {
    my ($self, @files) = @_;
    for my $f (@files) {
        mv($f, $self->temp_staging_directory);
    }
}

sub _verify_bwa_aln_did_happen {
    my $self = shift;
    my %p = @_;

    unless (-e $p{sai_file} && -s $p{sai_file}) {
        $self->error_message("Expected SAI file is $p{sai_file} nonexistent or zero length.");
        return;
    }

    my $count = $self->_inspect_log_file(
        log_file=>$p{log_file},
        log_regex=>'(\d+) sequences have been processed'
    );

    unless($count) {
        $self->error_message("Expected to see 'X sequences have been processed' in the log file where 'X' must be a nonzero number.");
        return;
    }

    my $expected_count;
    my $input_file = $self->temp_scratch_directory . "/" . $self->input_file;
    if($input_file =~ /\.bam/) {
        my $output_file;
        if(-s $self->flagstat_file) {
            $output_file = $self->flagstat_file;
        } else {
            $output_file = $self->temp_scratch_directory . "/input_bam.flagstat";
            my $cmd = Genome::Model::Tools::Sam::Flagstat->create(
                bam_file       => $input_file,
                output_file    => $output_file,
                use_version    => $self->samtools_version,
                include_stderr => 1,
            );
            unless ($cmd and $cmd->execute) {
                $self->error_message('Failed to create or execute flagstat command.');
                return;
            }
        }

        my $stats = Genome::Model::Tools::Sam::Flagstat->parse_file_into_hashref($output_file);
        unless($stats) {
            $self->status_message('Failed to get flagstat data  on input sequences from '.$output_file);
            return;
        }

        my $input_pass = $self->input_pass;
        if($input_pass eq 1) {
            $expected_count = $stats->{reads_marked_as_read1};
        } elsif($input_pass eq 2) {
            $expected_count = $stats->{reads_marked_as_read2};
        } elsif($input_pass eq 0) {
            $expected_count = $stats->{total_reads}; #fragment data
        }
    } else {
        my $wc_count_file = $self->temp_scratch_directory . '/wc.count';
        unless(Genome::Sys->shellcmd(
            cmd => 'wc -l ' . $input_file . ' > ' . $wc_count_file,
            input_files => [$input_file],
            output_files => [$wc_count_file],
        )) {
            $self->error_message('Failed to determine line count of FASTQ.');
            return;
        }

        my $wc_count = Genome::Sys->read_file($wc_count_file);
        ($expected_count) = $wc_count =~ /^(\d+)/;
        $expected_count /= 4;
    }

    if($count eq $expected_count) {
        $self->debug_message('Log reported expected count of ' . $expected_count . ' sequences processed.');
    } else {
        $self->error_message('Expected to process ' . $expected_count . ' sequences but processed ' . $count);
        return;
    }

    return 1;
}

sub _inspect_log_file {
    my $self = shift;
    my %p = @_;

    my $log_file = $p{log_file};
    unless ($log_file and -s $log_file) {
        $self->error_message("log file $log_file is not valid");
        return;
    }

    # XXX I am trying to avoid switching on $bwa_version in this module, but
    # I'm not entirely sure that storing 'log_format' in GMT::Bwa and defining
    # the logic for interpreting that property here is the correct separation
    # of functionality between GMT::Bwa and this module.
    my $last_line;
    my $log_format = Genome::Model::Tools::Bwa->log_format($self->aligner_version);
    if ($log_format eq 'new') {
        $last_line = `tail -4 $log_file | head -1`;
    } else {
        $last_line = `tail -1 $log_file`;
    }
    my $check_nonzero = 0;

    my $log_regex = $p{log_regex};
    if ($log_regex =~ m/\(\\d\+\)/) {
        $check_nonzero = 1;
    }

    if ($last_line =~ m/$log_regex/) {
        if ( !$check_nonzero || $1 > 0 ) {
            $self->debug_message('The last line of aligner.log matches the expected pattern');
            return $1;
        }
    }

    $self->error_message("The last line of $log_file is not valid: $last_line");
    return;
}


1;
