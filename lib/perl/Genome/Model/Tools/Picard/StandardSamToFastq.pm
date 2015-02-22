package Genome::Model::Tools::Picard::StandardSamToFastq;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Picard::StandardSamToFastq {
    is => 'Genome::Model::Tools::Picard',
    has_input => [
        input => {
            is  => 'String',
            doc => 'Input SAM/BAM file to extract reads from. Required.',
        },
    ],
    has_optional_input => [
        fastq => {
            is          => 'String',
            doc         => 'Output fastq file (single-end fastq or, if paired, first end of the pair fastq).',
        },
        second_end_fastq => {
            is          => 'String',
            doc         => 'Output fastq file (if paired, second end of the pair fastq).',
        },
        re_reverse => {
            is          => 'Boolean',
            default     => 1,
            doc         => 'Whether or not to re-reverse reads and qualities reported as negative strand.',
        },
        output_per_rg => {
            is          => 'Boolean',
            default     => 0,
            doc         => 'Output a fastq file per read group',
        },
        output_dir => {
            is          => 'String',
            doc         => 'Directory to write per-readgroup FASTQs if "output_per_rg" is set',
        },
        include_non_pf_reads => {
            is => 'Booelan',
            doc => '',
            default_value => 0,
        },
        include_non_primary_alignments => {
            is => 'Booelan',
            doc => '',
            default_value => 0,
        },
    ],
    doc => 'uses the standard picard "SamToFastq" command to convert a compliant BAM file to FASTQ(s).'
};

sub help_brief {
    'Tool to create FASTQ file from SAM/BAM using Picard';
}

sub help_detail {
    return <<EOS
    Tool to create FASTQ file from SAM/BAM using Picard.  For Picard documentation of this command see:
    http://picard.sourceforge.net/command-line-overview.shtml#SamToFastq
EOS
}

sub execute {
    my $self = shift;

    return unless $self->_validate_inputs();

    my $jar_path = $self->picard_path .'/SamToFastq.jar';
    unless (-e $jar_path) {
        die('Failed to find '. $jar_path .'!  This command may not be available in version '. $self->use_version);
    }

    my $make_arg = sub {
        my ($self, $arg, $is_boolean) = @_;
        my $val = $self->$arg;
        return unless defined $val;

        return unless $val; #Picard chokes when a conflicting boolean option is specified, even if it is "false" (indicating no conflict exists)

        if($is_boolean) {
            $val = $val? 'true' : 'false';
        }

        return join('=', uc($arg), $val);
    };

    my @args;

    push @args,
        map { $make_arg->($self, $_, 0) } ('input', 'fastq', 'second_end_fastq', 'output_dir','max_records_in_ram');
    push @args,
        map { $make_arg->($self, $_, 1) } ('re_reverse', 'output_per_rg', 'include_non_pf_reads', 'include_non_primary_alignments');

    my $convert_cmd = $jar_path .' net.sf.picard.sam.SamToFastq ' . join(' ', @args);
    $self->run_java_vm(
        cmd => $convert_cmd,
        input_files => [$self->input],
        skip_if_output_is_present => 0,
    );

    return 1;
}

sub _validate_inputs {
    my $self = shift;

    my @errors;
    if($self->output_per_rg) {
        unless($self->output_dir) {
            push @errors, 'output_dir must be specified when using "output_per_rg"';
        }
    } else {
        unless($self->fastq) {
            push @errors, 'fastq must be specified when not using "output_per_rg"';
        }
    }

    if(@errors) {
        map { $self->error_message($_) } @errors;
    }

    return !@errors;
}

1;
