package Genome::Model::Tools::Joinx::VcfMerge;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Joinx::VcfMerge {
    is => 'Genome::Model::Tools::Joinx',
    has_input => [
        input_files => {
            is => 'Text',
            is_many => 1,
            doc => 'List of vcf files to merge',
            shell_args_position => 1,
        },
    ],
    has_optional_input => [
        output_file => {
            is => 'Text',
            is_output => 1,
            doc => 'The output file (defaults to stdout)',
        },
        merge_strategy_file => {
            is => 'Text',
            doc => "The merge strategy file, provided to joinx to specify how to handle info field merging (-M option)",
        },
        sample_priority => {
            is => 'Text',
            doc => 'Sample priority. Order prefers data from file A over B. Unfiltered prefers data from the first filtered sample (intersect). Filters prefers data from the first unfiltered sample (union). Requires joinx1.6+',
            valid_values => ['order', 'unfiltered', 'filtered'],
        },
        clear_filters => {
            is => 'Boolean',
            default => 0,
            doc => 'Merged entries will have the FILTER column stripped out (-c option)',
        },
        merge_samples => {
            is => 'Boolean',
            default => 0,
            doc => 'Allow input files with overlapping samples (-s option)',
        },
        ratio_filter => {
            is => 'Text',
            doc => 'Require this ratio of inputs to agree on calls or else be filtered. Provide as a string: "ratio,filter name,filter description" ' .
                'I.E. "1.0,CNS,Sample consensus filter" will mark any variant where all input files did not agree (on position) as filtered with "CNS".',
        },
        use_bgzip => {
            is => 'Boolean',
            doc => 'zcats the input files into stdin, and bgzips the output',
            default => 0,
        },
        error_log => {
            is => 'Text',
            doc => 'where to redirect stderr from joinx, if desired',
        },
    ],
};

sub help_brief {
    "Merges one or more vcf files."
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
  gmt joinx vcf-merge a.vcf b.vcf ... --output-file merged.vcf
EOS
}

sub execute {
    my ($self) = @_;

    my @inputs = $self->_resolve_inputs();
    my $output = $self->_resolve_output();
    unless(@inputs) {
        if(defined($self->output_file)) {
            # make output file exist (for downstream tools)
            unless(system("touch $output") == 0) {
                die $self->error_message("Failed to touch $output");
            }
        }
        return 1;
    }

    my $flags = $self->_resolve_flags();
    my $joinx_bin_path = $self->joinx_path();
    my ($cmd, $params) = $self->_generate_joinx_command($joinx_bin_path,
            $flags, \@inputs, $output);
    my %params = %{$params};

    $self->status_message("Executing command: $cmd");
    Genome::Sys->shellcmd(%params);
    return 1;
}

sub _resolve_inputs {
    my ($self) = @_;

    # Grep out empty files
    my @non_empty_inputs;
    my @empty_inputs;
    for my $input ($self->input_files) {
        if(-s $input) {
            push(@non_empty_inputs, $input);
        } else {
            push(@empty_inputs, $input);
        }
    }
    if(@empty_inputs) {
        $self->warning_message("Found the following inputs were of zero size:\n\t".
            join("\n\t", @empty_inputs));
    }
    my @inputs = @non_empty_inputs;

    return @inputs;
}

sub _resolve_output {
    my ($self) = @_;

    if($self->use_bgzip && not $self->output_file){
       die $self->error_message("If use_bgzip is set, output_file must also be set, otherwise binary nonsense will spew forth.");
    }

    my $output = "-";
    $output = $self->output_file if defined($self->output_file);

    return $output;
}

sub _resolve_flags {
    my ($self) = @_;

    my $flags = "";
    if ($self->clear_filters) {
        $flags .= " -c";
    }
    if ($self->merge_strategy_file) {
        $flags .= " -M " . $self->merge_strategy_file;
    }
    if ($self->sample_priority) {
        if($self->sample_priority eq 'order') {
            $flags .= " -P o";
        } elsif ($self->sample_priority eq 'unfiltered') {
            $flags .= " -P u";
        } elsif ($self->sample_priority eq 'filtered') {
            $flags .= " -P f";
        } else {
            die $self->error_message("Invalid sample priority set: " .
                    $self->sample_priority);
        }
    }
    if ($self->merge_samples) {
        $flags .= " -s";
    }
    if ($self->ratio_filter) {
        $flags .= ' -R "' . $self->ratio_filter . '"';
    }

    return $flags;
}

# combine the parts of the joinx command together and return the command and
# params for Sys->shellcmd
sub _generate_joinx_command {
    my ($self, $joinx_bin_path, $flags, $inputs, $output) = @_;
    my @inputs = @{$inputs};

    @inputs = map {"<(zcat $_)"} @inputs if $self->use_bgzip;

    my $cmd = $joinx_bin_path . " vcf-merge $flags " . join(" ", @inputs);
    if($self->output_file) {
        my $log_part = '';
        $log_part = sprintf(' 2> %s', $self->error_log) if $self->error_log;

        if($self->use_bgzip) {
            $cmd .= sprintf('%s | bgzip -c > %s', $log_part, $output);
        } else {
            $cmd .= sprintf(' -o %s%s', $output, $log_part);
        }
    }

    my %params = (
        cmd => $cmd,
        input_files => $inputs, #without zcat wrapper
        allow_zero_size_output_files=>1,
    );
    $params{output_files} = [$output] if $self->output_file;
    return ($cmd, \%params);
}

1;
