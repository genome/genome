package Genome::Model::Build::Command::DetermineError;

use strict;
use warnings;

use Genome;

use Genome::Ptero::Utils;

use Ptero::HTTP;
use Ptero::Proxy::Workflow::Execution;

sub handle_failed_from_ptero {
    my $self = shift;
    my $ptero_proxy = shift;

    my @failed_steps = $self->_find_failed_ptero_steps($ptero_proxy);

    for my $step (@failed_steps) {
        my $found = $self->_set_status_from_ptero_step($step);
        last if $found;
    }
}

sub _find_failed_ptero_steps {
    my $self = shift;
    my $ptero_proxy = shift;

    my @all_steps = Genome::Ptero::Utils::get_all_executions_for_proxy($ptero_proxy);

    my @failed_executions = grep {
        Ptero::Statuses::is_abnormal($_->{status})
    } @all_steps;

    my @details = map {
        Ptero::Proxy::Workflow::Execution->new( url => $_->{details_url} )
    } @failed_executions;

    my @real_steps = grep { _is_lsf_job($_) } @details;

    return sort { _first_timestamp($b) cmp _first_timestamp($a) } @real_steps;
}

sub _set_status_from_ptero_step {
    my $self = shift;
    my $step = shift;

    my $stderr = _execution_stderr($step);
    return unless defined($stderr) && -s $stderr;

    my @error_data = parse_error_log($stderr);

    unless(grep { $_ } @error_data) {
        my $stdout = _job_stdout($step);
        @error_data = parse_output_log($stdout);
    }

    unless(grep { $_ } @error_data) {
        @error_data = find_die_or_warn_in_log($stderr);
    }

    $self->error_log($stderr);
    my %error_data;
    @error_data{qw(error_source_file error_source_line error_host error_date error_text)} = @error_data;
    for my $key (keys %error_data) {
        $self->$key($error_data{$key}) if defined $error_data{$key};
    }

    return 1;
}

sub _execution_stderr {
    my $execution_proxy = shift;

    return $execution_proxy->concrete_execution->{data}{stderr_log};
}

sub _job_stdout {
    my $execution_proxy = shift;

    my $lsf_info = Ptero::HTTP::make_request_and_decode_response(
        method => 'GET',
        url => $execution_proxy->concrete_execution->{data}{jobUrl},
    );

    my $temp_file = Genome::Sys->create_temp_file_path;
    return Genome::Sys->write_file($temp_file, $lsf_info->{stdout});
}

sub _is_lsf_job {
    my $execution_proxy = shift;

    my $url = $execution_proxy->concrete_execution->{data}{jobUrl};

    return $url && $url =~ m!^http://lsf!;
}

sub _first_timestamp {
    my $execution_proxy = shift;

    return $execution_proxy->concrete_execution->time_started;
}

sub _get_all_executions {
    my $workflow = shift;

    my @executions = @{$workflow->workflow_executions};
    my @child_executions;

    for my $e (@executions) {
        for my $child (@{ $e->child_workflow_proxies }) {
            push @child_executions, _get_all_executions($child);
        }
    }

    return (@executions, @child_executions);
}

1;
