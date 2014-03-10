#!/usr/bin/env perl

use Cwd qw(realpath);
use File::Basename;
require LWP::UserAgent;

print STDERR "Getting current snapshot version... ";
my $current_snapshot_version = current_snapshot_version();
print STDERR "$current_snapshot_version\n";

print STDERR "Checking current version...\n";
if (wait_for_valid_http_response($current_snapshot_version, 10)) {
    print STDERR "Web server is already running on $current_snapshot_version, skipping restart.\n";
} else {
    print STDERR "Restarting web server...\n";
    unless (restart_web_server($current_snapshot_version)) {
        exit 255;
    }
}

sub current_snapshot_version {
    my $snapshot_path = realpath('/gsc/scripts/opt/genome/current/web');
    return basename($snapshot_path);
}

sub restart_web_server {
    my $snapshot_version = shift;

    reload_genome_view();
    my $restart_succeeded = wait_for_valid_http_response($snapshot_version);
    unless ($restart_succeeded) {
        print STDERR "\n\nReload failed so attempting a stop and start.\n\n";
        return unless stop_genome_view();
        return unless start_genome_view();
        $restart_succeeded = wait_for_valid_http_response($snapshot_version);
    }

    return $restart_succeeded;
}

sub wait_for_valid_http_response {
    my $snapshot_version = shift;
    my $max_elapsed_seconds = shift || 120;
    my $start_time = time;

    my $url = 'https://127.0.0.1/view/genome/search/status.html';

    my $user_agent = LWP::UserAgent->new(timeout => 10);
    if ($user_agent->can('ssl_opts')) {
        $user_agent->ssl_opts(verify_hostname => 0);
    }

    my $response;
    do {
        $response = $user_agent->get($url);
        sleep(3);
    } until (valid_http_response($response, $snapshot_version) || (time - $start_time) >= $max_elapsed_seconds);

    return valid_http_response($response, $snapshot_version);
}

sub valid_http_response {
    my $response = shift || return;
    my $snapshot_version = shift || die;

    my $code = $response->code;

    if ($response->is_success) {
        my @content = split("\n", $response->decoded_content);
        if ( grep { $_ =~ /^<!-- Genome: .*$snapshot_version\/lib\/perl\/Genome.pm-->$/ } @content ) {
            print STDERR "OK: Web server is up (" . $response->status_line . ") on new version ($snapshot_version).\n";
            return 1;
        }
        else {
            print STDERR "Not OK: Web server is up (" . $response->status_line . ") but *not* on new version ($snapshot_version).\n";
            return;
        }
    }
    else {
        print STDERR "Not OK: Failed to get successful response (" . $response->status_line . ") from URL (". $response->filename . ").\n";
        return;
    }
}

sub reload_genome_view {
    my $reload_genome_view_cmd = "sudo /etc/init.d/genome_view reload";
    print STDERR "RUNNING ($reload_genome_view_cmd)...\n";
    my $exit_code = system("$reload_genome_view_cmd");
    if ($exit_code != 0) {
        print STDERR "Command ($reload_genome_view_cmd) did not succeed.\n";
        return;
    }
    return 1;
}

sub stop_genome_view {
    my $stop_genome_view_cmd = "sudo /etc/init.d/genome_view stop";
    print STDERR "RUNNING ($stop_genome_view_cmd)...\n";
    my $exit_code = system("$stop_genome_view_cmd");
    if ($exit_code != 0) {
        print STDERR "Command ($stop_genome_view_cmd) did not succeed.\n";
        return;
    }
    return 1;
}

sub start_genome_view {
    my $start_genome_view_cmd = "sudo /etc/init.d/genome_view start";
    print STDERR "RUNNING ($start_genome_view_cmd)...\n";
    my $exit_code = system("$start_genome_view_cmd");
    if ($exit_code != 0) {
        print STDERR "Command ($start_genome_view_cmd) did not succeed.\n";
        return;
    }
    return 1;
}
