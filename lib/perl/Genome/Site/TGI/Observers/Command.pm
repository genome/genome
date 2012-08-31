package Genome::Site::TGI::Observers::Command;

use strict;
use warnings;
use Command::V1;

Command::V1->add_observer(
    aspect => 'error_die',
    callback => \&command_death_handler,
);

sub command_death_handler {
    my $self = shift;
    my $aspect = shift;
    my %command_death_metrics = @_;

    #Set all defaults to 0 because postgres disallows empty string '' integers

    #parse bjobs -l for a 9 digit build id
    my $build_id = 0;
    $build_id = $1 if $ENV{LSB_JOBID} && `bjobs -l $ENV{LSB_JOBID}` =~ /build(\d{9})\D/;

    #The die message is parsed with a regex to glean extra information
    my $die_message = $command_death_metrics{die_message} || 0;
    $die_message =~ m{(.+?) at /.+?/Genome/(.+?) line (\d+)};
    my $inferred_message = $1 || 0;
    my $inferred_file = $2  || 0;
    my $inferred_line = $3  || 0;

    my $message = $command_death_metrics{error_message} || 0;
    my $package = $command_death_metrics{error_package} || 0;
    my $file = $command_death_metrics{error_file} || 0;
    my $subroutine = $command_death_metrics{error_subroutine} || 0;
    my $line = $command_death_metrics{error_line} || 0;

    my $username = Genome::Sys->username || 0;
    my $sudo_username = Genome::Sys->sudo_username || 0;

    return 1 unless $ENV{GENOME_LOG_COMMAND_ERROR};
    if ($ENV{GENOME_LOG_COMMAND_ERROR} eq 'default') {
        return 1 unless $build_id;
        return 1 unless $username eq 'apipe-builder';

        #Do not double commit error log entries
        my @errors = Genome::Model::Build::ErrorLogEntry->get(build_id=>$build_id);
        return 1 if scalar @errors;
    }

    #replace single quotes so they don't end the shell cmd
    $message =~ s/'/"/g;
    $inferred_message =~ s/'/"/g;

    my $includes = join(' ', map { '-I ' . $_ } UR::Util::used_libs);
    my $cmd = <<EOF
$^X $includes -e 'use Genome;
Genome::Model::Build::ErrorLogEntry->create(
message=>q{$message},
package=>q{$package},
file=>q{$file},
subroutine=>q{$subroutine},
line=>q{$line},
inferred_message=>q{$inferred_message},
inferred_file=>q{$inferred_file},
inferred_line=>q{$inferred_line},
build_id=>q{$build_id},
username=>q{$username},
sudo_username=>q{$sudo_username},
);
UR::Context->commit;'
EOF
    ;

    system($cmd);

    return 1;
}

1;
