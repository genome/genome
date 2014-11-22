package Genome::Model::Build::Command::Base;

class Genome::Model::Build::Command::Base {
    is          => 'Command::V2',
    is_abstract => 1,
    has         => [
        builds => {
            is                  => 'Genome::Model::Build',
            is_many             => 1,
            shell_args_position => 1,
            doc =>
              'Build(s) to use. Resolved from command line via text string.',
        },
    ],
};

sub _limit_results_for_builds {
    my ( $class, @builds ) = @_;

    my $user = getpwuid($<);
    unless ($user) {
        die "Could not get username from getpwuid.";
    }

    print STDERR "Screening builds that you are not able to modify... \n";

    my @apipe_builder_builds;
    my @run_by_builds;
    for my $build (@builds) {
        if (not $build->run_by) { # since we don't know who, we'll let them try
            push @run_by_builds, $build;
        }
        elsif ($user eq $build->run_by) {
            push @run_by_builds, $build;
        }
        elsif ($user eq 'apipe-builder') {
            if ( $class =~ /::Stop$/    || $class =~ /::Remove$/ ||
                 $class =~ /::Abandon$/ || $class =~ /::Start$/   ) {
                    push @apipe_builder_builds, $build;
            }
            else {
                next;
            }
        }
        else {
            next;
        }
    }

    if (@apipe_builder_builds) {
        print STDERR "\n* Allowing " . @apipe_builder_builds . " since you are apipe-builder.\n* ";
    }

    print STDERR "Found " . @run_by_builds . " builds out of " . @builds . " that can be modified by you.\n";

    return (@run_by_builds, @apipe_builder_builds);
}

1;

