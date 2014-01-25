package Genome::Env;

use strict;
use warnings;

use File::Basename;

# Make sure that all Genome environment variables that are
# set correspond to a file in Genome/Env/*.
check_genome_variables();

# Load all modules in Genome/Env/* and set their default value
# unless the environment variable is already defined
set_default_values();

sub check_genome_variables {
    for my $variable (keys %ENV) {
        next unless $variable =~ /^GENOME_/;
        my $package = __PACKAGE__ . "::$variable";
        eval "use $package";
        if ($@) {
            my @paths = get_all_submodule_paths();
            my @variables = get_variables_from_paths(@paths);
            print STDERR "Environment variable $variable set to " . $ENV{$variable} .
                " but there were errors using $package:\n" .
                "Available variables:\n\t" .
                join("\n\t", @variables) .
                "\n";
            exit 1;
        }
    }
    return 1;
}

sub set_default_values {
    my @paths = get_all_submodule_paths();
    my @variables = get_variables_from_paths(@paths);
    for my $variable (@variables) {
        next if exists $ENV{$variable};
        my $package = __PACKAGE__ . "::$variable";
        eval "use $package";
        die $@ if $@;

        no strict 'refs';
        no warnings 'once';
        my $default = ${$package . "::default_value"};
        use warnings 'once';
        use strict 'refs';
        next unless $default;
        $ENV{$variable} = $default;
    }
    return 1;
}

sub get_all_submodule_paths {
    my $base_path = __FILE__;
    $base_path =~ s/.pm$//;
    $base_path =~ s/ /\\ /g; # escape spaces for glob
    my @paths = glob("$base_path/*.pm");
    return @paths;
}

sub get_variables_from_paths {
    my @paths = @_;
    my @variables = map { basename($_) } @paths;
    map { $_ =~ s/.pm$// } @variables;
    return @variables;
}

sub NAME {
    my $class = shift;

    my $tail = ($class =~ /^Genome::Env::(.*)/)[0];
    my $NAME = uc($tail);
    unless ($NAME) {
        die 'Could not infer NAME from package';
    }
    return $NAME;
}

1;
