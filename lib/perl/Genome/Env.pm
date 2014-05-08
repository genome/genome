package Genome::Env;

use strict;
use warnings;

use Module::Find qw(findsubmod usesub);
use Set::Scalar qw();

sub used {
    return grep { /^GENOME_/ } keys %ENV;
}

sub allowed_modules {
    return usesub Genome::Env;
}

sub allowed {
    return map { NAME($_) } allowed_modules();
}

sub check_genome_variables {
    my @valid = allowed();
    my @used  = used();

    my $valid_set = Set::Scalar->new(@valid);
    my $used_set = Set::Scalar->new(@used);
    my $invalid_set = $used_set - $valid_set;

    unless ($invalid_set->is_empty) {
        print STDERR "Available environment variable(s):\n",
            '  ', join("\n  ", @valid), "\n\n";
        print STDERR "Unrecognized environment variable(s) found:\n",
            '  ', join("\n  ", map { join('=', $_, ($ENV{$_} || '')) } $invalid_set->members), "\n";
        return;
    }

    return 1;
}

sub set_default_values {
    my @valid = allowed_modules();
    for my $module (@valid) {
        my $name = NAME($module);
        if ($module->can('default_value')
            && !defined($ENV{$name})
        ) {
            $ENV{$name} = $module->default_value;
        }
    }
}

sub NAME {
    my $class = shift;

    my $tail = ($class =~ /^Genome::Env::(.*)/)[0];
    my $NAME = uc($tail);
    unless ($NAME) {
        die 'Could not infer NAME from module';
    }
    return $NAME;
}

# Load all modules in Genome/Env/* and set their default value
# unless the environment variable is already defined.
set_default_values();

# Make sure that all Genome environment variables that are
# set correspond to a file in Genome/Env/*.
check_genome_variables(); # instead of 1; this should be last statment
