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

1;
