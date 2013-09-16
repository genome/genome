use strict;
use warnings;

package Genome::Test::Factory::Test;
use base 'Test::Builder::Module';

use Exporter 'import';
our @EXPORT_OK = qw(test_setup_object);

use Params::Validate qw(validate :types);

sub test_setup_object {
    my $object_factory_class = shift;
    my %p = validate(@_, {
        test_name => {
            type => SCALAR,
            default => "setup_object created a valid object from $object_factory_class",
        },
        setup_object_args => {
            type => ARRAYREF,
            default => [],
        },
    });

    my @errors;
    my $test = sub {
        my @args = @{$p{setup_object_args}};
        my $object = $object_factory_class->setup_object(@args);
        unless ($object) {
            return 0;
        }

        @errors = $object->__errors__();
        if (@errors) {
            return 0;
        }

        return $object;
    };
    my $test_result = $test->();

    my $tb = __PACKAGE__->builder;
    $tb->ok($test_result, $p{test_name});

    unless ($test_result) {
        my $tab = ' ' x 2;
        my @lines;
        my $n_max = scalar(@errors);

        for (my $n = 1; $n <= $n_max; $n++) {
            my $e = $errors[$n - 1];
            push @lines, $tab x 2 . sprintf('Error %d of %d:', $n, $n_max);
            push @lines, $tab x 3 . 'Description: ' . $e->desc;
            push @lines, $tab x 3 . 'Properties: ' . join(', ', $e->properties);
        }
        $tb->diag(join("\n", @lines));
    }

    return $test_result;
}

1;
