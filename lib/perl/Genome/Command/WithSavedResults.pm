package Genome::Command::WithSavedResults;
use Genome::SoftwareResult::Default;

class Genome::Command::WithSavedResults {
    is => 'Command::V2',
    is_abstract => 1,
};

sub _init_subclass {
    my $subclass_name = shift;
    my $src = "package $subclass_name;\n" . <<EOS;
        use Moose; # inject into the namespace
        around execute => \\&Genome::SoftwareResult::Default::execute_wrapper;
EOS
    eval $src;
    if ($@) {
        die "error initializing $subclass_name from " . __PACKAGE_ . ": $@";
    }
    return 1;
}

1;

