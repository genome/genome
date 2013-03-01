package Genome::Command::WithSavedResults;
use Genome::SoftwareResult::Default;

class Genome::Command::WithSavedResults {
    is => 'Command::V2',
    type_has => [
        parallelize_by => { is => 'ARRAY', is_optional => 1,
                            doc => 'produce intermediate results and merge, grouping by this/these attributes' },     
    ],
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

sub _copyable_properties {
    # TODO: move this into a more central place.
    # The ::Default class needs it even for things which are not of this subclass.
    return Genome::SoftwareResult::Default::_copyable_properties(@_);
}

1;

