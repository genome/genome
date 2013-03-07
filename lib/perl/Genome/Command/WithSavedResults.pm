package Genome::Command::WithSavedResults;
use Genome::SoftwareResult::Default;

class Genome::Command::WithSavedResults {
    is => 'Command::V2',
    is_abstract => 1,
    type_has => [
        parallelize_by => { 
            is => 'ARRAY', is_optional => 1,
            doc => 'produce intermediate results and merge, grouping by this/these attributes' },     
    ],
    has_optional_input => [
        #output_dir  => { 
        #    is => 'FilesystemPath', 
        #    doc => 'override the output directory' 
        #},
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

    my $meta = $subclass_name->__meta__;
    unless ($meta->property("result_version")) {
        die "$subclass_name should implement result_version, typically with a default_value of '1'";
    }

    return 1;
}

sub _copyable_properties {
    # TODO: move this into a more central place.
    # The ::Default class needs it even for things which are not of this subclass.
    return Genome::SoftwareResult::Default::_copyable_properties(@_);
}

sub execute {
    my $self = shift;
    
    my $result_version = $self->result_version;
    my $method = "_execute_v$result_version";
    $method =~ s/\./_/g;
    unless ($self->can($method)) {
        die "no implementation ($method) for version $result_version!";
    }

    my @run_commands;
    my $parallelize_by = $self->__meta__->parallelize_by;
    if ($parallelize_by and @$parallelize_by) {
        if (@$parallelize_by > 1) {
            die "support for multiplexed paralleized_by not implimented!: @$parallelize_by";
        }
        my $prop = $parallelize_by->[0];
        my @values = $self->$prop;
        if (@values > 1) {
            my %props = $self->_copyable_properties($self->class);
            for my $value (@values) {
                print "breakdown for value " . $value->__display_name__ . "\n";
                $props{$prop} = [$value];
                my $partial = $self->class->create(%props);
                unless ($partial) {
                    die "failed to create partial for $prop " . $value->__display_name__;
                }
                push @run_commands, $partial;
            }
            
            print "run commands: @run_commands\n";
            for my $cmd (@run_commands) {
                $cmd->$method(@_);
            } 
            die "no merge logic yet!";
        }
    }
   
    return $self->$method(@_);
}

1;

