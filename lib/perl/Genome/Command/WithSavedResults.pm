package Genome::Command::WithSavedResults;
use Genome::SoftwareResult::Default;
use strict;
use warnings;

class Genome::Command::WithSavedResults {
    is => 'Command::V2',
    is_abstract => 1,
    type_has => [
        parallelize_by => { 
            is => 'ARRAY', is_optional => 1,
            doc => 'produce intermediate results and merge, grouping by this/these attributes' },     
    ],
    has => [
        result_version => {
            is_param => 1,
            is => 'Integer',
            is_abstract => 1,
            doc => 'the version of results, which may iterate as this logic iterates'
        },
        output_dir  => { 
            is_output => 1,
            is => 'FilesystemPath', 
            is_optional => 1,
            doc => 'override the output directory' 
        },
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
        die "error initializing $subclass_name from " . __PACKAGE__ . ": $@";
    }

    my @problems;
    my $meta = $subclass_name->__meta__;
    my $parallelize_by = $meta->parallelize_by;

    my $result_version_meta = $meta->property("result_version");
    unless ($result_version_meta) {
        push @problems, "$subclass_name should implement result_version, typically with a default_value of '1'";
    }

    my $versions = $result_version_meta->valid_values;
    for my $version (@$versions) {
        my $method1 = "_execute_v$version";
        unless ($subclass_name->can($method1)) {
            push @problems, "no method $method1 found, though a valid result_version includes '$version'";
        }
        # note this code doesn't work because the added parallelize_property isn't visible when this callback runs 
        if ($parallelize_by and @$parallelize_by) {
            my $method2 = "_merge_v$version";
            unless ($subclass_name->can($method2)) {
                push @problems, "no method $method2 found, though a valid result_version includes '$version'";
            }
        }
    }

    if (@problems) {
        for (@problems) { $subclass_name->error_message($_)  }
        die "error defining $subclass_name";
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
            delete $props{output_dir};
            delete $props{result};

            # TODO: compose a workflow here instead of a linear run
            my $n = 0; 
            for my $value (@values) {
                $n++;
                print "Breakdown $n: " . $value->__display_name__ . "\n";
                $props{$prop} = [$value];
                my $partial = $self->class->create(%props);
                unless ($partial) {
                    die "failed to create partial for $prop " . $value->__display_name__;
                }
                push @run_commands, $partial;
            }
            my @intermediate_results;
            for my $cmd (@run_commands) {
                my $retval = $cmd->execute(@_);
                my $result = $cmd->result();
                push @intermediate_results, $result;
            } 

            # TODO: this should be the final merging step in the workflow
            my $merge_method = "_merge_v$result_version";
            my $retval;
            if ($self->can($merge_method)) {
                $retval = $self->$merge_method(@intermediate_results);
            }
            else {
                $retval = $self->_default_merge(@intermediate_results);
            }
            my $merge_result = $self->result;
            for my $intermediate_result (@intermediate_results) {
                $intermediate_result->add_user(
                    label => 'composes',
                    user => $merge_result,
                );
            }
            return $retval;
        }
    }
   
    return $self->$method(@_);
}

sub _default_merge {
    my $self = shift;
    my @underlying = @_;
    
    my $subdir = $self->output_dir . '/underlying_results';
    Genome::Sys->create_directory($subdir);

    for my $r (@underlying) {
        my $sub_subdir = $r->id;
        my $path = $subdir . '/' . $sub_subdir;
        Genome::Sys->create_symlink($r->output_dir, $path);
    }

    return 1;
}

1;

