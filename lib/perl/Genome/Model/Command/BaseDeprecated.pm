package Genome::Model::Command::BaseDeprecated;

#:eclark 11/17/2009 Code review.

# get_model_class* methods at the bottom should be in Genome::Model, not here.
# create_directory and bsub_rusage probably don't even belong in this class

use strict;
use warnings;

use Genome;
use Carp;
use Regexp::Common;

class Genome::Model::Command::BaseDeprecated {
    is => 'Command::V1',
    is_abstract => 1,
    has => [
        model => {
            is => 'Genome::Model',
            id_by => 'model_id',
        },
        model_name => {
            is => 'Text',
            via => 'model',
            to => 'name',
        },
        name_pattern => {
            is => 'Text',
            shell_args_position => 99,
            is_optional => 1,
        },
    ],
    doc => "base class for commands that work with genome models",
};

sub command_name {
    my $class = ref($_[0]) || $_[0];
    return $class->SUPER::command_name unless $class eq __PACKAGE__;
    return 'genome model';
}

sub command_name_brief {
    my $class = ref($_[0]) || $_[0];
    return $class->SUPER::command_name_brief unless $class eq __PACKAGE__;
    return 'model';
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_)
        or return;

    if ( defined $self->model_id ) {
        unless ( $self->_verify_model ) {
            $self->delete;
            return;
        }
    }

    unless ($class->__meta__->property_meta_for_name("model")->is_optional or $self->model) {
        if ($self->name_pattern) {
            my $pattern = $self->name_pattern;
            if ($pattern) {
                my @models = Genome::Model->get(name => { operator => "like", value => '%' . $pattern . '%' });
                if (@models >1) {
                    $self->error_message(
                                         "No model specified, and multiple models match pattern \%${pattern}\%!\n"
                                         . join("\n", map { $_->name } @models)
                                         . "\n"
                                     );
                    $self->delete;
                    return;
                }
                elsif (@models == 1) {
                    $self->model($models[0]);
                }
            } else {
                # continue, the developer may set this value later...
            }
        } else {
            $self->error_message("No model or name_pattern exists");
            $self->delete;
            return;
        }
    }
    return $self;
}

sub _verify_model {
    my $self = shift;

    unless ( defined $self->model_id ) {
        $self->error_message("No model id given");
        return;
    }

    unless ( $self->model ) {
        $self->error_message( sprintf('Can\'t get model for id (%s) ', $self->model_id) );
        return;
    }

    return 1;
}

sub _sub_command_name_to_class_name_map{
    my $class = shift;
    return
        map { my ($type) = m/::(\w+)$/; $type => $_ }
                $class->sub_command_classes();
}

sub _get_sub_command_class_name{
    my $class = shift;
    my $sub_command_name = $class->sub_command_delegator(@_);
    unless ($sub_command_name) {
        # The subclassing column's value was probably undef, meaning this sub-command
        # should be skipped
        return;
    }

    # Does the sub-command exist?
    my %sub_command_types = $class->_sub_command_name_to_class_name_map();

    #this takes the db name of the sub class, 'foo bar' and turns it into a Class equivalent name FooBar
    my $key;
    my @words = split(/[-_\s]/, $sub_command_name);
    $key .= ucfirst $_ foreach @words;

    my $sub_command_type = $sub_command_types{$key};
    unless ($sub_command_type) {
        return;
    }

    return $sub_command_type;
}


sub sub_command_delegator {
    # This method is used by the mid-level (like ::Build::ReferenceAlignment::AlignReads modules
    # to return the right sub-sub-class like ::Build::ReferenceAlignment::AlignReads::Maq
    my($class,%params) = @_;

    if (not defined $params{'model_id'}) {
        return;
    }

    return unless $params{model_id} =~ /^$RE{num}{int}$/;

    my $model = Genome::Model->get(id => $params{'model_id'})
        or return;

    # Which property on the model will tell is the proper subclass to call?
    unless ($class->can('command_subclassing_model_property')) {
        #$class->error_message("class $class did not implement command_subclassing_model_property()");
        return;
    }
    my $subclassing_property = $class->command_subclassing_model_property();
    unless ($model->can($subclassing_property)) {
        $class->error_message("class $class command_subclassing_model_property() returned $subclassing_property, but that is not a property of a model");
        return;
    }

    my $value = $model->$subclassing_property;
    if ($value =~ m/^maq/) {
        return 'maq';
    } else {
        return $value;
    }

}

sub bsub_rusage {
    # override for tasks which require LSF resource requirements
    ''
}

sub create_temp_directory {
    my $self = shift;
    my $basename = shift;
    my $path = Genome::Sys->create_temp_directory($basename, @_)
        or die;
    $self->status_message("Created directory: $path");
    return $path;
}

sub create_directory {
    my ($self, $path) = @_;

    Genome::Sys->create_directory($path)
        or die;

    $self->status_message("Created directory: $path");

    return 1;
}

sub _ask_user_question {
    my $self = shift;
    my $question = shift;
    my $timeout = shift || 60;
    my $input;
    eval {
        local $SIG{ALRM} = sub { die "Failed to reply to question '$question' with in '$timeout' seconds\n" };
        $self->status_message($question);
        $self->status_message("Please reply: 'yes' or 'no'");
        alarm($timeout);
        chomp($input = <STDIN>);
        alarm(0);
    };
    if ($@) {
        $self->warning_message($@);
        return;
    }
    unless ($input =~ m/yes|no/) {
        $self->error_message("'$input' is an invalid answer to question '$question'");
        return;
    }
    return $input;
}

#< Models Classes and Subclasses >#
sub get_model_classes {
    my @classes = Genome::Sys::get_classes_in_subdirectory_that_isa(
        'Genome/Model',
        'Genome::Model',
    );

    unless ( @classes ) { # bad
        Carp::confess("No model subclasses found!");
    }

    return @classes;
}

sub get_model_subclasses {
    # should confess in get_model_classes
    return map { m#::([\w\d]+)$# } get_model_classes();
}

sub get_model_type_names {
    # should confess in get_model_classes
    return map { Genome::Utility::Text::camel_case_to_string($_, ' ') } get_model_subclasses();
}

1;
