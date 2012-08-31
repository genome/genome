package Genome::ProcessingProfile::Command::List;

use strict;
use warnings;
use Genome;
use Genome::ProcessingProfile;

class Genome::ProcessingProfile::Command::List {
    is => 'Command::SubCommandFactory',
    doc => 'list processing profiles by type'
};

sub _sub_commands_from { 'Genome::Model' }
sub _target_base_class { 'Genome::ProcessingProfile' }
sub _sub_commands_inherit_from { 'Genome::ProcessingProfile::Command::List::Base' }

sub _command_subclass {
    my ($self, $class_name) = @_;
    my ($subclass_name) = $class_name =~ /^Genome::ProcessingProfile::Command::List::(.*)$/;
    return $subclass_name;
}

sub _processing_profile_class {
    my ($self, $class_name) = @_;
    my $subclass_name = $self->_command_subclass($class_name);
    my $processing_profile_class = join('::', 'Genome::ProcessingProfile', $subclass_name);
    return $processing_profile_class;
}

sub _profile_string {
    my ($self, $class_name) = @_;
    my $subclass_name = $self->_command_subclass($class_name);
    my $profile_string = Genome::Utility::Text::camel_case_to_string($subclass_name);
    $profile_string =~ s/:://;
    return $profile_string;
}

sub _build_all_sub_commands {
    my $self = shift;

    # make the default commands per processing-profile
    my @subclasses = $self->SUPER::_build_all_sub_commands(@_);

    # add an extra one
    if ($self->class eq __PACKAGE__) {
        class Genome::ProcessingProfile::Command::List::All {
            is => 'Genome::ProcessingProfile::Command::List::Base',
            has => [
                subject_class_name  => {
                    is_constant => 1, 
                    value => 'Genome::ProcessingProfile' 
                },
                show => { default_value => 'id,type_name,name' },
            ],
            doc => 'list processing profiles'
        };

        # return the whole list
        push @subclasses, 'Genome::ProcessingProfile::Command::List::All';
    }

    return @subclasses;
}

sub _build_sub_command {
    my ($self,
        $class_name,            # Genome::ProcessingProfile::Command::List::Foo
        @inheritance,
    ) = @_;

    # the default _build_all_sub_commands() in the super class calls this

    my $pp_class_name = $self->_processing_profile_class($class_name);

    # params start with id and name, then sequencing_platform,
    # and then are in alpha order
    my $params_list = 'id,name' . 
        join('', 
            map { ",$_" } 
            sort { 
                ($a eq 'sequencing_platform' ? 1 : 2)
                cmp
                ($b eq 'sequencing_platform' ? 1 : 2)
            }
            $pp_class_name->params_for_class
        );

    # write the custom command for this processing profile
    class {$class_name} { 
        is => \@inheritance, 
        has => [
            subject_class_name => { 
                is_constant => 1, 
                value => $pp_class_name
            },
            show => { default_value => $params_list },
        ],
        doc => "list " . $self->_profile_string($class_name) . " processing profiles",
    };

    return $class_name;
}

1;

#$HeadURL$
#$Id$
