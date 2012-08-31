package Genome::ModelGroup::Command::Member::Tag;

use strict;
use warnings;

use Genome;

class Genome::ModelGroup::Command::Member::Tag {
    is => 'Genome::ModelGroup::Command::Member::Base', 
    has => [
        model_group => { 
            is => 'Genome::ModelGroup', 
            id_by => 'model_group_id',
            shell_args_position => 1,
            doc => 'model will only be tagged in this model-group',
        },
        model => {
            is => 'Genome::Model',
            id_by => 'genome_model_id',
            shell_args_position => 2,
            doc => 'model to be tagged',
        },
        tag_name => {
            is => 'Text',
            len => 255,
            shell_args_position => 3,
            is_optional => 1,
            doc => 'tags are key/value pairs- tag_name is the key',
        },
        tag_value => {
            is => 'Text',
            len => 255,
            shell_args_position => 4,
            default_value => 1,
            doc => 'tags are key/value pairs- tag_value is the value',
        }
    ],
};

sub help_synopsis {
    return <<"EOS"

Get list of tags for a model in a model group
genome model-group member tag --model-group 12 --model 34

Shortcut to do the command above:
genome model-group member tag 12 34

Set a tag named "mutant" to "1":
genome model-group member tag --model-group 12 --model 34 --tag-name mutant 

Set a tag named "user_liked" with a value of "jlolofie":
genome model-group member tag --model-group 56 --model 78 --tag-name user_liked --tag-value jlolofie

Shortcut to do the command above:
genome model-group member tag 56 78 user_liked jlolofie
EOS
}

sub help_brief {
    return <<EOS 
List or set tags for a model/model-group
EOS
}

sub help_detail {                           
    return <<EOS
You specify a model-group and a model in that group. If this is all your provide, a
list of existing tags is printed.

If you also supply a tag-name, that tag is set with default value "1".

Supply tag-value to explicitly set the value for that tag.

Tags refer to a model within a model-group. So, its possible to have a model tagged 
with  fav_color = yellow  in one model group but the same model tagged with
  fav_color = red  in a different model group.

[see also:  untag  ]
EOS
}

sub execute {

    my ($self) = @_;

    my $tag_name  = $self->tag_name();
    my $tag_value = $self->tag_value();
    my $model     = $self->model();
    my $model_id  = $model->genome_model_id();
    my $group     = $self->model_group();
    my $group_id  = $group->id();


    my $bridge = Genome::ModelGroupBridge->get(model_group_id => $group_id, model_id => $model_id);
    if (!$bridge) {
        die "ERROR: group $group_id doesnt contain model: $model_id";
    }

    if ($tag_name) {
        my $attr = Genome::MiscNote->create(
            subject     => $bridge,
            subject_id  => $bridge->id(),
            header_text => $tag_name,
            body_text   => $tag_value
        ) || die 'Error: couldnt create tag';

        print "Tag was set successfully:\ngroup: $group_id\nmodel: $model_id\n name: $tag_name\nvalue: $tag_value\n";
    } else {
        my $tags = $group->tags_for_model($model_id);
      
        if ($tags) {
            print "Existing tags for model-group $group_id, model $model_id:\n"; 
            for my $key (keys %$tags) {
                my $value = $tags->{$key};
                printf("%10s = %s\n", $key, $value);
            }
        } else {
            print "No tags exist for model-group $group_id, model $model_id:\n"; 
        }
    }

    return 1;
}

1;

