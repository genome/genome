package Genome::Config::Command::Grep;

use strict;
use warnings;

use Genome;

class Genome::Config::Command::Grep {
    is => 'Command::V2',
    has => [
        query => {
            is => 'Text',
            doc => 'The string to find in extant configurations',
            shell_args_position => 1,
        },
    ],
    has_optional => [
        show => {
            is => 'Text',
            doc => 'Properties of the matching config items to list',
            default_value => 'id,analysis_project,file_path',
        },
        menu_item => {
            is => 'Genome::Config::AnalysisMenu::Item',
            doc => 'restrict results to those matching this menu item',
        },
        only_concrete => {
            is => 'Boolean',
            doc => 'only include configurations that are either custom or have actually been used',
            default_value => 1,
        },
    ],
    has_transient_optional_output => [
        matching_items => {
            is => 'Genome::Config::Profile::Item',
            is_many => 1,
            doc => 'Those items that matched the search',
        },
    ],
};

sub execute {
    my $self = shift;

    my @params;

    if($self->menu_item) {
        push @params, 'analysismenu_item_id', $self->menu_item->id;
    }

    my @all_profile_items = Genome::Config::Profile::Item->get(@params);

    my @matching_profile_items;
    for my $item (@all_profile_items) {
        push @matching_profile_items, $item if $self->matches($item);
    }

    unless(@matching_profile_items) {
        $self->status_message("No matching items found.");
        return 1;
    }

    $self->matching_items(\@matching_profile_items);

    return UR::Object::Command::List->execute(
        subject_class_name => 'Genome::Config::Profile::Item',
        filter => 'id in [' . join(',', map $_->id, @matching_profile_items) . ']',
        show => $self->show,
    );
}

sub matches {
    my $self = shift;
    my $item = shift;

    return if $self->only_concrete and not $item->_is_concrete;

    my $query = $self->query;

    my $file_path = $item->file_path;

    return `grep $query $file_path`;
}

1;
