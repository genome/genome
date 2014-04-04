package Genome::Config::AnalysisProject::Command::UpdateConfigToCurrent;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::Command::UpdateConfigToCurrent {
    is => 'Command::V2',
    has_input => [
        profile_item => {
            is => 'Genome::Config::Profile::Item',
            doc => 'Specified config file will no longer be used',
            shell_args_position => 1,
        }
    ],
};

sub help_brief {
    return 'update an AnP config item to the current default.';
}

sub help_detail {
    return 'Using a given a profile config item for an analysis project [AnP], find the current default via the analysis menu item. The profile item is concrete if a config file has been copied to it. This happens when models are first created for an AnP, or if a custom config is added to an AnP. If the profile item is not concrete [no config file], this command is a no-op. If concrete, a new config item [non concrete] will be created.'
}

sub execute {
    my $self = shift;

    my $profile_item = $self->profile_item;
    if ( not $profile_item->is_concrete ) {
        $self->status_message('This profile item has not yet been concretized, replacing it with the current defaults is a no-op.');
        return 1;
    }

    my $menu_item = $profile_item->analysis_menu_item;
    if ( not $menu_item ) {
        $self->error_message('No analysis menu item associated with profile item! It is required to get the current default profile item!');
        return;
    }

    my $new_profile_item = Genome::Config::Profile::Item->create(
        analysis_project => $profile_item->analysis_project,
        status => $profile_item->status,
        analysis_menu_item => $menu_item,
    );

    $profile_item->status('disabled');

    $self->status_message('Successfully updated config ('.$profile_item->__display_name__.') to current ('.$new_profile_item->__display_name__.').');

    return 1;
}

1;
