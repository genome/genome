package Genome::SoftwareResult::WithNestedResults;

use strict;
use warnings;

use Genome;

class Genome::SoftwareResult::WithNestedResults {
    is => 'Genome::SoftwareResult',
    is_abstract => 1,
    has_transient_optional => [
        _user_data_for_nested_results => {
            is => 'HASH',
            doc => 'user information to pass along to any nested results',
        },
    ],
};

sub _preprocess_params_for_callback {
    my $class = shift;

    my %params = @_;

    unless(exists $params{_user_data_for_nested_results}) {
        if(exists $params{users} and ref($params{users}) eq 'HASH') {
            my %users = %{$params{users}};
            $params{_user_data_for_nested_results} = \%users;
        } else {
            $class->warning_message('No user information found to use for nested results.');
        }
    }

    return %params;
}

1;

