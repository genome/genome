package Genome::Test::Factory::InstrumentData::AlignmentResult;
use base qw(Genome::Test::Factory::Base);

use strict;
use warnings;

use Genome;
use Genome::Test::Factory::SoftwareResult::User;

sub generate_obj {
    my $self = shift;
    my %params = @_;

    my $build = delete($params{'build'});
    
    my $result_users;
    unless ($build) {
        $result_users = Genome::Test::Factory::SoftwareResult::User->setup_user_hash();
    } else {
        $result_users = Genome::Test::Factory::SoftwareResult::User->setup_user_hash_with_build($build);
    }
    return Genome::InstrumentData::AlignmentResult::Bwa->__define__(%params, _user_data_for_nested_results => $result_users);
}


1;
