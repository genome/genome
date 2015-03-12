package Genome::Test::Factory::InstrumentData::AlignmentResult;
use base qw(Genome::Test::Factory::Base);

use strict;
use warnings;

use Genome;
use Genome::Test::Factory::SoftwareResult::User;

sub generate_obj {
    my $self = shift;

    my $result_users = Genome::Test::Factory::SoftwareResult::User->setup_user_hash();
    return Genome::InstrumentData::AlignmentResult::Bwa->__define__(@_, _user_data_for_nested_results => $result_users);
}


1;
