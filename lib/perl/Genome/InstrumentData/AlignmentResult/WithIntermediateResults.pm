package Genome::InstrumentData::AlignmentResult::WithIntermediateResults;

use Genome;

use warnings;
use strict;

class Genome::InstrumentData::AlignmentResult::WithIntermediateResults {
    is => 'Genome::InstrumentData::AlignmentResult',
};

sub intermediate_result_user_label { return 'intermediate result'; }

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return if not $self;

    my $rm_iar_ok = $self->remove_intemediate_results;
    return if not $rm_iar_ok;

    return $self;
}

sub get_or_create_intermediate_result_for_params {
    my ($self, $params) = @_;

    Carp::confess('No params to generate intermediate result!') if not $params;
    Carp::confess('No params to generate intermediate result!') if not ref($params) eq 'HASH';
    Carp::confess('No params to generate intermediate result!') if not %$params;

    my $includes = join(' ', map { '-I ' . $_ } UR::Util::used_libs);
    my $parameters = join(', ', map($_ . ' => "' . (defined($params->{$_}) ? $params->{$_} : '') . '"', sort keys %$params));

    if(UR::DBI->no_commit()) {
        my $rv = eval "Genome::InstrumentData::IntermediateAlignmentResult->get_or_create($parameters);";
        if(!$rv or $@) {
            my $err = $@;
            die('Failed to generate intermediate result!' . ($err? $err : ' command returned false') );
        }
    } else {
        my $cmd = qq{$^X $includes -e 'use above "Genome"; Genome::InstrumentData::IntermediateAlignmentResult->get_or_create($parameters); UR::Context->commit;' };
        my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd) };
        if(!$rv or $@) {
            my $err = $@;
            die('Failed to generate intermediate result!' . ($err? $err : ' command returned false') );
        }
    }

    my $intermediate_result = Genome::InstrumentData::IntermediateAlignmentResult->get_with_lock(%$params);
    if ( not $intermediate_result ) {
        Carp::confess( $self->error_message("Failed to generate intermediate result!") );
    }

    if ( not $intermediate_result->users(user => $self, label => intermediate_result_user_label()) ) {
        $intermediate_result->add_user(user => $self, label => intermediate_result_user_label());
    }

    return $intermediate_result;
}

sub remove_intemediate_results {
    my $self = shift;

    my @iar_users = Genome::SoftwareResult::User->get(user => $self, label => intermediate_result_user_label);
    for my $iar_user ( @iar_users ) {
        my $iar = $iar_user->software_result;
        $iar_user->delete;
        $iar->delete;
    }

    return 1;
}

1;

