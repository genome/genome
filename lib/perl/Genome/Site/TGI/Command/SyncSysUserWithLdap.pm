package Genome::Site::TGI::Command::SyncSysUserWithLdap;

use strict;
use warnings;

use Genome;
use IO::File;
use Net::LDAP;

class Genome::Site::TGI::Command::SyncSysUserWithLdap{
    is => 'Command::V2',
    has => {
        max_changes_allowed => {
            is => 'Integer',
            default_value => 10,
            doc => 'The maximum number fo changes allowed to process the creates and deletes.',
        },
    },
    doc => 'Sync Genome sys users from LDAP users',
};

sub execute {
    my $self = shift;

    my $ldap_users = get_ldap_users();
    my @db_users = Genome::Sys::User->fix_params_and_get();

    my ($creates, $deletes) = get_changes($ldap_users,\@db_users);
    my $changes_count = @$creates + @$deletes;
    $self->status_message("Creates: %s", scalar(@$creates));
    $self->status_message("Deletes: %s", scalar(@$deletes));
    $self->status_message("Total: %s", $changes_count);
    $self->status_message('Max changes allowed: %s', $self->max_changes_allowed);
    if ($changes_count < 1) {
        $self->status_message("No differences found between database and ldap.");
        return 1;
    }
    elsif ($changes_count > $self->max_changes_allowed) {
        $self->status_message( "The number of changes exceeds the max number allowed. If this is expected, increase the --max-changes-allowed option to be higher than the number of changes to process.");
        return;
    }

    for my $u (@$creates) {
        my $email = $u->get_value('mail');
        $self->status_message("CREATE: $email\n");
        Genome::Sys::User->create(
            email => $email,
            name => $u->get_value('cn'),
            username => $u->get_value('uid'),
        );
    }

    for my $u (@$deletes) {
        $self->status_message("DELETE: " . $u->email . "\n");
        $u->delete();
    }

    $self->status_message('Done');
}

sub get_ldap_users {
    my $ldap = Net::LDAP->new('ipa1.gsc.wustl.edu', version => 3);
    my $mesg = $ldap->start_tls(verify => 'none');
    $mesg->code && die $mesg->error;

    $mesg = $ldap->bind;
    $mesg->code && die $mesg->error;

    $mesg = $ldap->search(
        base => "dc=gsc,dc=wustl,dc=edu",
        filter => "(objectClass=Person)",
        #filter => "(&(objectClass=Person)(uid=$username))",
    );
    $mesg->code && die $mesg->error;

    my %ldap_users;
    foreach my $ldap_user ($mesg->entries) {
        my $mail = $ldap_user->get_value('mail');
        next if not $mail; # skip if no email address
        $ldap_users{$mail} = $ldap_user;
        #$ldap_user->dump;
    }

    $mesg = $ldap->unbind;   # take down session
    $mesg->code && die $mesg->error;

    return \%ldap_users
}

sub get_changes {
    my ($ldap_users, $db_users) = @_;
    my @db_users = @$db_users;
    # email is called email in db, mail in ldap
    my $db_user = {};
    my (@creates, @deletes);

    # look for people in db but not ldap
    for my $u (@db_users) {

        my $email = $u->email();
        if (!$ldap_users->{$email}) {
            push @deletes, $u;
        } else {
            $db_user->{$email} = $u;
        }
    }

    # look for people in ldap but not db
    for my $mail (keys %$ldap_users) {
        my $u = $ldap_users->{$mail};

        if (!$db_user->{$mail}) {
            push @creates, $u;
        }
    }

    return ( \@creates, \@deletes );
}

1;

