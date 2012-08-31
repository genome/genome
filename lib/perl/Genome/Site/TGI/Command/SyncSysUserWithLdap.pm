package Genome::Site::TGI::Command::SyncSysUserWithLdap;

use strict;
use warnings;

use Genome;

class Genome::Site::TGI::Command::SyncSysUserWithLdap{
    is => 'Command::V2',
};

sub execute {
    my $self = shift;

    my $ldap_user = get_ldap_users();
    my @db_users = Genome::Sys::User->fix_params_and_get();

    my $changes = get_changes($ldap_user,\@db_users);
    my $creates = $changes->{'create'};  
    my $deletes = $changes->{'deletes'};  

    my $create_count = $creates ? scalar(@$creates) : 0;
    my $delete_count = $deletes ? scalar(@$deletes) : 0;
    my $changes_count = $create_count + $delete_count;

    if ($changes_count < 1) {
        $self->status("No differences found between database and ldap...exiting.\n");
        return 1;
    }

    if ($changes_count > 10) {
        print "Too many changes ($create_count creates, $delete_count deletes, $changes_count total). Sync manually if this is OK.";
        return;
    }

    for my $u (@{ $changes->{'create'} }) {
        $self->status("CREATE: " . $u->{'mail'} . "\n");
        Genome::Sys::User->create(
            email => $u->{'mail'},
            name => $u->{'cn'},
            username => $u->{'uid'},
        );
    } 

    for my $u (@{ $changes->{'delete'} }) {
        $self->status("DELETE: " . $u->email . "\n");
        $u->delete();
    } 

    $self->status("done- $create_count creates, $delete_count deletes, $changes_count total\n");
}

sub get_ldap_users {
    my $ldap_user = {};
    my @c = `ldapsearch -z 1000 -x`; # gets max 1000 records
    my @users;
    my $user;

    # go through each line of ldapsearch output
    for my $c (@c) {
        next if $c =~ /^\#/;
        chomp($c);

        if ($c =~ /^$/) {
            # process/destroy user
            push @users, $user;
            undef $user;
        } else {
            my ($key, $value) = split(/\:\s+/,$c);
            $user->{$key} = $value;       
        }
    }

    # filter out users that didnt have email address entry in ldap
    for my $u (@users) {
        next if !$u->{'mail'};
        $ldap_user->{$u->{'mail'}} = $u;
    }

    return $ldap_user;
}


sub get_changes {
    my ($ldap_user, $db_users) = @_;
    my @db_users = @$db_users;
    # email is called email in db, mail in ldap
    my $changes = {};
    my $db_user = {};

    # look for people in db but not ldap
    for my $u (@db_users) {

        my $email = $u->email();
        if (!$ldap_user->{$email}) {
            push @{$changes->{'delete'}}, $u;
        } else {
            $db_user->{$email} = $u;
        }
    }

    # look for people in ldap but not db
    for my $mail (keys %$ldap_user) {
        my $u = $ldap_user->{$mail};

        if (!$db_user->{$mail}) {
            push @{$changes->{'create'}}, $u;
        }
    }

    return $changes;
}

sub status {
    # Suppress status when running in cron.
    if ( -t STDOUT ) {
        print @_;
    }
}

1;
