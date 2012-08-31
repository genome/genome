
# sync genome_sys_user from ldap

#UR_DBI_NO_COMMIT=1

use Genome;
use Data::Dumper;

my @c = `ldapsearch -x`;

my @users;
my $user;

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

# called email in db, mail in ldap

my $ldap_user = {};
for my $u (@users) {
    next if !$u->{'mail'};
    $ldap_user->{$u->{'mail'}} = $u;
}

my @db_users = Genome::Sys::User->fix_params_and_get();

print scalar(@db_users);
print 

my @changes;

my $db_user = {};
for my $u (@db_users) {

    my $email = $u->email();
    if (!$ldap_user->{$email}) {
        push @changes, '- ' . $email;        
        $u->delete();
    } else {
        $db_user->{$email} = $u;
    }
}

for my $mail (keys %$ldap_user) {
    my $u = $ldap_user->{$mail};

    if (!$db_user->{$mail}) {
        Genome::Sys::User->create(  
            email => $u->{'mail'},
            name => $u->{'cn'}
        );
        push @changes, '+ ' . $u->{'mail'};
    }
}

print Dumper \@changes;

UR::Context->commit();



