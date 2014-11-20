package Genome::Env::GENOME_NO_REQUIRE_USER_VERIFY;
our $VERSION = $Genome::VERSION;

if(defined $ENV{GENOME_NO_REQUIRE_USER_VERIFY}) {
    warn 'The GENOME_NO_REQUIRE_USER_VERIFY environment variable no longer has any effect.  Use UR_NO_REQUIRE_USER_VERIFY instead.';
}

1;
