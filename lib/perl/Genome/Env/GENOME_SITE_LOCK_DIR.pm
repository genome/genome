package Genome::Env::GENOME_SITE_LOCK_DIR;

sub default_value { $ENV{GENOME_LOCK_DIR} || '/var/lock/genome/' }

1;
