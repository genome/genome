package Genome::Test::Factory::SoftwareResult::ImportedFile;
use base qw(Genome::Test::Factory::Base);

use strict;
use warnings;

use Genome;
use Genome::Test::Factory::SoftwareResult::User;
use Digest::MD5;

our @required_params = qw(file_content_hash);

sub generate_obj {
    my $self = shift;
    my %params = @_;

    return Genome::SoftwareResult::ImportedFile->__define__(%params);
}

sub create_file_content_hash {
    return Digest::MD5::md5_base64( rand );
}

1;
