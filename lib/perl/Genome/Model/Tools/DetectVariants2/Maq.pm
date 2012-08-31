package Genome::Model::Tools::DetectVariants2::Maq;

use strict;
use warnings;

use Genome;
use IO::File;

class Genome::Model::Tools::DetectVariants2::Maq {
    is => ['Genome::Model::Tools::DetectVariants2::Detector'],
};

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt detect-variants2 maq
EOS
}

sub help_detail {
    return <<EOS 
This tool runs maq for detection of SNVs and/or indels.
EOS
}

sub _detect_variants {
    my $self = shift;
    
    $self->error_message('MAQ support is not currently implemented in this version of the DetectVariants API.');
    die $self->error_message;
}



sub has_version {
    my $self = shift;
    my $version = shift;
    unless(defined($version)){
        $version = $self->version;
    }
    my @versions = Genome::Model::Tools::Maq->available_maq_versions;
    for my $v (@versions){
        if($v eq $version){
            return 1;
        }
    }
    return 0;
}

1;
