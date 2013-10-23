package Genome::Model::Tools::DetectVariants2::VarscanBase;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::DetectVariants2::VarscanBase {
    is_abstract => 1,
    is => ['Genome::Model::Tools::DetectVariants2::Detector'],
    has => [
    ],
};

# This method breaks the samtools params down into individual params for the varscan modules
sub _process_samtools_params {
    my ($self, $params) = @_;
    my $samtools_version;
    my $use_baq = 1;

    # Grab version if it exists
    if ($params =~ m/version/) {
        ($samtools_version) = ($params =~ m/--version\s*(\S+)/);
        $params =~ s/--version\s*(\S+)\s*//;
    }

    # Grab baq boolean
    if ($params =~ m/--nobaq/) {
        $use_baq = 0;
        $params =~ s/--nobaq\s*//;
    }

    return ($samtools_version, $use_baq, $params);
}

# Params should be set up as <samtools params>:<varscan params>
# If there is no : present, assume everything is varscan params (legacy processing profiles)
sub _split_params {
    my ($self, $params) = @_;
    my ($samtools_params, $varscan_params);
    if ($params =~ m/:/) {
        ($samtools_params, $varscan_params) = split ":", $params;
    } else {
        $varscan_params = $params;
    }
    return ($samtools_params, $varscan_params);
}

sub has_version {
    my $self = shift;
    my $version = shift;
    unless(defined($version)){
        $version = $self->version;
    }
    my @versions = Genome::Model::Tools::Varscan->available_varscan_versions;
    for my $v (@versions){
        if($v eq $version){
            return 1;
        }
    }
    return 0;  
}

1;
