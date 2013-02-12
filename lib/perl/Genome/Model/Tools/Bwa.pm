package Genome::Model::Tools::Bwa;

use strict;
use warnings;

use Genome;
use File::Basename;

my $DEFAULT = '0.5.9';

class Genome::Model::Tools::Bwa {
    is => 'Command',
    has => [
        use_version => { is => 'Version', is_optional => 1, default_value => $DEFAULT, doc => "Version of bwa to use, default is $DEFAULT" },
    ],
};

sub sub_command_sort_position { 12 }

sub help_brief {
    "Tools to run BWA or work with its output files.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
genome-model tools bwa ...    
EOS
}

sub help_detail {                           
    return <<EOS 
More information about the BWA suite of tools can be found at http://bwa.sourceforege.net.
EOS
}


my %BWA_VERSIONS = (
    '0.4.2'        => $ENV{GENOME_SW} . '/bwa/bwa-0.4.2-64/bwa',
    '0.4.9'        => $ENV{GENOME_SW} . '/bwa/bwa-0.4.9-64/bwa',
    '0.5.0'        => $ENV{GENOME_SW} . '/bwa/bwa-0.5.0-64/bwa',
    '0.5.1'        => $ENV{GENOME_SW} . '/bwa/bwa-0.5.1-64/bwa',
    '0.5.2'        => $ENV{GENOME_SW} . '/bwa/bwa-0.5.2-64/bwa',
    '0.5.3'        => $ENV{GENOME_SW} . '/bwa/bwa-0.5.3-64/bwa',
    '0.5.4'        => $ENV{GENOME_SW} . '/bwa/bwa-0.5.4-64/bwa',
    '0.5.5'        => $ENV{GENOME_SW} . '/bwa/bwa-0.5.5-64/bwa',
    '0.5.6'        => $ENV{GENOME_SW} . '/bwa/bwa-0.5.6-64/bwa',
    '0.5.7'        => $ENV{GENOME_SW} . '/bwa/bwa-0.5.7-64/bwa',
    '0.5.7-6'      => $ENV{GENOME_SW} . '/bwa/bwa-0.5.7-6-64/bwa',
    '0.5.8a'       => $ENV{GENOME_SW} . '/bwa/bwa-0.5.8a-64/bwa',
    '0.5.8c'       => $ENV{GENOME_SW} . '/bwa/bwa-0.5.8c-64/bwa',
    '0.5.9rc1'     => $ENV{GENOME_SW} . '/bwa/bwa-0.5.9rc1-64/bwa',
    '0.5.9'        => $ENV{GENOME_SW} . '/bwa/bwa-0.5.9-64/bwa',
    '0.5.9-pem0.1' => '/usr/bin/bwa-0.5.9-pem0.1',
    '0.5.9-i0.3'   => '/usr/bin/ibwa-0.5.9-0.3',
    '0.5.9-i0.4'   => '/gscuser/tabbott/bin/ibwa-0.5.9-0.4',
    '0.5.9-i0.5'   => '/gscuser/tabbott/bin/ibwa-0.5',
    '0.6.2'        => '/gscuser/iferguso/bin/bwa-0.6.2',
    'bwa'          => 'bwa',
);


sub bwa_path {
    my $self = $_[0];
    return $self->path_for_bwa_version($self->use_version);
}

sub available_bwa_versions {
    my $self = shift;
    return keys %BWA_VERSIONS;
}

sub path_for_bwa_version {
    my $class = shift;
    my $version = shift;

    if (defined $BWA_VERSIONS{$version}) {
        return $BWA_VERSIONS{$version};
    }
    die('No path for bwa version '. $version);
}

sub default_bwa_version {
    die "default samtools version: $DEFAULT is not valid" unless $BWA_VERSIONS{$DEFAULT};
    return $DEFAULT;
}
        
sub default_version { return default_bwa_version; }

sub supports_bam_input {
    my $class = shift;
    my $version = shift;

    my %ok_versions = (
       '0.5.9rc1'     => 1,
       '0.5.9'        => 1,
       '0.5.9-pem0.1' => 1,
       '0.5.9-i0.3'   => 1,
       '0.5.9-i0.4'   => 1,
    );

    return (exists $ok_versions{$version});

}

sub supports_multiple_reference {
    my $class = shift;
    my $version = shift;

    my %ok_versions = (
        '0.5.9-pem0.1' => 1,
        '0.5.9-i0.3'   => 1,
        '0.5.9-i0.4'   => 1,
    );

    return exists $ok_versions{$version};
}

sub supports_bwasw {
    my $class = shift;
    my $version = shift;

    # Although bwa 0.5.9 technically "supports" bwasw we don't want to use it
    # because it doesn't correctly set the flags and mate information in the
    # SAM output.
    my %ok_versions = (
        '0.6.1' => 1,
        '0.6.2' => 1,
    );

    return exists $ok_versions{$version};
}

1;

