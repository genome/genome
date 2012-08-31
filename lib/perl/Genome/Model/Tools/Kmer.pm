package Genome::Model::Tools::Kmer;

use strict;
use warnings;

use Genome;

my $DEFAULT = '1.3.2';

class Genome::Model::Tools::Kmer {
    is => 'Command',
    is_abstract => 1,
    has => [
        use_version => {
            is  => 'Version',
            doc => "genometools version to be used, default is $DEFAULT.",
            is_optional   => 1,
            default_value => $DEFAULT,
        },
    ],
};

my %GENOMETOOLS_VERSIONS = (
    '1.3.2'    => $ENV{GENOME_SW} . '/genometools/genometools-1.3.2',
);

sub path_for_genometools_version {
    my ($class, $version) = @_;
    $version ||= $DEFAULT;
    my $path = $GENOMETOOLS_VERSIONS{$version};
    if (Genome::Config->arch_os =~ /64/) {
        if ($path) {
            $path .= '-64';
        }
    }
    return $path if (defined $path && -d $path);
    die 'No path found for genometools version: '.$version;
}

sub default_genometools_version {
    die "default genometools version: $DEFAULT is not valid" unless $GENOMETOOLS_VERSIONS{$DEFAULT};
    return $DEFAULT;
}

sub genometools_path {
    my $self = shift;
    return $self->path_for_genometools_version($self->use_version) .'/bin/gt';
}
