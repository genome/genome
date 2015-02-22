package Genome::Model::Tools::BedTools;

use strict;
use warnings;

use Genome;

my $BEDTOOLS_DEFAULT = '2.9.0';

class Genome::Model::Tools::BedTools {
    is  => 'Command',
    is_abstract => 1,
    has_input => [
        use_version => {
            is  => 'Version',
            doc => 'BEDTools version to be used.  default_value='. $BEDTOOLS_DEFAULT,
            is_optional   => 1,
            default_value => $BEDTOOLS_DEFAULT,
        },
    ],
};


sub help_brief {
    "Tools to run BedTools.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt bed-tools ...
EOS
}

sub help_detail {
    return <<EOS
More information about the BedTools suite of tools can be found at http://code.google.com/p/bedtools/.
EOS
}

my %BEDTOOLS_VERSIONS = (
    '2.17.0' => $ENV{GENOME_SW} . '/bedtools/BEDTools-2.17.0',
    '2.16.2' => $ENV{GENOME_SW} . '/bedtools/BEDTools-2.16.2',
    '2.14.3' => $ENV{GENOME_SW} . '/bedtools/BEDTools-2.14.3',
    '2.9.0' => $ENV{GENOME_SW} . '/bedtools/BEDTools-2.9.0',
    '2.8.3' => $ENV{GENOME_SW} . '/bedtools/BEDTools-2.8.3',
    '2.6.1' => $ENV{GENOME_SW} . '/bedtools/BEDTools-2.6.1',
    '2.5.4' => $ENV{GENOME_SW} . '/bedtools/BEDTools-2.5.4',
    '2.3.2' => $ENV{GENOME_SW} . '/bedtools/BEDTools-2.3.2',
);

sub path_for_bedtools_version {
    my ($class, $version) = @_;
    $version ||= $BEDTOOLS_DEFAULT;
    my $path = $BEDTOOLS_VERSIONS{$version};
    if (Genome::Config->arch_os =~ /64/) {
        if ($path) {
            my $arch_path = $path .'-64';
            if (-d $arch_path) {
                $path = $arch_path;
            }
        }
    }
    return $path if (defined $path && -d $path);
    die 'No path found for bedtools version: '. $version;
}

sub default_bedtools_version {
    die "default bedtools version: $BEDTOOLS_DEFAULT is not valid" unless $BEDTOOLS_VERSIONS{$BEDTOOLS_DEFAULT};
    return $BEDTOOLS_DEFAULT;
}

sub bedtools_path {
    my $self = shift;
    return $self->path_for_bedtools_version($self->use_version);
}

1;
