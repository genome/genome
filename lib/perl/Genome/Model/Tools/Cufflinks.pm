package Genome::Model::Tools::Cufflinks;

use strict;
use warnings;

use Genome;
use File::Basename;
use IPC::Cmd;

my $DEFAULT = '1.3.0';

class Genome::Model::Tools::Cufflinks {
    is => 'Command',
    has => [
        use_version => { is => 'Version', is_optional => 1, default_value => $DEFAULT, doc => "Version of cufflinks to use, default is $DEFAULT" },
    ],
};


sub help_brief {
    "Tools to run Cufflinks or work with its output files.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
genome-model tools cufflinks ...    
EOS
}

sub help_detail {
    return <<EOS
More information about the Cufflinks aligner can be found at http://cufflinks.cbcb.umd.edu/.
EOS
}

# this is for old versions, from 1.3.0 onward we use a dpkg and heuristic for path resolution
my %CUFFLINKS_VERSIONS = (
    '0.7.0'  => $ENV{GENOME_SW} . '/cufflinks/cufflinks-0.7.0.Linux_x86_64',
    '0.8.0'  => $ENV{GENOME_SW} . '/cufflinks/cufflinks-0.8.0.Linux_x86_64',
    '0.8.2'  => $ENV{GENOME_SW} . '/cufflinks/cufflinks-0.8.2.Linux_x86_64',
    '0.8.3'  => $ENV{GENOME_SW} . '/cufflinks/cufflinks-0.8.3.Linux_x86_64',
    '0.9.0'  => $ENV{GENOME_SW} . '/cufflinks/cufflinks-0.9.0.Linux_x86_64',
    '0.9.1'  => $ENV{GENOME_SW} . '/cufflinks/cufflinks-0.9.1.Linux_x86_64',
    '0.9.2'  => $ENV{GENOME_SW} . '/cufflinks/cufflinks-0.9.2.Linux_x86_64',
    '0.9.3'  => $ENV{GENOME_SW} . '/cufflinks/cufflinks-0.9.3.Linux_x86_64',
    '1.0.0'  => $ENV{GENOME_SW} . '/cufflinks/cufflinks-1.0.0.Linux_x86_64',
    '1.0.1'  => $ENV{GENOME_SW} . '/cufflinks/cufflinks-1.0.1.Linux_x86_64',
    '1.0.3'  => $ENV{GENOME_SW} . '/cufflinks/cufflinks-1.0.3.Linux_x86_64',
    '1.1.0'  => $ENV{GENOME_SW} . '/cufflinks/cufflinks-1.1.0.Linux_x86_64',
    '1.2.1'  => $ENV{GENOME_SW} . '/cufflinks/cufflinks-1.2.1.Linux_x86_64',
);

our %MIN_VERSIONS = (
    gtf_to_sam => '1.0.1',
    gffread => '1.0.1',
    cuffmerge => '1.0.0',
);

sub _executable_names {
    qw/gtf_to_sam gffread cuffmerge cuffcompare cuffdiff cufflinks/
}

sub available_cufflinks_versions {
    my $self = shift;
    my @legacy_versions = keys %CUFFLINKS_VERSIONS;
    my @local_versions = qx/ update_alternatives --list cufflinks /;
    @local_versions = map { /^.*cufflinks(.*)/; $1; } @local_versions;
    return @legacy_versions, @local_versions;
}

for my $app (_executable_names) {
    my $code = sub {
        my $self = $_[0];
        my $version = $self->use_version;
        my $parse_version = version->parse($version);
        
        if (my $min_version = $MIN_VERSIONS{$app}) {
            unless ($parse_version >= version->parse($min_version)) {
                die("$app command not available with version: $version");
            }
        }

        my $path = `which ${app}${parse_version}`;
        chomp $path;
        if ($path) {
            # all new packaged versions will fall through here, starting with version 1.3.0
            return $path;
        }
        else {
            # try the legacy method
            return $self->path_for_cufflinks_version($self->use_version) ."/$app";
        }
    };

    Sub::Install::install_sub({
        code => $code,
        into => __PACKAGE__,
        as => $app . '_path'
    });
}

# NOTE: this is only used for the legacy path generator above
sub path_for_cufflinks_version {
    my $class = shift;
    my $version = shift;

    my $path = IPC::Cmd::can_run("cufflinks" . $version);
    return $path if ($path);

    if (defined $CUFFLINKS_VERSIONS{$version}) {
        return $CUFFLINKS_VERSIONS{$version};
    }
    die('No path for cufflinks version '. $version);
}

sub default_cufflinks_version {
    die "default cufflinks version: $DEFAULT is not valid" unless $CUFFLINKS_VERSIONS{$DEFAULT};
    return $DEFAULT;
}

1;

