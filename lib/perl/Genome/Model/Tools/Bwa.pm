package Genome::Model::Tools::Bwa;

use strict;
use warnings;

use Genome;
use File::Basename;

my $DEFAULT = '0.5.9';

class Genome::Model::Tools::Bwa {
    is => 'Command',
    has => [
        use_version => {
            is => 'Version',
            is_optional => 1,
            default_value => $DEFAULT,
            doc => "Version of bwa to use, default is $DEFAULT"
        },
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


sub _lookup_version {
    my $class = shift;
    my $version = shift;

    my %bwa_versions = (
        'bwa' => {
            path       => 'bwa',
            features   => [],
            log_format => 'old',
            index_type => 'includes_reverse',
        },
        '0.4.2' => {
            path       => $ENV{GENOME_SW} . '/bwa/bwa-0.4.2-64/bwa',
            features   => [],
            log_format => 'old',
            index_type => 'includes_reverse',
        },
        '0.4.9' => {
            path       => $ENV{GENOME_SW} . '/bwa/bwa-0.4.9-64/bwa',
            features   => [],
            log_format => 'old',
            index_type => 'includes_reverse',
        },
        '0.5.0' => {
            path       => $ENV{GENOME_SW} . '/bwa/bwa-0.5.0-64/bwa',
            features   => [],
            log_format => 'old',
            index_type => 'includes_reverse',
        },
        '0.5.1' => {
            path       => $ENV{GENOME_SW} . '/bwa/bwa-0.5.1-64/bwa',
            features   => [],
            log_format => 'old',
            index_type => 'includes_reverse',
        },
        '0.5.2' => {
            path       => $ENV{GENOME_SW} . '/bwa/bwa-0.5.2-64/bwa',
            features   => [],
            log_format => 'old',
            index_type => 'includes_reverse',
        },
        '0.5.3' => {
            path       => $ENV{GENOME_SW} . '/bwa/bwa-0.5.3-64/bwa',
            features   => [],
            log_format => 'old',
            index_type => 'includes_reverse',
        },
        '0.5.4' => {
            path       => $ENV{GENOME_SW} . '/bwa/bwa-0.5.4-64/bwa',
            features   => [],
            log_format => 'old',
            index_type => 'includes_reverse',
        },
        '0.5.5' => {
            path       => $ENV{GENOME_SW} . '/bwa/bwa-0.5.5-64/bwa',
            features   => [],
            log_format => 'old',
            index_type => 'includes_reverse',
        },
        '0.5.6' => {
            path       => $ENV{GENOME_SW} . '/bwa/bwa-0.5.6-64/bwa',
            features   => [],
            log_format => 'old',
            index_type => 'includes_reverse',
        },
        '0.5.7' => {
            path       => $ENV{GENOME_SW} . '/bwa/bwa-0.5.7-64/bwa',
            features   => [],
            log_format => 'old',
            index_type => 'includes_reverse',
        },
        '0.5.7-6' => {
            path       => $ENV{GENOME_SW} . '/bwa/bwa-0.5.7-6-64/bwa',
            features   => [],
            log_format => 'old',
            index_type => 'includes_reverse',
        },
        '0.5.8a' => {
            path       => $ENV{GENOME_SW} . '/bwa/bwa-0.5.8a-64/bwa',
            features   => [],
            log_format => 'old',
            index_type => 'includes_reverse',
        },
        '0.5.8c' => {
            path       => $ENV{GENOME_SW} . '/bwa/bwa-0.5.8c-64/bwa',
            features   => [],
            log_format => 'old',
            index_type => 'includes_reverse',
        },
        '0.5.9rc1' => {
            path       => $ENV{GENOME_SW} . '/bwa/bwa-0.5.9rc1-64/bwa',
            features   => ['bam_input'],
            log_format => 'old',
            index_type => 'includes_reverse',
        },
        '0.5.9' => {
            path       => $ENV{GENOME_SW} . '/bwa/bwa-0.5.9-64/bwa',
            # Although bwa 0.5.9 technically "supports" bwasw, we don't want to use
            # it because it doesn't correctly set the flags and mate information in
            # the SAM output.
            features   => ['bam_input'],
            log_format => 'old',
            index_type => 'includes_reverse',
        },
        '0.5.9-pem0.1' => {
            path       => '/usr/bin/bwa-0.5.9-pem0.1',
            features   => ['bam_input', 'multiple_reference_mode'],
            log_format => 'old',
            index_type => 'includes_reverse',
        },
        '0.5.9-i0.3' => {
            path       => '/usr/bin/ibwa-0.5.9-0.3',
            features   => ['bam_input', 'multiple_reference_mode'],
            log_format => 'old',
            index_type => 'includes_reverse',
        },
        '0.5.9-i0.4' => {
            path       => '/gscuser/tabbott/bin/ibwa-0.5.9-0.4',
            features   => ['bam_input', 'multiple_reference_mode'],
            log_format => 'old',
            index_type => 'includes_reverse',
        },
        '0.5.9-i0.5' => {
            path       => '/gscuser/tabbott/bin/ibwa-0.5',
            # TODO verify features
            features   => ['bam_input', 'multiple_reference_mode'],
            log_format => 'old',
            index_type => 'includes_reverse',
        },
        '0.6.0' => {
            path       => $ENV{GENOME_SW} . '/bwa/bwa-0.6.0/bwa',
            features   => ['bam_input'],
            log_format => 'old',
            index_type => 'includes_reverse',
        },
        '0.6.1' => {
            path       => $ENV{GENOME_SW} . '/bwa/bwa-0.6.1/bwa',
            features   => ['bam_input', 'bwasw'],
            log_format => 'new',
            index_type => 'skips_reverse',
        },
        '0.6.2' => {
            path       => '/usr/bin/bwa0.6.2',
            features   => ['bam_input', 'bwasw'],
            log_format => 'new',
            index_type => 'skips_reverse',
        },
        '0.7.2' => {
            path       => $ENV{GENOME_SW} . '/bwa/bwa-0.7.2/bwa',
            features   => ['bam_input', 'bwasw', 'mem'],
            log_format => 'new',
            index_type => 'skips_reverse',
        },
        '0.7.3a' => {
            path       => '/usr/bin/bwa0.7.3a',
            features   => ['bam_input', 'bwasw', 'mem'],
            log_format => 'new',
            index_type => 'skips_reverse',
        },
        '0.7.5' => {
            path       => '/usr/bin/bwa0.7.5',
            features   => ['bam_input', 'bwasw', 'mem', 'supplementary_alignment_flag'],
            log_format => 'new',
            index_type => 'skips_reverse',
        },
        '0.7.5a' => {
            path       => '/usr/bin/bwa0.7.5a',
            features   => ['bam_input', 'bwasw', 'mem', 'supplementary_alignment_flag'],
            log_format => 'new',
            index_type => 'skips_reverse',
        },
        # If you are adding a new version of bwa, it likely has the following features:
        #   features   => ['bam_input', 'bwasw', 'mem'],
        #   log_format => 'new',
        #   index_type => 'skips_reverse',
        # New versions of iBWA should include 'multiple_reference_mode' under features.
    );

    if ($version) {
        my $result = $bwa_versions{$version};
        return unless $result; # undef if version does not exist
        return wantarray ? %$result : $result;
    } else {
        return wantarray ? %bwa_versions : \%bwa_versions;
    }
}

sub _verify_bwa_version {
    my $class = shift;
    my $version = shift;

    die 'No knowledge of bwa version '. $version unless $class->_lookup_version($version);

    my %version_hash = $class->_lookup_version($version);
    die "Invalid version hash for bwa version $version" unless $class->_version_hash_keys_valid(keys %version_hash);

    # verify features
    my @features = @{$version_hash{features}};
    for my $feature (@features) {
        die "Invalid feature $feature for bwa version $version" unless $class->_is_valid_feature($feature);
    }

    # verify log format
    my $log_format = $version_hash{log_format};
    die "Invalid log format $log_format for bwa version $version" unless $class->_is_valid_log_format($log_format);

    # verify index type
    my $index_type = $version_hash{index_type};
    die "Invalid index type $index_type for bwa version $version" unless $class->_is_valid_index_type($index_type);
}

sub _version_hash_keys_valid {
    my $class = shift;
    my @given_keys = sort {$a cmp $b} @_;
    my @expected_keys = sort {$a cmp $b} qw(path features log_format index_type);

    # verify hash keys
    return 0 if @expected_keys != @given_keys;
    for (0 .. $#given_keys) {
        return 0 if $expected_keys[$_] ne $given_keys[$_];
    }

    return 1;
}

sub _is_valid_feature {
    my ($class, $feature) = @_;
    return grep { $feature eq $_ } qw(bam_input multiple_reference_mode mem bwasw supplementary_alignment_flag);
}

sub _is_valid_log_format {
    my ($class, $log_format) = @_;
    return grep { $log_format eq $_ } qw(new old);
}

sub _is_valid_index_type {
    my ($class, $index_type) = @_;
    return grep { $index_type eq $_ } qw(includes_reverse skips_reverse);
}

sub available_bwa_versions {
    my %available_versions = _lookup_version();
    return keys %available_versions;
}

sub default_bwa_version {
    die "default bwa version: $DEFAULT is not valid" unless _lookup_version($DEFAULT);
    return $DEFAULT;
}

sub default_version { return default_bwa_version; }

sub bwa_path {
    my $self = $_[0];
    return $self->path_for_bwa_version($self->use_version);
}

sub path_for_bwa_version {
    my ($class, $version) = @_;

    $class->_verify_bwa_version($version); # this will throw an exception if it fails

    return $class->_lookup_version($version)->{path};
}

sub supports {
    my ($class, $version, $feature) = @_;
    $class->_verify_bwa_version($version); # this will throw an exception if it fails
    my @features = @{$class->_lookup_version($version)->{features}};
    return grep { $_ eq $feature } @features;
}

sub supports_bam_input {
    my ($class, $version) = @_;
    return $class->supports($version, 'bam_input');
}

sub supports_multiple_reference {
    my ($class, $version) = @_;
    return $class->supports($version, 'multiple_reference_mode');
}

sub supports_mem {
    my ($class, $version) = @_;
    return $class->supports($version, 'mem');
}

sub supports_bwasw {
    my ($class, $version) = @_;
    return $class->supports($version, 'bwasw');
}

sub supports_supplementary_alignment_flag {
    my ($class, $version) = @_;
    return $class->supports($version, 'supplementary_alignment_flag');
}

sub index_extensions {
    my ($class, $version) = @_;
    # Newer versions of bwa include do not include reversed index files.

    $class->_verify_bwa_version($version); # this will throw an exception if it fails

    my @output_extensions = qw(amb ann bwt pac sa);
    my @reverse_output_extensions = qw(rbwt rpac rsa);

    push @output_extensions, @reverse_output_extensions
        if $class->_lookup_version($version)->{index_type} eq 'includes_reverse';
    return @output_extensions;
}

sub log_format {
    my ($class, $version) = @_;
    # Newer versions of bwa have a slightly different log format

    $class->_verify_bwa_version($version); # this will throw an exception if it fails
    return $class->_lookup_version($version)->{log_format};
}

1;

