package Genome::Model::Tools::Fastqc;

use strict;
use warnings;

use Genome;
use File::Basename;

my $FASTQC_DEFAULT = '0.10.0';
my $DEFAULT_MEMORY = 2;

class Genome::Model::Tools::Fastqc {
    is  => 'Command',
    has_input => [
        use_version => {
            is  => 'Version',
            doc => 'FastQC version to be used.',
            is_optional   => 1, 
            default_value => $FASTQC_DEFAULT,
        },
        maximum_memory => {
            is => 'Integer',
            doc => 'the maximum memory (Gb) to use when running Java VM. default_value='. $DEFAULT_MEMORY,
            is_optional => 1,
            default_value => $DEFAULT_MEMORY,
        },
    ],
};

sub help_brief {
    "Tools to run the Java toolkit FastQC and work with the output reports.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt fast-qc ...
EOS
}

sub help_detail {
    return <<EOS 
More information about the FastQC toolkit can be found http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/.
EOS
}

my %FASTQC_VERSIONS = (
    '0.10.0' => $ENV{GENOME_SW} . '/fastqc/FastQC-0.10.0',
    '0.6.1'  => $ENV{GENOME_SW} . '/fastqc/FastQC-0.6.1',
    '0.4.3'  => $ENV{GENOME_SW} . '/fastqc/FastQC-0.4.3',
    '0.4.1'  => $ENV{GENOME_SW} . '/fastqc/FastQC-0.4.1',
    '0.3'    => $ENV{GENOME_SW} . '/fastqc/FastQC-0.3',
    '0.2'    => $ENV{GENOME_SW} . '/fastqc/FastQC-0.2',
);

sub path_for_fastqc_version {
    my ($class, $version) = @_;
    $version ||= $FASTQC_DEFAULT;
    my $path = $FASTQC_VERSIONS{$version};
    return $path if defined $path;
    die 'No path found for fastqc version: '.$version;
}

sub default_fastqc_version {
    die "default fastqc version: $FASTQC_DEFAULT is not valid" unless $FASTQC_VERSIONS{$FASTQC_DEFAULT};
    return $FASTQC_DEFAULT;
}

sub fastqc_path {
    my $self = shift;
    return $self->path_for_fastqc_version($self->use_version);
}

sub run_java_vm {
    my $self = shift;
    my %params = @_;
    my $cmd = delete($params{'cmd'});
    unless ($cmd) {
        die('Must pass cmd to run_java_vm');
    }
    my $java_vm_cmd = 'java -Xmx'. $self->maximum_memory .'g '. $cmd;
    $params{'cmd'} = $java_vm_cmd;
    Genome::Sys->shellcmd(%params);
    return 1;
}



1;

