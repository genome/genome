package Genome::Model::Tools::DetectVariants2::BamWindow;

use strict;
use warnings;

use Cwd;
use Genome;

my $DEFAULT_VERSION = '0.4';


class Genome::Model::Tools::DetectVariants2::BamWindow{
    is => ['Genome::Model::Tools::DetectVariants2::Detector'],
    doc => "Produces a list of high confidence somatic snps and indels.",
# TODO ... make sure this works without old default snv and indel params default => '-q 1 -Q 15',
    # Make workflow choose 64 bit blades
    has_param => [
        lsf_resource => {
            default_value => 'rusage[mem=4000] select[type==LINUX64 && maxtmp>10000] span[hosts=1]',
        },
    ],

};

my %BAMWINDOW_VERSIONS = (
    '0.4'   => "/gscuser/bniu/gc6143/cnv/code/bamwindow-v0.4/bam-window"
    );

#---------------------------------------
sub _detect_variants {
    my $self = shift;
    my $params = $self->params;   

    my $cmd = $self->bamwindow_path; 
    if(defined($params)){
        $cmd .= " " . $params . " " . $self->aligned_reads_input . " >" . $self->_temp_staging_directory."/readcounts.wind";
    }

    print STDERR "Running bam-window command: $cmd\n";

    #the system expects a cnvs.hq file - ugh.
    my $cnvs = $self->_temp_staging_directory."/cnvs.hq";
    system("touch $cnvs");


    my $return = Genome::Sys->shellcmd(
        cmd => "$cmd",
    );
    unless($return) {
        $self->error_message("Failed to execute: Returned $return");
        die $self->error_message;
    }
    return $return;
}

#---------------------------------------

sub _sort_detector_output {
    return 1;
}

sub bamwindow_path {
    my $self = $_[0];
    return $self->path_for_bamwindow_version($self->version);
}

sub available_bamwindow_versions {
    my $self = shift;
    return keys %BAMWINDOW_VERSIONS;
}

sub path_for_bamwindow_version {
    my $class = shift;
    my $version = shift;

    if (defined $BAMWINDOW_VERSIONS{$version}) {
        return $BAMWINDOW_VERSIONS{$version};
    }
    die('No path for bamwindow version '. $version);
}

sub has_version {
    my $self = shift;
    my $version = shift;
    unless(defined($version)){
        $version = $self->version;
    }
    if(exists($BAMWINDOW_VERSIONS{$version})){
        return 1;
    }
    return 0;
}


1;
