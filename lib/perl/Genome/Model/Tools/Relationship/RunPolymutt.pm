package Genome::Model::Tools::Relationship::RunPolymutt;

use strict;
use warnings;
use Data::Dumper;
use Genome;           
use Genome::Info::IUB;
use POSIX;
our $DEFAULT_VERSION = '0.11';
use Cwd;
use File::Basename;
use File::Path;

class Genome::Model::Tools::Relationship::RunPolymutt {
    is => 'Command',
    has_input => [
        version => {
            is => 'Text',
            default => $DEFAULT_VERSION,
            doc => "Version to use",
        },
        denovo => {
            is=>'Text',
            is_optional=>1,
            default=>0,
        },
        output_vcf => {
            is=>'Text',
            is_output=>1,
        },
        glf_index => {
            is=>'Text',
        },
        dat_file => {
            is=>'Text',
        },
        ped_file => {
            is=>'Text',
        },
        threads => {
            is=>'Text',
            is_optional=>1,
            default=>4,
        },
        bgzip => {
            is_optional=>1,
            default=>1,
            doc=>'set this to 0 if you prefer uncompressed',
        },
        chr2process=> {
            is_optional=>1,
            default=>undef,
        },
       all_sites => {
            is => "Boolean",
            is_optional => 1,
            default => 0,
            doc => "Output calls on all sites (force genotype) in the glf files.",
       },
       roi_file => {
            is => "Path",
            is_optional => 1,
            doc => "Output calls on all sites (force genotype) in this roi file.",
       },
    ],
    has_param => [
        lsf_resource => {
            is => 'Text',
            default => "-R 'span[hosts=1] rusage[mem=1000] -n 4'",
        },
        lsf_queue => {
            is => 'Text',
            default => 'long',
        },
    ],
};

sub help_brief {
    "simulates reads and outputs a name sorted bam suitable for import into Genome::Model"
}

sub help_detail {
}

my %VERSIONS = (
    '0.02' => '/usr/bin/polymutt0.02',
    '0.10' => '/gscmnt/gc6126/info/medseq/launch_cleft_lip_phenotype_correlation/polymutt_binary/polymutt0.10',
    '0.11' => '/gscmnt/gc6126/info/medseq/launch_cleft_lip_phenotype_correlation/polymutt_binary/polymutt0.11',
    # New, experimental version that optionally runs on vcf and can force-genotype sites
    'vcf.0.01' => '/gscmnt/gc6126/info/medseq/launch_cleft_lip_phenotype_correlation/polymutt_binary/polymuttvcf.0.01',
);

sub path_for_version {
    my $class = shift;
    my $version = shift || $DEFAULT_VERSION;

    unless(exists $VERSIONS{$version}) {
        $class->error_message('No path found for polymutt version ' . $version);
        die $class->error_message;
    }

    return $VERSIONS{$version};
}

sub default_version {
    my $class = shift;

    unless(exists $VERSIONS{$DEFAULT_VERSION}) {
        $class->error_message('Default polymutt version (' . $DEFAULT_VERSION . ') is invalid.');
        die $class->error_message;
    }

    return $DEFAULT_VERSION;
}

sub available_versions {
    return keys(%VERSIONS);
}

#/gscuser/dlarson/src/polymutt.0.01/bin/polymutt -p 20000492.ped -d 20000492.dat -g 20000492.glfindex --minMapQuality 1 --nthreads 4 --vcf 20000492.standard.vcf
sub execute {
    $DB::single=1;
    my $self=shift;
    my $polymutt_cmd= $self->path_for_version($self->version);
    my $ped_file = $self->ped_file;
    my $dat_file = $self->dat_file;
    my $glf_index= $self->glf_index;
    my $threads = $self->threads;
    my $output_vcf = $self->output_vcf;
    my ($temp_output) = Genome::Sys->create_temp_file_path();
    my $cmd = $polymutt_cmd;
     $cmd .= " -p $ped_file";
     $cmd .= " -d $dat_file";
     $cmd .= " -g $glf_index";
     $cmd .= " --minMapQuality 1";
     $cmd .= " --nthreads $threads";

     # The output param name changed in the latest version
     if ($self->version eq "vcf.0.01") {
         $cmd .= " --out_vcf $temp_output";
     } else {
         $cmd .= " --vcf $temp_output";
     }

     if ($self->chr2process) {
         my $chrs = $self->chr2process;
         $cmd .= " --chr2process $chrs";
     }
    if($self->denovo) {
        $cmd .= " --denovo";
    }
    if($self->all_sites) {
        $cmd .= " --all_sites";
    }
    if($self->roi_file) {
        $cmd .= " --pos " . $self->roi_file;
    }

    print "Running: $cmd\n"; 
    my $rv = Genome::Sys->shellcmd(cmd=> $cmd);
    if($rv != 1) {
        return;
    }
    if($self->bgzip) {
        my $cmd = "bgzip -c $temp_output > $output_vcf";
        Genome::Sys->shellcmd(cmd=>$cmd); 
    }

    return 1;
}

1;
