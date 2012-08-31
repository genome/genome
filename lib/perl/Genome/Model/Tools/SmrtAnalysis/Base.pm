package Genome::Model::Tools::SmrtAnalysis::Base;

use strict;
use warnings;

use Genome;

my $DEFAULT_VERSION = '1.2.1';

my $DEFAULT_LSF_QUEUE = 'techd';
my $DEFAULT_LSF_RESOURCE = "-g /pacbio/smrtanalysis -M 4000000 -R 'select[type==LINUX64 && mem>=4000 && tmp>=40000] rusage[mem=4000,tmp=20000]'";

class Genome::Model::Tools::SmrtAnalysis::Base {
    is  => 'Command::V2',
    is_abstract => 1,
    has => [
        use_version => {
            default_value => Genome::Model::Tools::SmrtAnalysis::Base->default_version,
        },
        #seymour_home => {
        #    is => 'Text',
        #    doc => 'The base directory for the SMRT Analysis install',
        #    is_optional => 1,
        #    default_value => Genome::Model::Tools::SmrtAnalysis::Base->default_seymour_home,
        #},
        #sh_setup => {
        #    is => 'Text',
        #    doc => 'The sh script that is sourced before running each command',
        #    is_optional => 1,
        #    default_value => $SH_SETUP,
        #},
        #analysis_bin => {
        #    is_calculated => 1,
        #    calculate_from => ['seymour_home',],
        #    calculate => q{ $seymour_home .'/analysis/bin'; },
        #},
    ],
    has_optional_param => [
        lsf_queue => { default_value => $DEFAULT_LSF_QUEUE },
        lsf_resource => { default_value => $DEFAULT_LSF_RESOURCE },
    ],
};

sub default_version {
    return $DEFAULT_VERSION;
}

sub base_dir {
    return '/gscmnt/pacbio/production/git/pacbio';
}

sub r_lib {
    my $self = shift;
    my $r_lib = $self->base_dir .'/lib/R';
    unless (-d $r_lib) {
        die('Failed to find R lib: '. $r_lib);
    }
    return $r_lib;
}

sub seymour_home {
    my $self = shift;
    my $seymour_home = $self->base_dir .'/smrtanalysis-'. $self->use_version;
    unless (-d $seymour_home) {
        die('Failed to find seymour home: '. $seymour_home);
    }
    return $seymour_home;
}

sub sh_setup {
    my $self = shift;
    my $seymour_home = $self->seymour_home;
    my $sh_setup = $seymour_home .'/etc/setup.sh';
    unless (-f $sh_setup) {
        die('Failed to find sh setup: '. $sh_setup);
    }
    return $sh_setup;
}

sub analysis_bin {
    my $self = shift;
    my $seymour_home = $self->seymour_home;
    my $analysis_bin = $seymour_home .'/analysis/bin';
    unless (-d $analysis_bin) {
        die('Failed to find analysis/bin: '. $analysis_bin);
    }
    return $analysis_bin;
}

sub shellcmd {
    my $self = shift;
    my %params = @_;
    #unless ($params{cmd}) {
    #    die('Failed to provide command in cmd param!');
    #}
    #my $cmd = delete($params{cmd});
    #TODO: Figure out how this will work from tcsh, csh, etc.
    #my $new_cmd = '. '. $self->sh_setup .' && '. $cmd;
    #$params{cmd} = $new_cmd;
    Genome::Sys->shellcmd(%params);
    return 1;
}

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);

    # TODO: Remove this once tests have been created and run successfully as any UNIX user
    # There are problems with sourcing the right setup.sh for environment variables necessary for running smrtanalyis software using the correct login shell(bash only)
    my $user = Genome::Sys->username;
    unless ($user eq 'smrtanalysis') {
        die('Currently running the SMRT Analysis package is limited to the smrtanalysis user!');
    }

    #TODO: Check for proper environment variables
    #my $setup = $self->sh_setup;
    $ENV{SEYMOUR_HOME} = $self->seymour_home;
    $ENV{SEYMOUR_JAVA_CP} = '';

    return $self;
}
