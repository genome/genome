package Genome::Db::Ensembl::VepCache;

use strict;
use warnings;
use Genome;
use Sys::Hostname;

class Genome::Db::Ensembl::VepCache {
    is => "Genome::SoftwareResult::Stageable",
    has_param => [
        version => {
            is => 'Text',
            doc => 'Version of ensembl db',
        },
        species => {
            is => 'Text',
            doc => 'Species contained in cache',
        },
        sift => {
            is => 'Boolean',
            doc => 'Include sift, condel, and polyphen in cache',
        },
    ],
};

my %SPECIES_LOOKUP = (
    human => "homo_sapiens",
    mouse => "mus_musculus",
);

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);

    $self->_prepare_staging_directory;

    my $alternate = Genome::Db::Ensembl::VepCache->get_with_lock(version => $self->version,
                                                     species => $self->species,
                                                     sift => !($self->sift),
                                                     test_name => $self->test_name);
    if (defined $alternate) {
        my $cp_cmd = join(" ", "cp -r", $alternate->output_dir."/*", $self->temp_staging_directory);
        Genome::Sys->shellcmd(cmd => $cp_cmd);
    }

    else {
        my $version = $self->version;
        my $species = $self->species;
        if ($SPECIES_LOOKUP{$species}) {
            $species = $SPECIES_LOOKUP{$species};
        }
        print ("Download ensembl vep cache version $version of species $species");

        #ensembl-66 and later
        my $base_url = "ftp://ftp.ensembl.org/pub/release-VERSION/variation/VEP/SPECIES_vep_VERSION";
        #ensembl-64 and 65
        #my $base_url = "ftp://ftp.ensembl.org/pub/release-VERSION/variation/VEP/SPECIES/SPECIES_vep_VERSION";

        $base_url = $base_url.".tar.gz";

        my $temp_directory_path = $self->temp_staging_directory;

        my $tar_url = $base_url;
        $tar_url =~ s/VERSION/$version/g;
        $tar_url =~ s/SPECIES/$species/g;
        my $tar_file = join("/", $temp_directory_path, "$species"."_vep_$version.tar.gz");
        my $extracted_directory = join("/", $temp_directory_path, $species);
        my $wget_command = "wget $tar_url -O $tar_file";
        my $rv = Genome::Sys->shellcmd(cmd => $wget_command, output_files =>  [$tar_file]);
        unless($rv){
            $self->error_message("Failed to download vep cache version $version for $species");
            $self->delete;
            return;
        }

        my $extract_command = "tar -xzf $tar_file -C $temp_directory_path";
        $rv = Genome::Sys->shellcmd(cmd => $extract_command, input_files =>   [$tar_file], output_directories => [$extracted_directory]);
        unless($rv){
            $self->error_message("Failed to extract $tar_file");
            $self->delete;
            return;
        }

        $self->status_message("Finished downloading vep cache");
    }

    $self->_prepare_output_directory;
    $self->_promote_data;
    $self->_reallocate_disk_allocation;

    return $self;

}

sub resolve_allocation_subdirectory {
    my $self = shift;
    my $hostname = hostname;

    my $user = $ENV{'USER'};
    my $base_dir = sprintf("ensemblvepcache-%s-%s-%s-%s",           $hostname,       $user, $$, $self->id);
    my $directory = join('/', 'build_merged_alignments',$self->id,$base_dir);
    return $directory;
}

sub resolve_allocation_disk_group_name {
    $ENV{GENOME_DISK_GROUP_MODELS};
}

1;

