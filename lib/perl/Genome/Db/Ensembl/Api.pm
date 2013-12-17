package Genome::Db::Ensembl::Api;

use strict;
use warnings;
use Genome;
use Sys::Hostname;

class Genome::Db::Ensembl::Api {
    is => "Genome::SoftwareResult::Stageable",
    has_param => [
        version => {
            is => 'Text',
            doc => 'Version of ensembl API',
        },
    ],
};

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);

    $self->_prepare_staging_directory;

    my $version = $self->version;
    $self->status_message("Download ensembl API version $version");

    my @package_names = $self->package_names;
    my $base_url = "'http://cvs.sanger.ac.uk/cgi-bin/viewvc.cgi/PACKAGENAME.tar.gz?root=ensembl&only_with_tag=branch-ensembl-VERSION&view=tar'";

    my $temp_directory_path = $self->temp_staging_directory;

    for my $package_name (@package_names){
        my $tar_url = $base_url;
        $tar_url =~ s/PACKAGENAME/$package_name/;
        $tar_url =~ s/VERSION/$version/;
        my $tar_file = join("/", $temp_directory_path, "$package_name.tar.gz");
        my $extracted_directory = join("/", $temp_directory_path, $package_name);
        my $wget_command = "wget $tar_url -O $tar_file";
        my $rv = Genome::Sys->shellcmd(cmd => $wget_command, output_files =>  [$tar_file]);
        unless($rv){
            $self->error_message("Failed to download $package_name");
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

    }

    $self->status_message("Finished downloading ensembl API");
    
    $self->_prepare_output_directory;
    $self->_promote_data;
    $self->_reallocate_disk_allocation;

    return $self;

}

sub package_names {
    return qw/ ensembl ensembl-compara ensembl-variation ensembl-functgenomics /;
}

sub resolve_allocation_subdirectory {
    my $self = shift;
    my $hostname = hostname;

    my $user = $ENV{'USER'};
    my $base_dir = sprintf("ensemblapi-%s-%s-%s-%s",           $hostname,       $user, $$, $self->id);
    my $directory = join('/', 'build_merged_alignments',$self->id,$base_dir);
    return $directory;
}

sub resolve_allocation_disk_group_name {
    $ENV{GENOME_DISK_GROUP_MODELS};
}

sub prepend_api_path_and_execute {
    my $self = shift;
    my %shellcmd_params = @_;
    my @api_path = glob($self->output_dir."/ensembl*/modules");
    my $lib;
    if (@api_path){
        $lib = join(" ", $^X, '-S', map(join(" ", '-I', '"' . $_ . '"'), @api_path));
    }else{
        $self->error_message("No API path found for annotation build: " . $self->id);
        return;
    }

    $shellcmd_params{'cmd'} = join(" ", $lib, $shellcmd_params{'cmd'});
    my $rv = Genome::Sys->shellcmd(%shellcmd_params);

    return $rv;
}

1;

