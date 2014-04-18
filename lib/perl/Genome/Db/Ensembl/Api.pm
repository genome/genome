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
    #OLD CVS location
    #my $base_url = "'http://cvs.sanger.ac.uk/cgi-bin/viewvc.cgi/PACKAGENAME.tar.gz?root=ensembl&only_with_tag=branch-ensembl-VERSION&view=tar'";
    my $base_url = "https://github.com/Ensembl/PACKAGENAME/archive/release/VERSION.zip";

    my $temp_directory_path = $self->temp_staging_directory;

    #$self->download_vep($version, $temp_directory_path);

    for my $package_name (@package_names){
        my $tar_url = $base_url;
        $tar_url =~ s/PACKAGENAME/$package_name/;
        $tar_url =~ s/VERSION/$version/;
        my $tar_file = join("/", $temp_directory_path, "$package_name.zip");
        $self->download_and_extract($tar_url, $tar_file,
                        $temp_directory_path, $package_name."-release-".$version,
                        $package_name);
    }
    $self->move_vep($temp_directory_path."/ensembl-tools/scripts/variant_effect_predictor/",
        $temp_directory_path);

    $self->status_message("Finished downloading ensembl API");

    $self->_prepare_output_directory;
    $self->_promote_data;
    $self->_reallocate_disk_allocation;

    return $self;

}

sub move_vep {
    my $self = shift;
    my $from = shift;
    my $to = shift;
    my $mv_cmd = "mv $from $to";
    Genome::Sys->shellcmd(cmd => $mv_cmd);
}

sub download_and_extract {
    my $self = shift;
    my $tar_url = shift;
    my $tar_file = shift;
    my $extract_path = shift;
    my $extracted_directory_name = shift;
    my $final_location = shift;

    my $extracted_directory = join("/", $extract_path, $extracted_directory_name);

    my $wget_command = "wget $tar_url -O $tar_file";
    my $rv = Genome::Sys->shellcmd(cmd => $wget_command, output_files =>  [$tar_file]);
    unless($rv){
        $self->error_message("Failed to download $tar_file");
        $self->delete;
        return $rv;
    }

    my $extract_command = "unzip $tar_file -d $extract_path";
    $rv = Genome::Sys->shellcmd(cmd => $extract_command, input_files =>   [$tar_file], output_directories => [$extracted_directory]);
    unless($rv){
        $self->error_message("Failed to download and extract $tar_file");
        $self->delete;
        return;
    }

    my $final_directory = File::Spec->join($extract_path, $final_location);
    my $mv_cmd = "mv $extracted_directory $final_directory";
    Genome::Sys->shellcmd(cmd => $mv_cmd, input_files => [$extracted_directory], output_directories => [$final_directory]);
    return 1;
}

sub download_vep {
    my $self = shift;
    my $version = shift;
    my $download_path = shift;

    my $vep_name = "variant_effect_predictor";
    my $tar_url = "'http://cvs.sanger.ac.uk/cgi-bin/viewvc.cgi/ensembl-tools/scripts/$vep_name.tar.gz?root=ensembl&pathrev=branch-ensembl-VERSION&view=tar'";
    $tar_url =~ s/VERSION/$version/;

    return $self->download_and_extract($tar_url, $download_path."/$vep_name.tar.gz", $download_path, $vep_name);
}

sub package_names {
    return qw/ ensembl ensembl-compara ensembl-variation ensembl-funcgen ensembl-tools/;
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

sub vep_script {
    my $self = shift;
    my $script_name = shift;
    return join("/", $self->output_dir, "variant_effect_predictor", $script_name);
}
1;

