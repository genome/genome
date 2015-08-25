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
    my $base_url = "https://github.com/Ensembl/PACKAGENAME.git";

    my $temp_directory_path = $self->temp_staging_directory;

    for my $package_name (@package_names){
        my $git_url = $base_url;
        $git_url =~ s/PACKAGENAME/$package_name/;
        $git_url =~ s/VERSION/$version/;
        $self->download_and_extract($git_url,
                        $temp_directory_path,
                        "release/".$version,
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
    my $git_url = shift;
    my $temp_dir = shift;
    my $branch_name = shift;
    my $final_location = shift;

    my $download_dir = File::Spec->join($temp_dir, $final_location);
    my $clone_command = "git clone $git_url $download_dir";
    my $rv = Genome::Sys->shellcmd(cmd => $clone_command);
    unless($rv){
        $self->error_message("Failed to clone $git_url");
        $self->delete;
        return $rv;
    }

    my $current_location = `pwd`;
    my $git_command = "cd $download_dir; git checkout $branch_name; cd $current_location";
    $rv = Genome::Sys->shellcmd(cmd => $git_command);
    unless($rv){
        $self->error_message("Failed to checkout branch $branch_name");
        $self->delete;
        return $rv;
    }

    return 1;
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
    Genome::Config::get('disk_group_models');
}

sub prepend_api_path_and_execute {
    my $self = shift;
    my %shellcmd_params = @_;
    my @api_path = glob($self->output_dir."/ensembl*/modules");
    my $lib;
    if (@api_path){
        $lib = join(" ", $^X, '-S', map(join(" ", '-I', '"' . $_ . '"'), @api_path));
    }else{
        die $self->error_message("No API path found for annotation build: " . $self->id);
    }

    $shellcmd_params{'cmd'} = join(" ", $lib, $shellcmd_params{'cmd'});
    return Genome::Sys->shellcmd(%shellcmd_params);
}

sub vep_script {
    my $self = shift;
    my $script_name = shift;
    return join("/", $self->output_dir, "variant_effect_predictor", $script_name);
}
1;

