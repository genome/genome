package Genome::Model::RnaSeq::DetectFusionsResult::TophatFusionResult::Index;

use strict;
use warnings;

use Genome;
use LWP;
use LWP::UserAgent;
use File::Listing qw(parse_dir);

class Genome::Model::RnaSeq::DetectFusionsResult::TophatFusionResult::Index {
    is => 'Genome::Model::RnaSeq::DetectFusionsResult::Index::Base',
    has_param => [
        blast_ftp_url => {
            is_optional =>1,
            default => 'ftp://ftp.ncbi.nih.gov/blast/db/',
            doc => 'url for the ftp server to retrieve blast db files'
        },
    ],
    doc => 'create all the files needed for a successful run of TophatFusion'
};

#fusion_info.txt (optional) - omit for now

sub resolve_allocation_subdirectory {
    my $self = shift;
    return 'build_merged_alignments/tophatfusion-index/' . $self->id;
}

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_) or return;

    unless ($self->reference_build->species_name eq "human"){
        die($self->error_message("tophat-fusion only supports human fusion detection"));
    }

    $self->_prepare_staging_directory;

    #tophat fusion wants these specific filenames...
    for(qw(ensGene.txt.gz ensGtp.txt.gz refGene.txt.gz)){
        $self->_download_ucsc_table($_, $self->temp_staging_directory);
    }

    #yes..this path is hardcoded in tophat fusion
    my $dir = Genome::Sys->create_directory($self->temp_staging_directory . '/blast_human');
    for($self->_get_blast_files_to_download()){
        $self->_download_and_verify_blast_file($_, $dir);
    }

    $self->_prepare_output_directory;
    $self->_promote_data;
    $self->_reallocate_disk_allocation;

    return $self;
}

sub _get_blast_files_to_download {
    my $self = shift;

    my $ua = LWP::UserAgent->new;
    my $resp = $ua->get($self->blast_ftp_url);

    unless ($resp->is_success()){
        die($self->error_message($resp->status_line));
    }

    my @files = map { $_->[0] } parse_dir($resp->content);
    #human_genomic is hardcoded in tophat-fusion...
    @files = map { $self->blast_ftp_url . "/$_" } grep { $_ =~ /^human_genomic.+\.tar\.gz$/ || $_=~ /nt.+\.tar\.gz$/ } @files;
}

sub _download_and_verify_blast_file {
    my ($self, $remote_file, $dir) = @_;
    Genome::Sys->download_file_to_directory($remote_file, $dir);
    Genome::Sys->download_file_to_directory($remote_file . ".md5", $dir);

    my $filename = $dir . "/" . (split("/", $remote_file))[-1];

    my $file_digest = Genome::Sys->md5sum($filename);
    my $checksum = (split(/\s/,Genome::Sys->read_file($filename . ".md5")))[0];

    unless($file_digest eq $checksum){
        die($self->error_message("$filename appears to be corrupted, md5 checksum did not match!"));
    }

    Genome::Sys->shellcmd(
        cmd=>"tar -xzf $filename -C $dir",
        input_files => [$filename]
    ) or die ($self->error_message("Couldn't untar $filename!"));

    unlink($filename . ".md5") or die ($self->error_message("Couldn't remove $filename.md5"));
    unlink($filename) or die ($self->error_message("Couldn't remove $filename"));

    return 1;
}

sub _download_ucsc_table {
    my $self = shift;
    $self->SUPER::_download_ucsc_table(@_) or return;

    my ($remote_filename, $dir) = @_;
    if ($remote_filename =~ /refGene/){
        Genome::Sys->shellcmd(
            cmd => qq{cat $dir/refGene.txt | sort -V -k3,3 -k5,6  > $dir/refGene_sorted.txt},
            input_files => [$dir . "/refGene.txt"],
            output_files => [$dir . "/refGene_sorted.txt"]
        ) or die ($self->error_message("Couldn't sort refGene.txt!"));
        unlink($dir . "/refGene.txt");
    }

    return 1;

}

1;
