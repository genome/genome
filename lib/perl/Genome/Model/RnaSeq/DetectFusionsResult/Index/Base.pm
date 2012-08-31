package Genome::Model::RnaSeq::DetectFusionsResult::Index::Base;

use strict;
use warnings;

use Genome;
use LWP;
use LWP::UserAgent;
use File::Listing qw(parse_dir);


class Genome::Model::RnaSeq::DetectFusionsResult::Index::Base {
    is => 'Genome::SoftwareResult::Stageable',
    has_input => [
        reference_build => {
            is => "Genome::Model::Build::ReferenceSequence",
            doc => 'object representing the reference sequence version to use'
        },
    ],
    has_param => [
        ftp_server => {
            is_optional => 1,
            default => 'ftp://hgdownload.cse.ucsc.edu/goldenPath/[REF]/database/',
            doc => 'url for the ftp server to retrieve reference files, with [REF] specified in the url where the reference version goes',
        },
    ],
    doc => 'Base class for detect fusion result indexes',
};

sub resolve_allocation_disk_group_name {
    return 'info_genome_models';
}

sub _download_ucsc_table {
    my ($self, $remote_filename, $dir) = @_;

    Genome::Sys->download_file_to_directory($self->ftp_url . "/" . $remote_filename, $dir);

    my $local_file = $dir . "/" . $remote_filename;

    my $fh = Genome::Sys->open_gzip_file_for_reading($local_file);
    my $out_fh = Genome::Sys->open_file_for_writing( $dir . "/" . substr($remote_filename, 0, -3));

    while(my $line = <$fh>){
       $line =~ s/(^|\t)chr(.+)/$1$2/;
       print $out_fh $line;
    }

    $fh->close();
    $out_fh->close();
    unlink($local_file);

    return 1;
}

sub _ucsc_build_name {
    my $self = shift;
    my $version = $self->reference_build->version;

    #if we already have hgXX
    return $version if $version =~ /hg\d\d/;

    die($self->error_message("You have requested an unsupported reference sequence ($version)")) unless ($version - 18) >= 15;

    #if we have GRCh 36/37/etc
    return "hg" . ($version - 18);
}

sub ftp_url {
    my $self = shift;
    die ($self->error_message("This subclass does not implement this method")) unless $self->can("ftp_server");
    my $ref_sub_val = $self->_ucsc_build_name;
    my $ftp_server = $self->ftp_server;
    $ftp_server =~ s/\[REF\]/$ref_sub_val/;
    $ftp_server;
}

1;
