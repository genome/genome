package Genome::Model::Tools::Cmds::CreateMappedSnpArrayFile;

use strict;
use warnings;
use Genome;
use IO::File;
use Getopt::Long;

class Genome::Model::Tools::Cmds::CreateMappedSnpArrayFile {
    is => 'Command',
    has => [
        map_file => {
            type => 'String',
            is_optional => 0,
            doc => 'map.csv file giving coordinates of log2 snp array data'
        },
        snp_array_files => {
            type => 'String',
            is_optional => 0,
            doc => "A single-quoted string describing the input snp array data files, such as '/dir/*.log2' or '/dir/BRC*'.",
        },
        output_file => {
            type => 'String',
            is_optional => 0,
            doc => 'merged snp array data for a group of input files with column headers "CHR POS filename1 filename2..."'
        },
        map_headers => {
            type => 'Number',
            is_optional => 1,
            default => 1,
            doc => 'number of header lines to skip in the map file [1]'
        },
        snp_array_headers => {
            type => 'Number',
            is_optional => 1,
            default => 1,
            doc => 'number of header lines to skip in snp array files [1]'
        },
    ]
};

sub help_brief {
    'Create table of log2 snp array data [chr pos sample1 sample2...]'
}

sub help_detail {
    "This script creates a legible data file containing log2 snp array data with chromosome and position information. Can be used to create a mapped file for a single sample or for a group of samples which reference the same map.csv file."
}

sub execute {
    my $self = shift;

    #process input arguments
    my $outfile = $self->output_file;
    my $mapfile = $self->map_file;
    my $map_header_lines = $self->map_headers;
    my $snp_array_header_lines = $self->snp_array_headers;
    my @infiles = glob($self->snp_array_files);
    chomp @infiles;
    @infiles = sort @infiles; #so that files are always read & printed in same order
    
    #error checking:

    #check file lengths: first, find length of map file
    my $map_file_wc = `wc -l $mapfile`;
    $map_file_wc =~ s/^(\d+)\s+\w+$/$1/;

    #check for same length in snp-array-files
    for my $snp_file (@infiles) {
        my $snp_file_wc = `wc -l $snp_file`;
        $snp_file_wc =~ s/^(\d+)\s+\w+$/$1/;
        if ($snp_file_wc != $map_file_wc) {
            die "Found different lengths between map file $mapfile ($map_file_wc) and snp-array file $snp_file ($snp_file_wc).\n";
        }
    }

    #open filehandles
    my @filehandles;
    my $out_fh = new IO::File $outfile,"w";
    my $map_fh = new IO::File $mapfile,"r";
    for my $file (@infiles) {
        my $fh = new IO::File $file,"r";
        push @filehandles, $fh;
    }

    #print header line
    $out_fh->print("CHR\tPOS");
    for my $filename (@infiles) {
        $out_fh->print("\t$filename");
    }
    $out_fh->print("\n");

    #ignore header lines in input files
    for my $skip (1..$map_header_lines) {
        $map_fh->getline;
    }
    for my $file (@filehandles) {
        for my $skip (1..$snp_array_header_lines) {
            $file->getline;
        }
    }

    #write data to output file
    while (my $map_line = $map_fh->getline) {
        chomp $map_line;
        $out_fh->print($map_line);
        for my $file (@filehandles) {
            my $data = $file->getline;
            chomp $data;
            $out_fh->print("\t$data");
        }
        $out_fh->print("\n");
    }
    return 1;
}
1;
