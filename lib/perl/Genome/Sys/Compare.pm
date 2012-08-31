package Genome::Sys::Compare;

use strict;
use warnings;  

use Genome;
require Digest::MD5;
require IO::File;    
require Text::Table;
require File::Compare;


class Genome::Sys::Compare {
    is => 'UR::Singleton',
};


sub compare_dirs {
    my $class = shift;
    my @base_dirs = @_;

    my @file_basenames = $class->build_file_list_from_dirs(@base_dirs);
    #TODO put in an ignore intercept
    
    my @status;
    for my $file_basename (@file_basenames) {
        my @files = map { $_ . '/' . $file_basename } @base_dirs;
        my ($compare_result, $file_statuses) = $class->compare_files(@files);
        push @status, [$file_basename, $compare_result, $file_statuses];
    }

    my $compare_dirs_result = 0;
    $compare_dirs_result = 1 if (grep { $_->[1] == 1 } @status);
    $compare_dirs_result = -1 if (grep { $_->[1] == -1 } @status);

    my @headers = ('File', 'Different');
    if (@base_dirs > 2) {
        push @headers, @base_dirs;
    }
    @headers = map { "$_\n" . '-' x length($_) } @headers;

    return ($compare_dirs_result, \@headers, \@status);
}


sub compare_files {
    my $class = shift;
    my @files = @_;

    if (@files == 2) {
        return File::Compare::compare(@files) || '0';
    }
    elsif (@files > 2) {
        return $class->compare(@files);
    }
    else {
        die "Did not receive >= 2 files!\n";
    }
}


sub compare {
    my $class = shift;
    my @files = @_;

    my %hashes = ('0' => '0');
    my @file_status;
    for my $file (@files) {
        my $md5 = $class->md5($file);
        $hashes{$md5} = scalar keys %hashes unless defined $hashes{$md5};
        push @file_status, $hashes{$md5};
    }

    my $compare_result = 0;
    if (grep { $_ ne $file_status[0] } @file_status) {
        $compare_result = 1;
    }

    return ($compare_result, \@file_status);
}


sub print_table {
    my ($class, $headers, $data) = @_;
    my $tb = Text::Table->new(@$headers);
    $tb->load(@$data);
    print $tb;    
}


sub build_file_list_from_dirs {
    my $class = shift;
    my @dirs = @_;
    
    my %files;
    for my $dir (@dirs) {
        my @file_paths = qx(find $dir -type f);
        map { $_ =~ s/^$dir\/+// } @file_paths;
        chomp(@file_paths);
        @files{@file_paths} = @file_paths;
    }
    
    my @files = sort { $a cmp $b } keys %files;
    
    return @files;
}


sub md5 {
    my $class = shift;
    my $path = shift;
    return 0 unless (-s $path);
    my $file_handle = IO::File->new($path);
    return Digest::MD5->new->addfile($file_handle)->hexdigest;
}


1;
