package Genome::Model::Tools::Sam;

use strict;
use warnings;

use Genome; 
use File::Basename;
use POSIX;
use DateTime;
use IO::File;

my $DEFAULT = 'r963';
#3Gb
my $DEFAULT_MEMORY = 402653184;

class Genome::Model::Tools::Sam {
    is  => 'Command',
    has_input => [
        use_version => { 
            is  => 'Version', 
            doc => "samtools version to be used, default is $DEFAULT. ", 
            is_optional   => 1, 
            default_value => $DEFAULT,   
        },
        maximum_memory => {
            is  => 'Integer',
            doc => "the maximum memory available, default is $DEFAULT_MEMORY",
            is_optional   => 1,
            default_value => $DEFAULT_MEMORY,
        },
    ],
};

sub sub_command_sort_position { 12 }

sub help_brief {
    "Tools to run Sam or work with its output files.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
genome-model tools Sam ...    
EOS
}

sub help_detail {                           
    return <<EOS 
More information about the Sam suite of tools can be found at http://Samtools.sourceforege.net.
Everytime when we get a new version of samtools, we need update in this module and create new 
processing_profile/model for pipeline.
EOS
}


my %SAMTOOLS_VERSIONS = (
    '0.1.19' => $ENV{GENOME_SW} . '/samtools/samtools-0.1.19',
    r982    => $ENV{GENOME_SW} . '/samtools/samtools-0.1.18',
    r973    => $ENV{GENOME_SW} . '/samtools/samtools-0.1.17',
    r963    => $ENV{GENOME_SW} . '/samtools/samtools-0.1.16',
    r868    => $ENV{GENOME_SW} . '/samtools/samtools-0.1.12a',
    r783    => $ENV{GENOME_SW} . '/samtools/samtools-0.1.9',
    r613    => $ENV{GENOME_SW} . '/samtools/samtools-0.1.8',
    r599    => $ENV{GENOME_SW} . '/samtools/samtools-0.1.7ar599',
    r544    => $ENV{GENOME_SW} . '/samtools/samtools-0.1.7ar544',
    r510    => $ENV{GENOME_SW} . '/samtools/samtools-0.1.7a',
    r453    => $ENV{GENOME_SW} . '/samtools/samtools-0.1.6',
    r449    => $ENV{GENOME_SW} . '/samtools/samtools-0.1.5-32',
    r301wu1 => '/gscuser/dlarson/samtools/r301wu1',
    r320wu1 => '/gscuser/dlarson/samtools/r320wu1',
    r320wu2 => '/gscuser/dlarson/samtools/r320wu2',
    r350wu1 => '/gscuser/dlarson/samtools/r350wu1',
);


sub available_samtools_versions {
    my $self = shift;
    return keys(%SAMTOOLS_VERSIONS);
}

sub verify_sam_path {
    my ($proto, $version) = @_;

    unless ($version) {
        die 'No samtools version provided';
    }

    my $path = $SAMTOOLS_VERSIONS{$version};

    unless ($path and -d $path) {
        die 'No path found for samtools version: '. $version;
    }
    return $path;
}

sub path_for_samtools_version {
    my ($proto, $version) = @_;
    $version ||= $DEFAULT;
    my $path      = $proto->verify_sam_path($version);
    my $tool_path = $path . '/samtools';

    unless (-x $tool_path) {
        die 'No samtools path found for samtools version: '. $version;
    }
    return $tool_path;
}

#each samtools release version contains a bcftools
sub path_for_bcftools {
    my ($proto, $version) = @_;
    my $path      = $proto->verify_sam_path($version);
    my $tool_path = "$path/bcftools/bcftools";

    unless (-x $tool_path) {
        die $proto->error_message("bcftools: $tool_path is not executable");
    }
    return $tool_path;
}

sub path_for_vcfutils {
    my ($proto, $version) = @_;
    my $path      = $proto->verify_sam_path($version);
    my $tool_path = "$path/bcftools/vcfutils.pl";

    unless (-x $tool_path) {
        die $proto->error_message("vcfutils.pl: $tool_path is not executable");
    }
    return $tool_path;
}

sub default_samtools_version {
    die "default samtools version: $DEFAULT is not valid" unless $SAMTOOLS_VERSIONS{$DEFAULT};
    return $DEFAULT;
}

sub default_samtools_maximum_memory {
    return $DEFAULT_MEMORY;
}

sub samtools_path {
    my $self = shift;
    return $self->path_for_samtools_version($self->use_version);
}

sub samtools_pl_path {
    my $self = shift;
    my $dir  = dirname $self->samtools_path;
    my $path = "$dir/misc/samtools.pl";

    unless (-x $path) {
        $self->error_message("samtools.pl: $path is not executable");
        return;
    }
    return $path;
}

sub c_linkage_class {
    my $self = shift;

    $DB::single = $DB::stopper;
    my $version = $self->use_version;
    $version =~ s/\./_/g;

    my $class_to_use = __PACKAGE__ . "::CLinkage$version";

    #eval "use above '$class_to_use';";
    eval "use $class_to_use;";
    if ($@) {
        $self->error_message("Failed to use $class_to_use: $@");
        return undef;
    }

    return $class_to_use;
}

sub open_bamsam_in {
    my $self = shift;
    my $in_filename = shift;
    my ($type) = ($in_filename =~ /\.([^\.\s]+)\s*$/i);
    $type = uc($type);
    my $fh;
    if($type eq 'BAM') {
        $fh = new IO::File;
        #$fh->open('samtools view -h "' . $in_filename . '" | head -n 2000000 |');
        $fh->open('samtools view -h "' . $in_filename . '" |');
    }
    elsif($type eq 'SAM') {
        $fh = IO::File->new($in_filename);
    }
    else {
        die 'Unknown type specified for "' . $in_filename . "\".\n";
    }
    unless($fh) {
        die 'Failed to open "' . $in_filename . "\"\n.";
    }
    return $fh;
}

sub open_bamsam_out {
    my $self = shift;
    my $out_filename = shift;
    my ($type) = ($out_filename =~ /\.([^\.\s]+)\s*$/i);
    $type = uc($type);
    my $fh;
    if($type eq 'BAM') {
        $fh = new IO::File;
        $fh->open('| samtools view -S -b /dev/stdin > "' . $out_filename . '"');
    }
    elsif($type eq 'SAM') {
        $fh = IO::File->new($out_filename eq '-' ? stdout : '> ' . $out_filename);
    }
    else {
        die 'Unknown type specified for "' . $out_filename . "\".\n";
    }
    unless($fh) {
        die 'Failed to open "' . $out_filename . "\"\n.";
    }
    return $fh;
}

sub read_count {
    my $self = shift;
    my $filename = shift;

    # Check to see if there is a flagstat file beside the bam before line-counting the entire thing....
    if(-s $filename.".flagstat" ){
        my $flag_object = Genome::Model::Tools::Sam::Flagstat->create(bam_file => $filename, output_file => $filename.".flagstat");
        my $data = $flag_object->parse_file_into_hashref($filename.".flagstat");
        return $data->{'total_reads'};
    }

    my $samtools = $self->samtools_path;

    my ($type) = ($filename =~ /\.([^\.\s]+)\s*$/i);
    my $count_cmd;
    if ($type =~ /SAM/i) {
        $count_cmd = "grep -cv '^\@' $filename";
    } 
    elsif ($type =~ /BAM/i) {
        $count_cmd = "$samtools view $filename | wc -l";
    } 
    else {
        $self->error_message("Unknown type ($type) from filename ($filename).");
        return;
    }

    chomp(my $read_count = qx($count_cmd));
    ($read_count) = split(' ', $read_count);
    return $read_count;
}

sub read_length {
    my $self = shift;
    my $filename = shift;
    my $cmd = $self->samtools_path." view ".$filename."|";
    my $bam_fh = IO::File->new( $cmd);
    my $read = $bam_fh->getline;
    $bam_fh->close;
    my @read_fields = split /\t/, $read;
    my $read_length = length($read_fields[9]);
    return $read_length;
}

sub date {
    my $today = DateTime->today();
    $today =~ s/T.*//;
    return $today;
}

sub time {
    # Timestamps in SAM are to follow ISO8601, e.g. YYYY-MM-DD and YYYY-MM-DDThh:mmTZD
    # TZD = Z or +hh:mm or -hh:mm
    return DateTime->now() . "Z";
}

1;
