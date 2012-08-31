package Genome::Model::Tools::Bmr::WtfToBitmask;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;
use Bit::Vector;

class Genome::Model::Tools::Bmr::WtfToBitmask {
    is => 'Command',
    has => [
    wtf_file =>
    {
        type => 'String',
        is_optional => 0,
        doc => 'File to convert wiggle track format to whole genome bitmask. Expected to be 1-based coordinates and tab-separated',
    },
    output_file =>
    {
        type => 'String',
        is_optional => 1,
        doc => 'Name of the bitmask file to write to disk. If this option is not passed, then no output will be written',
    },
    exclude_regions => 
    {
        type => 'Boolean',
        is_optional => 1,
        default => 0,
        doc => 'Instead of flipping bits on from the wiggle file. Assume these are locations to flip off bits in the bitmask',

    },
    reference_index => 
    {
        type => 'String',
        is_optional => 1,
        default => '/gscmnt/839/info/medseq/reference_sequences/NCBI-human-build36/all_sequences.fa.fai',
        doc => 'samtools index of the reference genome for which you are creating a mask',

    },
    _bitmask => {
        type => 'HashRef',
        is_optional => 1,
    },
    ]
};


=cut
We should expect to see wiggle files like so indicating coverages
I suppose the assumption would be that unlisted bases were never targeted? In that case the coverage would be 0 or a third value? Hard to say.
fixedStep chrom=chr1 start=100090814 step=1
0
0
fixedStep chrom=chr1 start=100099616 step=1
1
1

=cut


sub execute {
    my $self=shift;
    #$DB::single = 1;
    #check inputs
    my $reference_index = $self->reference_index;
    if(-z $reference_index) {
        $self->error_message("$reference_index is of zero size or does not exist");
        return;
    }

    my $wtf_file = $self->wtf_file;
    if(-z $wtf_file) {
        $self->error_message("$wtf_file is of zero size or does not exist");
        return;
    }

    #read in reference sequences
    my $ref_fh = IO::File->new($reference_index);
    unless($ref_fh) {
        $self->error_message("Couldn't open $reference_index for parsing in chromosome lengths");
        return;
    }

    #ready to create the bitvectors
    my %genome;
    while(my $line = $ref_fh->getline) {
        chomp $line;
        my ($chr, $length) = split /\t/, $line;
        $genome{$chr} = Bit::Vector->new($length+1); #this throws an exception if it fails. Probably should be trapped at some point in the future. Adding 1 to allow for 1 based coordinates
        if($self->exclude_regions) {
            $genome{$chr}->Fill();  #flip to all on
        }

    }
    $ref_fh->close;

    my $wtf_fh = IO::File->new($wtf_file);
    unless($wtf_fh) {
        $self->error_message("Couldn't open $wtf_file wiggle file of regions");
        return;
    }
    my ($chromosome, $starting_origin, $block_start,$block_end, $span, $step, $type) = (0,0,0,0,1,0);
    while(my $line = $wtf_fh->getline) {
        chomp $line;
        next if $line =~ /^#/;  #support comments
        next if $line =~ /^$/;  #skip blank lines
        next if $line =~ /^track/; #Broad's wiggle header
        if($line =~ /fixedStep/) {
            ($chromosome,$starting_origin, $step) = $line =~ /^fixedStep\s+chrom=(\S+)\s+start=(\d+)\s+step=(\d+)/;
            unless($chromosome && $starting_origin && $step) {
                $self->error_message("Invalid fixedStep declaration.");
                return;
            }
            if ($chromosome =~ /^chr/) {
                $chromosome =~ s/^chr(.+)$/$1/;
            }
            ($span) = $line =~ /\s+span=(\d+)$/;
            $span = 1 unless $span;
            unless(exists($genome{$chromosome})) {
                $self->error_message("Chromosome $chromosome in line " . $wtf_fh->input_line_number . " does not exist in reference file.");
                return;
            }
            $type = 'fixed';
            $block_start = $starting_origin - $step;
        }
        elsif($line =~ /variableStep/) {
            $self->error_message("variableStep format unimplemented");
            return;
        }
        else {
            if($type eq 'fixed') {
                $block_start = $block_start + $step;
                $block_end = $block_start + $span - 1;
            }

            unless($block_start > 0 && $block_end > 0 && $block_start < $genome{$chromosome}->Size() && $block_end < $genome{$chromosome}->Size()) {
                $self->error_message("Coordinates on line " . $wtf_fh->input_line_number . " outside of the bounds of chromosome $chromosome. Skipping...");
                next;
            }
            #Assume non-zero equals on. Zero equals off
            my ($value) = $line =~ /^\s*(\d+)\s*$/;
            unless(defined $value) {
                $self->error_message("Non-numeric value on line " . $wtf_fh->input_line_number );
                return;
            }
            if($value != 0) {
                $genome{$chromosome}->Interval_Flip($block_start, $block_end);
            }
        }
    }
    $wtf_fh->close;

    $self->_bitmask(\%genome);

    if($self->output_file) {
        $self->write_genome_bitmask($self->output_file,$self->_bitmask);
    }
    return 1;
}

sub bitmask {
    my $self= shift;
    return $self->_bitmask;
}

sub write_genome_bitmask {
    my ($self,$filename,$genome_ref) = @_;
    unless($filename) {
        $self->error_message("No filename of file to write to");
        return;
    }
    unless($genome_ref) {
        $self->error_message("No bitmask to write to file");
        return;
    }
    #do some stuff to write this to a file without making it suck
    my $out_fh = IO::File->new($filename,">:raw");
    unless($out_fh) {
        $self->error_message("Unable to write to " . $filename);
        return;
    }
    my $header_string = join("\t", map {$_ => $genome_ref->{$_}->Size()} sort keys %$genome_ref);
    my $write_string = pack 'N/a*', $header_string;
    my $write_result = syswrite($out_fh,$write_string);
    unless(defined $write_result && $write_result == length($write_string)) {
        $self->error_message("Error writing the header");
        return;
    }
    for my $chr (sort keys %$genome_ref) {
        #first write the length in bytes 
        my $chr_write_string = $genome_ref->{$chr}->Block_Read();
        $write_result = syswrite $out_fh, pack("N",length($chr_write_string));
        unless(defined $write_result || $write_result != 4) {
            $self->error_message("Error writing the length of chromosome $chr");
            return;
        }
        $write_result = syswrite $out_fh, $genome_ref->{$chr}->Block_Read();
        unless(defined $write_result || $write_result != length($chr_write_string)) {
            $self->error_message("Error writing the header");
            return;
        }
    }
    $out_fh->close;
    return 1;
}

sub read_genome_bitmask {
    my ($self,$filename) = @_;
    unless($filename) {
        $self->error_message("No filename of file to write to");
        return;
    }
    #do some stuff to read this from a file without making it suck
    my $in_fh = IO::File->new($filename,"<:raw");
    unless($in_fh) {
        $self->error_message("Unable to read from " . $filename);
        return;
    }
    my $read_string;
    sysread $in_fh, $read_string, 4;
    my $header_length = unpack "N", $read_string;
    sysread $in_fh, $read_string, $header_length;
    my $header_string = unpack "a*",$read_string;
    my %genome = split /\t/, $header_string; #each key is the name, each value is the size in bits

    #now read in each one
    foreach my $chr (sort keys %genome) {
        $genome{$chr} = Bit::Vector->new($genome{$chr}); #this throws an exception if it fails. Probably should be trapped at some point in the future
        sysread $in_fh, $read_string, 4;
        my $chr_byte_length = unpack "N", $read_string;
        my $chr_read_string;
        sysread $in_fh, $chr_read_string, $chr_byte_length;
        $genome{$chr}->Block_Store($chr_read_string);
    }
    $in_fh->close;
    return \%genome;
}


1;

sub help_brief {
    "This reads in a WIG format file and converts it to a bitmask representation of the whole genome."
}

sub help_detail {
    <<'HELP';
    This script takes a file in WIG format, and uses it to set bits in a bit vector representation of the genome. This may be useful for quickly querying whether locations or regions overlap with others on a genome wide basis. By default, the bitmask is stored on the object. Thus calling this script with no --output-file option from the command line will result in no output. Very sad. This script was designed for coverage files, Thus non-zero values in the WIG file are flipped on in the bitmask and 0 values are left off. If you want the reverse behavior, use the --exclude-regions option. This script does not support track definitions nor the variableStep specification. See https://gscweb.gsc.wustl.edu/wiki/File_Format/Wiggle_Track_Format_(WIG_-_WTF) for more information on the wiggle format.
    The bitmask itself is a hash reference to a hash with chromosome names as keys and the value a Bit::Vector object where each position in the chromosome is represented by a bit. The length of each chromosome is determined from the provided index file. 
HELP
}


