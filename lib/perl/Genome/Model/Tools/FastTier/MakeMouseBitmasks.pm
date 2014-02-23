package Genome::Model::Tools::FastTier::MakeMouseBitmasks;

use strict;
use warnings;
use Bit::Vector;
use Genome;
use UR;
use IO::File;
use Data::Dumper;
use Set::IntSpan;

class Genome::Model::Tools::FastTier::MakeMouseBitmasks {
    is => 'Command',
    has => [
        output_directory => {
            type => 'Text',
            is_input => 1,
            doc => 'Location for tier bitmasks to be dropped',
        },
        reference_sequence => {
            type => 'Text',
            is_input => 1,
            doc => 'Reference sequence to use for tier mask creation',
        },
        tier1_bed => {
            is => 'Text',
            doc => 'Bed file containing all tier 1 regions',
        },
        tier1_output => {
            calculate_from => ['output_directory'],
            calculate => q{ "$output_directory/tier1.bitmask"; },
            is_output => 1,
        },
        tier2_output => {
            calculate_from => ['output_directory'],
            calculate => q{ "$output_directory/tier2.bitmask"; },
            is_output => 1,
        },
        tier3_output => {
            calculate_from => ['output_directory'],
            calculate => q{ "$output_directory/tier3.bitmask"; },
            is_output => 1,
        },
        tier4_output => {
            calculate_from => ['output_directory'],
            calculate => q{ "$output_directory/tier4.bitmask"; },
            is_output => 1,
        },
    ],
};


sub help_brief {
    "Used to generate bitmask files of teirs 1 through 4"
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
genome-model tools fast-tier make-tier-bitmask ...    
EOS
}

sub execute {
    my $self = shift;

    my $ref = $self->reference_sequence;
    my $ref_list_fh = Genome::Sys->open_file_for_reading("$ref.fai");
    my @chromosomes;

    while(my $line = $ref_list_fh->getline) {
        chomp $line;
        my ($chr) = split /\t/, $line;
        push @chromosomes, $chr; 
    }
    $ref_list_fh->close;

    my %genome;
    my $genome_size = 0;
    my $masked_genome_size = 0;
    for my $ref_chr (@chromosomes) {
        unless(open(FAIDX,"samtools faidx $ref $ref_chr |")) {
            die "Couldn't pipe samtools faidx\n";
        }
        my $header = <FAIDX>;
        my $chr = $ref_chr;
        my $chr_length = 0;
        my $cur_nstart = 0;
        my $cur_nstop = -1;
        my @n_blocks;
        while(my $line = <FAIDX>) {
            chomp $line;
            my @sequence = split //, $line;  
            while(my $base = shift @sequence) {
                $chr_length++;
                if($base !~ /[ACTG]/i) {
                    if($cur_nstart) {
                        $cur_nstop = $chr_length;
                    }
                    else {
                        $cur_nstart = $chr_length;
                        $cur_nstop = $chr_length;
                    }
                }
                else {
                    if($cur_nstart) {
                        push @n_blocks, [$cur_nstart, $cur_nstop];
                        $cur_nstart = 0;
                        $cur_nstop = -1;
                    }
                }
            }
        }
        $genome{$chr} = Bit::Vector->new($chr_length);
        while(my $interval = shift @n_blocks) {
            $genome{$chr}->Interval_Fill($interval->[0]-1,$interval->[1]-1);
        }
        $genome_size += $chr_length;
        $masked_genome_size += $chr_length - $genome{$chr}->Norm;
        unless(close(FAIDX)) {
            die "Error reading from samtools faidx pipe\n";
        }
    }


    #this makes UR retain fewer objects. Perhaps this will help with memory issues
    my $low=20000;
    my $high=100000;
    UR::Context->object_cache_size_lowwater($low);
    UR::Context->object_cache_size_highwater($high);

    #This script calculates the number of the bases in the genome covered by each current tier definition (as of 5/26/2009)

    #Tier 1 
    #Tier 1 contains all coding alterations and alterations to rna genes. To calculate its coverage we will scan through the transcript table and add coding exon bases and rna transcript bases to the set. RNA are tracked separately because Tier2 also contains the bases in coding transcripts.


    printf "Calculated genome size is %u\n", $genome_size;
    printf "Masked genome size is %u\n", $masked_genome_size;

    my $tier1 = $self->shadow_genome(\%genome);

    my $fh = Genome::Sys->open_file_for_reading($self->tier1_bed);
    while(my $line = $fh->getline){
        chomp $line;
        my ($chr,$start,$stop) = split /\s+/, $line;
        $self->add_region($tier1,$chr,$start,$stop);
    }
    $fh->close;


    my $tier2 = $self->shadow_genome(\%genome);
    $self->add_region($tier2,"1","0","1");
    
    my $tier3 = $self->complement_genome($tier1);
    my $tier4 = $self->shadow_genome(\%genome);
    $self->add_region($tier2,"1","0","1");

    print "Tier 1 contains ".$self->bases_covered($tier1)." bases.\n"; 
    print "Tier 2 contains ".$self->bases_covered($tier1)." bases.\n"; 
    print "Tier 3 contains ".$self->bases_covered($tier1)." bases.\n"; 
    print "Tier 4 contains ".$self->bases_covered($tier1)." bases.\n"; 

    $self->write_genome_bitmask($self->tier1_output, $tier1);
    $self->write_genome_bitmask($self->tier2_output, $tier2);
    $self->write_genome_bitmask($self->tier3_output, $tier3);
    $self->write_genome_bitmask($self->tier4_output, $tier4);

    return 1;
}

sub add_region {
    my $self = shift;
    my ($set,$chr, $start, $stop) = @_;
    if($start > $stop){
        $self->debug_message(" found a record with start > stop! ");
        ($stop,$start) = ($start,$stop);
    }
    unless(exists( $set->{$chr})){
        $self->debug_message(" found a chromosome that wasn't in the hash: ".$chr);
        return;
    }
    if($start >=0 && ($stop < $set->{$chr}->Size-1)){
        $set->{$chr}->Interval_Fill($start,$stop);    
    } else {
        $self->debug_message("Found a weird record off the edge of a chrom");
        return;
    }
}

sub write_genome_bitmask {
    my $self = shift;
    my ($filename,$genome_ref) = @_;
    unless($filename) {
        die("No filename of file to write to");
    }
    unless($genome_ref) {
        die("No bitmask to write to file");
    }
    #do some stuff to write this to a file without making it suck
    my $out_fh = IO::File->new($filename,">:raw");
    unless($out_fh) {
        die("Unable to write to " . $filename);
    }
    my $header_string = join("\t", map {$_ => $genome_ref->{$_}->Size()} sort keys %$genome_ref);
    my $write_string = pack 'N/a*', $header_string;
    my $write_result = syswrite($out_fh,$write_string);
    unless(defined $write_result && $write_result == length($write_string)) {
        die("Error writing the header");
    }
    for my $chr (sort keys %$genome_ref) {
        #first write the length in bytes 
        my $chr_write_string = $genome_ref->{$chr}->Block_Read();
        $write_result = syswrite $out_fh, pack("N",length($chr_write_string));
        unless(defined $write_result || $write_result != 4) {
            die("Error writing the length of chromosome $chr");
        }
        $write_result = syswrite $out_fh, $genome_ref->{$chr}->Block_Read();
        unless(defined $write_result || $write_result != length($chr_write_string)) {
            die("Error writing the header");
        }
    }
    $out_fh->close;
    return 1;
}

sub shadow_genome {
    my $self = shift;
    my $genome = shift;
    my %new;
    for my $chr (keys %$genome) {
        $new{$chr} = $genome->{$chr}->Shadow;
    }
    return \%new;
}

sub union_genomes {
    my $self = shift;
    my ($genome1, $genome2) = @_;
    my %union;
    for my $chr (keys %$genome1) {
        next unless defined $genome1->{$chr};
        $union{$chr} = $genome1->{$chr}->Clone;
        if(exists($genome2->{$chr})) {
            $union{$chr}->Union($genome2->{$chr},$union{$chr});  
        }
    }
    for my $chr (keys %$genome2) {
        next unless defined $genome2->{$chr};
        if(!exists($genome1->{$chr})) {
            $union{$chr} = $genome2->{$chr}->Clone;  
        }
    }
    return \%union;
}

sub in_place_union_genomes {
    my $self = shift;
    my ($genome1, $genome2) = @_;
    for my $chr (keys %$genome1) {
        next unless defined $genome1->{$chr};
        if(exists($genome2->{$chr})) {
            $genome1->{$chr}->Union($genome1->{$chr},$genome2->{$chr});  
        }
    }
    for my $chr (keys %$genome2) {
        next unless defined $genome2->{$chr};
        if(!exists($genome1->{$chr})) {
            $genome1->{$chr} = $genome2->{$chr}->Clone;  
        }
    }
}

sub difference_genomes {
    my $self = shift;
    my ($genome1, $genome2) = @_;
    my %difference;
    for my $chr (keys %$genome1) {
        next unless defined $genome1->{$chr};
        $difference{$chr} = $genome1->{$chr}->Clone;
        if(exists($genome2->{$chr})) {
            $difference{$chr}->Difference($difference{$chr},$genome2->{$chr});  
        }
    }
    return \%difference;
}
sub in_place_difference_genomes {
    my $self = shift;
    my ($genome1, $genome2) = @_;
    for my $chr (keys %$genome1) {
        next unless defined $genome1->{$chr};
        if(exists($genome2->{$chr})) {
            $genome1->{$chr}->Difference($genome1->{$chr},$genome2->{$chr});  
        }
    }
}
sub complement_genome {
    my $self = shift;
    my ($genome) = @_;
    my %result;
    for my $chr (keys %$genome) {
        next unless defined $genome->{$chr};
        $result{$chr} = $genome->{$chr}->Clone;
        $result{$chr}->Complement($result{$chr});#in-place calc. Perhaps more mem efficient
    }
    return \%result;
}
sub in_place_complement_genome {
    my $self = shift;
    my ($genome) = @_;
    for my $chr (keys %$genome) {
        next unless defined $genome->{$chr};
        $genome->{$chr}->Complement($genome->{$chr});#in-place calc. Perhaps more mem efficient
    }
}


sub bases_covered {
    my $self = shift;
    my ($genome) = @_;
    my $total = 0;
    for my $chr (keys %$genome) {
        next unless defined $genome->{$chr};
        $total += $genome->{$chr}->Norm();
    }
    return $total;
}
