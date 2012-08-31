package Genome::Model::ReferenceSequence::Command::GenerateBitmask; 

use warnings;
use strict;
use Genome;
use IO::File;
use Bit::Vector;

class Genome::Model::ReferenceSequence::Command::GenerateBitmask {
    is => 'Genome::Command::Base',
    has_input => [
    refseq => { 
        is => 'Genome::Model::Build::ImportedReferenceSequence', # we're going to remove "Imported" soon
        id_by => 'refseq_id',
        shell_args_position => 1, 
        doc => 'the reference build, specified by name (like NCBI-human-build36)'
    },
    bases => {
        is => 'Text',
        shell_args_position => 2,
        doc => 'List of bases to bitmask.',
        valid_values => ['AT','CG','CpG','GC'],
    },
    ],
    doc => 'Generate a set of bitmask data for a given specific reference sequence.'
};

sub help_synopsis {
    my $class = shift;
    return <<EOS;
genome model reference-sequence generate-bitmask NCBI-human-build36 [AT|CG|CpG|GC]

genome model reference-sequence generate-bitmask -r NCBI-human-build36 [AT|CG|CpG|GC]

genome model reference-sequence generate-bitmask MYMODELNAME-buildMYBUILDVERSION [AT|CG|CpG|GC]

EOS
}

sub help_detail {

    my $class = shift;
    return <<'EOS';
Creates a bitmask file for a given reference sequence and stores it with the data for that reference.
Description of each "bases" option:
AT = All A & T nucleotides
CG = All C & G nucleotides minus CpG
CpG = All CpG islands
GC = All G & C nucleotides
EOS
}

sub execute {
    my $self = shift;
    my $refseq = $self->refseq;
    my $bases = $self->bases;
    #unless ($bases =~ /^(AT|CG|CpG)$/) {
    #    $self->error_message("Current options for bitmask files are AT, CG, or CpG.");
    #    return;
    #}
    $self->status_message("Building a $bases bitmask files for refseq " . $refseq->__display_name__ . "...");

    #Calculate bitmask filename and check for current existence
    my $data_directory = $refseq->data_directory;
    my $fasta_file = $data_directory . "/all_sequences.fa";
    my $index_file = $data_directory . "/all_sequences.fa.fai";
    print "data_dir: $data_directory\n"; #REMOVE LATER
    print "fasta_path: $fasta_file\n"; #REMOVE LATER
    print "index_path: $index_file\n"; #REMOVE LATER
    my $filename = 'all_sequences.' . $bases . '_bitmask';
    my $full_path_filename = $data_directory . '/' . $filename;
    if (-e $full_path_filename) {
        $self->error_message("$full_path_filename already exists!");
        return;
    }
    $self->status_message("Final file will be at $full_path_filename.");

    #Create list of chromosomes in reference sequence
    my @chrs;
    my $index_fh = new IO::File $index_file,"r";
    while (my $line = $index_fh->getline) {
        my ($chr,$length) = split /\t/,$line;
        push @chrs,$chr;
    }

    #Perform work in temp directory
    my $temp_dir = Genome::Sys->create_temp_directory();
    my $temp_path_filename = $temp_dir . '/' . $filename;

    #set up hash to record bitmasks
    my %genome;
    
    #grab reference sequence for each chromosome using samtools faidx
    for my $chr (@chrs) {
        my $cmd = "samtools faidx $fasta_file $chr |";
        print $cmd . "\n";
        open(REF,$cmd);
        my $header = <REF>;
        unless ($header =~ /^>/) {
            $self->error_message("Header line did not start with '>' as expected. Header: $header");
            return;
        }
        my @sequence_array = <REF>;
        my $seq = join("",@sequence_array);
        $seq =~ s/\n//g;
        my $length = length($seq);
        print $length . "\n";

        #apply the bitmask
        my $mask;
        if ($bases eq "AT") {
            ($mask = $seq) =~ s/[AT]/1/g;
        }
        elsif ($bases eq "CG") {
            ($mask = $seq) =~ s/CG/00/g;
            $mask =~ s/[CG]/1/g;
        }
        elsif ($bases eq "CpG") {
            ($mask = $seq) =~ s/CG/11/g;
        }
        elsif ($bases eq 'GC') {
            ($mask = $seq) =~ s/[GC]/1/g;
        }
        else {
            $self->error_message("How in the heck did you get down to here - this bitmask is not one of the current options.");
            return;
        }
        $mask =~ s/[^1]/0/g;

        #create the appropriate bit vector
        $mask = "0" . $mask; #add zero in first element for 1-based indexing
        my $revmask = reverse($mask); #Bit::Vector->new_Bin reverses order for significance purposes

        $length = length($revmask);
        print "$length\n";
        $genome{$chr} = Bit::Vector->new_Bin(length($mask),$revmask);

    }#end, @chrs loop

    #write bitmask to file
    $self->write_genome_bitmask($temp_path_filename,\%genome);

    #Copy the results to real disk after completion and reallocate
    Genome::Sys->copy_file($temp_path_filename,$full_path_filename);
    $self->status_message("Resizing the build disk allocation...");
    $refseq->reallocate();

    #Finished!
    $self->status_message("Tasks Completed.");

    return 1;
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

1;
