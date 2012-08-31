package Genome::Model::Tools::Bmr::CompareWiggles;

use strict;
use warnings;

use Genome;
use IO::File;
use Bit::Vector;

class Genome::Model::Tools::Bmr::CompareWiggles {
    is => 'Genome::Command::Base',
    has_input => [
    refseq_build_name => {
        is => 'String',
        is_optional => 1,
        default => 'NCBI-human-build36',
        doc => 'The reference sequence build, used to gather and generate bitmask files for base masking and sample coverage.',
    },
    roi_bedfile => {
        type => 'String',
        is_optional => 0,
        doc => 'BED file used to limit background regions of interest when calculating background mutation rate',
    },
    wiggle_file_list => {
        type => 'Csv',
        is_optional => 0,
        doc => 'Wiggle files detailing genome-wide coverage of each sample in dataset, separated by commas',
    },
    output_wiggle => {
        type => 'String',
        is_optional => 0,
        doc => 'In this case, a wiggle file which shows what sites were "covered" by at least one of the input wiggle files',
    },
#    rejected_mutations => {
#        type => 'String',
#        is_optional => 1,
#        doc => 'File to catch mutations that fall in the ROI list location-wise, but have a gene name which does not match any of the genes in the ROI list. Default operation is to print to STDOUT.',
#    },
    ]
};

sub help_brief {
    "Find positions covered by at least one of the input wiggle files"
}

sub help_detail {
    "Find positions covered by at least one of the input wiggle files"
}

sub execute {
    my $self = shift;

    #resolve refseq
    my $ref_build_name = $self->refseq_build_name;
    my ($ref_model_name,$ref_build_version) = $ref_build_name =~ /^(\S+)-build(\S*)$/;
    my $ref_model = Genome::Model->get(name=>$ref_model_name);
    my $ref_build = $ref_model->build_by_version($ref_build_version);

    #create an empty genome bitmask to record results of comparison
    my $roi_bitmask = $self->create_empty_genome_bitmask($ref_build);
    
    #load ROIs into a hash %ROIs -> chr -> gene -> start = stop;
    my %ROIs;
    my $roi_bedfile = $self->roi_bedfile;
    my $bed_fh = new IO::File $roi_bedfile,"r";
    while (my $line = $bed_fh->getline) {
        chomp $line;
        my ($chr,$start,$stop,$exon_id) = split /\t/,$line;
        (my $gene = $exon_id) =~ s/^([^\.]+)\..+$/$1/;
        if ($chr eq "M") { $chr = "MT"; } #for broad roi lists
        if (exists $ROIs{$chr}{$gene}{$start}) {
            next if $stop < $ROIs{$chr}{$gene}{$start};
        }
        $ROIs{$chr}{$gene}{$start} = $stop;
    }
    $bed_fh->close;

    #Parse wiggle files
    my $wiggle_list = $self->wiggle_file_list;
    my @wiggle_files = split ",",$wiggle_list;

    #loop through samples to calculate class coverages
    my $bitmask_conversion = Genome::Model::Tools::Bmr::WtfToBitmask->create(
        reference_index => $ref_build->full_consensus_sam_index_path,
    );

    for my $file (@wiggle_files) {

        #create sample coverage bitmask;
        $bitmask_conversion->wtf_file($file);
        if ($bitmask_conversion->is_executed) {
            $bitmask_conversion->is_executed('0');
        }
        $bitmask_conversion->execute;

        unless ($bitmask_conversion) {
            $self->error_message("Not seeing a bitmask object for file $file.");
            return;
        }
        my $cov_bitmask = $bitmask_conversion->bitmask;
        unless ($cov_bitmask) {
            $self->error_message("Not seeing cov_bitmask from bitmask_cov->bitmask.");
            return;
        }

    }#end, for my wiggle file

    #clean up the object memory
    $bitmask_conversion->delete;
    undef $bitmask_conversion;

    return 1;
}

sub create_empty_genome_bitmask {
    my $self = shift;
    my $ref_build = shift;
    my %genome;
    my $ref_index_file = $ref_build->full_consensus_sam_index_path;
    my $ref_fh = new IO::File $ref_index_file,"r";
    while (my $line = $ref_fh->getline) {
        chomp $line;
        my ($chr,$length) = split /\t/,$line;
        $genome{$chr} = Bit::Vector->new($length + 1); #adding 1 for 1-based coordinates
    }
    $ref_fh->close;
    return \%genome;
}

sub read_genome_bitmask {
    my ($self,$filename) = @_;
    unless($filename) {
        $self->error_message("File $filename not found.");
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

sub count_bits {
    my ($self,$vector,$start,$stop) = @_;
    my $count = 0;
    for my $pos ($start..$stop) {
        if ($vector->bit_test($pos)) {
            $count++;
        }
    }
    return $count;
}

1;



#more specific hash structure
#%COVMUTS->gene->class->(coverage,#mutations)
#
#classes
#
#CG.C.transit.T   CG.G.transit.A
#CG.C.transver.A CG.G.transver.T
#CG.C.transver.G CG.G.transver.C
#
#CpG.C.transit.T   CpG.G.transit.A
#CpG.C.transver.A CpG.G.transver.T
#CpG.C.transver.G CpG.G.transver.C
#
#AT.A.transit.G   AT.T.transit.C
#AT.A.transver.T AT.T.transver.A
#AT.A.transver.C AT.T.transver.G
#
#and Indels
