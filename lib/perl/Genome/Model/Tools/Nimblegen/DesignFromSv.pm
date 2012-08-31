package Genome::Model::Tools::Nimblegen::DesignFromSv;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;
use Sort::Naturally qw( nsort );
use Genome::Sys;

class Genome::Model::Tools::Nimblegen::DesignFromSv {
    is => 'Command',
    has => [
    sv_file => {
        type => 'String', is_optional => 0,
        doc => "A HQfiltered formatted file of SV sites to generate probe regions for. Assumes STDIN if not specified",
    },
    output_file => {
        type => 'String', is_optional => 1,
        doc => "Output file. Assumes STDOUT if not specified",
    },
    assembly_format => {
        type => 'Boolean', is_optional => 1, default => 0,
        doc => "input file is assembly format",
    },
    span => {
        type => 'Integer', is_optional => 1, default => 200,
        doc => "The region to be spanned",
    },
    include_mitochondrial_sites => {
        type => 'Boolean', is_optional => 1, default => 0,
        doc => "Whether or not to remove sites on the mitochondria or non-chromosomal contigs",
    },
    include_unplaced_contig_sites => {
        type => 'Boolean', is_optional => 1, default => 1,
        doc => "Whether or not to remove sites on the unplaced contigs of the chromosome",
    },
    include_y_chrom_sites => {
        type => 'Boolean', is_optional => 1, default => 1,
        doc => "Whether or not to include sites on the Y chromosome in the output (if cases are all female)",
    },
    reference_index => {
        type => 'String', is_optional => 0,
        default => "/gscmnt/839/info/medseq/reference_sequences/NCBI-human-build36/all_sequences.fa.fai",
        doc => "samtools index of the reference sequence",
    },
    resolution => {
        type => 'Integer', is_optional => 1, default => 10000,
        doc => "Filter out the resolution > this number and not output it to nimblegen list."
    },
    count_file => {
        type => 'String', is_optional => 1,
        doc => "Count the whole bases to be covered."
    },
    filtered_out_file => {
        type => 'String',
        is_optional => 1,
        doc => "Save those in .capture but not in .nimblegen. Writes to STDERR if undefined"
    },
    ]
};

sub execute {
    my $self=shift;
    my $reference_index = $self->reference_index;
    my $include_mitochondrial_sites = $self->include_mitochondrial_sites;
    my $include_unplaced_contig_sites = $self->include_unplaced_contig_sites;
    my $include_y_chrom_sites = $self->include_y_chrom_sites;

    # Hash to help convert hg19 unplaced contig names to the UCSC equivalent
    my %ucsc_unplaced_contigs = (( 'GL000207.1', 'chr18_gl000207_random' ), ( 'GL000226.1', 'chrUn_gl000226' ), ( 'GL000229.1', 'chrUn_gl000229' ), ( 'GL000231.1', 'chrUn_gl000231' ), ( 'GL000210.1', 'chr21_gl000210_random' ), ( 'GL000239.1', 'chrUn_gl000239' ), ( 'GL000235.1', 'chrUn_gl000235' ), ( 'GL000201.1', 'chr9_gl000201_random' ), ( 'GL000247.1', 'chrUn_gl000247' ), ( 'GL000245.1', 'chrUn_gl000245' ), ( 'GL000197.1', 'chr8_gl000197_random' ), ( 'GL000203.1', 'chr17_gl000203_random' ), ( 'GL000246.1', 'chrUn_gl000246' ), ( 'GL000249.1', 'chrUn_gl000249' ), ( 'GL000196.1', 'chr8_gl000196_random' ), ( 'GL000248.1', 'chrUn_gl000248' ), ( 'GL000244.1', 'chrUn_gl000244' ), ( 'GL000238.1', 'chrUn_gl000238' ), ( 'GL000202.1', 'chr11_gl000202_random' ), ( 'GL000234.1', 'chrUn_gl000234' ), ( 'GL000232.1', 'chrUn_gl000232' ), ( 'GL000206.1', 'chr17_gl000206_random' ), ( 'GL000240.1', 'chrUn_gl000240' ), ( 'GL000236.1', 'chrUn_gl000236' ), ( 'GL000241.1', 'chrUn_gl000241' ), ( 'GL000243.1', 'chrUn_gl000243' ), ( 'GL000242.1', 'chrUn_gl000242' ), ( 'GL000230.1', 'chrUn_gl000230' ), ( 'GL000237.1', 'chrUn_gl000237' ), ( 'GL000233.1', 'chrUn_gl000233' ), ( 'GL000204.1', 'chr17_gl000204_random' ), ( 'GL000198.1', 'chr9_gl000198_random' ), ( 'GL000208.1', 'chr19_gl000208_random' ), ( 'GL000191.1', 'chr1_gl000191_random' ), ( 'GL000227.1', 'chrUn_gl000227' ), ( 'GL000228.1', 'chrUn_gl000228' ), ( 'GL000214.1', 'chrUn_gl000214' ), ( 'GL000221.1', 'chrUn_gl000221' ), ( 'GL000209.1', 'chr19_gl000209_random' ), ( 'GL000218.1', 'chrUn_gl000218' ), ( 'GL000220.1', 'chrUn_gl000220' ), ( 'GL000213.1', 'chrUn_gl000213' ), ( 'GL000211.1', 'chrUn_gl000211' ), ( 'GL000199.1', 'chr9_gl000199_random' ), ( 'GL000217.1', 'chrUn_gl000217' ), ( 'GL000216.1', 'chrUn_gl000216' ), ( 'GL000215.1', 'chrUn_gl000215' ), ( 'GL000205.1', 'chr17_gl000205_random' ), ( 'GL000219.1', 'chrUn_gl000219' ), ( 'GL000224.1', 'chrUn_gl000224' ), ( 'GL000223.1', 'chrUn_gl000223' ), ( 'GL000195.1', 'chr7_gl000195_random' ), ( 'GL000212.1', 'chrUn_gl000212' ), ( 'GL000222.1', 'chrUn_gl000222' ), ( 'GL000200.1', 'chr9_gl000200_random' ), ( 'GL000193.1', 'chr4_gl000193_random' ), ( 'GL000194.1', 'chr4_gl000194_random' ), ( 'GL000225.1', 'chrUn_gl000225' ), ( 'GL000192.1', 'chr1_gl000192_random' ));

    # Depending on the bool flags the user sets, sites on some chromosomes will be excluded
    my %valid_chrs = map{ chomp; $_ => 1 } `cut -f 1 $reference_index`; # Used to check input files for valid ref names
    my %include_chrs = ( $include_unplaced_contig_sites ? %valid_chrs : ( map{ $_ => 1 } ( 1..22, qw( X Y MT ))));
    delete $include_chrs{MT} unless( $include_mitochondrial_sites );
    delete $include_chrs{Y} unless( $include_y_chrom_sites );

    my $fh = IO::File->new($reference_index,"r");
    unless($fh) {
        $self->error_message("Unable to open the reference sequence index: $reference_index");
        return;
    }

    #read in the index to get the chromosome lengths
    my %chromosome_lengths;
    while(my $line = $fh->getline) {
        chomp $line;
        my ($chr, $length) = split /\t/, $line;
        $chromosome_lengths{$chr} = $length;
    }
    $fh->close;

    my $output_fh;
    if(defined $self->output_file) {
        $output_fh = IO::File->new($self->output_file,"w+");
        unless($output_fh) {
            $self->error_message("Unable to open file " . $self->output_file . " for writing.");
            return;
        }
    }
    else {
        $output_fh = IO::File->new_from_fd(fileno(STDOUT),"w");
        unless($output_fh) {
            $self->error_message("Unable to open STDOUT for writing.");
            return;
        }
    }

    my $count_fh;
    if(defined $self->count_file) {
        $count_fh = IO::File->new($self->count_file,"a+");
        unless($count_fh) {
            $self->error_message("Unable to open file " . $self->count_file . " for writing.");
            return;
        }
    }
    else {
        $count_fh = IO::File->new_from_fd(fileno(STDOUT), "w");
        unless($count_fh){
            $self->error_message("Unable to open STDOUT for writing.");
            return;
        }
    }

    my $filtered_out_fh;
    if(defined $self->filtered_out_file) {
        $filtered_out_fh = IO::File->new($self->filtered_out_file,"w");
        unless($filtered_out_fh) {
            $self->error_message("Unable to open file ". $self->filtered_out_file . " for writing.");
            return;
        }
    }
    else {
        $filtered_out_fh = IO::File->new_from_fd(fileno(STDERR), "w");
        unless($filtered_out_fh){
            $self->error_message("Unable to open STDERR for writing.");
            return;
        }
    }

    my $input_fh;
    if(defined $self->sv_file) {
        $input_fh = IO::File->new($self->sv_file,"r");
        unless($input_fh) {
            $self->error_message("Unable to open file ". $self->sv_file . " for reading.");
            return;
        }
    }
    else {
        $input_fh = IO::File->new_from_fd(fileno(STDIN),"r");
        unless($input_fh) {
            $self->error_message("Unable to open STDIN for reading.");
            return;
        }
    }
    my %cover = ();
    while(my $line = $input_fh->getline) {
        next if $line =~ /^#/;  #skip comments
        chomp $line;

        if($self->assembly_format) {
            my @list = split(/\t/,$line);
            my @sub_list = @list[1 .. 6];
            my $new_col = join(".",@sub_list);
            $line = "$new_col\t".$line;
        }

        #my ($id,$chr1,$outer_start,$inner_start,$chr2,$inner_end,$outer_end,$type,$orient, $minsize) = split /\s+/, $line;
        my ($id, )=split("\t", $line);
        my ($chr1,$outer_start,$inner_start,$chr2,$inner_end,$outer_end) = ($id =~ /(\S+)\.(-*\d+)\.(-*\d+)\.(\S+)\.(-*\d+)\.(-*\d+)/);
        if(!defined $chr1 || !defined $chr2) {
            print "$line\n";
        }
        unless( defined $include_chrs{$chr1} && defined $include_chrs{$chr2} ) {
            printf $filtered_out_fh "Chr to ignore: %s\n", $line;
            next;
        }

        my $outer_start_ = $outer_start - $self->span;
        my $inner_start_ = $inner_start + $self->span;
        my $inner_end_ = $inner_end - $self->span;
        my $outer_end_ = $outer_end + $self->span;

        if($outer_start - $self->span < 1) {
            $outer_start_ = 1;
        }
        if($outer_start - $self->span > $chromosome_lengths{$chr1}) {
            printf $filtered_out_fh "Out of bounds: %s\n", $line;
            next;
        }
        if($inner_start + $self->span < 1) {
            $inner_start_ = 1;
        }
        if($inner_start + $self->span > $chromosome_lengths{$chr1}) {
            $inner_start_ = $chromosome_lengths{$chr1};
        }
        if($outer_end + $self->span < 1){
            printf $filtered_out_fh "Out of bounds: %s\n", $line;
            next;
        }
        if($outer_end + $self->span > $chromosome_lengths{$chr2}) {
            $outer_end_ = $chromosome_lengths{$chr2};
        }
        if($inner_end - $self->span < 1){
            $inner_end_ = 1;
        }
        if($inner_end - $self->span > $chromosome_lengths{$chr2}) {
            $inner_end_ = $chromosome_lengths{$chr2};
        }

        # filter out those resolution > 2k
        if($inner_start_ - $outer_start_ > $self->resolution || $outer_end_ - $inner_end_ > $self->resolution) {
            printf $filtered_out_fh "Resolution too high: %s\n", $line;
            next;
        }

        # record how many base pair has been covered
        for(my $i = $outer_start_; $i <= $inner_start_; $i++) {
            ${$cover{$chr1}}{$i} = 1 if(! defined $cover{$chr1}{$i});
        }
        for(my $i = $inner_end_; $i <= $outer_end_; $i++) {
            ${$cover{$chr1}}{$i} = 1 if(! defined $cover{$chr1}{$i});
        }

        #if($type !~ /INS/) {#|| ($type =~ /DEL/ && $minsize > 1000)) {
        $chr1 = $ucsc_unplaced_contigs{$chr1} if( defined $ucsc_unplaced_contigs{$chr1} );
        $chr2 = $ucsc_unplaced_contigs{$chr2} if( defined $ucsc_unplaced_contigs{$chr2} );
        $chr1 = "chr$chr1" unless( $chr1 =~ m/^chr/ );
        $chr2 = "chr$chr2" unless( $chr2 =~ m/^chr/ );
        printf $output_fh "%s\t%d\t%d\t%d\t%s\n",$chr1,$outer_start_, $inner_start_, (($inner_start_) - ($outer_start_)), $line;
        printf $output_fh "%s\t%d\t%d\t%d\t%s\n",$chr2,$inner_end_, $outer_end_, (($outer_end_) - ($inner_end_)), $line;
        #}
        #else {
        #    printf $output_fh "chr%s\t%d\t%d\t%d\t%s\n",$chr1,$outer_start - 100, $outer_end + 100, (($outer_end + 100) - ($outer_start - 100)), $line;
        #}
    }

    my $inall = 0;
    my $chr;
    my $base;
    foreach $chr (keys %cover) {
        foreach $base (keys %{$cover{$chr}}) {
            $inall++;
        }
    }
    printf $count_fh "%s\t%d bps covered\n", $self->output_file, $inall;

    return 1;
}

sub help_brief {
    "Takes an Sv file and produces a list of regions to target for validation.";
}

sub help_detail {
    return <<EOS
Takes an SV file and produces a list of regions to target for validation via Nimblegen Solid Phase Capture Array.
EOS
}

1;
