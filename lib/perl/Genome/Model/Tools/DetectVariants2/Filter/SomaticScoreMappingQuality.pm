package Genome::Model::Tools::DetectVariants2::Filter::SomaticScoreMappingQuality;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;
use Readonly;
use Genome::Info::IUB;

my $DEFAULT_VERSION = '0.2';
my $READCOUNT_COMMAND = 'bam-readcount';

class Genome::Model::Tools::DetectVariants2::Filter::SomaticScoreMappingQuality{
    is => 'Genome::Model::Tools::DetectVariants2::Filter',
    has => [
        bam_readcount_version => {
            is => 'Version',
            is_input=>1, 
            is_optional => 1,
            default_value => $DEFAULT_VERSION,
            doc => "Version of bam-readcount to use, default is $DEFAULT_VERSION"
        },
        bam_readcount_params => {
            is => 'String',
            is_optional => 1,
            is_input=>1, 
            default_value => "-q 1",
            doc => "Parameters to pass to bam-readcount"
        },
        min_mapping_quality => {
            type => 'String',
            default => '40',
            is_optional => 1,
            is_input => 1,
            doc => 'minimum average mapping quality threshold for high confidence call',
        },
        min_somatic_score => {
            type => 'String',
            default => '40',
            is_optional => 1,
            is_input => 1,
            doc => 'minimum somatic quality threshold for high confidence call',
        },
    ],
};

my %READCOUNT_VERSIONS = (
    '0.2' => $ENV{GENOME_SW} . '/samtools/readcount/readcount-v0.2/' . $READCOUNT_COMMAND,
);

sub help_brief {
    return "This module takes in somatic sniper output and filters it to high confidence variants";
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
    gmt filter-variants high-confidence --sniper-file sniper,out --tumor-bam somefile.bam 
EOS
}

sub help_detail {                           
    return <<EOS 
This module takes in somatic sniper output and filters it to high confidence variants
EOS
}

sub _variant_type { 'snvs' };

sub _filter_variants {
    my $self = shift;

    my $tumor_bam_file = $self->aligned_reads_input;
    my $input_file = $self->input_directory."/snvs.hq.bed";

    my $output_file = $self->_temp_staging_directory."/snvs.hq.bed";
    # Writes into temp, then sorts into where it should be
    my $lq_output_file = $self->_temp_scratch_directory."/snvs.lq.bed";
    my $sorted_lq_output_file = $self->_temp_staging_directory."/snvs.lq.bed";

    #test architecture to make sure we can run read count program
    unless (`uname -a` =~ /x86_64/) {
       $self->error_message("Must run on a 64 bit machine");
       die;
    }

    #check on BAM file
    unless(-e $tumor_bam_file) {
        $self->error_message("Tumor bam file: " . $tumor_bam_file . " does not exist");
        die;
    }

    my $fh = Genome::Sys->open_file_for_reading($input_file);
    my $ofh = Genome::Sys->open_file_for_writing($output_file);
    my $lq_output_fh = Genome::Sys->open_file_for_writing($lq_output_file);

    my ($tfh,$temp_path) = Genome::Sys->create_temp_file;
    $temp_path =~ s/\:/\\\:/g;

    #read sniper and skip indels
    my $somatic_threshold = $self->min_somatic_score;
    my @sniper_lines;
    while(my $line = $fh->getline) {
        chomp $line;
        my ($chr, $start, $stop, $ref_iub, $somatic_score, $tumor_read_depth) = split /\t/, $line;
        my ($ref, $iub) = split "/", $ref_iub;
        #filtering out thing with a ref = *
        if( $ref eq "*"){
            print $lq_output_fh $line."\n";;
            next;
        }
        # here we are filtering out lines which do not meet the somatic-score threshold
        if($somatic_score >= $somatic_threshold) {
            # Positions are 0 based (input from bed), but the bam file will have 1 based positions. Correct this offset so that we look at the correct position when examining the bam.
            print $tfh join ("\t", ($chr, $start+1, $stop)) . "\n";
            push @sniper_lines, $line;
        }
        else {
            print $lq_output_fh $line."\n";
        }
    }
    $tfh->close;

    #If nothing passed the somatic score threshold, don't bother looking at mapping scores
    unless(@sniper_lines) {
        $lq_output_fh->close;
        $ofh->close;
        $self->sort_lq_output($lq_output_file, $sorted_lq_output_file);
        return 1;
    }
    #Run readcount program 
    my $readcount_command = sprintf("%s %s -l %s %s |",$self->readcount_path, $self->bam_readcount_params, $temp_path, $tumor_bam_file);
    $self->debug_message("Running: $readcount_command");

    my $readcounts = IO::File->new($readcount_command);

    while(my $count_line = $readcounts->getline) {
        chomp $count_line;
        my ($chr, $pos, $ref, $depth, @base_stats) = split /\t/, $count_line;
        
        my $current_variant = shift @sniper_lines;
        last unless $current_variant;
        my ($vchr, $vstart, $vstop, $ref_iub, $somatic_score, $tumor_read_depth) = split /\t/, $current_variant;
        # Positions are 0 based (input from bed), but the bam file will have 1 based positions. Correct this offset.
        $vstart++;
        my ($vref, $viub) = split "/", $ref_iub;

        #check if the sniper line was present in the readcount output
        while($vchr ne $chr && $vstart != $pos && @sniper_lines) {
            $self->debug_message("Skipped $current_variant");
            
            print $lq_output_fh $current_variant ."\n";;
            $current_variant = shift @sniper_lines;
            last unless $current_variant;
            ($vchr, $vstart, $vstop, $ref_iub, $somatic_score, $tumor_read_depth) = split /\t/, $current_variant;
            ($vref, $viub) = split "/", $ref_iub;
            # Positions are 0 based (input from bed), but the bam file will have 1 based positions. Correct this offset.
            $vstart++;
        }
        last unless $current_variant;
        
        my %bases;
        for my $base_stat (@base_stats) {
            my ($base,$reads,$avg_mq, $avg_bq) = split /:/, $base_stat;
            #Leaving bases coded as '=' unhandled
            next if($base eq "=");
            $bases{$base} = $avg_mq;
        }

        my @vars = Genome::Info::IUB->variant_alleles_for_iub($vref,$viub);
        my $variant_is_hq = 0;
        foreach my $var (@vars) {
            if(exists($bases{$var}) && $bases{$var} >= $self->min_mapping_quality) {
                print $ofh $current_variant, "\n";
                $variant_is_hq =1;
                last;
            }
        }
        # Unless the IUB had one possible allele that was HQ, make it LQ
        unless ($variant_is_hq) {
            print $lq_output_fh $current_variant ."\n";
        }
    }

    # If we ran out of readcount lines before we ran out of input variant lines, print the rest to LQ
    for my $line (@sniper_lines) { 
        print $lq_output_fh $line ."\n";
    }

    unless($readcounts->close()) {
        $self->error_message("Error running " . $self->readcount_path);
        die;
    }
    $lq_output_fh->close;
    $ofh->close;
    $fh->close;

    # Sort the LQ output file (This is necessary since we filter things out twice, once by somatic score, once by mapping quality...they will be out of order)
    $self->sort_lq_output($lq_output_file, $sorted_lq_output_file);

    return 1;
}

sub sort_lq_output {
    my $self = shift;
    my $lq_output_file = shift;
    my $sorted_lq_output_file = shift;
    my @sort_input = ($lq_output_file);
    unless ( Genome::Model::Tools::Joinx::Sort->execute(input_files => \@sort_input, output_file => $sorted_lq_output_file) ) {
        $self->error_message("Failed to sort the LQ output $lq_output_file into $sorted_lq_output_file using Joinx::Sort");
        die $self->error_message;
    }
    return 1;
}

sub readcount_path {
    my $self = $_[0];
    return $self->path_for_readcount_version($self->bam_readcount_version);
}

sub available_readcount_versions {
    my $self = shift;
    return keys %READCOUNT_VERSIONS;
}

sub path_for_readcount_version {
    my $class = shift;
    my $version = shift;

    if (defined $READCOUNT_VERSIONS{$version}) {
        return $READCOUNT_VERSIONS{$version};
    }
    die('No path for bam-readcount version '. $version);
}

sub default_readcount_version {
    die "default bam-readcount version: $DEFAULT_VERSION is not valid" unless $READCOUNT_VERSIONS{$DEFAULT_VERSION};
    return $DEFAULT_VERSION;
}

1;
