package Genome::Model::Tools::Analysis::TelomereLength;


use warnings;
use strict;
use Genome;
use Workflow;
use FileHandle;
use IO::File;
use Math::Round;
use String::Approx 'amatch','adist';
use Genome::Info::IUB;
use Cwd ('getcwd','abs_path');

class Genome::Model::Tools::Analysis::TelomereLength {
    is => ['Command'],
    has => [
    roi_file => { 
        is  => 'String',
        is_input=>1, 
	is_optional => 1,
        doc => 'name of the file containing regions to examine telomeric reads; tab-delimited (chr start stop) 1/line',
    },
    roi => { 
        is  => 'String',
        is_input=>1, 
	is_optional => 1,
        doc => 'chr:start-stop (comma separated for multiple regions) ',
    },
    somatic_id => {
	is  => 'String',
        is_input=>1, 
	is_optional => 1,
        doc => 'somatic variation id. If >1 id, comma-separate',
    },
    group_id => {
	is  => 'String',
        is_input=>1, 
	is_optional => 1,
        doc => 'model group id containing many somatic variation id.  Will supercede --somatic-id ',
    },    
    edit_dist => { 
        is  => 'String',
        is_input=>1, 
        is_optional => 1,
        default_value => '15',
        doc => 'The maximum Levenshtein distance allowed for a 100bp read to be classified as telomeric',
    },
    verbose => { 
        is  => 'Boolean',
        is_input=>1, 
        is_optional => 1,
        #default_value => 'output.pdf',
        doc => 'Returns more detailed output of telomeric read count in reach ROI and whether it is mapped or unmapped',
	default => 0,
    },
    intersect_region_file => {
	is_input => 1,
	is_optional => 1,
	doc => 'tab delimited region file similar to --roi-file.  Used only in debugging to see how many reads from your ROI intersect with user defined ROIs',
    },

    ],
};

sub help_brief {
    "Given an annotation file, gives an output of transition/transversion, cpg islands, and cpgs within cpg islands.";
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt analysis telomere-length  --roi-file=SJMEL001003-0260.alltier.snv --output-file=SJMEL.pdf --plot-title="SJMEL001003-0260 --window-size=10"

EOS
}

sub help_detail {                           
    return <<EOS 
    This tool summarizes the mutation spectrum sequence context.  It produces a stacked barplot for each mutation cateogry showing the proportion of bases around the point of interest (position 0) +/- window_size basepairs.
EOS
    }

sub execute {
    my $self = shift;
    $DB::single = 1;

    ##define the region-of-interest by either command line or via a file##
    my $ROI_file = $self->roi_file; #specify ROIs in a file
    my $ROI = $self->roi;  #specify ROI on commandline
    #RegionOfInterest contain all the regions that we are going to examine for telomeric reads
    my @RegionOfInterest;
    if($ROI_file) {
	@RegionOfInterest = getROI($ROI_file); 
    }elsif($ROI) {
	@RegionOfInterest = split(/,/,$ROI);
    }else { #default to region defined for build37
	@RegionOfInterest = qw(1:10001-10467 12:95159-95735 12:133841522-133841895 2:243152477-243152612 5:10001-11813 X:155260017-155260464);
	my $str = join(";",@RegionOfInterest);
	print STDERR "User not define ROI via --roi or --roi-file.  Default to ROI for build37! ";
	print STDERR "$str\n";
    }
    #####################################################################

    my $somatic_id = $self->somatic_id;
    my $group_id = $self->group_id;
    my $verbose = $self->verbose;
    my $edit_dist = $self->edit_dist;
    my $intersect_file = $self->intersect_region_file;

    my @models=();
    if($group_id) {
	my $model_group = Genome::ModelGroup->get($group_id);
	@models = $model_group->models;
    }elsif($somatic_id) {
	my @modelIDs = split(/,/,$somatic_id);
	#@models = map{ Genome::Model->get($_) } @modelIDs;
	@models = Genome::Model->get(\@modelIDs);
    }
    

    foreach my $somatic_model(@models) {
	
	my $sample = $somatic_model->subject->source_common_name;
	#my $sample_label = $somatic_model->last_succeeded_build->tumor_build->model->subject->source_common_name;
	$sample = $somatic_model->id if(!$sample); #label defaults to model ID if common name cannot be found.
	
	my $normal_build = $somatic_model->normal_model->last_succeeded_build;
	my $tumor_build = $somatic_model->tumor_model->last_succeeded_build;
	my $normal_telomere = get_telomere_readcount4build($normal_build,\@RegionOfInterest,$edit_dist,$intersect_file);
	my $tumor_telomere = get_telomere_readcount4build($tumor_build,\@RegionOfInterest,$edit_dist,$intersect_file);


	my $N_bam_stat = Genome::Model::Tools::Sam::Flagstat->parse_file_into_hashref($normal_build->whole_rmdup_bam_flagstat_file);
	my $T_bam_stat = Genome::Model::Tools::Sam::Flagstat->parse_file_into_hashref($tumor_build->whole_rmdup_bam_flagstat_file);
        #my $N_map_undup_reads = $N_bam_stat->{'reads_mapped'} - $N_bam_stat->{'reads_marked_duplicates'};
        my $N_total_reads = $N_bam_stat->{'total_reads'};
	my $T_total_reads = $T_bam_stat->{'total_reads'};
	my $n_t_ratio = $N_total_reads/$T_total_reads; #factor to shift the tumor read depth to
   
	$tumor_telomere = normalize_read_depth($tumor_telomere,$n_t_ratio);

	my $x=1;
	if(!$verbose) {
	    print_summary($normal_telomere,$tumor_telomere,$sample,$edit_dist);
	}
	else {
	    print_detail($normal_telomere,$tumor_telomere,$sample,\@RegionOfInterest,$edit_dist);
	}
    }
    
    #process_somatic_variation_models(\@models,\@RegionOfInterest,$edit_dist);

    return 1;                               
}



sub get_telomere_readcount4build {
#takes a normal/turmo ref-align build, return its telomeric read count

    my $build = shift;
    my $ROI = shift;
    my $edit_dist = shift;
    my $region_file = shift;

    my $BAM = $build->whole_rmdup_bam_file;
    my $read_count = get_telomere_count($BAM,$ROI,$edit_dist,$region_file);
    #normalize the raw read count
    #$read_count = normalize_read_depth($read_count,$bam_total_reads);

    return($read_count);
}

sub normalize_read_depth {
#normalize the raw telomere readcount by the ratio of tumor and normal total reads
    my $read_count = shift;
    my $n_t_ratio = shift;

    foreach my $roi (keys %$read_count) {
	foreach my $cat( keys %{ $read_count->{$roi} }) {
	    $read_count->{$roi}->{$cat} =round($read_count->{$roi}->{$cat} * $n_t_ratio);
	}
    }

    return $read_count;
}


sub getROI {

    my $file = shift;

    my @ROIs=();
    open(FILE, $file) or die "Unable to open the file $file due to $!";
    while(<FILE>) {
	next if(/\#/);#skip comment lines
	chomp;
	my($chr,$start,$stop) = split(/\s+/,$_);
	next if(!$chr && !$start && !$stop);
	my $roi = "$chr:${start}-${stop}";
	push(@ROIs,$roi);
    }
    close FILE;

    return @ROIs;
}

sub print_summary {

    my $normal = shift;
    my $tumor = shift;
    my $sample = shift;
    my $edit_dist = shift;

    my $filter_status;
    $filter_status = (0) ? "phred_filter" : "no_phred_filter";
    my $read_status;
    $read_status = (0) ? "all_reads" : "unique_reads";

    my $normal_summary = generate_total_counts($normal); 
    my $tumor_summary = generate_total_counts($tumor);

    my $telomere_reads_normal = $normal_summary->{'TELOMERE_READS'};
    my $telomere_reads_tumor  = $tumor_summary->{'TELOMERE_READS'};
    my $total_ratio = $telomere_reads_tumor/$telomere_reads_normal;
    print "$sample\t$edit_dist\t$telomere_reads_normal\t$telomere_reads_tumor\t$total_ratio\t$filter_status\t$read_status\n";



}

sub print_detail {
#for debugging purpose...usually not needed
    my $normal = shift;
    my $tumor = shift;
    my $sample = shift;
    my $ROI = shift;
    my $edit_dist = shift;

    foreach my $roi (keys %$normal) {
	my $N_map = $normal->{$roi}->{'mapped'} || 0;
	my $N_mapROI = $normal->{$roi}->{'mapped_intersect'} || 0;
	my $N_mapPercent =  ($N_map == 0) ? 'NA' : nearest(0.001,$N_mapROI/$N_map);
	my $N_unmap = $normal->{$roi}->{'unmapped'} || 0;
	my $N_unmapROI = $normal->{$roi}->{'unmapped_intersect'} || 0;
	my $N_unmapPercent = ($N_unmap == 0) ? 'NA' : nearest(0.001,$N_unmapROI/$N_unmap);

	my $T_map = $tumor->{$roi}->{'mapped'} || 0;
	my $T_mapROI = $tumor->{$roi}->{'mapped_intersect'} || 0;
	my $T_mapPercent =  ($T_map == 0) ? 'NA' : nearest(0.001,$T_mapROI/$T_map);
	my $T_unmap = $tumor->{$roi}->{'unmapped'} || 0;
	my $T_unmapROI = $tumor->{$roi}->{'unmapped_intersect'} || 0;
	my $T_unmapPercent = ($T_unmap == 0) ? 'NA' : nearest(0.001,$T_unmapROI/$T_unmap);
	print "$sample\t$edit_dist\t$roi\t$N_map\t$N_mapROI\t$N_mapPercent\t$N_unmap\t$N_unmapROI\t$N_unmapPercent\t$T_map\t$T_mapROI\t$T_mapPercent\t$T_unmap\t$T_unmapROI\t$T_unmapPercent\n";

    }
   
}

sub get_telomere_count {
#this version does not require all telomeric repeats to be perfectly consecutive
#use approximate match based on Levenshtein distance

    my $bamfile = shift;  #path to the bam file
    my $ROI = shift; #ref to list of ROIs
    my $max_edit = shift; #number of "edits" allowed and still qualify for a match
    my $region_file = shift; 

    #telomere pattern on both strands up to 100bp in all 6 frameshifts
    my @patternA = generate_all_6_frames('CCCTAA');
    my @patternB = generate_all_6_frames('TTAGGG');

    #my $patternA = 'CCCTAA'x16 . 'CCCT';
    #my $patternB = 'TTAGGG'x16 . 'TTAG'; 

    my $intersect_region = get_intersect_region_file($region_file);
    #my $intersect_region = get_region_file($region_file);

    my $failedmatch=0;
    my $match_unmapped_reads=0;
    my $match_mapped_reads=0;
    my $telomere_count={};
    #my $total_reads=0;

    my $all_reads=0;
    foreach my $roi (@$ROI) {
	my $cmd;
	if($all_reads) {
	    $cmd = "samtools view $bamfile $roi"; #extract all reads in a region
	}else {
	    $cmd = "samtools view -F 0x400 $bamfile $roi"; #extract all unique reads (defined by Picard) in a region
	}
	my $total_reads=0;
	open(SAMTOOL_VIEW, "$cmd|") or die "Can't open the file handle due to $!";
	while(<SAMTOOL_VIEW>) {
	    $total_reads++;
	    chomp;
	    my @list = split(/\t/,$_);
	    my $readStart = $list[3];
	    my $readSeq = $list[9];
	    my $readQUAL = $list[10];
	    my $CIGAR = $list[5];
	    my $chr = $list[2];
	    my $bitscore = $list[1];
	    
	    #replace any base with N for bases with very low phred score
	    #if($phred_filter) {
		#$readSeq = process_seq_read_by_base_quality($readSeq,$readQUAL,20);
	    #}

	    my @matchedA = amatch($readSeq,["i $max_edit"],@patternA);
	    my @matchedB = amatch($readSeq,["i $max_edit"],@patternB);

	    #my $matchedA = amatch($patternA, ["i $max_edit"], $readSeq);
	    #my $matchedB = amatch($patternB, ["i $max_edit"], $readSeq);

	    if(!@matchedA && !@matchedB) {  #skip if a read fail to approximate match on both strands
		$failedmatch++;
		#my $distA = adist($patternA, $readSeq);
		#my $distB = adist($patternB, $readSeq);
		$telomere_count->{$roi}->{'failed_similarity'}++;
		next;
	    }
	    if(@matchedA && @matchedB) {
		print STDERR "Weird case.  Cannot match to both strand but it did\n";
	    }

	    if(4 & $bitscore) {
		#if($CIGAR eq '*') {
		$match_unmapped_reads++;
		$telomere_count->{$roi}->{'unmapped'}++;
		if(overlaps_with_ROI($chr,$readStart,$intersect_region) ) {
		    $telomere_count->{$roi}->{'unmapped_intersect'}++; 
		}
	    }else {
		$match_mapped_reads++;	
		$telomere_count->{$roi}->{'mapped'}++;
		if(overlaps_with_ROI($chr,$readStart,$intersect_region) ) {
		    $telomere_count->{$roi}->{'mapped_intersect'}++;
		}
	    }
	}
	close SAMTOOL_VIEW;
	$telomere_count->{$roi}->{'Total_reads'}=$total_reads;
	$telomere_count->{$roi}->{'unmapped_intersect'} = 0 if(!exists($telomere_count->{$roi}->{'unmapped_intersect'}));
	$telomere_count->{$roi}->{'mapped_intersect'} = 0 if(!exists($telomere_count->{$roi}->{'mapped_intersect'}));

    }
    return  $telomere_count;

}

sub generate_all_6_frames {
#takes a Telomere 6peats pattern and generate all possible readingframe 
#sequence reads 100bp long.

    my $pattern = shift;  #base pattern

    my $ref_seq = $pattern x 20;
    my @all_frames;
    for(0 .. 5) {
	my $seq = substr($ref_seq,$_,100);
	push(@all_frames,$seq);
    }

    return @all_frames;


}

sub get_intersect_region_file {

    my $file = shift;

    my $region={};
    if(!$file || !-e $file) { #if not defined/exists, return empty hash ref
	return $region;
    }
    open(FILE,$file) or die "Can't open the file $file due to $!";
    while(<FILE>) {
	next if(/\#/);
	chomp;
	my($chr,$start,$stop) = split(/\s+/,$_);
	$region->{$chr}->{"${start}_${stop}"}=1;
    }
    close FILE;
    
    return $region;

}

sub overlaps_with_ROI {
#check to see if a bam read overlaps with any regions in the ROI file.
#this is used to find out what % of MAPPED reads that have telomere repeats that also
#overlaps with our "defined" telomere-region

    my $chr = shift;
    my $bam_read_start = shift;
    my $intersect_regions = shift;
    my $bam_read_end = $bam_read_start+99;

    my $overlap_status=0;
    if(!exists($intersect_regions->{$chr})) {
	$overlap_status=0;
    }else {
	foreach my $roi (keys %{$intersect_regions->{$chr}}) {
	    my ($start,$stop) = split(/\_/,$roi);
	    if(($bam_read_start >=$start and $bam_read_start <= $stop) or ($bam_read_end >= $start and $bam_read_end <= $stop)  ) {
		$overlap_status=1;
		last;
	    }
	
	}
    }

    return $overlap_status;

}

sub generate_total_counts {
#sum up read counts for all ROI over all categories

    my $read_count = shift;

    my $summary={'TOTAL_READS' => 0,   
                 'FAILED2MATCH'=> 0,
                 'MAP' => 0,
		 'MAP_ROI' => 0,
		 'UNMAP' => 0,
		 'UNMAP_ROI' => 0
                };
    foreach my $roi(keys %$read_count) {
	$summary->{'TOTAL_READS'} += $read_count->{$roi}->{'Total_reads'};
	$summary->{'FAILED2MATCH'} += $read_count->{$roi}->{'failed_similarity'};
	$summary->{'MAP'} += $read_count->{$roi}->{'mapped'};
	$summary->{'MAP_ROI'} += $read_count->{$roi}->{'mapped_intersect'};
	$summary->{'UNMAP'} += $read_count->{$roi}->{'unmapped'};
	$summary->{'UNMAP_ROI'} += $read_count->{$roi}->{'unmapped_intersect'};
    }

    $summary->{'TELOMERE_READS'}=$summary->{'MAP'}+$summary->{'UNMAP'};

    return $summary;

}


sub rem_white_space {

    my $string = shift;

    $string =~ s/^\s+//;
    $string =~ s/\s+$//;

    return $string;
}
