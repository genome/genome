package Genome::Model::Tools::Analysis::TelomereLength;


use warnings;
use strict;
use Genome;
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
    build => {
       is  => 'String',
       is_input=>1, 
       is_optional => 1,
       doc => 'specify the alignment reference build.  Acceptable value= "36" or "37" ',
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
    bam_list_file => {
        is  => 'String',
        is_input=>1, 
        is_optional => 1,
        doc => 'tab-delimited file of lable, NormalBAM_path, TumorBAM_path',
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
gmt analysis telomere-length --build=37 --group-id=12345 --edit-dist=10 
gmt analysis telomere-length --build=37 --somatic-id=456789,44556677 --edit-dist=2
gmt analysis telomere-length --bam-list-file SJMEL.bam.list --roi 12:133841453-133841911 --edit 2
EOS
}

sub help_detail {                           
    return <<EOS 
    This tool estimate the telomere length between normal and tumor samples by counting the number of telomeric reads in both samples.  There are 3 ways to specify which samples to run (choose only 1) --somatic-id,--group-id,-bam-list-file.  This tool needs a list of ROIs to search for telomeric reads.  There are 3 ways to specify ROIs.  1.) --build=36/37.  This will use a pre-generated ROI file which contain > 95% of telomeric reads.  This is the recommended way to run this script.  2.)--roi=12:1000-4000,21:50-400  You can also specify individual regions (comma separated).  3.)--roi-file=telomeric.roi.file  If you have your own ROI you want to search for telomeric reads, make a tab-delimited file containing all ROIs (chr start stop).  You can also set a similarity threshhold for what is a telomeric reads by --edit-dist=5.  Edit distance=5 means that in 100bp reads, you are willing to tolerate 5 mismatches.  If --verbose is NOT enabled, the output will be to STDOUT and consists of tab-delimited 5 columns: sample_label  edit_dist  Num_reads_Normal  Num_reads_Tumor  Tumor/Normal.  If --verbose is ON, it will print out much detailed information on number of each for each of the ROI.  For most common usage, you do not need --verbose on.

EOS
    }

sub execute {
    my $self = shift;
    $DB::single = 1;

    my $ref_seq_build = $self->build;

    ##define the region-of-interest by either command line or via a file##
    my $ROI_file = $self->roi_file; #specify ROIs in a file
    my $ROI = $self->roi;  #specify ROI on commandline
    #RegionOfInterest contain all the regions that we are going to examine for telomeric reads
    my @RegionOfInterest;
    if($ref_seq_build) {
	if($ref_seq_build ==36) {
	    @RegionOfInterest = getROI("/gscmnt/gc6111/info/medseq/telomere/build36.telomere.merged.bed");
	}elsif($ref_seq_build ==37) {
	    @RegionOfInterest = getROI("/gscmnt/gc6111/info/medseq/telomere/build37.telomere.merged.bed");
	}else {
	    print STDERR "Illegal build ref version.  Aborting...\n";
	    return 0;
	}
    }
    elsif($ROI_file) {
	@RegionOfInterest = getROI($ROI_file); 
    }elsif($ROI) {
	@RegionOfInterest = split(/,/,$ROI);
    }else { #default to region defined for build37
	print STDERR "Need to specify --build OR --roi-file OR --roi. Aborting...\n";
	return 0;
    }
    #####################################################################

    my $bamlist = $self->bam_list_file;
    my $somatic_id = $self->somatic_id;
    my $group_id = $self->group_id;
    my $verbose = $self->verbose;
    my $edit_dist = $self->edit_dist;
    my $intersect_file = $self->intersect_region_file;

    unless($bamlist) {
    #user specified model_id, group_id
	my @models=();
	if($group_id) {
	    my $model_group = Genome::ModelGroup->get($group_id);
	    @models = $model_group->models;
	}elsif($somatic_id) {
	    my @modelIDs = split(/,/,$somatic_id);
	    #@models = map{ Genome::Model->get($_) } @modelIDs;
	    @models = Genome::Model->get(\@modelIDs);
	}
	process_somatic_models(\@models,\@RegionOfInterest,$edit_dist,$intersect_file,$verbose);
    }
    else {
    #user specified bam list file
	my $valid_bamfiles = process_bam_list($bamlist);
	process_via_bam_list_file($valid_bamfiles,\@RegionOfInterest,$edit_dist,$intersect_file,$verbose);
    }



    return 1;                               
}

sub process_via_bam_list_file {

    my $bamfiles = shift;
    my $ROIs = shift;       #array ref of ROIs
    my $edit_dist = shift; 
    my $intersect_file =shift;
    my $verbose = shift;
    

    foreach my $label (keys %$bamfiles) {
	my $normalBAM = $bamfiles->{$label}->{'normal'};
	my $tumorBAM = $bamfiles->{$label}->{'tumor'};
	my $N_bam_flagstat = "${normalBAM}.flagstat";
	my $T_bam_flagstat = "${tumorBAM}.flagstat";

	#check to see if all the BAM files are there
	my $N_bam_stat = Genome::Model::Tools::Sam::Flagstat->parse_file_into_hashref($N_bam_flagstat);
	my $T_bam_stat = Genome::Model::Tools::Sam::Flagstat->parse_file_into_hashref($T_bam_flagstat);
	my $N_total_reads = $N_bam_stat->{'total_reads'};
	my $T_total_reads = $T_bam_stat->{'total_reads'};
	if($T_total_reads ==0) {print STDERR "Error with Tumor Total Reads: skipping.\n"; next;}
	my $n_t_ratio = $N_total_reads/$T_total_reads; #factor to shift the tumor read depth to
	
	my $normal_telomere = get_telomere_count2($normalBAM,$ROIs,$edit_dist,$intersect_file);
	my $tumor_telomere  = get_telomere_count2($tumorBAM,$ROIs,$edit_dist,$intersect_file);
	$tumor_telomere = normalize_read_depth($tumor_telomere,$n_t_ratio);

	print_result($label,$normal_telomere,$tumor_telomere,$edit_dist,$ROIs,$verbose);

    }
    my $x=1;
    


}


sub process_somatic_models {

    my $models = shift; #array ref of somatic models
    my $ROIs = shift;    #array ref of ROIs
    my $edit_dist = shift; 
    my $intersect_file =shift;
    my $verbose = shift;

    foreach my $somatic_model(@$models) {
	my $sample;
	if($somatic_model->subject->isa('Genome::Individual') ) {
	    $sample = $somatic_model->subject->common_name; 
	}else {
	    $sample = $somatic_model->subject->source_common_name;
	    #$sample = $somatic_model->subject->name;
	}
	$sample = $somatic_model->id if(!$sample); #label defaults to model ID if common name cannot be found.

	my $normal_build = $somatic_model->normal_model->last_succeeded_build;
	my $tumor_build = $somatic_model->tumor_model->last_succeeded_build;
	my $normal_telomere = get_telomere_readcount4build($normal_build,$ROIs,$edit_dist,$intersect_file);
	my $tumor_telomere = get_telomere_readcount4build($tumor_build,$ROIs,$edit_dist,$intersect_file);

	if(keys %$normal_telomere ==0 || keys %$tumor_telomere==0) {
	    print STDERR "Error encountered for $sample.  Skipping...\n";
	    next;
	}
	my $N_bam_stat = Genome::Model::Tools::Sam::Flagstat->parse_file_into_hashref($normal_build->whole_rmdup_bam_flagstat_file);
	my $T_bam_stat = Genome::Model::Tools::Sam::Flagstat->parse_file_into_hashref($tumor_build->whole_rmdup_bam_flagstat_file);
        #my $N_map_undup_reads = $N_bam_stat->{'reads_mapped'} - $N_bam_stat->{'reads_marked_duplicates'};
        my $N_total_reads = $N_bam_stat->{'total_reads'};
	my $T_total_reads = $T_bam_stat->{'total_reads'};
	my $n_t_ratio = $N_total_reads/$T_total_reads; #factor to shift the tumor read depth to
   
	$tumor_telomere = normalize_read_depth($tumor_telomere,$n_t_ratio);

	print_result($sample,$normal_telomere,$tumor_telomere,$edit_dist,$ROIs,$verbose);

    }

}

sub print_result {

    my $label = shift;
    my $normal_telomere = shift;
    my $tumor_telomere = shift;
    my $edit_dist = shift;
    my $ROIs = shift;
    my $verbose = shift;

    my $x=1;
    if(!$verbose) {
	print_summary($normal_telomere,$tumor_telomere,$label,$edit_dist);
    }
    else {
	print_detail($normal_telomere,$tumor_telomere,$label,$ROIs,$edit_dist);
    }

}

sub get_telomere_readcount4build {
#takes a normal/turmo ref-align build, return its telomeric read count

    my $build = shift;
    my $ROI = shift;
    my $edit_dist = shift;
    my $region_file = shift;

    my $BAM = $build->whole_rmdup_bam_file;
    my $read_count = get_telomere_count2($BAM,$ROI,$edit_dist,$region_file);
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

    #my $filter_status;
    #$filter_status = ($phred_filter) ? "phred_filter" : "no_phred_filter";
    #my $read_status;
    #$read_status = ($all_reads) ? "all_reads" : "unique_reads";

    my $normal_summary = generate_total_counts($normal); 
    my $tumor_summary = generate_total_counts($tumor);

    my $telomere_reads_normal = $normal_summary->{'TELOMERE_READS'};
    my $telomere_reads_tumor  = $tumor_summary->{'TELOMERE_READS'};
    my $total_ratio = ($telomere_reads_normal == 0) ? 'NA' : nearest(0.001,$telomere_reads_tumor/$telomere_reads_normal);
    #print "$sample\t$edit_dist\t$telomere_reads_normal\t$telomere_reads_tumor\t$total_ratio\t$filter_status\t$read_status\n";
    print "$sample\tED_$edit_dist\t$telomere_reads_normal\t$telomere_reads_tumor\t$total_ratio\n";



}

sub print_detail {
#for debugging purpose...usually not needed
    my $normal = shift;
    my $tumor = shift;
    my $sample = shift;
    my $ROI = shift;
    my $edit_dist = shift;

    foreach my $roi (keys %$normal) {
	my $N_total = $normal->{$roi}->{'Total_reads'};
	my $N_map = $normal->{$roi}->{'mapped'};
	my $N_mapROI = $normal->{$roi}->{'mapped_intersect'};
	my $N_mapPercent =  ($N_map == 0) ? 'NA' : nearest(0.001,$N_mapROI/$N_map);
	my $N_unmap = $normal->{$roi}->{'unmapped'};
	my $N_unmapROI = $normal->{$roi}->{'unmapped_intersect'};
	my $N_unmapPercent = ($N_unmap == 0) ? 'NA' : nearest(0.001,$N_unmapROI/$N_unmap);

	my $T_total = $tumor->{$roi}->{'Total_reads'};
	my $T_map = $tumor->{$roi}->{'mapped'};
	my $T_mapROI = $tumor->{$roi}->{'mapped_intersect'};
	my $T_mapPercent =  ($T_map == 0) ? 'NA' : nearest(0.001,$T_mapROI/$T_map);
	my $T_unmap = $tumor->{$roi}->{'unmapped'};
	my $T_unmapROI = $tumor->{$roi}->{'unmapped_intersect'};
	my $T_unmapPercent = ($T_unmap == 0) ? 'NA' : nearest(0.001,$T_unmapROI/$T_unmap);
	print "$sample\tED_${edit_dist}\t$roi\t$N_total\t$N_map\t$N_mapROI\t$N_mapPercent\t$N_unmap\t$N_unmapROI\t$N_unmapPercent\tnormal\n";
        print "$sample\tED_${edit_dist}\t$roi\t$T_total\t$T_map\t$T_mapROI\t$T_mapPercent\t$T_unmap\t$T_unmapROI\t$T_unmapPercent\ttumor\n";

    }
   
}

sub get_telomere_count2 {
#this version does not require all telomeric repeats to be perfectly consecutive
#use approximate match

    my $bamfile = shift;  #path to the bam file
    my $ROI = shift; #ref to list of ROIs
    my $max_edit = shift; #number of "edits" allowed and still qualify for a match
    my $region_file = shift; 

    #telomere pattern on both strands up to 100bp
    my @patternA = generate_all_6_frames('CCCTAA');
    my @patternB = generate_all_6_frames('TTAGGG');

    #####
    my $patternA = 'CCCTAA';
    my $patternB = 'TTAGGG'; 
    #####

    my $intersect_region = get_intersect_region_file($region_file);
    
    my $failedmatch=0;
    my $match_unmapped_reads=0;
    my $match_mapped_reads=0;
    my $telomere_count={};
    my $weird=0;

    return $telomere_count if(!-e $bamfile);

    my $all_reads=0;
    foreach my $roi (@$ROI) {
	my $cmd;
	#print STDERR "processing $roi\n";
	if($all_reads) {
	    $cmd = "samtools view $bamfile $roi"; #extract all reads in a region
	}else {
	    $cmd = "samtools view -F 0x400 $bamfile $roi"; #extract all unique reads in a region
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

            ##############################################
	    #approximate match in 6 reading frames for each read is very slow and time consuming
	    #do a perfect match first, if a read does not have at least 1 telomere repeats (6bp)
            #no futher processing will be done;  This should reduce processing of a lot of reads.

	    my $telomere_peatsA = () = $readSeq =~ /$patternA/gio;
	    my $telomere_peatsB = () = $readSeq =~ /$patternB/gio;
	    my $telomere_repeats = ($telomere_peatsA > $telomere_peatsB) ? $telomere_peatsA : $telomere_peatsB;
	    if($telomere_repeats < 1) { #skip approximate match for any reads with Telomere peats less than specified
		$telomere_count->{$roi}->{'skipped'}++;
		$telomere_count->{$roi}->{'failed_similarity'}++;
		next;
	    }


            #############################################
	    my @matchedA = amatch($readSeq,["i $max_edit"],@patternA);
	    my @matchedB = amatch($readSeq,["i $max_edit"],@patternB);

	    #my $matchedA = amatch($patternA, ["i $max_edit"], $readSeq);
	    #my $matchedB = amatch($patternB, ["i $max_edit"], $readSeq);

	    if(!@matchedA && !@matchedB) {  #skip if a read fail to approximate match on both strands
		$failedmatch++;
		$telomere_count->{$roi}->{'failed_similarity'}++;
		next;
	    }
	    if(@matchedA && @matchedB) {
		$weird++;
		#this section of code deals with a bug in amatch.  Sometimes string with greater than 65 edit distance
                #will still result in a match.  This code will check the actual edit distance.
                #if edit distance from both strand is > the cutoff, the read is discarded.
		#my $min_distA = min(adist($readSeq,@patternA));
		#my $min_distB = min(adist($readSeq,@patternB));
		#my $x=1;
		#if($min_distA > $max_edit && $min_distB > $max_edit) {
		#    print "Bad read discarded: $readSeq\n";
		#    next;
		#}
		#print STDERR "Weird case.  Cannot match to both strand but it did\n";
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
		#print TEMP "$chr\t$readStart\t",$readStart+100,"\n" if(defined($location_file));
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
	$telomere_count->{$roi}->{'unmapped'} = 0 if(!exists($telomere_count->{$roi}->{'unmapped'}));
	$telomere_count->{$roi}->{'mapped'} = 0 if(!exists($telomere_count->{$roi}->{'mapped'}));
	$telomere_count->{$roi}->{'failed_similarity'} = 0 if(!exists($telomere_count->{$roi}->{'failed_similarity'}));
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


sub process_bam_list {

    my $file = shift;
    my $check_only = shift;

    open(FILE, $file) or die "Can't open the file $file due to $!";
    my $BAM_files={};
    while(<FILE>) {
	chomp;
	my($label,$normal_bam,$tumor_bam) = split(/\t/,$_);
	my $N_bam_flagstat = "${normal_bam}.flagstat";
	my $T_bam_flagstat = "${tumor_bam}.flagstat";

	if(all_files_exist($label,$tumor_bam,$normal_bam,$T_bam_flagstat,$N_bam_flagstat)) {    
	    $BAM_files->{$label} = { 'normal'=> $normal_bam,
				     'tumor' => $tumor_bam,
	                           };
	}
    }
    close FILE;    
    return $BAM_files;

}

sub all_files_exist {

    my $label = shift;
    my $tumorBAM = shift;
    my $normalBAM = shift;
    my $tumorBAM_flagstat = shift;
    my $normalBAM_flagstat = shift;

    my $status=1;
    if(!-e $tumorBAM) {
	print STDERR  "$label Tumor BAM cannot be found\n";
	$status=0;
    }
    if(!-e $normalBAM) {
	print STDERR  "$label Normal BAM cannot be found\n";
	$status=0;
    }
    if(!-e $tumorBAM_flagstat) {
	print STDERR  "$label Tumor BAM flagstat cannot be found\n";
	$status=0;
    }
    if(!-e $normalBAM_flagstat) {
	print STDERR  "$label Normal BAM flagstat cannot be found\n";
	$status=0;
    }

    print STDERR "All required files (BAM, flagstat) are found for $label\n" if($status);

    return $status;
}


sub rem_white_space {

    my $string = shift;

    $string =~ s/^\s+//;
    $string =~ s/\s+$//;

    return $string;
}
