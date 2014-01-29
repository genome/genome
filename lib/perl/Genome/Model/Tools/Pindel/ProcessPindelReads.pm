package Genome::Model::Tools::Pindel::ProcessPindelReads;

use warnings;
use strict;

use Genome;
use Workflow;
use Genome::Statistics;

class Genome::Model::Tools::Pindel::ProcessPindelReads {
    is => [ 'Genome::Command::Base' ],
    has => [
        input_file => {
            is => 'String',
            doc => 'Fully qualified path to raw pindel output to be processed.',
        },
        output_file => {
            is => 'String',
            doc => 'Fully qualified path to desired output of processed pindel reads',
        },
        reference_build_id => {
            is => 'Text',
            doc => 'The build-id of a reference sequence build',
            is_optional => 0,
            is_input => 1,
        },
        reference_sequence_input => {
            calculate_from => ['reference_build_id'],
            calculate => q{ Genome::Model::Build->get($reference_build_id)->cached_full_consensus_path('fa') },
            doc => 'Location of the reference sequence file',
        },
        mode => {
            is => 'String',
            default => 'to_bed',
            valid_values => [
                'to_bed',
                'somatic_filter',
                'read_support',
                'vaf_filter',
            ],
            doc => 'What to do with raw pindel reads',
        },
        sort_output => {
            is => 'Boolean',
            doc => 'This flag is set by default, unset to prevent sorting of the output bed file.',
            default => 0,
        },
        variant_freq_cutoff => {
            is => 'Text',
            doc => " This is the minimum variant freq for read-support",
            default => "0.0",
        },
        capture_data => {
            is => 'Boolean',
            doc => "Set this flag to dump samtools view output to disk, which is slower, but avoids 'out of memory' errors.",
            default => 0,
        },
    ],
    has_optional => [
        create_hq_raw_reads => {
            is => 'Boolean',
            doc => 'Set this if you want to generate an HQ raw reads file',
            default => 0,
        },
        aligned_reads_input => {
            is => 'String',
            doc => 'Full path to aligned reads bam file',
        },
        control_aligned_reads_input => {    
            is => 'String',
            doc => 'Full path to control aligned reads bam file',
        },
        lq_output_file => {
            is => 'String',
            doc => 'Full path to lq_output_file',
        },
        hq_raw_output_file => {
            is => 'String',
            doc => 'Full path to new indels.hq file',
        },
        lq_raw_output_file => {
            is => 'String',
            doc => 'Full path to new indels.lq file',
        },
        read_support_output_file => {
            is => 'String',
            doc => 'Full path to the read support output',
        },
        big_output_file => {
            is => 'String',
            doc => 'Indels over size 100 go in here',
        },
    ],
    has_transient_optional => [
        _input_fh => {
            is => 'IO::File',
            doc => 'This is a file handle to the input file',
        },
        _output_fh => {
            is => 'IO::File',
            doc => 'This is a file handle to the output file',
        },
        _lq_output_fh => {
            is => 'IO::File',
            doc => 'This is a file handle to the lq output file',
        },
        _hq_raw_output_fh => {
            is => 'IO::File',
            doc => 'This is a file handle to the hq raw output file',
        },
        _lq_raw_output_fh => {
            is => 'IO::File',
            doc => 'This is a file handle to the lq raw output file',
        },
        _support_field_index => {
            is => 'String',
            doc => 'store the index of the Supports field in pindel input',
        },
        _read_support_fh => {
            is => 'IO::File',
            doc => 'This is a file handle to the read support output file',
        },
        _big_output_fh => {
            is => 'IO::File',
            doc => 'Indels over size 100 go here',
        },
        _refseq => {
            is => 'Text',
            doc => 'This is used to store the cached location of the reference, to avoid containtly using the accessor on the refseq build.',
        },
    ],
};

sub help_synopsis {
    return <<"EOS"
gmt pindel process-pindel-reads --input-file pindel.outfile --output-file pindel.adapted 
EOS
}

sub help_detail {                           
    return <<EOS 
Transforms a file from pindel output format to bed format
EOS
}

sub execute {
    my $self = shift;
    $self->_refseq($self->reference_sequence_input);
    $self->debug_message("Dumping reads from samtools view to temp files due to excessive read depth.") if $self->capture_data;
    print "Dumping reads from samtools view to temp files due to excessive read depth.\n" if  $self->capture_data;
    my $big_output_file = $self->big_output_file;
    my $hq_raw_output_file = $self->hq_raw_output_file;
    my $lq_raw_output_file = $self->lq_raw_output_file;
    my $read_support_output_file = $self->read_support_output_file;

    #open the big_output file-handle unless running "somatic_filter"
    if($self->mode ne 'somatic_filter'){
        $self->_big_output_fh(Genome::Sys->open_file_for_writing($big_output_file));
    }

    #open the raw reads file handles
    if($self->create_hq_raw_reads){
        $self->_hq_raw_output_fh(Genome::Sys->open_file_for_writing($hq_raw_output_file));
        $self->_lq_raw_output_fh(Genome::Sys->open_file_for_writing($lq_raw_output_file));
    }

    #open read-support file handle if needed
    if(defined($read_support_output_file)){
        $self->_read_support_fh(Genome::Sys->open_file_for_writing($read_support_output_file));
    }
    
    #set output to a temp file if sorting output, otherwise dump it right into the final destination
    my $sort_output_flag = $self->sort_output;
    my $bed_mode =  $self->mode eq 'to_bed';
    my $sort_output = ($sort_output_flag || $bed_mode);
    my $output = $sort_output ? $self->output_file.".temp" : $self->output_file;

    #process the raw pindel reads, calling $self->$mode once each read has been read into memory
    unless($self->process_source($self->input_file,$output,$self->_refseq)){
        die $self->error_message("Failed to get a return value from process_source.");
    }

    #close the main file handles
    $self->_input_fh->close;
    $self->_output_fh->close;

    print $output . "\n";

    #sort bed output if sort_output flag is set
    if($sort_output){
        my $real_output = $self->output_file;
        my @output_list = ($output);
        my $sort_cmd = Genome::Model::Tools::Joinx::Sort->create(input_files => \@output_list, output_file => $real_output);
        unless($sort_cmd->execute){
            die $self->error_message("Could not complete sort of output file.");
        };
        unlink($output);
    }
    #shut everything down as needed
    if($self->mode ne 'somatic_filter'){
        $self->_big_output_fh->close;
    }
    if($self->create_hq_raw_reads){
        $self->_hq_raw_output_fh->close;
        $self->_lq_raw_output_fh->close;
    }
    if(defined($read_support_output_file)){
        $self->_read_support_fh->close;
    }
    return 1;
}

sub process_source { 
    my $self = shift;
    my ($input, $output, $refseq) = @_;
    my %events;
    my $mode = $self->mode;
    my $input_fh = Genome::Sys->open_file_for_reading($input);
    my $output_fh = Genome::Sys->open_file_for_writing($output);
    $self->_input_fh($input_fh);
    $self->_output_fh($output_fh);

    my ($chrom,$pos,$size,$type);
    while(my $line = $input_fh->getline){
        my $normal_support=0;
        my $read = 0;
        if($line =~ m/^#+$/){
            my @event = (); 
            my $call = $input_fh->getline;
            chomp $call;
            push @event, $call;
            my $reference = $input_fh->getline;
            chomp $reference;
            push @event, $reference;
            my @call_fields = split /\s/, $call;
            my $support = $call_fields[15];

            unless(defined($support)){
                print "No support. Call was:   ".$call."\n";
                die;
            }

            for (1..$support){
                $line = $input_fh->getline;
                chomp $line;
                push @event, $line;
                $read=$line;
            }

            $self->$mode(\@event);
        }
    }
    return 1;
}

sub to_bed {
    my $self = shift;
    my $event = shift;
    my @event = @{$event};
    my $output_fh = $self->_output_fh;
    my $bed = $self->get_bed_line(\@event);
    unless( defined($bed) ){
        return undef;
    }
    print $output_fh $bed."\n";
    return $bed;
}

sub somatic_filter {
    my $self = shift;
    my $event = shift;
    my $hq_raw = $self->_hq_raw_output_fh;
    my $lq_raw = $self->_lq_raw_output_fh;
    my @event = @{ $event };
    my ($call, $reference, @support) = @event;
    my $normal_support = 0;
    for my $support (@support){
        if($support =~ m/normal/){
            $normal_support++;
        } 
    }
    if($self->create_hq_raw_reads && ($normal_support==0)){
        $self->print_raw_read(\@event,$hq_raw);
    }
    elsif ($normal_support && $self->create_hq_raw_reads){
        $self->print_raw_read(\@event,$lq_raw);
    }
    else {
        die $self->error_message("Nothing to do!");
    }
    return 1;
}

sub read_support {
    my $self = shift;
    my $event = shift;
    my $hq_raw = $self->_hq_raw_output_fh;
    my $lq_raw = $self->_lq_raw_output_fh;
    my @event = @{ $event };
    my ($call, $reference, @support) = @event;

    my $tumor_bam = $self->aligned_reads_input;
    my $normal_bam = $self->control_aligned_reads_input;
    unless(-s $tumor_bam){
        die $self->error_message("Could not locate tumor_bam at: ".$tumor_bam);
    }
    unless(-s $normal_bam){
        die $self->error_message("Could not locate normal_bam at: ".$tumor_bam);
    }

    my @call_fields = split /\s/, $call;

    my $type = $call_fields[1];
    my $size = $call_fields[2];
    my $chr = $call_fields[7];
    my $start = $call_fields[9];
    my $stop = $call_fields[10];
    my $lower_range = $call_fields[12];
    my $upper_range = $call_fields[13];
    my $pos_support = $call_fields[18];
    my $neg_support = $call_fields[21];

    my $pos_percent=0;
    if($neg_support==0){
        $pos_percent = 1.0;
    } else {
        $pos_percent = sprintf("%.2f", $pos_support / ($pos_support + $neg_support));
    }

    my $reads = scalar(@support);

    #my ($chr,$start,$stop)

    # Call samtools over the variant start-stop to get overlapping reads
    #my @results = `samtools view $tumor_bam $chr:$stop-$stop`;

    my $tumor_read_support=0;
    my $tumor_read_sw_support=0;
    my $normal_read_support=0;
    my $normal_read_sw_support=0;

    my $temp = Genome::Sys->create_temp_file_path;

    my $cap = $self->capture_data;

    if($cap){
        my $tsam_cmd = "samtools view $tumor_bam $chr:$stop-$stop > $temp";
        #Genome::Sys->shellcmd( cmd => $tsam_cmd);
        if(system($tsam_cmd)){
            die $self->error_message("Failed to run the command: $tsam_cmd");
        }

        my $tfh = Genome::Sys->open_file_for_reading($temp);

        while(my $result = $tfh->getline){
            chomp $result;
            my @details = split /\t/, $result;
            # Parse overlapping reads for cigar strings containing insertions or deletions
            if($details[5] =~ m/[ID]/){
                $tumor_read_sw_support++;
            }
            else {
                $tumor_read_support++;
            }
        }
        $tfh->close;

        # Call samtools over the variant start-stop in the normal bam to get overlapping reads
        my $nsam_cmd = "samtools view $normal_bam $chr:$stop-$stop > $temp";
        #Genome::Sys->shellcmd( cmd => $nsam_cmd);
        if(system($nsam_cmd)){
            die $self->error_message("Failed to run the command: $nsam_cmd");
        }
        my $nfh = Genome::Sys->open_file_for_reading($temp);

        while(my $result = $nfh->getline){
            chomp $result;
            my @details = split /\t/, $result;
            # Parse overlapping reads for insertions or deletions
            if($details[5] =~ m/[ID]/){
                $normal_read_sw_support++;
            }
            else {
                $normal_read_support++;
            }

        }
        $nfh->close;
    }
    else {
        my @results = `samtools view $tumor_bam $chr:$stop-$stop`;
        for my $result (@results){
            chomp $result;
            my @details = split /\t/, $result;
            # Parse overlapping reads for cigar strings containing insertions or deletions
            if($details[5] =~ m/[ID]/){
                $tumor_read_sw_support++;
            }
            else {
                $tumor_read_support++;
            }
        }

        @results = `samtools view $normal_bam $chr:$stop-$stop`;

        for my $result (@results){
            chomp $result;
            my @details = split /\t/, $result;
            # Parse overlapping reads for insertions or deletions
            if($details[5] =~ m/[ID]/){
                $normal_read_sw_support++;
            }
            else {
                $normal_read_support++;
            }

        }
    }

    my $is_lq = 1;

    my $p_value = Genome::Statistics::calculate_p_value(
                    $normal_read_support, 
                    $normal_read_sw_support, 
                    $tumor_read_support, 
                    $tumor_read_sw_support);
    if($p_value == 1) {
        $p_value = Genome::Statistics::calculate_p_value(
                    $normal_read_sw_support, 
                    $normal_read_support, 
                    $tumor_read_sw_support, 
                    $tumor_read_support);
    }

    my $read_support = join("\t",($reads,$tumor_read_support,$tumor_read_sw_support,$normal_read_support,$normal_read_sw_support,$pos_percent,$p_value))."\n";
    #my ($chr,$start,$stop,$refvar,$pindel_reads,$t_reads,$t_sw_reads,$n_reads,$n_sw_reads,$ps, $p_value) = split "\t", $line;

    if($p_value <= .15) { #assuming significant smith waterman support, trust the fishers exact test to make a germline determination
        $is_lq=0;
    }
    if(($tumor_read_sw_support + $tumor_read_support < 10) && ($reads > $tumor_read_sw_support)) { #low coverage area, and pindel found more reads than were smith waterman mapped available-- rescue this from pvalue filter
        $is_lq=0;
    }
    my $bed = $self->get_bed_line(\@event);
    if(defined($self->_read_support_fh)){
        my $read_support_fh = $self->_read_support_fh;
        print $read_support_fh $bed."\t".$read_support;
    }

    if($self->create_hq_raw_reads && ($is_lq==0)){
        $self->print_raw_read(\@event,$hq_raw);
    }
    elsif($is_lq && $self->create_hq_raw_reads){
        $self->print_raw_read(\@event,$lq_raw);
    }
    else {
        die $self->error_message("Nothing to do!");
    } 
    if(system("rm -f $temp")){
        die $self->error_message("Failed to complete the command: rm -f $temp");
    }
    
    return 1;
}

sub vaf_filter {
    my $self = shift;
    my $event = shift;
    my $hq_raw = $self->_hq_raw_output_fh;
    my $lq_raw = $self->_lq_raw_output_fh;
    my @event = @{ $event };
    my ($call, $reference, @support) = @event;
    my $tumor_bam = $self->aligned_reads_input;

    unless(-s $tumor_bam){
        die $self->error_message("Could not locate tumor_bam at: ".$tumor_bam);
    }

    my @call_fields = split /\s/, $call;

    my $type = $call_fields[1];
    my $size = $call_fields[2];
    my $chr = $call_fields[7];
    my $start = $call_fields[9];
    my $stop = $call_fields[10];

    my $reads = scalar(@support);
    my $temp = Genome::Sys->create_temp_file_path;
    my $tumor_read_support=0;
    my $tumor_read_sw_support=0;

    my $cap = $self->capture_data;

    if($cap){
        my $tsam_cmd = "samtools view $tumor_bam $chr:$stop-$stop > $temp";
        if(system($tsam_cmd)){
            die $self->error_message("Failed to run the command: $tsam_cmd");
        }
        $tumor_read_support = $self->line_count($temp);        
    }
    else {
        my @results = `samtools view $tumor_bam $chr:$stop-$stop`;
        $tumor_read_support = scalar(@results);
    }



    # Call samtools over the variant start-stop to get overlapping reads from the tumor bam
    #my @results = `samtools view $tumor_bam $chr:$stop-$stop`;
    #my $tumor_read_support=scalar(@results);

    # Currently, there are some instances when samtools view runs for along time, and either does not return,
    # or returns nothing, resulting in zero tumor_read_support, which itself causes a divide by zero error.
    # SO, in the interests of testing this filter on real data, I have added an override which sets 
    # tumor_read_support to 1 and fails the site, but notifies the user.  rlong 8/24/2011


    my $override=0;
    unless($tumor_read_support > 0){
        $self->debug_message("Found ".$tumor_read_support." reads in the tumor, but pindel found: " . $reads . " supporting the variant at ".join("\t",($chr,$start,$stop))."\n");
        $tumor_read_support =1;
        $override = 1;
    }

    my $variant_calls = scalar( grep{ m/tumor/} @support);
    my $vaf = $variant_calls / $tumor_read_support;
    my $vaf_cutoff_met = 0;
    my $cutoff = $self->variant_freq_cutoff;
    if ($vaf >= $cutoff) {
        $vaf_cutoff_met = 1;
    } 

    if($override){
        $vaf_cutoff_met=0;
    }

    my $bed = $self->get_bed_line(\@event);
    chomp $bed if $bed;

    if($self->create_hq_raw_reads && ($vaf_cutoff_met)){
        $self->print_raw_read(\@event,$hq_raw);
    }
    elsif((!$vaf_cutoff_met) && $self->create_hq_raw_reads){
        $self->print_raw_read(\@event,$lq_raw);
    }
    else {
        die $self->error_message("Nothing to do!");
    } 
    return 1;
}

sub line_count {
    my $self = shift;
    my $input = shift;
    unless( -e $input ) {
        die $self->error_message("Could not locate file for line count: $input");
    }
    my $result = `wc -l $input`;
    my ($answer)  = split /\s/,$result;
    return $answer
}

sub print_raw_read {
    my $self = shift;
    my $event = shift;
    my $raw_fh = shift;
    my @event = @{ $event };
    print $raw_fh "####################################################################################################\n";
    for my $line (@event){
        print $raw_fh $line."\n";
    }
    return 1;
}

sub get_bed_line {
    my $self = shift;
    my $event = shift;
    my @event = @{$event};
    my ($call,$reference,$read) = @event;
    my @bed_line = $self->parse($call, $reference, $read);
    unless((@bed_line)&& scalar(@bed_line)==5){
        return;
    }
    my $bed = join("\t",@bed_line);
    return $bed;
}

sub parse {
    my $self=shift;
    my ($call, $reference, $first_read) = @_;
    my $refseq = $self->_refseq;
    my @call_fields = split /\s+/, $call;
    my $type = $call_fields[1];
    my $size = $call_fields[2];
    my ($chr, $start, $stop);
    $chr = $call_fields[7];
    $start = $call_fields[9];
    $stop = $call_fields[10];    
    my $support = $call_fields[15];
    my ($ref, $var);
    if($type =~ m/D/) {
        $var =0;
        ###Make pindels coordinates(which seem to be last undeleted base and first undeleted base) 
        ###conform to our annotators requirements
        $stop = $stop -1;
        ###also deletions which don't contain their full sequence should be dumped to separate file
        my $allele_string;
        my $start_for_faidx = $start+1; 
        my $sam_default = Genome::Model::Tools::Sam->path_for_samtools_version;
        my $faidx_cmd = "$sam_default faidx " . $self->_refseq . " $chr:$start_for_faidx-$stop"; 
        my @faidx_return= `$faidx_cmd`;
        shift(@faidx_return);
        chomp @faidx_return;
        $allele_string = join("",@faidx_return);

        $ref = $allele_string;
    }
    elsif($type =~ m/I/) {
        #misunderstanding of bed format
        #0 based numbers teh gaps so an insertion of any number of bases between base 10 and 11 in 1base
        #is 10 10 in bed format
        #$start = $start - 1;
        $ref=0;
        my ($letters_until_space) =   ($reference =~ m/^([ACGTN]+) /);
        my $offset_into_first_read = length($letters_until_space);
        $var = substr($first_read, $offset_into_first_read, $size);
        $stop = $stop - 1;
    }
    if($size >= 100) {
        my $big_fh = $self->_big_output_fh;
        print $big_fh join("\t",($chr,$start,$stop,$size,$support))."\n";
        return undef;
    }
    return ($chr,$start,$stop,$ref,$var);
}
