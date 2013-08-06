package Genome::Model::Tools::Validation::BuildRemappingContigs;

use strict;
use warnings;

use Genome;
use Command;
use Genome::Sys;    #for file parsing etc
use POSIX; #for rounding
use Bio::PrimarySeq;    #necessary to do pairwise alignment with Bio::dpAlign
use Bio::Tools::dpAlign;    #for pairwise alignment of the contigs
use Bio::SimpleAlign;   #dpAlign returns alignments in this format
use Bio::AlignIO;
use Error qw(:try); #for BioPerl exception handling


class Genome::Model::Tools::Validation::BuildRemappingContigs {
    is => 'Command',
    has => [
        input_file => {
            type => 'String',
            is_optional => 0,
            doc => 'The input file used for assembly. This is used to track which variants assembled',
        },
        tumor_assembly_file => {
            type => 'String',
            doc => 'The tumor assembly output file',
        },
        tumor_assembly_breakpoints_file => {
            type => 'String',
            is_optional => 0,
            doc => 'The file of tumor breakpoints picked by assembly in fasta format',
        },
        normal_assembly_file => {
            type => 'String',
            doc => 'The normal assembly output file',
        },
        normal_assembly_breakpoints_file => {
            type => 'String',
            is_optional => 0,
            doc => 'The file of normal breakpoints picked by assembly in fasta format',
        },
        relapse_assembly_breakpoints_file => {
            type => 'String',
            is_optional => 1,
            doc => 'The file of relapse breakpoints picked by assembly in fasta format',
        },
        reference_sequence => {
            type => 'String',
            is_optional => 0,
            example_values => ['\'/gscmnt/839/info/medseq/reference_sequences/NCBI-human-build36/all_sequences.fa\' (GC standard Build36)'],
            doc => 'The samtools indexed fasta sequence of the reference the indels were predicted against.'
        },
        output_file => {
            type => 'String',
            is_optional => 0,
            doc => 'File to dump all of the contigs to',
        },
        contig_size => {
            type => 'Integer',
            is_optional => 0,
            default => 500, #force 100bp reads to align across the variant
            doc => 'The intended size of the contigs. If contigs overlap then they may be merged',
        },
        append_indel_alleles => {
            type => 'Boolean',
            is_optional => 1,
            default => 0,
            doc => 'append the indel alleles to the contig name. Requires that the input file was generated with indel alleles appended.'
        },
        minimum_local_overlap_size => {
            type => 'Integer',
            is_optional => 0,
            default => 150, 
            doc => 'During contig merging, the minimum length of the region around the indel to compare across contigs for similarity. Leave this alone unless read lengths are not 100',
        },
        samtools_version => {
            type => 'String',
            is_optional => 1,
            default => 'r783',
            doc => "The gsc version string for an installed version of samtools. You probably don't want to change this.",
        },
        _samtools_exec => {
            type => 'String',
            is_optional => 1,
        },
        _alignment_factory_object => {
            type => 'Object',
            is_optional => 1,
        },
        _generation_stats => {
            type => 'Hashref',
            is_optional => 1,
            default => {},
        },        
    ]
};


sub execute {
    my $self=shift;


    #we will need to handle the following
    #Reading in of contigs
    #   * make sure to account for strand
    #normal assembly contigs - what will this mean if   1) they are the same
    #                                                   2) they are different
    #                                                   3) they don't exist or exist in a different order    
    #uniform contig size
    #detection and merging of overlapping contigs
    #then spit out contigs for matching
    #   * contig names should contain all info necessary to count the variant in question

    #check that we are on a 64bit system and can run samtools
    
    unless (POSIX::uname =~ /64/) {
        $self->error_message("This script requires a 64bit system");
        return;
    }

   
    #check that the reference exists before doing anything
    Genome::Sys->validate_file_for_reading($self->reference_sequence); #this should croak if the file is invalid

    #set the same executable path on the object
    $self->_samtools_exec(Genome::Model::Tools::Sam->path_for_samtools_version($self->samtools_version));

    #open the output file for writing
    my $output_fh = IO::File->new($self->output_file,"w");
    unless($output_fh) {
        $self->error_message("Unable to open " . $self->output_file . " for writing.");
        return;
    }

    #parse the assembly input file
    my $input_fh = Genome::Sys->open_file_for_reading($self->input_file); #this should die if it fails
    my %expected_contigs;
    my %alleles;
    while(my $line = $input_fh->getline) {
        chomp $line;
        my @fields = split /\t/, $line;
        next if $fields[0] =~ /^#/; #skip header or comments
        next if $fields[7] < 3; #skip small indels
        my $id = join(".",@fields[0,1,3,4,6,7],"+-");
        $expected_contigs{$id} = {};
        $alleles{$id} = join("_",@fields[8..9]);
    }
    

    #read in files
    my $tumor_breakpoints = $self->tumor_assembly_breakpoints_file;
    my $normal_breakpoints = $self->normal_assembly_breakpoints_file;
    my $tumor_assembly_file = $self->tumor_assembly_file;
    my $normal_assembly_file = $self->normal_assembly_file;


    my $tumor_contigs = $self->read_in_breakpoints($tumor_breakpoints, $tumor_assembly_file, 'tumor');
    unless($tumor_contigs) {
        $self->error_message("Unable to parse $tumor_breakpoints");
        return;
    }
    my $normal_contigs = $self->read_in_breakpoints($normal_breakpoints, $normal_assembly_file, 'normal');
    unless($normal_contigs) {
        $self->error_message("Unable to parse $normal_breakpoints");
        return;
    }


    #pool all the contigs
    my @contigs = (@$tumor_contigs,@$normal_contigs);
    my @resized_contigs = ();

    #add in the relapse if provided

    if($self->relapse_assembly_breakpoints_file) {
        my $relapse_breakpoints = $self->relapse_assembly_breakpoints_file;
        my $relapse_contigs = $self->read_in_breakpoints($relapse_breakpoints, 'relapse');
        unless($relapse_contigs) {
            $self->error_message("Unable to parse $relapse_breakpoints");
            return;
        }
        push @contigs, @$relapse_contigs;
    }

    #resize the contigs
    for my $contig (@contigs) {
        $self->resize_contig($contig); #this edits the contig in place to reach a desired size
        #print STDERR "> resized contig\n",$contig->{'sequence'},"\n";

        unless($contig->{'contig_start'} < $contig->{'contig_stop'}) {
            #we will skip those until they are fixed
            $self->status_message("Trimmed genomic coordinates' start is greater than stop for variant starting at " . $contig->{'pred_pos1'} . "\n");
        }
        else {
            push @resized_contigs, $contig;
        }
    }

    #for testing purposes now want to report if any of the contigs overlap at all
    #what we would really want to do is sort all contigs tumor/normal and report overlap
    my $stats_hash = $self->_generation_stats;
    $stats_hash->{number_of_overlapping_contigs} = 0;
    $stats_hash->{number_of_contigs_total} = 0;

    my $current_chr = '';
    my $current_start = 0;
    my $current_stop = 0;
    my @overlapping_contigs = ();
    for my $contig (sort {_sort_contigs($a,$b)} @resized_contigs) {
        if(exists($expected_contigs{$contig->{id}})) {
            $expected_contigs{$contig->{id}}{$contig->{source}} = 1;
        }
        else {
            $self->status_message("Unexpected contig id: " . $contig->{id});
        }

        $stats_hash->{number_of_contigs_total}++;
        #these should be sorted by chromosome, start, stop now
        #check for overlap with the current region

        if($current_chr eq $contig->{pred_chr1} && $current_start <= $contig->{contig_start} && $current_stop >= $contig->{contig_start}) {
            $self->status_message("Overlap detected for " . join(".",$contig->{pred_chr1},$contig->{pred_pos1},$contig->{pred_pos2},$contig->{pred_type},$contig->{pred_size}));
            $stats_hash->{number_of_overlapping_contigs}++;
            if($current_stop < $contig->{contig_stop}) {
                #roll into region we are intersecting with
                $current_stop = $contig->{contig_stop};
            }
        }
        else {
            #no overlap, handle last region's set of contigs
            if(@overlapping_contigs > 1) {
                @overlapping_contigs = $self->handle_overlap(@overlapping_contigs);
            }
            if(@overlapping_contigs) {
                #Here we will print out a contig
                #make sure that the info to count the contig is present and make sure that lines are shortish
                foreach my $unique_contig (@overlapping_contigs) {                    
                    if($self->append_indel_alleles){
                        print $output_fh ">",join("_",(@$unique_contig{qw( pred_chr1 pred_pos1 pred_pos2 pred_type source)},$alleles{$unique_contig->{id}}));
                    } else {
                        print $output_fh ">",join("_",@$unique_contig{qw( pred_chr1 pred_pos1 pred_pos2 pred_type source) });
                    }

                    #need to code in the range on the contig for the variant as well as the range on the reference to count. For non-overlapping contigs this is simple. Let's also code overlap status
                    printf $output_fh " Overlap:%d",@overlapping_contigs - 1;   #this should code the number of other contigs overlapping the contig
                    printf $output_fh " Ref:%s.%d.%d Con:%d.%d",@$unique_contig{qw( assem_chr1 assem_pos1 assem_pos2 contig_location_of_variant microhomology_contig_endpoint )};

                    #Here we print out annotation information
                    #swap the coords if they are reversed gah!
                    if($unique_contig->{assem_pos1} > $unique_contig->{assem_pos2}) {
                        $self->status_message("Genomic coordinates of indel make no sense for variant starting at " . $unique_contig->{'pred_pos1'} . "\n");
                        ($unique_contig->{assem_pos1},$unique_contig->{assem_pos2}) = ($unique_contig->{assem_pos2},$unique_contig->{assem_pos1});
                    }

                    my ($annotation_start,$annotation_end,$ref,$var);
                    if($unique_contig->{assem_type} eq 'INS') {
                        $annotation_start = $unique_contig->{assem_pos1} - 1; #base before the event
                        $annotation_end = $annotation_start + 1;
                        $ref = 0;
                        $var = uc(substr($unique_contig->{sequence},$unique_contig->{contig_location_of_variant} - 1, $unique_contig->{assem_size})); #subtracting one in order to change to index as contig_location is the first base of the variant
                    }
                    else {
                        $annotation_start = $unique_contig->{assem_pos1};
                        $annotation_end = $annotation_start + $unique_contig->{assem_size} - 1;
                        $ref = uc($self->fetch_flanking_sequence($unique_contig->{assem_chr1},$annotation_start,$annotation_end));
                        $var = 0;
                    }
                    
                    printf $output_fh " Anno:%s.%d.%d.%s.%s.%s", $unique_contig->{assem_chr1},$annotation_start,$annotation_end,$ref,$var,$unique_contig->{assem_type};
                    printf $output_fh "\n";

                    #print sequence with each line containing 80bp
                    my $sequence = $unique_contig->{sequence};
                    while($sequence) {
                        print $output_fh substr($sequence,0,80,""),"\n";
                    }
                }
            }
            @overlapping_contigs = ();
            $current_chr = $contig->{pred_chr1};
            $current_start = $contig->{contig_start};
            $current_stop = $contig->{contig_stop};
        }
        push @overlapping_contigs, $contig;
    }


    #----------------
    #flush the buffer and output the last contig (same code as above, should probably be dropped into a function later)
    if(@overlapping_contigs){
        if(@overlapping_contigs > 1) {
            @overlapping_contigs = $self->handle_overlap(@overlapping_contigs);
        }
        foreach my $unique_contig (@overlapping_contigs) {                    
            print STDERR "dumping $unique_contig->{id}\n";
            if($self->append_indel_alleles){
                print $output_fh ">",join("_",(@$unique_contig{qw( pred_chr1 pred_pos1 pred_pos2 pred_type source)},$alleles{$unique_contig->{id}}));
            } else {
                print $output_fh ">",join("_",@$unique_contig{qw( pred_chr1 pred_pos1 pred_pos2 pred_type source) });
            }

            #need to code in the range on the contig for the variant as well as the range on the reference to count. For non-overlapping contigs this is simple. Let's also code overlap status
            printf $output_fh " Overlap:%d",@overlapping_contigs - 1;   #this should code the number of other contigs overlapping the contig
            printf $output_fh " Ref:%s.%d.%d Con:%d.%d",@$unique_contig{qw( assem_chr1 assem_pos1 assem_pos2 contig_location_of_variant microhomology_contig_endpoint )};

            #Here we print out annotation information
            #swap the coords if they are reversed gah!
            if($unique_contig->{assem_pos1} > $unique_contig->{assem_pos2}) {

                $self->status_message("Genomic coordinates of indel make no sense for variant starting at " . $unique_contig->{'pred_pos1'} . "\n");
                ($unique_contig->{assem_pos1},$unique_contig->{assem_pos2}) = ($unique_contig->{assem_pos2},$unique_contig->{assem_pos1});
            }

            my ($annotation_start,$annotation_end,$ref,$var);
            if($unique_contig->{assem_type} eq 'INS') {
                $annotation_start = $unique_contig->{assem_pos1} - 1; #base before the event
                $annotation_end = $annotation_start + 1;
                $ref = 0;
                $var = uc(substr($unique_contig->{sequence},$unique_contig->{contig_location_of_variant} - 1, $unique_contig->{assem_size})); #subtracting one in order to change to index as contig_location is the first base of the variant
            }
            else {
                $annotation_start = $unique_contig->{assem_pos1};
                $annotation_end = $annotation_start + $unique_contig->{assem_size} - 1;
                $ref = uc($self->fetch_flanking_sequence($unique_contig->{assem_chr1},$annotation_start,$annotation_end));
                $var = 0;
            }
            
            printf $output_fh " Anno:%s.%d.%d.%s.%s.%s", $unique_contig->{assem_chr1},$annotation_start,$annotation_end,$ref,$var,$unique_contig->{assem_type};
            printf $output_fh "\n";

            #print sequence with each line containing 80bp
            my $sequence = $unique_contig->{sequence};
            while($sequence) {
                print $output_fh substr($sequence,0,80,""),"\n";
            }
        }
    }
    #----------------




    #figure out how many indels have no contig at all
    #and how many have a contig in normal and how many from normal
    foreach my $id (keys %expected_contigs) {
        $stats_hash->{number_attempted}++;
        if(!keys %{$expected_contigs{$id}}) {
            $stats_hash->{number_failed_assembly}++;
        }
        else {
            $stats_hash->{number_passed_assembly}++;
        }
        if(exists($expected_contigs{$id}{tumor}) && $expected_contigs{$id}{tumor} == 1) {
            $stats_hash->{number_with_tumor}++;
        }
        if(exists($expected_contigs{$id}{normal}) && $expected_contigs{$id}{normal} == 1) {
            $stats_hash->{number_with_normal}++;
        }
        if($self->relapse_assembly_breakpoints_file && exists($expected_contigs{$id}{relapse}) && $expected_contigs{$id}{relapse} == 1) {
            $stats_hash->{number_with_relapse}++;
        }
        if(exists($expected_contigs{$id}{tumor}) && exists($expected_contigs{$id}{normal}) && $expected_contigs{$id}{normal} == 1 && $expected_contigs{$id}{tumor} == 1) {
            $stats_hash->{number_with_both}++;
        }
        if(exists($expected_contigs{$id}{tumor}) || exists($expected_contigs{$id}{relapse})) {
            $stats_hash->{number_with_contig_in_either_tumor_or_relapse}++;
        }
    }

            

    #print stats
    foreach my $key (keys %$stats_hash) {
        $self->status_message("$key: " . $stats_hash->{$key});
    }
        
    return 1;

}

sub read_in_breakpoints {
    my ($self, $breakpoint_file, $assembly_file, $source) = @_;
    my $assembly_info_array = Genome::Model::Tools::TigraSv->parse_assembly_file($assembly_file);
    my %assembly_info = map {$_->{prefix} => $_} @$assembly_info_array;

    my $fh = Genome::Sys->open_file_for_reading($breakpoint_file);
    if($fh) {
        my @contigs;
        my $current_contig = {};
        my $current_contig_sequence = "";
        while(my $line = $fh->getline) {
            chomp $line;
            if($line =~ /^>/) { #it's a header line
                #resolve last contig
                if(%$current_contig) {
                    $current_contig->{'sequence'} = $current_contig_sequence;
                    if(defined $source) {
                        $current_contig->{'source'} = $source;
                    }
                    #adjust strand
                    if($current_contig->{'strand'} ne '+') {
                        #print STDERR "> original contig\n$current_contig_sequence\n";
                        $current_contig->{'sequence'} =~ tr/ACGTacgt/TGCAtgca/;
                        $current_contig->{'sequence'} = reverse $current_contig->{'sequence'};
                        $current_contig->{'strand'} = '+';

                        #also need to swap the genomic coordinates
                        ($current_contig->{'contig_start'},$current_contig->{'contig_stop'}) = ($current_contig->{'contig_stop'},$current_contig->{'contig_start'});

                        #lastly need to adjust the position of the indel
                        #roughly this is: the length of the contig - (the old position - 1) equals the new position in the reversed contig. The last base of the indel is: contig_location + size - 1 if an insertion and contig_location - 1 if deletion
                        my $type_size_toggle_var = $current_contig->{'assem_type'} eq 'DEL' ? 0 : 1;
                        my $temp_var_location = $current_contig->{'contig_location_of_variant'};
                        $current_contig->{'contig_location_of_variant'} = $current_contig->{'length'} - $current_contig->{'microhomology_contig_endpoint'} + 1;
                        $current_contig->{'microhomology_contig_endpoint'} = $current_contig->{'length'} - $temp_var_location + 1;
                        #print STDERR "> reverse complemented contig\n",$current_contig->{'sequence'},"\n";
                    }

                    #check to make sure we didn't get back something crazy
                    if($current_contig->{'assem_type'} !~ /INS|DEL|ITX/i) {
use Data::Dumper;
print Dumper($current_contig);
                        $self->statusmessage("Skipping contig that assembled as a type other than insertion, tandem duplication (ITX) or deletion with variant starting at " . $current_contig->{'pred_pos1'});
                    }
                    else {
                        #check that Ins field makes sense
                        unless($current_contig->{'contig_location_of_variant'} <= $current_contig->{'microhomology_contig_endpoint'}) {
                            $self->status_message("Microhomology makes no sense for variant starting at " . $current_contig->{'pred_pos1'} . "\n");
                        }

                        #check that coordinates make sense
                        unless($current_contig->{'contig_start'} < $current_contig->{'contig_stop'}) {
                            #we will skip those until they are fixed
                            $self->status_message("Genomic coordinates make no sense for variant starting at " . $current_contig->{'pred_pos1'} . "\n");
                            $self->status_message(join(" - ",($current_contig->{'contig_start'},$current_contig->{'contig_stop'})) . "\n");
                        }
                        else {
                            push @contigs, $current_contig;
                        }
                    }

                }
                $current_contig = $self->parse_breakpoint_contig_header($line, \%assembly_info);
                $current_contig_sequence = "";
            }
            else {
                $current_contig_sequence .= $line;
            }
        }
        #resolve last indel of file
        if(%$current_contig) {
            $current_contig->{'sequence'} = $current_contig_sequence;
            if(defined $source) {
                $current_contig->{'source'} = $source;
            }
            #adjust strand
            if($current_contig->{'strand'} ne '+') {
                #print STDERR "> original contig\n$current_contig_sequence\n";
                $current_contig->{'sequence'} =~ tr/ACGTacgt/TGCAtgca/;
                $current_contig->{'sequence'} = reverse $current_contig->{'sequence'};
                $current_contig->{'strand'} = '+';

                #also need to swap the genomic coordinates
                ($current_contig->{'contig_start'},$current_contig->{'contig_stop'}) = ($current_contig->{'contig_stop'},$current_contig->{'contig_start'});

                #lastly need to adjust the position of the indel
                #roughly this is: the length of the contig - (the old position - 1) equals the new position in the reversed contig. The last base of the indel is: contig_location + size - 1 if an insertion and contig_location - 1 if deletion
                my $type_size_toggle_var = $current_contig->{'assem_type'} eq 'DEL' ? 0 : 1;
                my $temp_var_location = $current_contig->{'contig_location_of_variant'};
                $current_contig->{'contig_location_of_variant'} = $current_contig->{'length'} - $current_contig->{'microhomology_contig_endpoint'} + 1;
                $current_contig->{'microhomology_contig_endpoint'} = $current_contig->{'length'} - $temp_var_location + 1;
                #print STDERR "> reverse complemented contig\n",$current_contig->{'sequence'},"\n";
            }

            #check to make sure we didn't get back something crazy
            if($current_contig->{'assem_type'} !~ /INS|DEL|ITX/i) {
                $self->status_message("Skipping contig that assembled as a type other than insertion, tandem duplication (ITX) or deletion with variant starting at " . $current_contig->{'pred_pos1'});
            }
            else {
                #check that Ins field makes sense
                unless($current_contig->{'contig_location_of_variant'} <= $current_contig->{'microhomology_contig_endpoint'}) {
                    $self->status_message("Microhomology makes no sense for variant starting at " . $current_contig->{'pred_pos1'} . "\n");
                }

                #check that coordinates make sense
                unless($current_contig->{'contig_start'} < $current_contig->{'contig_stop'}) {
                    #we will skip those until they are fixed
                    $self->status_message("Genomic coordinates make no sense for variant starting at " . $current_contig->{'pred_pos1'} . "\n");
                    $self->status_message(join(" - ",($current_contig->{'contig_start'},$current_contig->{'contig_stop'})) . "\n");

                }
                else {
                    push @contigs, $current_contig;
                }
            }

        }
        return \@contigs;
    }
    else {
        #propagate to caller
        return;
    }
}

#this parses the header into a hash containing the relevant info
#contig "objects" should be better defined and there needs to be some sort of determination that the header is actually generating a valid object

sub parse_breakpoint_contig_header {
    my ($self, $header_line, $assembly_info) = @_;
    my %contig;
    $header_line =~ s/^.//;    #remove the caret
    my @header_fields = split ",", $header_line;
    for my $field (@header_fields) {
        my ($tag,$entry) = $field =~ /(\w+):*(\S*)/;
        if($tag =~ /^ID/) {
            die "Sequence prefix $entry not found in assembly info!" unless exists $assembly_info->{$entry};
            my $e = $assembly_info->{$entry};
            #this contains info about the original call
            @contig{ qw( pred_chr1 pred_pos1 pred_chr2 pred_pos2 pred_type pred_size pred_orientation) } = @$e{ qw(chr1 pos1 chr2 pos2 type ori) };

            #also store the id for identification later
            $contig{'id'} = $entry;
        }
        elsif($tag =~ /^CrossMatch/) {
            #this contains info about what the crossmatch parsing thing thought was the variant
            #remember that the reported coordinates include any microhomology
            #In the future Ken says variant supporting reads would be reads that overlap this region with at least N bases of overlap where N is the microhomology size
            @contig{ qw( assem_chr1 assem_pos1 assem_chr2 assem_pos2 assem_type assem_size assem_orientation) } = split('\|', $entry);

            #convert ITX into INS
            if($contig{assem_type} eq 'ITX') {
                $contig{assem_type} = 'INS';
            }
        }
        elsif($tag =~ /^Strand/) {
            #the strand of the contig. If minus then it needs to be reverse complemented
            $contig{'strand'} = $entry;
        }
        elsif($tag =~ /^Length/) {
            #the reported length of the contig sequence
            $contig{'length'} = $entry;
        }
        elsif($tag =~ /^Ins/) {
            #the microhomology of the breakpoint
            #eg 23-24 means 2 bp of microhomology in the contig itself. 
            #   200-- means no microhomology
            my ($start, $stop) = $entry =~ /(\d+)\-(\d*)/;
            if($stop) {
                #then there IS microhomology
                $contig{'microhomology_contig_endpoint'} = $stop;  #the last base of microhomology on the contig
            }
            else {
                #this is the first base of the variant
                #this will ONLY ever happen for deletions
                $contig{'microhomology_contig_endpoint'} = $start;
            }

            #set the location of the variant in the contig
            $contig{'contig_location_of_variant'} = $start;
        }
        elsif($tag =~ /Ref_start/) {
            #this is the position of the contig on the reference
            #note that this might not correspond to the first base of the contig. For now we're going to include mismatching bases as insertions/mismatches
            $contig{'contig_start'} = $entry;
        }
        elsif($tag =~ /Ref_end/) {
            #this is the end position of the contig on the reference
            #note that this might not correspond to the last base of the contig. For now we're going to include mismatching bases as insertions/mismatches
            $contig{'contig_stop'} = $entry;
        }
        elsif($tag =~ /Contig_start/) {
            #this is the position of the first aligned base of the contig
            #note that this might not correspond to the first base of the contig. For now we're going to include mismatching bases as insertions/mismatches
            $contig{'align_start'} = $entry;
        }
        elsif($tag =~ /Contig_end/) {
            #this is the position of the last aligned base of the contig
            #note that this might not correspond to the last base of the contig. For now we're going to include mismatching bases as insertions/mismatches
            $contig{'align_stop'} = $entry;
        }
        else {
            #we'll ignore all other fields for now
        }    
    }
    return \%contig;
}

sub resize_contig {
    my ($self, $contig) = @_;

    my $desired_size = $self->contig_size;

    #the leftmost position of the variant should be available as well as the length. For contigs where the ratio of variant sequence to contig is high, this MAY be biased towards having variants extend into the rightmost flank so we will need to do the math and appropriately trim
    #calculations that will be needed for padding or trimming
    my $contig_length = $contig->{'length'};
    my $indel_size = $contig->{'assem_size'};
    my $type_size_toggle_var = $contig->{'assem_type'} eq 'DEL' ? 0 : 1;
    my $left_flank_size = $contig->{'contig_location_of_variant'} - 1; #the number of bases of the contig preceeding the variant
    my $right_flank_size = $contig_length - $contig->{'microhomology_contig_endpoint'}; #the number of bases of the contig succeeding the variant
    
    my $chr = $contig->{'assem_chr1'};
    my $contig_start = $contig->{'contig_start'};
    my $contig_stop = $contig->{'contig_stop'};

    my $desired_left_flank_size = ceil(($desired_size - ($contig->{'microhomology_contig_endpoint'} - $contig->{'contig_location_of_variant'} + 1)) / 2);   #round up to preferentially add to the left flank
    my $desired_right_flank_size = floor(($desired_size - ($contig->{'microhomology_contig_endpoint'} - $contig->{'contig_location_of_variant'} + 1)) / 2); #round down to preferentially shorten at the right flank

    my $change_needed_to_left_flank_size = $desired_left_flank_size - $left_flank_size;
    my $change_needed_to_right_flank_size = $desired_right_flank_size - $right_flank_size;

    if($change_needed_to_left_flank_size < 0) {
        #trim the existing contig
        substr($contig->{'sequence'},0,abs($change_needed_to_left_flank_size),"");
    }
    elsif($change_needed_to_left_flank_size > 0) {
        my $lstart = $contig_start - $change_needed_to_left_flank_size;
        my $lend = $contig_start - 1;
        my $additional_lseq = $self->fetch_flanking_sequence($chr,$lstart,$lend);
        unless(defined $additional_lseq) {
            $self->status_message("Unable to fetch additional sequence for padding the left flanking sequence");
            return;
        }
        $contig->{'sequence'} = join("",$additional_lseq, $contig->{'sequence'});
    }

    if($change_needed_to_right_flank_size < 0) {
        #trim the existing contig
        substr($contig->{'sequence'},$change_needed_to_right_flank_size, abs($change_needed_to_right_flank_size),"");
    }
    elsif($change_needed_to_right_flank_size > 0) {
        my $rstart = $contig_stop + 1;
        my $rend = $contig_stop + $change_needed_to_right_flank_size;
        my $additional_rseq = $self->fetch_flanking_sequence($chr,$rstart,$rend);
        unless(defined $additional_rseq) {
            $self->status_message("Unable to fetch additional sequence for padding the right flanking sequence");
            return;
        }

        #otherwise, pad the contig
        $contig->{'sequence'} = join("", $contig->{'sequence'}, $additional_rseq);
    }

    #generally update all coordinates here
    #on the left side
    #basically need to alter all contig relevant coordinates
    #and set the leftmost ref coordinate to the trimmed coordinate
    $contig->{'contig_start'} -= $change_needed_to_left_flank_size;
    $contig->{'contig_location_of_variant'} += $change_needed_to_left_flank_size;
    $contig->{'microhomology_contig_endpoint'} += $change_needed_to_left_flank_size;

    #on the right side, should only have to update the genomic coord
    $contig->{'contig_start'} += $change_needed_to_right_flank_size;

    #lastly set length to the correct number
    $contig->{'length'} = length $contig->{'sequence'};
    return 1;
}

sub fetch_flanking_sequence {
    my ($self,$chr,$start,$stop) = @_;

    my $sam_executable_path = $self->_samtools_exec;
    my $refseq = $self->reference_sequence;
    my $cmd = "$sam_executable_path faidx $refseq $chr:$start-$stop";

    my ($header,@seq) = `$cmd`;
    unless(defined $header) {
        $self->error_message("Error fetching sequence for $chr:$start-$stop");
        return;
    }
    else {
        return join("",map { chomp; lc($_); } @seq); #this should change the sequence into a single string regardless of length
    }
}

sub _sort_contigs {
    my ($a, $b) = @_;
    
    if($a->{pred_chr1} eq $b->{pred_chr1}) { 
        if($a->{contig_start} == $b->{contig_start}) {
            return $a->{contig_stop} <=> $b->{contig_stop};
        }
        else {
            return $a->{contig_start} <=> $b->{contig_start};
        }
    } 
    else { 
        return $a->{pred_chr1} cmp $b->{pred_chr1}; 
    } 
}

sub handle_overlap {
    my ($self, @contigs) = @_;

    my $stats_hash = $self->_generation_stats;

    #set up the pairwise alignment object if necessary
    my $alignment_factory;
    if($self->_alignment_factory_object) {
        $alignment_factory = $self->_alignment_factory_object;
    }
    else {
        #don't really know what appropriate parameters are so just use defaults
        $alignment_factory = new Bio::Tools::dpAlign(-alg => Bio::Tools::dpAlign::DPALIGN_ENDSFREE_MILLER_MYERS); 
        unless($alignment_factory) {
            $self->error_message("Unable to create a Bio::Tools::dpAlign object for pairwise alignment");
            die;
        }

        $self->_alignment_factory_object($alignment_factory);   #store on object so we only need one of these things
    }
    
    #this should handle the case where we want to merge multiple contigs
    my @contigs_considering = @contigs;
    my $last_contig_number = 0;

    while($last_contig_number != @contigs_considering) {
        $last_contig_number = @contigs_considering;

        my @contigs_considered = ();
        while(my $contig1 = pop @contigs_considering) {

            my $contig1_seq =  Bio::PrimarySeq->new( -seq => $contig1->{'sequence'}, -id  => join("-",$contig1->{id},$contig1->{source}),);
            my @contigs_to_continue_examining = ();

            while(my $contig2 = pop @contigs_considering) {
                my $contig2_seq = Bio::PrimarySeq->new( -seq => $contig2->{'sequence'}, -id => join("-",$contig2->{id},$contig2->{source}),);

                #do alignment
                my $alignment = $alignment_factory->pairwise_alignment($contig1_seq, $contig2_seq);
                my $alnout = new Bio::AlignIO(-format => 'pfam', -fh => \*STDERR);
                
                #experimentally grab the key region and only assess the mismatches there...
                #
                #contig_location_of_variant microhomology_contig_endpoint
                my ($contig1_var_start, $contig1_var_end) = @{$contig1}{qw(contig_location_of_variant microhomology_contig_endpoint)};
                my ($contig2_var_start, $contig2_var_end) = @{$contig2}{qw(contig_location_of_variant microhomology_contig_endpoint)};
                try {
                    my $contig1_alignment_start = $alignment->column_from_residue_number(join("-",$contig1->{id},$contig1->{source}),$contig1_var_start);
                    my $contig1_alignment_end = $alignment->column_from_residue_number(join("-",$contig1->{id},$contig1->{source}),$contig1_var_end);
                    my $contig2_alignment_start = $alignment->column_from_residue_number(join("-",$contig2->{id},$contig2->{source}),$contig2_var_start);
                    my $contig2_alignment_end = $alignment->column_from_residue_number(join("-",$contig2->{id},$contig2->{source}),$contig2_var_end);
                    my $msa_start = $contig2_alignment_start > $contig1_alignment_start ? $contig1_alignment_start : $contig2_alignment_start;
                    my $msa_end = $contig2_alignment_end > $contig1_alignment_end ? $contig2_alignment_end : $contig1_alignment_end;
                    my $sub_alignment_length = $msa_end - $msa_start + 1;
                    my $critical_region_size = $self->minimum_local_overlap_size;
                    my $size_difference = $critical_region_size - $sub_alignment_length;
                    my $pad = ceil($size_difference / 2);
                    $msa_start -= $pad;
                    $msa_start = 0 if $msa_start < 0;
                    my $sub_alignment = $alignment->slice($msa_start, $msa_end + $pad); #going past the end seems ok. Just need to make sure we're not negative


                    if($sub_alignment->percentage_identity >= 98 && $sub_alignment->gap_line !~ /[-]/ && $alignment->length >= ($self->contig_size - 5)) {  
                        #this is the simple case, the contigs are nearly 100% identical and we do nothing. This is essentially disappearing contig2
                        $stats_hash->{near_perfect_merged_overlaps}++;
                    }
                    else {
                        push @contigs_to_continue_examining, $contig2;
                        #write the alignment to stderr for debugging
                        $self->status_message(join("",">",$contig1->{id},"\n",$contig1->{sequence}));
                        $self->status_message(join("",">",$contig2->{id},"\n",$contig2->{sequence}));

                        $alnout->write_aln($alignment);
                    }
                }
                otherwise {
                    #some exception happened with the alignments.
                    #we'll just refuse to merge the contigs
                    push @contigs_to_continue_examining, $contig2;
                    #write the alignment to stderr for debugging
                    $self->status_message(join("",">",$contig1->{id},"\n",$contig1->{sequence}));
                    $self->status_message(join("",">",$contig2->{id},"\n",$contig2->{sequence}));

                    $alnout->write_aln($alignment);
                };


            }
            @contigs_considering = @contigs_to_continue_examining;
            push @contigs_considered, $contig1; #this contig is representative of the set
        }
        @contigs_considering = @contigs_considered;
    }
    return @contigs_considering;
}
        



1;

sub help_brief {
    "generates contigs from assembly results for remapping of reads"
}

sub help_detail {
    <<'HELP';
This commmand attempts to generate a reference sequence containing variant contigs appropriate for read remapping to determine validation status
HELP
}
