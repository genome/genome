package Genome::InstrumentData::AlignmentResult::Crossmatch;

use strict;
use warnings;
use File::Basename;

use Bio::SeqIO;
use Carp;
 
use Genome;

use Genome::InstrumentData::AlignmentResult::Crossmatch::NISC::Assembly::SeqCM;
use Genome::InstrumentData::AlignmentResult::Crossmatch::GTB::Sequencing::ReadIO;
use Genome::InstrumentData::AlignmentResult::Crossmatch::NISC::Assembly::Alignment;

class Genome::InstrumentData::AlignmentResult::Crossmatch {
    is => 'Genome::InstrumentData::AlignmentResult',
    has_constant => [
        aligner_name => { value => 'Crossmatch', is_param=>1 },
    ],
    has_transient_optional => [
         _Crossmatch_sam_cmd => { is=>'Text' }
    ]
};

sub required_arch_os { 'x86_64' }

# fill me in here with what compute resources you need.
sub required_rusage { 
    my $self = shift;
    my %params = $self->decomposed_aligner_params;
    my $cores_to_use = $params{'cores_to_use'};
    if ($cores_to_use =~ /(\d+)/) { 
        return "-R 'select[model!=Opteron250 && type==LINUX64 && tmp>90000 && mem>16000] span[hosts=1] rusage[tmp=90000, mem=16000]' -M 16000000 -n $1";
    } else {
        print "Could not determine number of cores to use from params! Defaulting to 4.";
        return "-R 'select[model!=Opteron250 && type==LINUX64 && tmp>90000 && mem>16000] span[hosts=1] rusage[tmp=90000, mem=16000]' -M 16000000 -n 4";
    }
}

sub _run_aligner {
    my $self = shift;
    my @input_pathnames = @_;
    
    my $tmp_dir = $self->temp_scratch_directory;
    my $fq_file = "$tmp_dir/merged.fq";
    my $staging_sam_file = "$tmp_dir/all_sequences.sam";
    
    # get refseq info
    my $ref_file = $self->get_reference_sequence_index->full_consensus_path('fa');

    my $crossmatch_path = Genome::Model::Tools::Crossmatch->path_for_crossmatch_version($self->aligner_version);
    # TODO (iferguso) figure out how to get nonmatching to work correctly. it currently seems to cause crashes.
    my $nonmatching = 0;

    my %aligner_params = $self->decomposed_aligner_params;
    my $cores_to_use = $aligner_params{'cores_to_use'};
    if ($cores_to_use =~ /(\d+)/) {
        $cores_to_use = $1;
    } else {
        print "Could not determine number of cores to use from params! Defaulting to 4.";
    }


    # TODO (iferguso) i don't believe cross_match handles paired end reads. if it does then this needs to be changed:
    # TODO (iferguso) should do something to handle the log files
    foreach my $input_file (@input_pathnames) {
        my $in_seq_obj = Bio::SeqIO->new(
            -file   => $input_file,
            -format => 'fastq',
        );

        # reached end indicates whether we have anything left in the input fq
        my $reached_end = 0;
        # these two make sure we've finished every job before moving on
        my $started_jobs = 0;
        my $finished_jobs = 0;
        # this is a counter of how many fa/qual files (chunks) we've split off of the input fq
        my $minor = 0;
        # just pushing running chunks onto an array, that get removed after completion
        my @cvs;


        # yes "jobs still executing" is the same as "not every job finished", but i want
        #   to be absolutely sure that everything is done and finished before we move on
        # as long as we still have fq to read in OR jobs still executing OR not every job finished
        while ($reached_end == 0 || @cvs > 0 || $started_jobs != $finished_jobs) {
            # now, as long as there are open jobs slots AND we have fq to read in...
            while (scalar @cvs < $cores_to_use && $reached_end == 0) {
                my %chunk;

                #### STEP 1: Convert fastq files into fasta+qual files
                # TODO (iferguso) instead of testing the filetype by looking at the extension, i'm just assuming it is a .fq file
                #if ($input_file =~ /(.+)\.fq$/) {
                if ($input_file =~ /(.+)/) {
                    $chunk{'fastq_infile'} = $input_file;
                    $chunk{'fasta_infile'} = "$1.$minor.fa";
                    $chunk{'qual_infile'}  = "$1.$minor.fa.qual";
                    $chunk{'align_output'} = "$1.$minor.cm.aligned.out";
                    $chunk{'sam_output'}   = "$1.$minor.sam.aligned.out";
                    $chunk{'nonmatching_fa_output'}   = "$1.$minor.fa.nonmatching";
                    $chunk{'nonmatching_qual_output'} = "$1.$minor.fa.nonmatching.qual";
                    $chunk{'log_output'}  = "$1.$minor.fa.log";
                    $chunk{'stdout_file'} = "$1.$minor.cm.stdout";
                    $chunk{'stderr_file'} = "$1.$minor.cm.stderr";
                } else {
                    $self->error_message(
                        "Problem parsing filename of input fq $input_file.\n".
                        "Possibly didn't end with .fq so couldn't extract name?");
                    $self->die($self->error_message);
                }
                $self->debug_message("Iteration $minor of splitting fq '$input_file' into '".$chunk{'qual_infile'}."' and '".$chunk{'fasta_infile'}."'.");

                my $out_seq_obj = Bio::SeqIO->new(
                    -file   => ">".$chunk{'fasta_infile'},
                    -format => 'fasta'
                );
                 
                my $out_qual_obj = Bio::SeqIO->new(
                    -file   => ">".$chunk{'qual_infile'},
                    -format => 'qual',
                );
                         
                my $counter = 0;
                while ($counter < 25000 && $reached_end == 0) {
                    my $fastq_obj = $in_seq_obj->next_seq;
                    if ($fastq_obj) {
                        $out_seq_obj->write_seq($fastq_obj);
                        $out_qual_obj->write_seq($fastq_obj);
                        $counter++;
                    } else { $reached_end = 1; }
                }

                my $cv; 

                # start our job
                if ($nonmatching) {
                    my $cmdline = $crossmatch_path . sprintf(' %s %s %s -output_nonmatching_queries > %s',
                        $chunk{'fasta_infile'}, $ref_file, $aligner_params{align_params}, $chunk{'align_output'});

                    $cv = Genome::Utility::AsyncFileSystem->shellcmd(
                        cmd             => $cmdline,
                        '>'             => $chunk{'stdout_file'},
                        '2>'            => $chunk{'stderr_file'},
                        input_files     => [ $chunk{'fasta_infile'}, $chunk{'qual_infile'}, $ref_file ],
                        output_files    => [ $chunk{'align_output'}, $chunk{'nonmatching_fa_output'}, $chunk{'nonmatching_qual_output'}, $chunk{'log_output'} ],
                        skip_if_output_is_present => 0,
                    );
                } else {
                    my $cmdline = $crossmatch_path . sprintf(' %s %s %s > %s',
                        $chunk{'fasta_infile'}, $ref_file, $aligner_params{align_params}, $chunk{'align_output'});

                    $cv = Genome::Utility::AsyncFileSystem->shellcmd(
                        cmd             => $cmdline,
                        '>'             => $chunk{'stdout_file'},
                        '2>'            => $chunk{'stderr_file'},
                        input_files     => [ $chunk{'fasta_infile'}, $chunk{'qual_infile'}, $ref_file ],
                        output_files    => [ $chunk{'align_output'}, $chunk{'log_output'} ],
                        skip_if_output_is_present => 0,
                    );
                }
                
                # save the actual async shellcmd object in the hash as well
                $chunk{'cv'} = $cv;
                # push this chunk's hash to the big list
                push @cvs, \%chunk;
                print "\n$cv\n";
                # increase our counter by one for naming purposes
                $minor++;
                # increase the number of jobs we've started
                $started_jobs++;
            }
            my @finished = _select_cv(@cvs);
            # remove finished from @cvs
            foreach my $f (@finished) {
                @cvs = grep { $_ != $f } @cvs;
            }

            foreach my $cv (@finished) {
                #### STEP 3: Parse outputs in SAM
                {
                    my $sam_output_fh = IO::File->new(">".$cv->{'sam_output'});
                    if ( !$sam_output_fh ) {
                        $self->error_message("Error opening sam output file for writing $!");
                        return;
                    }
                    $self->debug_message("For converting, opened ".$cv->{'sam_output'});
                    ### The following was adapted from code by Nancy Hansen <nhansen@mail.nih.gov>

                    ##########################################################
                    # Author:	Nancy Hansen
                    # Date:		3/27/09
                    # Function:	Parse a cross_match file, and output the
                    #               alignments in SAM format.
                    ##########################################################

                    my $paired; # option to specify that reads are paired
                    my $fastq_offset = 33; # offset for fastq read objects
                    my $hard_clip; # option to hard clip read sequences (i.e., not include unaligned ends)
                    my $norm_score = 40; # maximum score to normalize to
                    my $sq_file; # optional file containing the "SQ" sequence dictionary for the alignments
                    my $no_header; # option to not write a header to the file
                    my $max_bases; # option to use only first max_bases bases of each read

                    # TODO could be useful. might want to integrate into params from pipeline
                    #GetOptions("no_header" => \$no_header, "paired" => \$paired, 
                    #           "offset=i" => \$fastq_offset, "hard_clip" => \$hard_clip,
                    #           "sq_file=s" => \$sq_file, "max_bases=i" => \$max_bases );

                    #my $fastq_file = $ARGV[0];
                    #my $cons_fasta_offset = $ARGV[1];
                    # TODO this is obtuse and not logical. it's $ma because that's the "$input_file" (fq) that the finished job came from
                    my $fastq_file = $cv->{'fastq_infile'};
                    my $cons_fasta_offset = 0;

                    my $rh_read_seqs = {};
                    if ($fastq_file){ # read in quality values
                        my $read_io = Genome::InstrumentData::AlignmentResult::Crossmatch::GTB::Sequencing::ReadIO->new(
                            -file => $fastq_file,
                            -format => 'fastq', 
                            -qual_offset => $fastq_offset);

                        while (my $read_obj = $read_io->next_read()){
                            my $name = $read_obj->name();
                            $rh_read_seqs->{$name} = $read_obj;
                        }
                    }

                    my $cm_obj = Genome::InstrumentData::AlignmentResult::Crossmatch::NISC::Assembly::SeqCM->new(-file => $cv->{'align_output'});
                    my @sam_strings = ();
                    # divide into different hits:
                    my $rh_hit_length = {}; # store names of hits and their lengths

                    while (my $align_obj = $cm_obj->next_alignment()){
                        my $hit = $align_obj->hit();
                        my $hit_length = $align_obj->hit_length();
                        $rh_hit_length->{$hit} = $hit_length;
                        push @sam_strings, _sam_string($align_obj, $cons_fasta_offset, $paired, $hard_clip, $max_bases, $rh_read_seqs);
                    }

                    # write header:

                    $sam_output_fh->print(qq!\@HD\tVN:1.0\n!) if (!$no_header);

                    if ($sq_file){
                        my $sq_lines = `cat $sq_file`;
                        foreach my $sq_line (split /\n/, $sq_lines){
                            if ($sq_line =~ /^\@SQ/){
                                $sam_output_fh->print("$sq_line\n");
                            } else {
                                croak "Invalid SQ line in $sq_file:\n$sq_line\n";
                            }
                        }
                    }

                    my $max_score = 0;
                    foreach my $hit (keys %{$rh_hit_length}){
                        my $hit_length = $rh_hit_length->{$hit};
                        $sam_output_fh->print(qq!\@SQ\tSN:$hit\tLN:$hit_length\n!) if ((!$sq_file) && (!$no_header));
                    }

                    my $score_factor = ($max_score) ? $norm_score/$max_score : 'NA';

                    $sam_output_fh->print(qq!\@PG\tID:cross_match\tVN:1.080812\n!) if (!$no_header);

                    foreach my $sam_string (@sam_strings){
                        $sam_output_fh->print($sam_string);
                    }

                }

                #### STEP 4: Filter alignments
                {
                    unless (-s $cv->{'sam_output'}) {
                        $self->error_message("Unable to filter alignments. Sam file ".$cv->{'sam_output'}." is zero length, so something went wrong.");
                        $self->die($self->error_message);
                    }

                    # put your output file here, append to this file!
                        #my $output_file = $self->temp_staging_directory . "/all_sequences.sam"
                    die "Failed to process sam command line, error_message is ".$self->error_message unless $self->_filter_sam_output($cv->{'sam_output'}, $staging_sam_file);
                }

                $finished_jobs++;
            }
        }
    }

    return 1;
}

sub fillmd_for_sam {
    # It appears that pulse_cm2sam does not do md strings
    return 1;
}

sub _filter_sam_output {
    my ($self, $aligned_sam_file, $all_sequences_sam_file) = @_;
    # ccarey: bfast?! I assume this is a copypasta
    # some aligners output unaligned reads to a separate file. bfast appears to put them in the same file

    my $aligned_fh = IO::File->new( $aligned_sam_file );
    if ( !$aligned_fh ) {
            $self->error_message("Error opening output sam file for reading: $!");
            return;
    }
    $self->debug_message("Opened $aligned_sam_file");

    my $all_seq_fh = IO::File->new(">>$all_sequences_sam_file");
    if ( !$all_seq_fh ) {
        $self->error_message("Error opening all seq sam file for writing: $!");
        return;
    }
    $self->debug_message("Opened $all_sequences_sam_file");
    
    while (<$aligned_fh>) {
        #write out the aligned map, excluding the default header- all lines starting with @.
        $all_seq_fh->print($_) unless $_ =~ /^@/;
    }

    $aligned_fh->close;
    $all_seq_fh->close;
    return 1;
}

sub decomposed_aligner_params {
    my $self = shift;
    my $params = $self->aligner_params || "-raw -tags -bandwidth 3 -penalty -1 -gap_init -1 -gap_ext -1 -masklevel 0 -minscore 70";
    
    my @spar = split /\:/, $params;

    # TODO (iferguso) this could be changed: we currently specify number of cores to use right here:
    return ('align_params' => $spar[0], 'cores_to_use' => 4);
}

sub aligner_params_for_sam_header {
    my $self = shift;
    
    my %params = $self->decomposed_aligner_params;
    
    return 'cross_match' . $params{align_params};
}

sub _select_cv {
    while (1) {
        my @done = grep { $_->{'cv'}->ready } @_;

        if (@done) {
            return @done;
        } else {
            AnyEvent->one_event;
        }
    }
} 

### Following adapted from code by Nancy Hansen <nhansen@mail.nih.gov>
sub _sam_string {
    my $align = shift;

    #added
    my $cons_fasta_offset = shift;
    my $paired = shift;
    my $hard_clip = shift;
    my $max_bases = shift;
    my $rh_read_seqs = shift;

    my $query = $align->query();
    my $reference = $align->hit();
    my $query_string = $align->query_string();
    my $hit_left_end = $align->hit_left_end()+$cons_fasta_offset;
    my $comp = $align->comp();
    my $score = $align->score();
    my $flag = 0;
    $flag |= ($paired) ? 1 : 0; # paired reads?
    $flag |= ($comp eq 'C') ? 16 : 0;

    # create cigar string:
    my $cigar_string = $align->cigar_string(-hard_clip => $hard_clip, -max_bases => $max_bases);

    my $qual_string = ($rh_read_seqs->{$query}) ? $rh_read_seqs->{$query}->qual_string() : '';
    my $sequence = ($rh_read_seqs->{$query}) ? $rh_read_seqs->{$query}->seq() : '';
    if ($comp eq 'C') {
        $sequence = reverse $sequence;
        $sequence =~ tr/ATGCatgc/TACGtacg/;
        $qual_string = reverse $qual_string;
    }

    if ($max_bases){ # clip cigar sequence, qual string, and adjust left end
        if ($comp eq 'U'){ # clip from end
            $sequence =~ s:^(.{$max_bases}).*$:$1:;
            $qual_string =~ s:^(.{$max_bases}).*$:$1:;
        } else {
            my $orig_length = length $sequence; 
            my $query_remaining = $align->query_remaining();
            my $number_clipped = ($orig_length - $query_remaining - $max_bases > 0) ? $orig_length - $query_remaining - $max_bases : 0;
            $sequence =~ s:^.*(.{$max_bases})$:$1:;
            $qual_string =~ s:^.*(.{$max_bases})$:$1:;
            $hit_left_end += $number_clipped;

            # need to adjust for insertions/deletions from reference in clipped portion:
            my $unclipped_cigar = $align->cigar_string(-hard_clip => $hard_clip);
            my $clipped_cigar_length = length $cigar_string;
            $unclipped_cigar =~ s:^(.*)(.{$clipped_cigar_length})$:$1:; # take only the cigar for the portion clipped off

            my ($clipped_insertions, $clipped_deletions) = (0,0);
            while (my $new_op = ($unclipped_cigar =~ s:^(\d+[MIDNSHP])::) ? $1 : undef){
                my ($no_bases, $op) = ($new_op =~ /^(\d+)([MIDNSHP])$/) ? ($1, $2) : (0,0);
                $clipped_insertions += $no_bases if ($op eq 'I');
                $clipped_deletions += $no_bases if ($op eq 'D');
            }

            $hit_left_end += $clipped_deletions - $clipped_insertions;
        }
    }

    return qq!$query\t$flag\t$reference\t$hit_left_end\t$score\t$cigar_string\t*\t0\t0\t$sequence\t$qual_string\n!;
}

sub prepare_reference_sequence_index {
    my $class = shift;

    $class->debug_message("Crossmatch doesn't need any index made, doing nothing.");

}
