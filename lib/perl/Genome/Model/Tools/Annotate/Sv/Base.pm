package Genome::Model::Tools::Annotate::Sv::Base;

use strict;
use warnings;
use Genome;
use File::Basename;
use Sort::Naturally;
use List::Util qw(min max);


class Genome::Model::Tools::Annotate::Sv::Base{
    is => "Command::V2",
    has => [
        input_file => {
            is => 'String',
        },
        output_file => {
            is => 'String',
        },
        annotation_build_id => {
            is => 'String',
        },
        annotation_build => {
            is => 'Genome::Model::Build::ImportedAnnotation',
            id_by => 'annotation_build_id',
        },
        flanking_distance => {
            is => 'Integer',
            is_input => 1,
            default => 50000,
            doc => 'Distance, in bp, for genes to be considered "flanking"',
        },
        chrA_column => {
            is => 'Integer',
            default => 2,
            doc => "1-based index of column that contains the chromosome of breakpoint A",
        },
        bpA_column => {
            is => 'Integer',
            default => 3,
            doc => "1-based index of column that contains the position of breakpoint A",
        },
        chrB_column => {
            is => 'Integer',
            default => 5,
            doc => "1-based index of column that contains the chromosome of breakpoint B",
        },
        bpB_column => {
            is => 'Integer',
            default => 6,
            doc => "1-based index of column that contains the position of breakpoint B",
        },
        event_type_column => {
            is => 'Integer',
            default => 8,
            doc => '1-based index of column that contains the event type',
        },
        orient_column => {
            is => 'Integer',
            default => 9,
            doc => '1-based index of column that contains the orientation',
        },
        score_column => {
            is => 'Integer',
            default => 13,
            doc => '1-based index of column that contains the assembly score',
        },

    ],
};

sub execute {
    my $self = shift;
    my $output_file    = $self->output_file;

    my $out = Genome::Sys->open_file_for_writing($output_file) 
        or die "failed to open $output_file for writing\n";

    my $infile = $self->input_file;
    my $build  = $self->annotation_build
        or die "failed to get annotation build ".$self->annotation_build_id;

    my $in = Genome::Sys->open_file_for_reading($infile) or die "failed to open $infile for reading\n";

    my $breakpoints_list;

    while (my $line = $in->getline) {
        chomp $line;
        next if $line =~ /^\s*$/ || $line =~ /^#/;
        last if $line =~ /-------------------/; #to accomodate some of our squaredancer files with LSF output in them

        #parse either sd or bd file
        my @fields = split /\t/, $line;

        my $chrA = $fields[$self->chrA_column -1];
        my $bpA  = $fields[$self->bpA_column -1];
        my $chrB = $fields[$self->chrB_column -1];
        my $bpB  = $fields[$self->bpB_column -1];
        my $event  = $fields[$self->event_type_column -1];
        my $orient = $fields[$self->orient_column -1 ];
        my $score  = $fields[$self->score_column -1 ];

        #FusionTranscript is picky on orientation
        if ($orient =~ /,/) {
            unless ($score =~ /,/) {
                die $self->error_message("Orientation does not match score in line: $line");
            }
            $orient = $self->pick_orient($orient, $score);
        }

        $breakpoints_list = $self->add_breakpoints_to_chromosome($line, $chrA, $bpA, $chrB, $bpB, $event, $orient, $breakpoints_list);
    }
    
    $in->close;

    $breakpoints_list = $self->fill_in_transcripts($breakpoints_list, $build);
    
    my @transcripts_to_cache;
    for my $chr (nsort keys %{$breakpoints_list}) {
        my @transcripts;
        for my $item (@{$breakpoints_list->{$chr}}) {
            if ($item->{breakpoint_link}) {
                next;
            }
            #push @transcripts, map{$_->id} @{$item->{transcripts_between_breakpoints}};
            if ($item->{chrA} eq $chr) {
                push @transcripts, map{$_->id} @{$item->{transcripts_crossing_breakpoint_a}};
                #push @transcripts, map{$_->id} @{$item->{transcripts_flanking_breakpoint_a}};
            }
            if ($item->{chrB} eq $chr) {
                push @transcripts, map{$_->id} @{$item->{transcripts_crossing_breakpoint_b}};
                #push @transcripts, map{$_->id} @{$item->{transcripts_flanking_breakpoint_b}};
            }
        }
        my @cached = Genome::TranscriptStructure->get(
            chrom_name => $chr,
            data_directory => $build->data_directory."/annotation_data",
            transcript_id => \@transcripts,
        );
    }

    my @column_names = qw(chrA bpA chrB bpB event);

    my $content = $self->process_breakpoint_list($breakpoints_list);
        #my @analysis_column_names = $instance->column_names;
    my @analysis_column_names = $self->column_names;
    unless (@analysis_column_names) {
        die $self->error_message("There is no valid column names for ".__PACKAGE__);
    }
    @column_names = (@column_names, @analysis_column_names);

    $out->print(join "\t", @column_names);
    $out->print("\n");
    
    for my $chr (nsort keys %{$breakpoints_list}) {
        for my $item (@{$breakpoints_list->{$chr}}) {
            if ($item->{breakpoint_link}) {
                next;
            }
            my $key = $self->get_key_from_item($item);
            $out->print(join("\t", $item->{chrA}, $item->{bpA}, $item->{chrB}, $item->{bpB}, $item->{event},  @{$content->{$key}})."\n"); 
        }
    }
    $out->close;
    return 1;

}

#assume only two items in each array for now
sub pick_orient {
    my ($self, $orient, $score) = @_;

    my @orients = split /,/, $orient;
    my @scores  = split /,/, $score;

    if (@orients > 2 or @scores > 2) {
        die $self->error_message("Assume only two items in the array for now");
    }

    #How to handle conflict orietnation like +-,-+ ?
    if ($scores[1] > $scores[0]) {
        return $orients[1];
    }
    else {
        return $orients[0];
    }
}


#TODO NEEDS to be rewritten - I don't think it is right.
#First of all, it stops at one end position rather than getting all intervals that cross the breakpoint.
#Second of all, it only considers the 2nd breakpoint of the SV, not the first breakpoint.  I think
#annotations crossing either breakpoint should probably be considered.
sub annotate_interval_matches {
    #both breakpoints need to match within some wiggle room
    my $self = shift;
    my $positions = shift;
    my $annotation = shift;
    my $annot_length = shift;
    my $tag = shift;
    my $breakpoint_key = shift;

    for my $chr (keys %$positions) {
        my @sorted_items = sort {$a->{$breakpoint_key}<=>$b->{$breakpoint_key}} (@{$positions->{$chr}});
        my @sorted_positions = map{$_->{$breakpoint_key}} @sorted_items;
        my %annotated_output;

        my @chromEnds = sort {$a<=>$b} keys %{$annotation->{$chr}};
        for my $pos (@sorted_positions) {
            while (@chromEnds>0 && $pos>$chromEnds[0]+$annot_length) {
                shift @chromEnds;
            }
            next unless @chromEnds>0;
            for my $start (keys %{$$annotation{$chr}{$chromEnds[0]}}) {
                if ($pos>=$start-$annot_length) {
                    for my $var (@{$$annotation{$chr}{$chromEnds[0]}{$start}}){
                        for my $position_item (@{$positions->{$chr}}) {
                            if ($position_item->{$breakpoint_key} eq $pos) {
                                push @{$position_item->{$tag}->{$breakpoint_key}}, $var;
                            }
                        }
                    }
                }
            }
        }
    }
    return 1;
}

sub get_var_annotation {
    my ($self, $item, $annotation_ref) = @_;
    
    my $varreport = "N/A";
    my @vars;
    my $frac = $self->overlap_fraction;
    
    if (defined $annotation_ref) {
        for my $var (@$annotation_ref) {
            my $pos1 = min($item->{bpB}, $var->{chromEnd});
            my $pos2 = max($item->{bpA}, $var->{chromStart});
            my $overlap = $pos1-$pos2+1;
            my $ratio1 = $overlap/(abs($item->{bpB}-$item->{bpA})+1);
            my $ratio2 = $overlap/(abs($var->{chromEnd}-$var->{chromStart})+1);
            if ($ratio1 >= $frac && $ratio2 >= $frac ) {
                push @vars, $var;
            }
        }

        if (@vars) {
            $varreport = join(",",map{$_->{name}} @vars);
        }
    }
    return $varreport;
}

sub read_ucsc_annotation{
    my ($self, $file) = @_;
    my %annotation;
    my $in = Genome::Sys->open_file_for_reading($file) || die "Unable to open annotation: $file\n";
    while (<$in>) {
        chomp;
        next if /^\#/;
        my $p;
        my @extra;
        ($p->{bin},$p->{chrom},$p->{chromStart},$p->{chromEnd},$p->{name},@extra) = split /\t+/;
        $p->{chrom} =~ s/chr//;
        $p->{extra} = \@extra;
        push @{$annotation{$p->{chrom}}{$p->{chromEnd}}{$p->{chromStart}}}, $p;
    }
    $in->close;
    return \%annotation;
}

sub is_between_breakpoints {
    # Used to get genes that are flanked by deletion breakpoints
    # The entire transcript has to be within the two breakpoints
    # Annotation of the individual breakpoints will give the transcripts interrupted by breakpoints
    my ($self, $start, $stop, $transcript) = @_;
    if ( $start > $stop ) { ($start, $stop) = ($stop, $start); }


    if ( $transcript->transcript_start() >= $start and $transcript->transcript_start() <= $stop and
            $transcript->transcript_stop()  >= $start and $transcript->transcript_stop() <= $stop ) {
            return 1;
    }
    return 0;
}

sub fill_in_transcripts {
    my $self = shift;
    my $breakpoints_list = shift;
    my $build = shift;
    for my $chr (keys %$breakpoints_list) {
        my $chr_breakpoint_list = $breakpoints_list->{$chr};
        my $transcript_iterator = $build->transcript_iterator(chrom_name => $chr);
        die "transcript iterator not defined for chr $chr" unless ($transcript_iterator);
        while (my $transcript = $transcript_iterator->next) {
            for my $item (@$chr_breakpoint_list) {
                if ($item->{event} eq "DEL" and $self->is_between_breakpoints($item->{bpA}, $item->{bpB}, $transcript)) {
                    push (@{$item->{transcripts_between_breakpoints}}, $transcript);
                }
                if ($item->{chrA} eq $chr){
                    if ($transcript->within_transcript($item->{bpA})) {
                        push (@{$item->{transcripts_crossing_breakpoint_a}}, $transcript);
                    }
                    elsif ($transcript->within_transcript_with_flanks($item->{bpA}, $self->flanking_distance)) {
                        push (@{$item->{transcripts_flanking_breakpoint_a}}, $transcript);
                    }
                }
                if ($item->{chrB} eq $chr) {
                    if ($transcript->within_transcript($item->{bpB})) {
                        if (defined $item->{breakpoint_link}) {
                            my $hash = $item->{breakpoint_link};
                            push (@{$hash->{transcripts_crossing_breakpoint_b}}, $transcript);
                        }
                        else {
                            push (@{$item->{transcripts_crossing_breakpoint_b}}, $transcript);
                        }
                    }
                    elsif ($transcript->within_transcript_with_flanks($item->{bpB}, $self->flanking_distance)) {
                        if (defined $item->{breakpoint_link}) {
                            my $hash = $item->{breakpoint_link};
                            push (@{$hash->{transcripts_flanking_breakpoint_b}}, $transcript);
                        }
                        else {
                            push (@{$item->{transcripts_flanking_breakpoint_b}}, $transcript);
                        }
                    }
                }
            }
        }
    }
    return $breakpoints_list;
}

sub add_breakpoints_to_chromosome {
    my ($self, $line, $chrA, $bpA, $chrB, $bpB, $event, $orient, $breakpoints_list) = @_;

    for my $var ($chrA,$bpA,$chrB,$bpB,$event,$orient) {
        unless (defined $var) { die "DID not define necessary variables for call:\n$line\n"; }
    }
    my %hash = (
        chrA  => $chrA,  bpA    => $bpA, 
        chrB  => $chrB,  bpB    => $bpB, 
        event => $event, orient => $orient,
    );
    push (@{$breakpoints_list->{$chrA}}, \%hash);

    unless ($chrA eq $chrB) {
        push (@{$breakpoints_list->{$chrB}}, {%hash, breakpoint_link => \%hash});
    }
    return $breakpoints_list;
}

sub process_breakpoint_list {
    #override in subclass;
    #interface for sv annotators
};

sub column_names {
    #override in subclass;
    #interface for sv annotator
};

sub get_key_from_item {
    my $class = shift;
    my $item = shift;
    return join("--", $item->{chrA}, $item->{bpA}, $item->{chrB}, $item->{bpB}, $item->{event});
}


1;
