package Genome::Model::Tools::Annotate::Sv;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Annotate::Sv {
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
        annotators_to_run => {
            is => 'String',
            is_many => 1,
            default => ['Transcripts'],
        },
    ],
};

sub execute {
    my $self = shift;

    my $out = Genome::Sys->open_file_for_writing($self->output_file);

    my @column_names = qw(chrA bpA chrB bpB event size score);
    for my $analysis_name ($self->annotators_to_run) {
        my $class_name = 'Genome::Model::Tools::Annotate::Sv::'.$analysis_name;
        my @analyais_column_names = $class_name->column_names;
        unless (@analyais_column_names) {
            die $self->error_message("There is no valid column names for $class_name");
        }
        push @column_names, @analyais_column_names;
    }

    $out->print(join "\t", @column_names);
    $out->print("\n");

    my $infile = $self->input_file;
    my $build  = $self->annotation_build;

    my $in = Genome::Sys->open_file_for_reading($infile);

    my $breakpoints_list;

    while (my $line = $in->getline) {
        chomp $line;
        next if $line =~ /^\s*$/ || $line =~ /^#/;
        last if $line =~ /-------------------/; #to accomodate some of our squaredancer files with LSF output in them

        #parse either sd or bd file
        my ($id,$chrA,$bpA,undef,$chrB,$bpB,undef,$event,$orient,undef,$size,$samples,$score) = split /\t/,$line;
        $orient = $1 if $line =~ /([\+\-]*)\,Ins\:/;
        next if ($samples =~ /normal/);

        $breakpoints_list = $self->add_breakpoints_to_chromosome($line, $chrA, $bpA, $chrB, $bpB, $event, $orient, $size, $score, $breakpoints_list);
    }
    
    $in->close;

    $breakpoints_list = $self->fill_in_transcripts($breakpoints_list, $build);
    my %annotation_content;

    foreach my $analysis ($self->annotators_to_run) {
        my $class_name = "Genome::Model::Tools::Annotate::Sv::$analysis";
        my $content = $class_name->process_breakpoint_list($breakpoints_list);
        $annotation_content{$analysis} = $content;
    }
    
    foreach my $chr (sort {$a<=>$b} keys %{$breakpoints_list}) {
        foreach my $item (@{$breakpoints_list->{$chr}}) {
            my @all_content;
            if ($item->{breakpoint_link}) {
                next;
            }
            foreach my $analysis (keys %annotation_content) {
                my $content = $annotation_content{$analysis}->{join("--", $item->{chrA}, $item->{bpA}, $item->{chrB}, $item->{bpB}, $item->{event})};
                @all_content = (@all_content, @$content);
            }
            $out->print(join("\t", $item->{chrA}, $item->{bpA}, $item->{chrB}, $item->{bpB}, $item->{event}, $item->{size}, $item->{score}, @all_content)."\n"); 
        }
    }
    $out->close;
    return 1;
}

sub crosses_breakpoint {
    # Return all transcripts spanning the given position
    my ( $self, $transcript, $position ) = @_;

    if ( $position >= $transcript->transcript_start() and $position <= $transcript->transcript_stop() ) {
        return 1;
    }

    return 0;
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
    foreach my $chr (keys %$breakpoints_list) {
        my $chr_breakpoint_list = $breakpoints_list->{$chr};
        my $transcript_iterator = $build->transcript_iterator(chrom_name => $chr);
        die "transcript iterator not defined for chr $chr" unless ($transcript_iterator);
        while (my $transcript = $transcript_iterator->next) {
            foreach my $item (@$chr_breakpoint_list) {
                if ($item->{event} eq "DEL" and $self->is_between_breakpoints($item->{bpA}, $item->{bpB}, $transcript)) {
                    push (@{$item->{transcripts_between_breakpoints}}, $transcript);
                }
                if ($item->{chrA} eq $chr and $self->crosses_breakpoint($transcript, $item->{bpA})) {
                    push (@{$item->{transcripts_crossing_breakpoint_a}}, $transcript);
                }
                if ($item->{chrB} eq $chr and $self->crosses_breakpoint($transcript, $item->{bpB})) {
                    if (defined $item->{breakpoint_link}) {
                        $DB::single=1;
                        my $hash = $item->{breakpoint_link};
                        push (@{$hash->{transcripts_crossing_breakpoint_b}}, $transcript);
                    }
                    else {
                        push (@{$item->{transcripts_crossing_breakpoint_b}}, $transcript);
                    }
                }
            }
        }
    }
    return $breakpoints_list;
}

sub add_breakpoints_to_chromosome {
    my ($self, $line, $chrA, $bpA, $chrB, $bpB, $event, $orient, $size, $score, $breakpoints_list) = @_;

    for my $var ($chrA,$bpA,$chrB,$bpB,$event,$orient,$size,$score) {
        unless (defined $var) { die "DID not define necessary variables for call:\n$line\n"; }
    }
    my %hash = (
        chrA  => $chrA,  bpA    => $bpA, 
        chrB  => $chrB,  bpB    => $bpB, 
        event => $event, orient => $orient,
        size  => $size,  score  => $score
    );
    push (@{$breakpoints_list->{$chrA}}, \%hash);

    unless ($chrA eq $chrB) {
        push (@{$breakpoints_list->{$chrB}}, {%hash, breakpoint_link => \%hash});
    }
    return $breakpoints_list;
}

1;

