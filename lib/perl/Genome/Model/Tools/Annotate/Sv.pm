package Genome::Model::Tools::Annotate::Sv;

use strict;
use warnings;

use Genome;
use File::Basename;

my @has_param;

BEGIN{
    my @annotators = (
        'Genome::Model::Tools::Annotate::Sv::Transcripts',
        'Genome::Model::Tools::Annotate::Sv::FusionTranscripts',
        'Genome::Model::Tools::Annotate::Sv::Dbsnp',
        'Genome::Model::Tools::Annotate::Sv::Segdup',
        'Genome::Model::Tools::Annotate::Sv::RepeatMasker',
    );
    foreach my $module (@annotators) {
        my $module_meta = UR::Object::Type->get($module);
        my @module_path = split /::/, $module;
        my $module_short_name = Genome::Utility::Text::camel_case_to_string($module_path[-1]);
        $module_short_name =~ s/\s+/_/g;
        my @p = $module_meta->properties;
        foreach my $p (@p) {
            if ($p->can("is_input") and $p->is_input) {
                my $name = $p->property_name;
                $name = $module_short_name."_".$name;
                my %data = %{$p};
                for my $key (keys %data) {
                    delete $data{$key} if $key =~ /^_/;
                }
                delete $data{id};
                delete $data{db_committed};
                delete $data{class_name};
                delete $data{is_input};
                $data{is_optional} = 1;
                $data{is_param} = 1;
                $data{doc} .= " -- for the $module_short_name annotator";
                $data{property_name} = $name;
                push @has_param, $name, \%data;
            }
        }
    }
}

class Genome::Model::Tools::Annotate::Sv {
    is => 'Command::V2',
    has_param => \@has_param,
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
        annotator_list  => {
            is => 'String',
            is_many => 1,
            default => ['Transcripts', 'FusionTranscripts', 'Dbsnp'],
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
    ],
};

sub execute {
    my $self = shift;
    my $output_file    = $self->output_file;
    my @annotator_list = $self->annotator_list;

    my $out = Genome::Sys->open_file_for_writing($output_file) 
        or die "failed to open $output_file for writing\n";

    my @column_names = qw(chrA bpA chrB bpB event);
    for my $analysis_name (@annotator_list) {
        my $class_name = 'Genome::Model::Tools::Annotate::Sv::'.$analysis_name;
        my @analysis_column_names = $class_name->column_names;
        unless (@analysis_column_names) {
            die $self->error_message("There is no valid column names for $class_name");
        }
        @column_names = (@column_names, @analysis_column_names);
    }

    $out->print(join "\t", @column_names);
    $out->print("\n");

    my $infile = $self->input_file;
    my $build  = $self->annotation_build;

    my $in = Genome::Sys->open_file_for_reading($infile) or die "failed to open $infile for reading\n";

    my $breakpoints_list;

    while (my $line = $in->getline) {
        chomp $line;
        next if $line =~ /^\s*$/ || $line =~ /^#/;
        last if $line =~ /-------------------/; #to accomodate some of our squaredancer files with LSF output in them

        #parse either sd or bd file
        my @fields = split /\t/, $line;
        my $chrA = $fields[$self->chrA_column -1];
        my $bpA = $fields[$self->bpA_column -1];
        my $chrB = $fields[$self->chrB_column -1];
        my $bpB = $fields[$self->bpB_column -1];
        my $event = $fields[$self->event_type_column -1];
        my $orient = $fields[$self->orient_column -1 ];

        $breakpoints_list = $self->add_breakpoints_to_chromosome($line, $chrA, $bpA, $chrB, $bpB, $event, $orient, $breakpoints_list);
    }
    
    $in->close;

    $breakpoints_list = $self->fill_in_transcripts($breakpoints_list, $build);
    my %annotation_content;

    for my $type (@annotator_list) {
        my $class_name = "Genome::Model::Tools::Annotate::Sv::$type";
        my $module_meta = UR::Object::Type->get($class_name);
        my $module_short_name = Genome::Utility::Text::camel_case_to_string($type);
        $module_short_name =~ s/\s+/_/g;
        my @p = $module_meta->properties;
        my %params;
        foreach my $property (@p) {
            if ($property->can("is_input") and $property->is_input) {
                my $property_name = $module_short_name."_".$property->property_name;
                if ($self->$property_name) {
                    $params{$property->property_name} = $self->$property_name;
                }
            }
        }
        my $instance = $class_name->create(%params);
        my $content = $instance->process_breakpoint_list($breakpoints_list);
        $annotation_content{$type} = $content;
    }
   
    for my $chr (sort {$a<=>$b} keys %{$breakpoints_list}) {
        for my $item (@{$breakpoints_list->{$chr}}) {
            my @all_content;
            if ($item->{breakpoint_link}) {
                next;
            }
            my $key = Genome::Model::Tools::Annotate::Sv::Base->get_key_from_item($item);
            for my $analysis (@annotator_list) {
                my $content = $annotation_content{$analysis}->{$key};
                @all_content = (@all_content, @$content);
            }
            $out->print(join("\t", $item->{chrA}, $item->{bpA}, $item->{chrB}, $item->{bpB}, $item->{event},  @all_content)."\n"); 
        }
    }
    $out->close;
    return 1;
}

sub crosses_breakpoint {
    # Return all transcripts spanning the given position
    my ( $self, $transcript, $position ) = @_;

    return $transcript->within_transcript($position);
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
                if ($item->{chrA} eq $chr and $self->crosses_breakpoint($transcript, $item->{bpA})) {
                    push (@{$item->{transcripts_crossing_breakpoint_a}}, $transcript);
                }
                if ($item->{chrB} eq $chr and $self->crosses_breakpoint($transcript, $item->{bpB})) {
                    if (defined $item->{breakpoint_link}) {
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

1;

