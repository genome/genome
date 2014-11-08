package Genome::Individual::Command::View;

use strict;
use warnings;

use Genome;
use List::MoreUtils "uniq";
use Genome::Utility::Text qw(justify side_by_side);
use Genome::Utility::List qw(in);

class Genome::Individual::Command::View {
    doc => "Display basic information about an individual.",
    is => ['Genome::Command::Viewer', 'Genome::Command::WithColor'],
    has => [
        individual => {
            is => 'Genome::Individual',
            shell_args_position => 1,
            doc => 'Genome::Individual',
        },
        filter => {
            is => 'Boolean',
            is_optional => 1,
            default_value => 1,
            doc => 'Only show models that have specified processing-profiles',
        },
        processing_profiles => {
            is => 'Genome::ProcessingProfile',
            is_many => 1,
            is_optional => 1,
            doc => 'If --filter-by-processing-profiles, then only show models with these processing-profiles',
        },
    ],
};

sub help_synopsis {
    my $self = shift;
    my $result .= <<EOP;
    Displays basic information about an individual.
EOP
    return $result;
}

sub help_detail {
    my $self = shift;
    my $result .= <<EOP;
Displays information about the individual and the associated samples.
EOP
    return $result;
}

sub write_report {
    my ($self, $width, $handle) = @_;

    $self->_resolve_default_processing_profiles();
    $self->_write_basic_info($width, $handle);
    $self->_write_sample_info($width, $handle);
    $self->_write_model_info($width, $handle);
}

sub _resolve_default_processing_profiles {
    my ($self) = @_;
    return if defined($self->processing_profiles);

    my @default_pp_ids = (
            2635769, # reference-alignment
            2756469, # whole genome somatic variation
            2756470, # exome somatic variation
            2754795, # rna-seq
            2649924, # clin-seq
            );
    my @pps;
    for my $pp_id (@default_pp_ids) {
        push(@pps, Genome::ProcessingProfile->get($pp_id));
    }
    $self->processing_profiles(\@pps);
}

sub _write_basic_info {
    my ($self, $width, $handle) = @_;
    my $ind = $self->individual;

    my $aka = $ind->common_name ?
            sprintf("a.k.a. %s", $ind->common_name) : '';
    printf $handle "\nIndividual: %s (%s) %s\n",
            $self->_color($ind->name, 'bold'),
            $ind->id,
            $self->_color($aka, 'cyan');
    printf $handle "Taxon: %s (%s)\n",
            $self->_color($ind->taxon->name, 'bold'),
            $ind->taxon->id;

    if(defined($ind->father) or defined($ind->mother)) {
        my $father = $ind->father;
        my $father_string = 'undef';
        if(defined($father)) {
            my $father_name = $self->_color($father->name, 'bold');
            my $father_id = $father->id;
            $father_string = sprintf("%s (%s)", $father_name, $father_id);
        }

        my $mother = $ind->mother;
        my $mother_string = 'undef';
        if(defined($mother)) {
            my $mother_name = $self->_color($mother->name, 'bold');
            my $mother_id = $mother->id;
            $mother_string = sprintf("%s (%s)", $mother_name, $mother_id);
        }
        printf $handle "Father: %s   Mother: %s\n",
                    $father_string, $mother_string;
    }

    my ($names, $values) = $self->_get_attributes($ind);
    if($names and $values) {
        print $handle side_by_side(
                [$names, $values],
                fill => ['.', ' '],
                max_width => $width,
        );
    }
}

sub _get_attributes {
    my ($self, $ind) = @_;

    my @attributes_blacklist =
            qw(id mother_id father_id taxon_id common_name);
    my $labels = '';
    my $values = '';
    for my $attr ($ind->attributes) {
        unless(in($attr->attribute_label, @attributes_blacklist)) {
            $labels .= $attr->attribute_label . "\n";
            $values .= sprintf("%s\n",
                    $self->_color($attr->attribute_value, 'bold'));
        }
    }
    return $labels, $values;
}

sub _write_sample_info {
    my ($self, $width, $handle) = @_;

    my $ind = $self->individual;
    my @samples = $ind->samples;

    print $handle "\n\n";
    unless(scalar(@samples)) {
        print $handle $self->_color("No Samples", 'red') . "\n";
        return;
    }

    my $hformat = $self->_color("%6s %-12s %-27s  %-14s %-12s", 'bold') . "\n";
    my $centered_str = justify(" -- Instrument Data --", 'center', 27);
    my $header = sprintf($hformat, " ", " ", $centered_str, "Extraction", "Common");
    $header .= sprintf($hformat, " ", 'Sample ID', 'wgs exome rna other unknown',
            'Type', 'Name');
    printf $handle $header;

    my $i = 0;
    for my $sample (@samples) {
        $i += 1;
        my $line = justify("$i)", 'right', 6);
        $line .= " ";

        $line .= justify($sample->id, 'left', 12);
        $line .= " ";

        my @counts = _get_instrument_data_counts($sample);
        $line .= justify($counts[0], 'right', 3, '.', '');
        $line .= " ";
        $line .= justify($counts[1], 'right', 5, '.', '');
        $line .= " ";
        $line .= justify($counts[2], 'right', 3, '.', '');
        $line .= " ";
        $line .= justify($counts[3], 'right', 5, '.', '');
        $line .= " ";
        $line .= justify($counts[4], 'right', 7, '.', '');

        my $extraction_type = $sample->extraction_type ||
                $self->_color("Unknown", 'red');
        $line .= "  ". justify($extraction_type, 'left', 15);

        my $common_name = $sample->common_name || '';
        if($common_name eq 'tumor') {
            $common_name = $self->_color('tumor ', 'magenta');
        } elsif($common_name eq 'normal') {
            $common_name = $self->_color('normal', 'cyan');
        } elsif($common_name) {
            $common_name = $self->_color($common_name, 'bold');
        }
        $line .= justify($common_name, 'left', 12);
        print $handle "$line\n";
    }
    return;
}

sub _get_instrument_data_counts {
    my ($sample) = @_;

    my @ids = $sample->instrument_data;
    my %instrument_data_index;
    for my $id (@ids) {
        my $category = _determine_instrument_data_type($id);
        if(in($category, keys %instrument_data_index)) {
            $instrument_data_index{$category} += 1;
        } else {
            $instrument_data_index{$category} = 1;
        }
    }
    my $wgs     = $instrument_data_index{wgs}     || '';
    my $exome   = $instrument_data_index{exome}   || '';
    my $rna     = $instrument_data_index{rna}     || '';
    my $other   = $instrument_data_index{other}   || '';
    my $unknown = $instrument_data_index{unknown} || '';

    return $wgs, $exome, $rna, $other, $unknown;
}

sub _determine_instrument_data_type {
    my ($instrument_data) = @_;

    my $sample_type = $instrument_data->sample_type;
    if($sample_type and $sample_type =~ m/rna/) {
        return 'rna';
    }

    if($instrument_data->can('target_region_set_name')) {
        my $trsn = $instrument_data->target_region_set_name;
        return 'wgs' unless($trsn);

        my $fl = Genome::FeatureList->get(name => $trsn);
        return 'exome' if($fl->content_type eq 'exome');
        return 'other' if(not $fl or not $fl->content_type);
        return 'unknown';
    } else {
        return 'unknown';
    }
}

sub _write_model_info {
    my ($self, $width, $handle) = @_;

    my $ind = $self->individual;
    my @samples = $ind->samples;
    return unless scalar(@samples);

    printf $handle "\n=== %s ===\n",
            $self->_color("Model Information", 'bold');
    if($self->filter) {
        my @pp_ids = map {$_->id} $self->processing_profiles;
        printf $handle "* Only showing models with the following ".
                "processing-profiles:\n    %s\n",
                join(", ", @pp_ids);
    }

    my $i = 0;
    for my $sample (@samples) {
        $i += 1;
        my @models = $sample->models;
        next unless scalar(@models);

        my @filtered_models = @models;
        if($self->filter) {
            my @pp_ids = map {$_->id} $self->processing_profiles;
            @filtered_models = grep
                    {in($_->processing_profile->id, \@pp_ids)} @models;
        }

        next unless scalar(@filtered_models);
        printf $handle "For sample %s)\n", $i;

        my $header = sprintf("    %s %s %s %s %s",
                justify("Model ID", 'left', 12),
                justify("Class", 'left', 22),
                justify("PP ID", 'left', 9),
                justify("Status", 'left', 14),
                justify("Latest Build ID", 'left', 12));
        print $handle $self->_color($header, 'bold') . "\n";

        for my $model (@filtered_models) {
            my $model_class = $model->class;
            $model_class =~ s/Genome::Model:://;

            my $status = $model->status;
            if($status eq 'Succeeded') {
                $status = $self->_color($status, 'green');
            } else {
                $status = $self->_color($status, 'red');
            }
            my $latest_build = $model->latest_build;
            my $lb_id = $latest_build ? $latest_build->id : '';

            my $line = sprintf("    %s %s %s %s %s\n",
                    justify($model->id, 'left', 12),
                    justify($model_class, 'left', 22),
                    justify($model->processing_profile->id, 'left', 9),
                    justify($status, 'left', 14),
                    justify($lb_id, 'left', 12));
            print $handle $line;
        }

        my $num_filtered_out = scalar(@models) - scalar(@filtered_models);
        if($num_filtered_out > 0) {
            printf $handle "    + %s models with other processing-profiles.\n",
                    $num_filtered_out;
        }
    }
}

1;
