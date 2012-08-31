package Genome::Individual::Command::CpuDiskSummary;

use strict;
use warnings;

use Genome;

class Genome::Individual::Command::CpuDiskSummary {
    is => 'Command::V2',
    has => [
        patients    => { is => 'Genome::Individual', 
                        shell_args_position => 1,
                        require_user_verify => 0, 
                        is_input => 1,
                        is_many => 1,
                        doc => 'patients on which to report', 
                    },
        output_file => { is => 'FilePath',
                        is_optional => 1,
                        default_value => '-',
                        doc => 'write to this file instead of STDOUT',
                    },
        format      => { is => 'Text',
                        valid_values => ['table','raw'],
                        default_value => 'table',
                        doc => 'set to raw for a raw dump of the data'
                    },
        builds      => {
                        is => 'Genome::Model::Build',
                        is_optional => 1,
                        is_many => 1,
                        require_user_verify => 0,
                        doc => 'limit report to these somatic builds'
                    },
        after       => {
                        is => 'Date',
                        is_optional => 1,
                        doc => 'ignore patients with builds which started before this date (omit obsolete processing)'
                    },
    ],
    doc => 'summarize resources used for a particular patient'
};

sub sub_command_sort_position { 2 }

sub execute {
    my $self = shift;

    #

    my $after = $self->after;

    my @patients = $self->patients;
    my @patient_ids = map { $_->id } @patients; 
    
    my @samples = Genome::Sample->get(source_id => \@patient_ids);
    my @sample_ids = map { $_->id } @samples;

    my @all_models = Genome::Model->get(subject_id => [@patient_ids, @sample_ids]);
    my @all_builds = Genome::Model::Build->get(model_id => [ map { $_->id } @all_models ], -hints => ['disk_allocation','newest_workflow_instance']);
    my @all_events = Genome::Model::Event->get(build_id => [ map { $_->id } @all_builds ]);
    my @all_sr = Genome::SoftwareResult->get("users.user_id" => [ map { $_->id } @all_builds ]);
  
    my @da = Genome::Disk::Allocation->get(
      owner_id => [
        map { $_->id } (@all_sr, @all_events, @all_builds)
      ]
    );

    if (1) {
        $DB::single = 1;
        my @all_wf_names = map { $_->workflow_name } @all_builds;
        my @all_workflows = Workflow::Model::Instance->get(name => \@all_wf_names, -hint => ['cache_workflow']);
        my @children = Workflow::Operation::Instance->get(parent_instance_id => [map { $_->id } @all_workflows]);
        my @execs = Workflow::Operation::InstanceExecution->get(instance_id => [ map { $_->id } @children ]);
        while (@execs) {
            my @child_instances = Workflow::Operation::Instance->get(parent_execution_id => [ map { $_->id } @execs ]);
            @execs = Workflow::Operation::InstanceExecution->get(instance_id => [ map { $_->id } @child_instances ]);
        }
        
        my @force_build_ids;
        if (my @builds = $self->builds) {
            @force_build_ids = map { $_->id } @builds;
        }
    }

    my @row_names = (
        'Individual',
        'Tumor Reference alignment Gbases', 
        'Normal Reference alignment Gbases',
        'Tumor Reference alignment (GB)', 
        'Normal Reference alignment (GB)',
        'Tumor Reference alignment (slot*hr)', 
        'Normal Reference alignment (slot*hr)',
        'Tumor Somatic variation (GB)', 
        'Tumor Somatic variation (slot*hr)',
        '---',
    );
    my %expected_row_names;
    @expected_row_names{@row_names} = @row_names;
   
    my @columns; 
    for my $patient (sort { $a->common_name cmp $b->common_name } @patients) {
        # one column per individual, with hash of row/key => value pairs
        my %f;
        push @columns, \%f;

        $f{"Individual"} = $patient->common_name || $patient->name;

        my @patient_models = Genome::Model->get(
            subject_id => [$patient->id, map { $_->id } $patient->samples],     # old/new
        );

        my @patient_builds = Genome::Model::Build->get(
            model_id => [ map { $_->id } @patient_models ]
        );

        my $first_build_date = UR::Context->now;
        for my $build (@patient_builds) {
            if ($build->date_scheduled lt $first_build_date) {
                $first_build_date = $build->date_scheduled;
            }
        }

        if ($after and $first_build_date lt $after) {
            warn "skipping individual " . $patient->__display_name__ . " because its first build date of $first_build_date is before $after";
            next;
        }

        my %models_by_type;
        for my $model (@patient_models) {
            my $list = $models_by_type{$model->type_name} ||= [];
            push @$list, $model;
        }

        # use the somatic models to identify whether a sample is tumor/normal/unknown
        my %tumor_normal_status_for_sample_id;
        my @somatic_models = (
            map { @$_ }
            grep { defined $_ }
            map { $models_by_type{$_} }
            ('somatic','somatic variation')
        );
        if (not @somatic_models) {
            $self->warning_message("no somatic models for " . $patient->__display_name__ . "\n");
        }
        for my $somatic_model (@somatic_models) {
            my $tumor_id = $somatic_model->tumor_model->subject_id;
            my $normal_id = $somatic_model->normal_model->subject_id;
            $tumor_normal_status_for_sample_id{$tumor_id} = 'Tumor';
            $tumor_normal_status_for_sample_id{$normal_id} ||= 'Normal';  # normal is lower precedence than tumor
        }

        my %seen_builds;
        my %seen_dirs;
        my %subject_gbases;

        #warn "**** LOOPING OVER BUILDS ***********************************************\n";

        for my $build (@patient_builds) {
            next if $seen_builds{$build->id};
            $seen_builds{$build->id} = $build;

            $self->status_message(
                "patient " . $patient->__display_name__ . ' using build ' . $build->__display_name__ . "\n"
            );

            my $model = $build->model;
            my $subject = $model->subject;
            my $type = ucfirst($model->type_name);
            
            my $tn_status = $tumor_normal_status_for_sample_id{$subject->id};
            if ($tn_status) {
                $type = $tn_status . ' ' . $type;
            }

            $f{"$type build"} .= $build->__display_name__ . " ";
            $f{"$type (slot*hr)"} += (eval { $build->cpu_slot_hours } || 0);
            
            my $gb_regular;
            if (my $alloc1 = $build->disk_allocation) {
                $gb_regular = ($alloc1->kilobytes_requested / (1024**2));
            }
            else {
                $gb_regular = 0;
                warn "no allocation for build " . $build->__display_name__;
            }

            $f{"$type GB regular"} += $gb_regular; 
            $f{"$type (GB)"} += $gb_regular;
            
            my @sr = Genome::SoftwareResult->get("users.user_id" => $build->id);
            for my $sr (@sr) {
                my $sr_class = $sr->class;
                my $alloc = $sr->disk_allocation;
                unless ($alloc) {
                    warn "missing disk allocation for software result $sr " . $sr->id;
                    next;
                }
                my $path = $alloc->absolute_path;
                if ($seen_dirs{$path}) {
                    next;
                }
                else {
                    $seen_dirs{$path} = 1;
                }
                my $size = $alloc->kilobytes_requested / (1024**2);
                $f{"$type GB $sr_class"} += $size;
                $f{"$type (GB)"} += $size;
            }
            
            my $gbases = eval { $build->metric(name => 'instrument data total kb')->value/1_000_000 } || 0;
            $subject_gbases{$subject->id} ||= 0;
            if ($gbases > $subject_gbases{$subject->id}) {
                $subject_gbases{$subject->id} = $gbases;
                $f{"$type Gbases"} = eval { $build->metric(name => 'instrument data total kb')->value/1_000_000 };
            }

            # note the widest value so we can build a table below
            my $max_length = 0;
            for (values %f) {
                my $length = length($_);
                $max_length = $length if $max_length < $length;
            }
            $f{max_length} = $max_length;

            # note any new hash keys which appear in %f so we can append them to the report 
            for my $key (keys %f) {
                unless ($expected_row_names{$key}) {
                    $expected_row_names{$key} = 1;
                    push @row_names, $key;
                }
            }
        
        } #next build

    } # next individual

    my $outfh = Genome::Sys->open_file_for_writing($self->output_file);
    if ($self->format eq 'table') {
        no warnings;
        my $max_row_name_length = 0;
        for (@row_names) {
            my $length = length($_);
            $max_row_name_length = $length if $max_row_name_length < $length;
        }

        for my $row_name (@row_names) {
            #my $value = ' ' x ($max_row_name_length - length($row_name)) . $row_name;

            $outfh->print($row_name,"\t");

            for my $col (@columns) {
                my $value = $col->{$row_name};
                #$value = ' ' x ($col->{max_length} - length($value)) . $value;
                $outfh->print($value);
            }
            continue {
                $outfh->print("\t");
            }
            $outfh->print("\n");
        }
    }
    elsif($self->format eq 'raw') {
        $outfh->print(Data::Dumper::Dumper(\@columns));
    }
    else {
        die "unknown output format " . $self->format . "???";
    }

    my @o = $UR::Context::current->all_objects_loaded('UR::Object');
    my @c = grep { $_->__changes__ } grep { $_->class !~ /^UR::/ and $_->class !~ /Command/  and $_->class eq /Workflow::Operation(::Instance|)$/ } @o;
    if (@c) {
        print "*************** CHANGES *********\n";
        print scalar(@o),"\n";
        print scalar(@c),"\n";
        print join("\n",map { "$_" } @c),"\n";
        print Data::Dumper::Dumper($c[0], [ $c[0]->__changes__ ]);
        IO::File->new(">/tmp/obj.dat")->print(Data::Dumper::Dumper($c[0]));
        warn "exiting because therea are changes which would be committed to the database...\n";
        exit;
    }

    return 1;
}

1;

