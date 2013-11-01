package Genome::Model::Command::Export::Metadata;

# TODO: this is a legacy script that embeds a bunch of conditional logic
# These should be broken-out into methods on the objects being exported.

use Genome;
use strict;
use warnings;

class Genome::Model::Command::Export::Metadata {
    is => 'Command::V2',
    has_input => [
        models => { 
            is => 'Genome::Model', 
            is_many => 1,
            shell_args_position => 1, 
            doc => 'models to dump (recurse inputs)', 
        },
        output_path => {
            is => 'FilesystemPath',
            default_value => '-',
            doc => 'the path to which to stream output (default to stdout)',
        },
    ],
    has_param => [
        refdata_only => { 
            is => 'Boolean', 
            doc => 'exclude data related to the sequencing subject, and get only reference data, processing profiles, etc.', 
            is_optional => 1,
        },
        verbose => {
            is => 'Boolean',
            default_value => 0,
            is_optional => 1,
            doc => 'verbose output',
        }
    ],
    has_optional_transient => [
        _terminal_width => { is => 'Number', },
    ],
    doc => 'serialize models for later import into other GMS instances',
};

sub help_synopsis {
    return <<EOS
genome model export metadata id=2891454740 > 2891454740.dat

genome model export metadata id=2891454740 --refdata > 2891454740.refdata.dat

genome model export metadata "name like myproject%" > myproject.dat

# untested
GENOME_SYS_ID=gms100 genome model export data "name like myproject%" | GENOME_SYS_ID=gms200 genome model import data --update -
EOS
}

sub help_detail {
    return <<EOS
Stream models from the system in serialized form, prepared for import into another GMS instance.
EOS
}

sub execute {
    my $self = shift;
    my @models = $self->models;
   
    my $tput_cols = `tput cols`;
    chomp $tput_cols;
    $tput_cols ||= 100;
    $self->_terminal_width($tput_cols);

    my $out_fh = Genome::Sys->open_file_for_writing($self->output_path);
    
    unless (@models) {
        $self->error_message("no models specified!");
        return;
    }

    $DB::single = 1;
    my $refdata_only = $self->refdata_only;

    # When just dumping refdata_only, don't dump subject-centric objects, but still recurse through them
    
    my %exclude;
    if ($refdata_only) {
        %exclude = map { $_ => 1 } qw/
            Genome::Model::GenotypeMicroarray 
            Genome::Model::ReferenceAlignment 
            Genome::Model::SomaticVariation 
            Genome::Model::RnaSeq 
            Genome::Model::ClinSeq
            Genome::InstrumenaData
            Genome::InstrumentDataAttribute 
            Genome::Library 
            Genome::Sample 
            Genome::Individual 
            Genome::PopulationGroup 
        /;
    }
    %exclude = (%exclude, 'Genome::Disk::Group' => 1);
    %exclude = (%exclude, 'Genome::Model::Metric' => 1);

    # Build the "sanitize" hash, to replace internal names with dummy names.
    # This is outside the repo intentionally.
    # This translates local values to those which are distributable outside of TGI.

    my $sanitize_file = $ENV{GENOME_HOME} . "/export/sanitize.csv";
    print "\nExpected sanitize file: $sanitize_file\n\n";
    unless (-e $sanitize_file) {
        die "Expected external file $sanitize_file to exist to sanitize data.  Disable this if you are sure you can dump data unsanitized.";
    }

    my %sanitize_map = (
        # TGI has two IDs for build 37: merge into one.
        # The one we merge away is only around in TGI for legacy reasons.
        'sata420/info/model_data/2857786885/build102671028' => 'ams1102/info/model_data/2869585698/build106942997',
        'sata420/info/model_data/2857786885/build106942997' => 'ams1102/info/model_data/2869585698/build106942997',
        102671028 => 106942997,
        '/gscmnt/' => '/opt/gms/GMS1/fs/',
        'WUGC' => 'GMS1',
    );

    my @rows = IO::File->new($sanitize_file)->getlines;
    chomp @rows;
    for my $row (@rows) {
        next unless $row;
        next if $row =~ /^#/;
        my ($old,$new) = split(/,/,$row);
        $sanitize_map{$old} = $new;
    }

    # Queue each of the specified objects, and recurse through related data.

    my %queue;
    for my $obj (@models) {
        $self->add_to_dump_queue($obj, \%queue, \%exclude, \%sanitize_map);
    }

    # Get the disk groups

    my @group_names = qw/info_alignments info_apipe info_genome_models info_apipe_ref/;
    my @groups = Genome::Disk::Group->get(group_name => \@group_names);
    for my $group (@groups) {
        $self->add_to_dump_queue($group, \%queue, \%exclude, \%sanitize_map);
    }

    # Dump everything queued.

    my %defaults = (
        user_name => "genome",
        run_by => "genome"
    );

    my $depth = 0;
    for my $cls (sort keys %queue) {
        $out_fh->print("# $cls\n"); 
        my $obj_hash = $queue{$cls};
        for my $id (sort keys %$obj_hash) {
            my $obj = $obj_hash->{$id};
            delete $obj->{__get_serial};
            delete $obj->{__load};
            for my $field (keys %defaults) {
                if (exists $obj->{$field}) {
                    $obj->{$field} = $defaults{$field};
                }
                if (exists $obj->{db_committed}{$field}) {
                    $obj->{db_committed}{$field} = $defaults{$field};
                }
            }
            my $txt = UR::Util::d($obj);
            for my $old (keys %sanitize_map) {
                my $new = $sanitize_map{$old};
                if ($txt =~ s/$old/$new/g) {
                    $self->p("replaced $old with $new");
                }
            }
            $out_fh->print($txt,"\n");
        }
    }

    if ($self->verbose) {
        $self->p("export complete");
    }
    else {
        $self->p("");
        $self->p("export complete");
        print STDERR "\n";
    }
    return 1;
}


# TODO: instead of conditional logic, move this method
# to each of the types of data we intend to dump.

my $depth = 0;

sub add_to_dump_queue {
    my $self = shift;
    my $obj = shift;
    my $queue = shift;
    my $exclude = shift;
    my $sanitize_map = shift;

    Carp::confess("bad params!") unless $sanitize_map;
    
    my $id = $obj->id;
    my $final_class = $obj->class;

    my $base_class = $final_class;
    for my $base (
        qw/
            Genome::Model 
            Genome::Model::Build 
            Genome::Subject 
            Genome::Model::Event 
            Genome::InstrumentData 
            Genome::ProcessingProfile 
            Genome::SoftwareResult
        /
    ) {
        if ($final_class->isa($base)) {
            $base_class = $base;
            last;
        }
    }

    if ($queue->{$base_class}{$id}) {
        return;
    }

    if ($base_class->isa("UR::Value") and not $base_class->isa("Genome::Db")) {
        return;
    }

    $depth++;
    
    if ($exclude->{$final_class} or $exclude->{$base_class}){
        $self->p("skip $final_class " . $obj->__display_name__);
    }
    else {
        $self->p("export $final_class ($base_class): " . $obj->__display_name__);
        $queue->{$base_class}{$id} = $obj;
        
        my @sr_assoc = Genome::SoftwareResult::User->get(user_id => $id);
        for my $sr_assoc (@sr_assoc) {
            $self->add_to_dump_queue($sr_assoc, $queue, $exclude, $sanitize_map);
        }

        my @sr = Genome::SoftwareResult->get(id => [ map { $_->software_result_id } @sr_assoc ]);
        for my $sr (@sr) {
            $self->add_to_dump_queue($sr, $queue, $exclude, $sanitize_map);
        }

        my @disks = grep { $obj->isa($_->owner_class_name) } Genome::Disk::Allocation->get(owner_id => $id);
        if ($exclude->{$base_class}) {
            # whent excluding a disk, ensure its group is still included
            for my $disk (@disks) {
                my $group = $disk->group;
                $self->add_to_dump_queue($group, $queue, $exclude, $sanitize_map);
            }    
        }
        else {
            for my $d (@disks) {
                $self->add_to_dump_queue($d, $queue, $exclude, $sanitize_map)
            }
        }
    }

    # TODO: this conditional logic stuff should really be broken out into methods on the class,
    # with a standard API, instead of conditional logic.  The idea is still maturing and is 
    # pretty green, though, so currently it helps to have it all in one place, where relationships
    # can be seen easily.

    if ($obj->isa("Genome::Disk::Allocation")) {
        my $group = $obj->group;
        $self->add_to_dump_queue($group, $queue, $exclude, $sanitize_map) if $group;
        my $volume = $obj->volume;
        if ($volume) {
            $self->add_to_dump_queue($volume, $queue, $exclude, $sanitize_map);
            my @assignments = $volume->assignments;
            for my $a (@assignments) {
                $self->add_to_dump_queue($a, $queue, $exclude, $sanitize_map);
            }
        }
    }

    for my $ext (qw/Input Param Metric/) {
        my $related_class = $base_class . "::$ext";
        if (UR::Object::Type->get($related_class)) {
            my $owner_method;
            my $value_method;
            my $value_method2;
            if ($obj->isa("Genome::Model")) {
                $owner_method = "model_id";
                $value_method = "value";
            }
            elsif ($obj->isa("Genome::Model::Build")) {
                $owner_method = "build_id";
                $value_method = "value";
            }   
            elsif ($obj->isa("Genome::SoftwareResult")) {
                $owner_method = "software_result_id";
                if ($ext eq 'Metric') {
                    $value_method = "metric_value";
                }
                else {
                    $value_method = "value_obj";
                    $value_method2 = "value_id";
                }
            }
            else {
                next;
            }
            my @assoc = $related_class->get($owner_method => $obj->id);
            for my $a (@assoc) {
                my $v = $a->$value_method;
                if (not defined $v and $value_method2) {
                    my $id = $a->$value_method2;
                    die if not defined $id;
                    $v = UR::Value::Text->get($id);
                }
                my $vid = (ref($v) ? $v->id : $v);

                unless ($sanitize_map->{$vid} and $sanitize_map->{$vid} eq $obj->id) {
                    $self->add_to_dump_queue($a, $queue, $exclude, $sanitize_map) unless $exclude->{$final_class};
                    $self->add_to_dump_queue($v, $queue, $exclude, $sanitize_map) if ref $v;
                }
            }
        }
    }

    if ($obj->isa("Genome::Model::Build")) {
        $obj->status("Dummy");
        my $e = $obj->the_master_event;
        die unless $e->event_status eq "Dummy";
        $e->{db_committed}{event_status} = "Dummy";
        my @e = $obj->events();
        for my $e (@e) {
            $self->add_to_dump_queue($e,$queue, $exclude, $sanitize_map);
        }

        $self->add_to_dump_queue($obj->model, $queue, $exclude, $sanitize_map);

        if ($obj->isa("Genome::Model::Build::ReferenceSequence")) {
            my @i = Genome::Model::Build::ReferenceSequence::AlignerIndex->get(reference_build_id => $obj->id, test_name => undef);
            for my $i (@i) {
                my $dir = $i->output_dir;
                next if $dir and $dir =~ /gscarchive/;
                next unless $i->id eq '117803766';        # TODO: make this smarter
                $self->add_to_dump_queue($i, $queue, $exclude, $sanitize_map);
            }
            my @prev_builds = grep { $_->isa("Genome::Model::Build::ReferenceSequence") } values %{ $queue->{"Genome::Model::Build"} };
            if (@prev_builds) {
                $DB::single = 1;
                my @converters1 = map { Genome::Model::Build::ReferenceSequence::Converter->get(source_reference_build => $obj, destination_reference_build => $_) } @prev_builds;
                my @converters2 = map { Genome::Model::Build::ReferenceSequence::Converter->get(destination_reference_build => $obj, source_reference_build => $_) } @prev_builds;
                for my $converter (@converters1, @converters2) {
                    $self->add_to_dump_queue($converter, $queue, $exclude, $sanitize_map);
                }
            }
        }
    }

    if ($obj->isa("Genome::Model")) {
        $self->add_to_dump_queue($obj->subject, $queue, $exclude, $sanitize_map);
        $self->add_to_dump_queue($obj->processing_profile, $queue, $exclude, $sanitize_map);
    }

    if ($obj->isa("Genome::ProcessingProfile")) {
        my @pp = Genome::ProcessingProfile::Param->get("processing_profile_id" => $obj->id);
        for my $pp (@pp) {
            $self->add_to_dump_queue($pp, $queue, $exclude, $sanitize_map);
        }
    }

    if ($obj->isa("Genome::SubjectAttribute")) {
        my $nomenclature_ambiguous_value = $obj->nomenclature;
        my $n = Genome::Nomenclature->get(name => $nomenclature_ambiguous_value);
        my $f = Genome::Nomenclature::Field->get(id => $nomenclature_ambiguous_value);
        if ($n and not $f) {
            # old attributes have a text nomenclature we got by name
            $f = Genome::Nomenclature::Field->get(nomenclature_id => $n->id, name => $obj->attribute_label);
        }
        elsif ($f and not $n) {
            # for new ones, the nomenclature is the fk to the nomenclature_field
            $n = $f->nomenclature;
        }
        elsif (not $n and not $f) {
            # continue
        }
        else {
            die "odd nomenclature field value: $nomenclature_ambiguous_value";
        }
        $self->add_to_dump_queue($n, $queue, $exclude, $sanitize_map) if ($n);
        $self->add_to_dump_queue($f, $queue, $exclude, $sanitize_map) if ($f);
    }

    if ($obj->isa("Genome::FeatureList")) {
        my $reference_build = $obj->reference;
        if ($reference_build) {
            $self->add_to_dump_queue($reference_build, $queue, $exclude, $sanitize_map);
        }
    }

    my $parent;
    if (!$exclude->{$final_class} and ($obj->isa("Genome::InstrumentData") or $obj->isa("Genome::Subject"))) {
        my @a = $obj->attributes;
        for $a (@a) {
            $self->add_to_dump_queue($a, $queue, $exclude, $sanitize_map);
        }

        if ($obj->isa("Genome::InstrumentData")) {
            $parent = $obj->library;
            if (my $target_region_set_name = $obj->target_region_set_name) {
                my @feature_lists = Genome::FeatureList->get(name => $target_region_set_name);
                for my $f (@feature_lists) {
                    $self->add_to_dump_queue($f, $queue, $exclude, $sanitize_map);
                }
            }
        }
        elsif ($obj->isa("Genome::Library")) {
            $parent = $obj->sample;
        }
        elsif ($obj->isa("Genome::Sample")) {
            $parent = $obj->source;
        }
        elsif ($obj->isa("Genome::Individual")) {
            $parent = $obj->taxon;
        }
        elsif ($obj->isa("Genome::PopulationGroup")) {
            $parent = $obj->taxon;
        }
        
        if ($parent) {
            $self->add_to_dump_queue($parent, $queue, $exclude, $sanitize_map);
        }

        my @alloc = 
          grep { $obj->isa($_->owner_class_name) }
          Genome::Disk::Allocation->get(owner_id => $obj->id);

        if (@alloc == 0 and $obj->isa("Genome::InstrumentData")) {
            my $path = $obj->bam_path;
            unless ($path) {
                $path = $obj->archive_path;
            }
            unless ($path) {
                print(Data::Dumper::Dumper($obj));
                die "NO PATH on $obj $obj->{id}\n";
            }

            my $dir = File::Basename::dirname($path);
            my ($genome_sys_id, $disk, $group, $data);
            if ($dir =~ m|/gscmnt/(.+?)/(.+?)/(.+)$|) {
                $genome_sys_id = 'GMS1';
                $disk = $1;
                $group = $2;
                $data = $3;
            }
            elsif ($dir =~ m|/opt/gms/(.+?)/fs/(.+?)/(.+?)/(.+)$|) {
                $genome_sys_id = $1;
                $disk = $2;
                $group = $3;
                $data = $4;
            }
            else {
                die "cannot parse path $dir!!!";
            }

            #print("### DISK: sys:$genome_sys_id disk:$disk group:$group subdir:$data\n");
            my $alloc = UR::Object::create("Genome::Disk::Allocation",
                disk_group_name => 'reads',
                reallocation_time => UR::Time->now,
                mount_path => "/opt/gms/$genome_sys_id/fs/$disk",
                group_subdirectory => $group,
                allocation_path => $data,
                owner_class_name => $obj->class,
                owner_id => $obj->id,
                kilobytes_requested => 1,
            );
            $self->add_to_dump_queue($alloc,$queue,$exclude,$sanitize_map); 
        }
    }

    $depth--;
}

# a print function to log output

my $last_len = 0;
sub p {
    my $self = shift;
    STDERR->autoflush(1);
    my $max_width = $self->_terminal_width;
    for my $m (@_) {
        my $msg;
        if ($self->verbose) {
            $msg = (" " x $depth) . $m . "\n";
        }
        else {
            $msg = $m;
            $msg = substr($msg,0,$max_width) if length($msg) > $max_width;
            my $len = length($msg);
            $msg = "\r$msg";
            my $delta = $last_len - $len;
            if ($delta > 0) {
                $msg .= (' ' x $delta);
            }
            $last_len = $len;
        }
        print STDERR $msg;
    }
}

1;

