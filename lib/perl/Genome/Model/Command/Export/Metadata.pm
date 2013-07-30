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

    # Build the "sanitize" hash, to replace internal names with dummy names.
    # This is outside the repo intentionally.
    # This translates local values to those which are distributable outside of TGI.

    my $sanitize_file = $ENV{GENOME_HOME} . "/export/sanitize.csv";
    unless (-e $sanitize_file) {
        die "Expected external file $sanitize_file to exist to sanitize data.  Disable this if you are sure you can dump data unsanitized.";
    }

    my %sanitize_map = (
        # TGI has two IDs for build 37: merge into one.
        # The one we merge away is only around in TGI for legacy reasons.
        'sata420/info/model_data/2857786885/build102671028' => 'ams1102/info/model_data/2869585698/build106942997',
        'sata420/info/model_data/2857786885/build106942997' => 'ams1102/info/model_data/2869585698/build106942997',
        102671028 => 106942997,
        '/gscmnt/' => '/opt/gms/fs/',
        'WUGC' => 'GMS',
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
        add_to_dump_queue($obj, \%queue, \%exclude, \%sanitize_map);
    }

    # Get the disk groups

    my @group_names = qw/info_alignments info_apipe info_genome_models info_apipe_ref/;
    my @groups = Genome::Disk::Group->get(group_name => \@group_names);
    for my $group (@groups) {
        add_to_dump_queue($group, \%queue, \%exclude, \%sanitize_map);
    }

    # Dump everything queued.

    my %defaults = (
        user_name => "genome",
        run_by => "genome"
    );

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
                    print STDERR "replaced $old with $new\n";
                }
            }
            $out_fh->print($txt,"\n");
        }
    }

    return 1;
}


# TODO: instead of conditional logic, move this method
# to each of the types of data we intend to dump.

my $depth = 0;

sub add_to_dump_queue {
    my $obj = shift;
    my $queue = shift;
    my $exclude = shift;
    my $sanitize_map = shift;

    Carp::confess("bad params!") unless $sanitize_map;
    
    my $id = $obj->id;
    my $final_class = $obj->class;

    my $base_class = $final_class;
    for my $base (qw/
        Genome::Model 
        Genome::Model::Build 
        Genome::Subject 
        Genome::Model::Event 
        Genome::InstrumentData 
        Genome::ProcessingProfile 
        Genome::SoftwareResult/
    ) {
        if ($final_class->isa($base)) {
            $base_class = $base;
            last;
        }
    }

    if ($queue->{$base_class}{$id}) {
        return;
    }

    if ($base_class->isa("UR::Value")) {
        return;
    }

    $depth++;
    
    if ($exclude->{$final_class} or $exclude->{$base_class}) {
        p("skip $final_class " . $obj->__display_name__);
    }
    else {
        p("got $final_class ($base_class) " . $obj->__display_name__);
        $queue->{$base_class}{$id} = $obj;
        
        my @sr_assoc = Genome::SoftwareResult::User->get(user_id => $id);
        for my $sr_assoc (@sr_assoc) {
            add_to_dump_queue($sr_assoc, $queue, $exclude, $sanitize_map);
        }

        my @sr = Genome::SoftwareResult->get(id => [ map { $_->software_result_id } @sr_assoc ]);
        for my $sr (@sr) {
            add_to_dump_queue($sr, $queue, $exclude, $sanitize_map);
        }

        my @disks = grep { $obj->isa($_->owner_class_name) } Genome::Disk::Allocation->get(owner_id => $id);
        if ($exclude->{$base_class}) {
            # whent excluding a disk, ensure its group is still included
            for my $disk (@disks) {
                my $group = $disk->group;
                add_to_dump_queue($group, $queue, $exclude, $sanitize_map);
            }    
        }
        else {
            for my $d (@disks) {
                add_to_dump_queue($d, $queue, $exclude, $sanitize_map)
            }
        }
    }

    if ($obj->isa("Genome::Disk::Allocation")) {
        my $group = $obj->group;
        my $volume = $obj->volume;
        my @assignments = $volume->assignments;
        add_to_dump_queue($group, $queue, $exclude, $sanitize_map);
        add_to_dump_queue($volume, $queue, $exclude, $sanitize_map);
        for my $a (@assignments) {
            add_to_dump_queue($a, $queue, $exclude, $sanitize_map);
        }
    }

    for my $ext (qw/Input Param/) {
        my $related_class = $base_class . "::$ext";
        if (UR::Object::Type->get($related_class)) {
            my $owner_method;
            my $value_method;
            if ($obj->isa("Genome::Model")) {
                $owner_method = "model_id";
                $value_method = "value";
            }
            elsif ($obj->isa("Genome::Model::Build")) {
                $owner_method = "build_id";
                $value_method = "value";
            }   
            elsif ($obj->isa("SoftwareResult")) {
                $owner_method = "software_result_id";
                $value_method = "value_obj";
            }
            else {
                next;
            }
            my @assoc = $related_class->get($owner_method => $obj->id);
            for my $a (@assoc) {
                my $v = $a->$value_method;
                unless ($sanitize_map->{$v->id} and $sanitize_map->{$v->id} == $obj->id) {
                    add_to_dump_queue($a, $queue, $exclude, $sanitize_map) unless $exclude->{$final_class};
                    if ($v->isa("Genome::InstrumentData")) {
                        $DB::single = 1;
                    }
                    add_to_dump_queue($v, $queue, $exclude, $sanitize_map);
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
            add_to_dump_queue($e,$queue, $exclude, $sanitize_map);
        }

        add_to_dump_queue($obj->model, $queue, $exclude, $sanitize_map);
    }

    if ($obj->isa("Genome::Model")) {
        add_to_dump_queue($obj->subject, $queue, $exclude, $sanitize_map);
        add_to_dump_queue($obj->processing_profile, $queue, $exclude, $sanitize_map);
    }

    if ($obj->isa("Genome::ProcessingProfile")) {
        my @pp = Genome::ProcessingProfile::Param->get("processing_profile_id" => $obj->id);
        for my $pp (@pp) {
            add_to_dump_queue($pp, $queue, $exclude, $sanitize_map);
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
            next;
        }
        else {
            die "odd nomenclature field value: $nomenclature_ambiguous_value";
        }
        add_to_dump_queue($n, $queue, $exclude, $sanitize_map) if ($n);
        add_to_dump_queue($f, $queue, $exclude, $sanitize_map) if ($f);
    }

    my $parent;
    if (!$exclude->{$final_class} and ($obj->isa("Genome::InstrumentData") or $obj->isa("Genome::Subject"))) {
        my @a = $obj->attributes;
        for $a (@a) {
            add_to_dump_queue($a, $queue, $exclude, $sanitize_map);
        }

        if ($obj->isa("Genome::InstrumentData")) {
            $parent = $obj->library;
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
            add_to_dump_queue($parent, $queue, $exclude, $sanitize_map);
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
            $dir =~ m|/gscmnt/(.*?)/(.*?)/(.*)$| or die "unusual path!";
            my $disk = $1;
            my $group = $2;
            my $data = $3;

            print("### DISK: /opt/gms/fs $disk $group $data\n");
            #my $alloc = Genome::Disk::Allocation->create(
            #    disk_group_name => 'reads',
            #    reallocation_time => UR::Time->now,
            #    mount_path => "/opt/gms/fs/$disk",
            #    group_subdirectory => $group,
            #    allocation_path => $data,
            #);
            #add_to_dump_queue($alloc,$queue,$exclude,$sanitize_map); 
        }
    }

    $depth--;
}

# a print function to log output

sub p {
    for my $m (@_) {
        print STDERR ((" " x $depth),$m, "\n");
    }
}

1;

