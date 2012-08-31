#!/usr/bin/env genome-perl

use Genome;

#UR::Context->object_cache_size_highwater(200_000);
#UR::Context->object_cache_size_lowwater(20_000);



use warnings;
use strict;

use Getopt::Long;

my ($rebuild, $help, $add_all, $add, $lock_name, $chunk);

GetOptions(
    "rebuild" => \$rebuild,
    "help" => \$help,
    "add-all" => \$add_all,
    "add=s" => \$add,
    "lock=s" => \$lock_name,
    "chunk=s" => \$chunk
);

my @types = qw(
user
processing_profile
work_order
model
model_group
individual
flowcell
sample
population_group
taxon
library
disk_group
disk_volume
imported_instrument_data
);

my @types_to_add;

if ($help) {
    my $type_str = join("\n",@types);
    print "\nUSAGE: gcsearch_load_from_db --add [type1,type2, ...]\n\ntypes:\n$type_str\n\n";
    exit 1;
}

if (!$add) {
    if ($add_all) {
        @types_to_add = @types;
        print "trying to add all types of objects to search index\n\n";
    } else {
        print "\nError: must specify --add [type1,type2, ...] or --add-all\n\n";
        exit 1;
    }
} else {
    @types_to_add = split(/,/, $add);

    for my $t (@types_to_add) {
        if (! grep /$t/, @types) {
            my $type_str = join("\n",@types);
            print "\nError: \'$t\' is not a valid type:\n$type_str\n\n";
            exit 1;
        }
    }
}

#Just clear the cache for each entry instead of building all the search result views for now
Genome::Search->get()->refresh_cache_on_add(1);

my $time;

my $lock_resource = $ENV{GENOME_LOCK_DIR} . '/gcsearch/db_loader';

if (Genome::Config->dev_mode()) {
    $lock_resource .= '_dev';
}

if ($lock_name) {
    $lock_resource .= "_$lock_name";
}

my $lock = Genome::Sys->lock_resource(resource_lock=>$lock_resource, max_try=>0);
unless ($lock) {
    die "could not lock, another instance must be running.";
}

if($rebuild) {
    Genome::Search->clear;
}

my $f = get_functions();

for my $t (@types_to_add) {

    if (!defined($f->{$t})) {
        die "oh noes no function for $t\n";
    }
    $f->{$t}->();
}

Genome::Sys->unlock_resource(resource_lock=>$lock);

exit;





sub status($) {
    my $msg = shift;
    print STDERR $msg;
}

sub loading($) {
    my $type = shift;
    status 'Loading ' . $type . '...';

    $time = time;
}

sub done() {
    my $diff = time - $time;
    status "done. ($diff seconds)\n";
}

sub get_functions {

    my $f = {};

    Genome::InstrumentData::Solexa->get(); #preload

    $f->{'user'} = sub {
        loading "users";
        my @u = Genome::Sys::User->get();
        Genome::Search->add(@u);
        done;
    };

    $f->{'processing_profile'} = sub {
        loading "processing profiles";
        my @processing_profiles = Genome::ProcessingProfile->get(-hint => ['params']); #models would've done this later, but we need "them all" so do it in advance
        Genome::Search->add(@processing_profiles);
        done;
    };

    $f->{'work_order'} = sub {
        loading "work orders";
        my @wo = Genome::WorkOrder->get();
        Genome::Search->add(@wo);
        done;
    };

    $f->{'model_group'} = sub {
        loading "model groups";
        my @model_groups = Genome::ModelGroup->get();
        Genome::Search->add(@model_groups);
        done;
    };

    $f->{'individual'} = sub {
        loading "individuals";
        my @individuals = Genome::Individual->get( common_name => {operator => 'ne', value => undef } );
        print 'adding ' . scalar(@individuals) . ' individuals';
        Genome::Search->add(@individuals);
        done;
    };

    $f->{'model'} = sub {
        loading "models";
        my @models;

        eval {
            @models = Genome::Model->get(type_name =>
            {operator => 'IN',
                value => [  'rna seq',
                            'reference alignment',
                            'somatic',
                            'convergence',
                            'metagenomic composition shotgun',
                            'assembly',
                            'amplicon assembly',
                            'de novo assembly',
                            'genotype microarray' ]
            });

        };

        if($@) {
            #Somebody currently has a model in the DB without having deployed the module.  Fall back on a safe set of model types
            @models = Genome::Model->get(type_name => {operator => 'IN', value => ['reference alignment', 'somatic', 'assembly', 'amplicon assembly', 'de novo assembly', 'genotype microarray', 'virome screen', 'convergence'] });
        }

#        Genome::Model::Build->get(-hint => ['events', 'status']); #Load these all at once rather than letting each model query for its own in turn (but do this after the models are loaded)

        my $num_chunks = 21;
        my $chunk_size = int(scalar(@models) / $num_chunks);
        my ($i, $j);

        die 'need --chunk for models' if !defined($chunk);

        print "objects loaded: " . $UR::Context::all_objects_cache_size . "\n";

        $i = $chunk * $chunk_size;
        warn "\nmodel chunk $chunk (offset: $i length: $chunk_size)\n";
        Genome::Search->add(splice(@models, $i, $chunk_size));

        open(my $fh, ">>/gscuser/jlolofie/tmp/search_transition/models_added." . $ENV{'LSB_JOBID'});
        for my $m (splice(@models, $i, $chunk_size)) {
            print $fh $m->genome_model_id() . "\n";
        }
        close($fh);

        done;
    };

    $f->{'flowcell'} = sub {
        loading "flow cells";
        my @flow_cells = Genome::InstrumentData::FlowCell->get();
        my $num_chunks = 21;
        my $chunk_size = int(scalar(@flow_cells) / $num_chunks);
        my ($i, $j);

        die 'need --chunk for flowcell cuz it doesnt fit in memory ;)' if !defined($chunk);

        $i = $chunk * $chunk_size;
        $j = $i + $chunk_size;
        print "\nflowcell chunk $chunk ($i to $j) size: $chunk_size\n";

        my @fc = splice(@flow_cells, $i, $chunk_size);
        Genome::Search->add(@fc);

        open(my $fh, ">>/gscuser/jlolofie/tmp/search_transition/flowcells_added." . $ENV{'LSB_JOBID'});
        for my $fc (@fc) {
            print $fh $fc->flow_cell_id . "\n";
        }
        close($fh);

        done;
    };

    $f->{'sample'} = sub {
        loading "samples";
        my @samples = Genome::Sample->get( name => {operator => 'ne', value => undef } );
        Genome::Search->add(@samples);
        done;
    };

    $f->{'population_group'} = sub {
        loading "population groups";
        my @population_groups = Genome::PopulationGroup->get();
        Genome::Search->add(@population_groups);
        done;
    };

    $f->{'taxon'} = sub {
        loading "taxons";
        my @taxons = Genome::Taxon->get();
        Genome::Search->add(@taxons);
        done;
    };

    $f->{'library'} = sub {
        $DB::single = 1;
        loading "libraries";
        my @libraries = Genome::Library->get(); #Would hint samples and taxons, but were already loaded above
        Genome::Search->add(@libraries);
        done;
    };

    $f->{'disk_group'} = sub {
        loading "disk groups";
        my @disk_groups = Genome::Disk::Group->get();
        Genome::Search->add(@disk_groups);
        done;
    };

    $f->{'disk_volume'} = sub {
        loading "disk volumes";
        my @disk_volumes = Genome::Disk::Volume->get( -hint => ['assignments']);
        Genome::Search->add(@disk_volumes);
        done;
    };

    $f->{'imported_instrument_data'} = sub {
        loading "imported_instrument_data";
        my @iid = Genome::InstrumentData::Imported->get();
        print 'adding ' . scalar(@iid) . ' imported_instrument_data';
        Genome::Search->add(@iid);
        done;
    };

    return $f;
}




