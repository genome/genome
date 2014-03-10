package Genome::Model::Convergence::Command::FindSignificantAttributes;

use strict;
use warnings;

use Genome;

class Genome::Model::Convergence::Command::FindSignificantAttributes {
    is => 'Genome::Command::Base',
    has => [
        build => {
            shell_args_position => 1,
            is => 'Genome::Model::Build::Convergence',
            id_by => 'build_id',
            doc => 'the build on which to run the analysis'
        },
        build_id => { is => 'Text', is_input => 1, is_output => 1 },
        data_file => {
            is => 'String',
            doc => 'Location of the ARFF file used as input to Weka.  (If the provided file exists it will be used; if not it will be generated.)',
        },
        output_file => {
            is => 'String',
            doc => 'Location to write the output from Weka',
        },
    ],
    has_optional => [
        annotation_build => {
            is => 'Genome::Model::Build::ImportedAnnotation',
            id_by => 'annotation_build_id',
            doc => 'the annotation build used to annotate the variants of subbuilds in the convergence build--if not provided, will attempt to be taken from a subbuild',
        },
        annotation_build_id => { is => 'Text', is_input => 1, is_output => 1},

        mutation_types => {
            is => 'ARRAY',
            doc => 'The types of mutations to consider significant',
            default_value => [qw(
                        nonsense
                        missense
                        splice_site
                        splice_region
                        nonstop
                        cryptic_splice_site
                        splice_site_del
                        splice_site_ins
                        splice_region_del
                        splice_region_ins
                        frame_shift_del
                        in_frame_del
                        frame_shift_ins
                        in_frame_ins
                    )],
            is_constant => 1, #if this is adapted to a command line parameter, will have to be something other than an ARRAY
        },
    ],
};

sub help_brief {
    "Perform cross-genome analysis on the builds input to the convergence build",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
genome model convergence find-significant-attributes --build 102074749
EOS
}

sub help_detail {
    return <<EOS
Uses Weka to generate a decision tree using the presence of a mutation in transcripts
to classify a genome as "normal" or "tumor".
EOS
}

sub execute {
    my $self = shift;

    if(Genome::Sys->check_for_path_existence($self->output_file)) {
        $self->error_message('Output file already exists. Will not continue.');
        return;
    }

    my $build = $self->build;

    my @all_subbuilds = $build->all_subbuilds_closure;

    unless(scalar @all_subbuilds) {
        $self->status_message('No subbuilds found in model.  Skipping execution.');
        return 1;
    }

    #TODO Replace this with a check for any build able to produce annotation data
    my @applicable_subbuilds = grep($_->type_name eq 'reference alignment', @all_subbuilds);
    unless(scalar @applicable_subbuilds) {
        $self->status_message('No applicable subbuilds found in model.  Skipping execution.');
        return 1;
    }

    unless($self->annotation_build) {
        #TODO Be smarter about finding this; possibly verify all subbuilds used same one or support multiples, etc.
        my $annotation_build = $applicable_subbuilds[0]->model->annotation_reference_build;

        unless($annotation_build) {
            $self->error_message('No annotation build provided and one could not be determined from the subbuilds.');
            return;
        } else {
            $self->annotation_build($annotation_build);
        }
    }

    my $arff_file = $self->data_file;
    if(Genome::Sys->check_for_path_existence($arff_file)) {
        #TODO How do we know this data file has anything to do with the build in question?
        $self->debug_message('Found data file ' . $arff_file . ' for use in analysis.');
    } else {
        $self->debug_message('Generating data file ' . $arff_file . ' for use in analysis.');

        $self->generate_data_file(@applicable_subbuilds) or return;
    }

    return $self->analyze_data;
}

sub analyze_data {
    my $self = shift;

    my $arff_file = $self->data_file;
    my $output_file = $self->output_file;

    my $cmd = 'java -Xmx2048m -cp $ENV{GENOME_SW_LEGACY_JAVA}/weka.jar ' . 'weka.classifiers.trees.J48 -x 10 -B -C 0.25 -M 2 -t ' . $arff_file . ' > ' . $output_file;
    #my $cmd = 'java -Xmx2048m -cp $ENV{GENOME_SW_LEGACY_JAVA}/weka.jar ' . 'weka.classifiers.rules.DecisionTable -x 10 -X 1 -S "weka.attributeSelection.BestFirst -D 2 -N 30 -S 2" -t ' . $arff_file . ' -R > ' . $output_file;

    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$arff_file],
        output_files => [$output_file],
    );

    return 1;
}

sub generate_data_file {
    my $self = shift;
    my @builds = @_;

    my $data_fh = Genome::Sys->open_file_for_writing($self->data_file);

    $self->_generate_header($data_fh) or return;
    $self->_generate_data_section($data_fh, @builds) or return;

    close($data_fh);

    return 1;
}

sub _generate_header {
    my $self = shift;
    my $data_fh = shift;

    my $tx_iterator = $self->annotation_build->transcript_iterator;

    print $data_fh
        '@RELATION genomes', "\n", #TODO a better relation name--maybe the convergence model name?
        "\n";

    my $num_transcripts = 0;
    for(my $tx = $tx_iterator->next; $tx; $tx = $tx_iterator->next, $num_transcripts++) {
        print $data_fh
           '@ATTRIBUTE ', $self->_sanitize_name($tx->transcript_name), ' {0,1}', "\n";
    }

    print $data_fh
        '@ATTRIBUTE sample_type {tumor,normal}', "\n",
        "\n";

    return 1;
}

sub _generate_data_section {
    my $self = shift;
    my $data_fh = shift;
    my @builds = @_;

    print $data_fh
        '@DATA', "\n";

    BUILD: for my $build (@builds) {

        my $subject = $build->model->subject;
        unless($subject->class eq 'Genome::Sample') {
            #TODO Get the sample from whatever the subject is if possible
            $self->warning_message('Skipping build #' . $build->id . '. Unhandled subject type ' . $subject->class . ' encountered.');
            next BUILD;
        }

        #Check before looking up the transcripts, but will be printed at end of line
        my $subject_cn = $subject->common_name;
        unless ($subject_cn eq 'tumor' or $subject_cn eq 'normal') {
            $self->warning_message('Skipping build #' . $build->id . '. Common name not present or unhandled: ' . ($subject_cn? $subject_cn : '<NULL>'));
            next BUILD;
        }

        my $annotation_build = $self->annotation_build;
        my $annotation_path = $self->_resolve_annotation_file($build);

        unless($annotation_path) {
            $self->warning_message('Skipping build #' . $build->id . '. No annotation file found.');
            next BUILD;
        }

        my $annotation_fh = Genome::Sys->open_file_for_reading($annotation_path);

        #TODO Super-ugly hash... transcript names may not be "in order" in the annotated file due to overlappingness, etc.
        my %transcript_names_present;

        while(my $tx_name = $self->_next_annotation_tx_name($annotation_fh)) {
            $transcript_names_present{$tx_name} = 1;
        }

        $self->_unload_transcripts; #Lame hack--see sub for details
        my $tx_iterator = $annotation_build->transcript_iterator;

        my $present_count = 0;
        for(my $tx = $tx_iterator->next; $tx; $tx = $tx_iterator->next) {
            my $tx_present;

            if(exists $transcript_names_present{$tx->transcript_name}) {
                $tx_present = "1";
                $present_count += 1;
            } else {
                $tx_present = "0";
            }

            print $data_fh $tx_present . ',';
        }

        $self->status_message('Found ' . $present_count . ' mutated transcripts for build #' . $build->id);

        print $data_fh "$subject_cn\n";
    }

    return 1;
}

#TODO There should really be an accessor somewhere on builds that produce this file
sub _resolve_annotation_file {
    my $self = shift;
    my $build = shift;

    my $variant_directory = $build->variants_directory;

    for my $filename ('filtered.indelpe.snps.post_annotation', 'filtered.variants.post_annotation') {
        my $potential_annotation_path = join('/', $variant_directory, $filename);

        if(Genome::Sys->check_for_path_existence($potential_annotation_path)) {
            return $potential_annotation_path;
        }
    }

    return;
}

#Transform the name to meet Weka's specifications
sub _sanitize_name {
    my $self = shift;
    my $name = shift;
    my $type = 'tx'; #attribute name must begin with a letter, so make one up

    $name =~ s/[ ]/_/g; #We don't like spaces (alternatively could quote names everywhere...)
    $name = join('_', $type, $name);
    return $name;
}

sub _next_annotation_tx_name {
    my $self = shift;
    my $annotation_fh = shift;

    my $result;

    do {
        my $line = <$annotation_fh>;

        return unless $line;

        #For now we only care about $tx_name and $trv_type, but here are all the other columns
        my ($chr, $start, $stop, $ref, $var, $type, $gene_name, $tx_name, $tx_species, $tx_source,
            $tx_version, $strand, $tx_status, $trv_type, $c_pos, $amino_acid_change, $ucsc_cons,
            $domain, $all_domains)
            = split("\t", $line);

        $result = $tx_name;

        unless(grep($trv_type eq $_, @{ $self->mutation_types })) {
            $result = '-'; #match the unannotated output--undef would terminate the loop
        }

    } while ($result eq '-');

    return $result;
}

#FIXME Lame hack!... The iterator follows the order of transcripts in the file-based datasource only if UR doesn't think
#all objects have already been loaded--this code borrowed from UR::Context->clear_cache
sub _unload_transcripts {
    my $self = shift;

    my $class_name = 'Genome::Transcript';
    for my $obj ($UR::Context::current->all_objects_loaded_unsubclassed($class_name)) {
        # Check the type against %local_dont_unload again, because all_objects_loaded()
        # will return child class objects, as well as the class you asked for.
        my $obj_type = ref $obj;
#        next if ($local_dont_unload{$obj_type} ||
#                 grep {$obj_type->isa($_) } @{$args{'dont_unload'}});
        $obj->unload;
    }
    my @obj = grep { defined($_) } values %{ $UR::Context::all_objects_loaded->{$class_name} };
    if (@obj) {
        $self->warning_message("Skipped unload of $class_name objects during clear_cache: "
            . join(",",map { $_->id } @obj )
            . "\n"
        );
        if (my @changed = grep { $_->__changes__ } @obj) {
            require YAML;
            $self->error_message(
                "The following objects have changes:\n"
                . Data::Dumper::Dumper(\@changed)
                . "The clear_cache method cannot be called with unsaved changes on objects.\n"
                . "Use reverse_all_changes() first to really undo everything, then clear_cache(),"
                . " or call sync_database() and clear_cache() if you want to just lighten memory but keep your changes.\n"
                . "Clearing the cache with active changes will be supported after we're sure all code like this is gone. :)\n"
            );
            exit 1;
        }
    }
    no warnings qw(once);
    delete $UR::Context::all_objects_loaded->{$class_name};
    delete $UR::Context::all_objects_are_loaded->{$class_name};
    delete $UR::Context::all_params_loaded->{$class_name};

    return 1;
}

1;
