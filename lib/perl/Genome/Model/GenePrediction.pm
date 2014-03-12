package Genome::Model::GenePrediction;

use strict;
use warnings;
use Genome;
use Carp 'confess';

class Genome::Model::GenePrediction {
    is => 'Genome::ModelDeprecated',
    is_abstract => 1,
    subclassify_by => 'subclass_name',
    has => [
        subclass_name => {
            calculate_from => 'processing_profile_id',
            calculate => sub {
                my $pp_id = shift;
                confess "No processing profile id given to gene prediction model!" unless defined $pp_id;
                my $pp = Genome::ProcessingProfile::GenePrediction->get($pp_id);
                confess "Could not get processing profile $pp_id to create gene prediction model!" unless $pp;
                my $subclass = $pp->subclass_name;
                $subclass =~ s/ProcessingProfile/Model/;
                return $subclass;
            },
        },
        domain => {
            is => 'Text',
            valid_values => ['bacterial', 'archaeal', 'eukaryotic'],
            doc => 'Domain of organism',
            is_param => 1,
        },
        subject => {
            is => 'Genome::Taxon',
            id_by => 'subject_id',
        },
    ],
    has_optional => [
        assembly_model_links => {
            is => 'Genome::Model::Link',
            is_many => 1,
            reverse_as => 'to_model',
            where => [ role => 'gene_prediction' ],
        },
        assembly_model => {
            is => 'Genome::Model',
            via => 'assembly_model_links',
            to => 'from_model',
        },
        gram_stain => {
            via => 'subject',
            to => 'gram_stain_category',
        },
        # The species latin name on some taxons includes a strain name, which needs to be removed
        organism_name => {
            calculate_from => ['subject'],
            calculate => q( 
                my $latin_name = $subject->species_latin_name;
                $latin_name =~ s/\s+/_/g; 
                my $first = index($latin_name, "_");
                my $second = index($latin_name, "_", $first + 1);
                return $latin_name if $second == -1;
                return substr($latin_name, 0, $second);
            ),
        },
        ncbi_taxonomy_id => {
            via => 'subject',
            to => 'ncbi_taxon_id',
        },
        assembly_contigs_file => {
            is => 'FilePath',
            is_input => 1,
        },
    ],
};

sub create {
    my $class = shift;
    my %params = @_;

    if ($class eq __PACKAGE__ or $class->__meta__->is_abstract) {
        my $self = $class->SUPER::create(%params);
        return $self;
    }

    my $create_assembly_model = delete $params{create_assembly_model} || 0;
    my $assembly_processing_profile_name = delete $params{assembly_processing_profile_name};
    my $start_assembly_build = delete $params{start_assembly_build} || 0;
    my $assembly_contigs_file = delete $params{assembly_contigs_file};

    my $self = $class->SUPER::create(%params) or return;

    # Taxon must have certain fields defined
    my $taxon = $self->subject;
    $self->validate_taxon($taxon);

    # If given a path to a contigs file, don't go through the trouble of getting/creating an assembly model
    # TODO Eventually, this behavior should be phased out so we are completely reliant on assembly builds, but
    # right now there aren't any assembly models for eukaryotic organisms.
    if (defined $assembly_contigs_file) {
        unless (-e $assembly_contigs_file and -s $assembly_contigs_file) {
            confess "No file or zero-sized file at $assembly_contigs_file, cannot create gene prediction model!";
        }

        # It's expected that this is a fasta file, make sure it is
        my $format_guesser = Bio::Tools::GuessSeqFormat->new(-file => $assembly_contigs_file);
        my $format = $format_guesser->guess;
        unless ($format eq 'fasta') {
            confess "Assembly contigs file is not in fasta format (format seems to be $format)! Currently only fasta format is supported!";
        }

        $self->add_input(
            value_class_name => 'UR::Value::FilePath',
            value_id => $assembly_contigs_file,
            name => 'assembly_contigs_file',
        );

        return $self;
    }
    
    # Either get a denovo assembly model or create one if the create_assembly_model flag is defined
    my $assembly_model = $self->find_or_create_assembly_model(
        create_assembly_model => $create_assembly_model,
        assembly_processing_profile_name => $assembly_processing_profile_name,
    );
    
    # Make a link from the assembly model to this model. This link will be used to kick off a new build
    # of the gene prediction model every time its linked assembly model completes a build.
    my $to_rv = $assembly_model->add_to_model(
        to_model => $self,
        role => 'gene_prediction',
    );
    unless ($to_rv) {
        $self->error_message("Could not create a link from the assembly model to the gene prediction model!");
        confess;
    }

    # If there's not a successful build on the assembly model, start one if the proper flag is set
    my $build = $assembly_model->last_succeeded_build;
    $build = $assembly_model->current_running_build unless $build;
    unless ($build) {
        if ($start_assembly_build) {
            $self->status_message("Starting build of assembly model " . $assembly_model->id);
            $build = $self->start_assembly_build($assembly_model->id);
            $self->status_message("Build started!");
            # TODO Either wait for the assembly to finish or make sure the assembly model kicks
            # off the gene prediction build when it's done.
        }
        else {
            $self->status_message(
                "The --start-assembly-build option is not set, so automatic build of the assembly model will not occur.\n" .
                "Either rerun this command with the --start-assembly-build flag or manually kick off a build by running:\n" .
                "genome model build start --model " . $assembly_model->id);
            confess "No assembly build found and assemble flag not set";
        }
    }

    $self->add_input(
        value_class_name => 'UR::Value::FilePath',
        value_id => $self->find_assembly_contigs_file(),
        name => 'assembly_contigs_file',
    );

    return $self;
}

# Starts a build of the given assembly model
sub start_assembly_build {
    my ($self, $assembly_model) = @_;
    my $assembly_model_id = $assembly_model->id;

    unless ($assembly_model->isa('Genome::Model::DeNovoAssembly') or 
            $assembly_model->isa('Genome::Model::ImportedAssembly')) {
        $self->error_message("Model $assembly_model_id is not an assembly model! Not starting build!");
        confess;
    }

    my $start_command = Genome::Model::Build::Command::Start->create(
        model_identifier => $assembly_model_id,
    );
    unless ($start_command) {
        $self->error_message("Could not create the build start command object!");
        confess;
    }

    my $rv = $start_command->execute;
    unless (defined $rv and $rv == 1) {
        $self->error_message("Could not start build of assembly model!");
        confess;
    }

    my $build = $start_command->build;
    unless ($build) {
        $self->error_message("Could not get build from create command object!");
        confess;
    }

    return $build;
}

sub get_most_recent_model {
    my ($self, $assembly_models) = @_;
    my @sorted_models = sort {$b->creation_date cmp $a->genome_model_id} @$assembly_models;
    return shift @sorted_models;
}

# Checks that the necessary fields on the taxon are defined and confesss if any aren't
sub validate_taxon {
    my ($self, $taxon) = @_;

    my @missing_fields;
    # Gram stain doesn't apply to some (all?) eukaryotes
    push @missing_fields, "gram stain" if $self->domain =~ /Bacterial/ and not defined $taxon->gram_stain_category;
    push @missing_fields, "domain" unless defined $taxon->domain;

    if (@missing_fields) {
        $self->error_message("These fields on taxon " . $taxon->id .
            " need to be defined!\n" . join("\n", @missing_fields));
        confess;
    }

    return 1;
}

# Returns the path to the assembly build's contigs.bases file
sub find_assembly_contigs_file {
    my $self = shift;
    my $assembly_model = $self->assembly_model;
    my $assembly_build = $assembly_model->last_succeeded_build;
    my $path = $assembly_build->data_directory . "/edit_dir/contigs.bases";
    return $path;
}

# Find an assembly model with the same subject as this model or create one
# if the proper flag is set and a processing profile name is provided
sub find_or_create_assembly_model {
    my $self = shift;
    my %params = @_;
    my $create_assembly_model = $params{create_assembly_model};
    my $assembly_processing_profile_name = $params{assembly_processing_profile_name};
    my $taxon = $self->subject;

    my @denovo = Genome::Model::DeNovoAssembly->get(
        subject_name => $taxon->name,
        subject_class_name => 'Genome::Taxon',
    );
    my @imported = Genome::Model::ImportedAssembly->get(
        subject_name => $taxon->name,
        subject_class_name => 'Genome::Taxon',
    );

    my $assembly_model;
    if (@denovo) {
        $assembly_model = $self->get_most_recent_model(\@denovo);
        $self->status_message("Using de novo assembly model " . $assembly_model->id);
    }
    elsif (@imported) {
        $assembly_model = $self->get_most_recent_model(\@imported);
        $self->status_message("Using imported assembly model " . $assembly_model->id);
    }
    else {
        if ($create_assembly_model) {
            $self->status_message("Did not find any de novo or imported assemblies and create " .
                "assembly model flag is set, creating a new assembly model using processing " .
                "profile with name " . $assembly_processing_profile_name);

            unless (defined $assembly_processing_profile_name) {
                $self->error_message("Not provided with assembly model processing profile name! Can't create model!");
                confess;
            }

            my $assembly_define_obj = Genome::Model::Command::Define::DeNovoAssembly->create(
                processing_profile_name => $assembly_processing_profile_name,
                subject_name => $taxon->name,
            );
            unless ($assembly_define_obj) {
                $self->error_message("Could not create assembly model define object!");
                confess;
            }

            my $define_rv = $assembly_define_obj->execute;
            unless (defined $define_rv and $define_rv == 1) {
                $self->error_message("Trouble while attempting to define new assembly model!");
                confess;
            }

            my $model_id = $assembly_define_obj->result_model_id;
            $assembly_model = Genome::Model::DeNovoAssembly->get($model_id);
            unless ($assembly_model) {
                $self->error_message("Could not get newly created assembly model with ID $model_id!");
                confess;
            }

            $self->status_message("Successfully created assembly model with ID $model_id, now assigning data!");

            my $assembly_assign_obj = Genome::Model::Command::InstrumentData::Assign::AllCompatible->create(
                model => $assembly_model,
            );
            unless ($assembly_assign_obj) {
                $self->error_message("Could not create instrument data assignment object!");
                confess;
            }

            my $assign_rv = $assembly_assign_obj->execute;
            unless (defined $assign_rv and $assign_rv == 1) {
                $self->error_message("Trouble while attempting to assign instrument data to model!");
                confess;
            }

            $self->status_message("Instrument data has been successfully assigned to model!");
        }
        else {
            $self->status_message(
                "Could not find any assembly models with the taxon you provided. If you would like to create an " .
                "assembly model for use with this gene prediction model, run the following command\n" .
                "genome model define de-novo-assembly --processing-profile-name " . 
                $assembly_processing_profile_name . " --subject-name " . $taxon->name . "\n\n" .
                "That command should give you a model ID. Use it to assign instrument data to the assembly model:\n" .
                "genome instrument-data assign --model-id <MODEL_ID> --all\n\n" .
                "Now you have an assembly model with instrument data! Rerun this define command with the " .
                "--start-assembly-build option, which will start a build of that assembly model and kick off gene prediction " .
                "once that build has completed!"
            );
            confess "Could not find any assemblies for taxon " . $taxon->id;
        }
    }

    return $assembly_model;
}

sub _resolve_type_name_for_class {
    return 'gene prediction';
}

1;

