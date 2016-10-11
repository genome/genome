#!/usr/bin/perl

use strict;
use warnings;

use File::Basename;
use Genome;
use Data::Dumper;

my $usage = "$0 <WORK_ORDER_ID> <DIR_TO_COPY_TO> <?LIST_OF_MODEL_SUBJECT_NAMES_TO_COPY .. only copy these models?>"; # TODO <LIST_OF_SAMPLE_NAMES>
print "USAGE: $usage\n" and exit unless @ARGV && @ARGV >= 2;

my $wo_id = $ARGV[0];
my $gp = Genome::Project->get( $wo_id );
die "Can't find genome::project for work-order, $wo_id\n" unless $gp;

my $copy_dir = $ARGV[1];
die "Invalid directory to copy builds to, $copy_dir\n" unless -d $copy_dir;

my %subjects_to_copy;
if( $ARGV[2] ) {
    die "Can't file of list of samples to copy or file is empty, $ARGV[2]\n" unless -s $ARGV[2];
    my $fh = Genome::Sys->open_file_for_reading( $ARGV[2] );
    while( my $line = $fh->getline ) {
        chomp $line;
        $subjects_to_copy{ $line } = 1;
    }
    $fh->close;
}

# project parts
my @project_parts = $gp->parts;
die "Exiting .. can't find any project parts for work order, $wo_id\n" unless @project_parts;
#print Dumper \@project_parts; exit;

# display list of samples and inst_data/models available for work-order
print "Checking workorder for instrument-data and models\n";
&map_inst_data_and_models_to_samples;

my @models = grep{ $_->entity_class_name eq 'Genome::Model::DeNovoAssembly' } @project_parts;
die "Exiting .. didn't find any de-novo assembly project parts\n" unless @models;

print "Finding succeeded builds\n";
# find succeeded builds
my @succeeded_builds;
for my $model_part( @models ) {
    my $model = Genome::Model->get( $model_part->entity_id );
    die "Can't find genome::model for model id ".$model_part->entity_id."\n" unless $model;

    next if $ARGV[2] && ! exists $subjects_to_copy{ $model->subject->name };
    #print "Skipping model, ".$model->name.", which is not found in the input list of subjects to copy\n" and next
    #if $ARGV[2] && ! exists $subjects_to_copy{ $model->subject->name };

    next if model_is_newbler_assembly_with_multiple_inst_data( $model );

    my $build = $model->last_succeeded_build;
    print "Skipping model, ".$model->name.", sample ".$model->subject->name.", no succeeded build found\n" and next
    unless $build;

    push @succeeded_builds, $build;
}

print "Found ".scalar( @succeeded_builds )." succeeded builds to copy\n";#.join( "\n", map{ $_->data_directory }  @succeeded_builds ),"\n";
exit unless @succeeded_builds > 0;
for( @succeeded_builds ) {
    print "\t".$_->model->subject->name.' '.$_->data_directory."\n";
}
print "Continue and copy these assemblies? yes/no?\n";
my $ans = <STDIN>;
exit unless $ans =~ /yes/i;

# copy succeeded builds;
for my $build( @succeeded_builds ) {
    my $build_dir = $build->data_directory;
    die "No data directory found, $build_dir\n" unless -d $build_dir;

    # project dir
    my $output_dir = $copy_dir.'/'.$build->model->subject->name;
    $output_dir =~ s/\s+/_/g;
    Genome::Sys->create_directory( $output_dir ) unless -d $output_dir;

    # create build id empty file
    my $build_id_file = $output_dir."/build_id=".$build->id;
    Genome::Sys->touch( $build_id_file ) if not -e $build_id_file;

    my $assembler = $build->model->processing_profile->assembler_name;
    # create subdirs
    my @sub_dirs = subdirs_for_assembler( $assembler );
    die "Can't determine which subdirs to create for assembler, $assembler\n" unless @sub_dirs;
    for my $dir_name( @sub_dirs ) {
        my $sub_dir = $output_dir."/$dir_name";
        Genome::Sys->create_directory( $sub_dir ) unless -d $sub_dir;
    }

    # edit_dir files to copy
    my $builds_edit_dir = assemblers_edit_dir( $assembler );
    die "Can't determine locations of builds edit_dir for assembler, $assembler\n" unless $builds_edit_dir;
    for my $file( glob("$build_dir/$builds_edit_dir/*") ) {
        my $file_name = basename( $file );
        my $to = $output_dir.'/edit_dir/'.$file_name;
        print "skipping, $file already exists\n" and next if -s $to;
        print "Copying $file to $to\n";
        Genome::Sys->copy_file( $file, $to );
    }

    my @addl_files = additional_files_for_assembler( $assembler, $build );
    if( @addl_files ) {
        while( @addl_files ) {
            my( $from, $to ) = split(/\s+/, shift @addl_files);
            $from = $build_dir.'/'.$from;
            $to   = $output_dir.'/'.$to;
            die "Can't find file, $from to link\n" unless -s $from;
            print "Linking $from to $to\n";
            Genome::Sys->create_symlink( $from, $to ) if not -l $to;
        }
    }
}

exit;

sub model_is_newbler_assembly_with_multiple_inst_data {
    my $model = shift;
    my @inst_data = $model->instrument_data;
    my $assembler = $model->processing_profile->assembler_name;
    return 1 if $assembler =~ /newbler/i && @inst_data > 1;
    return;
}

sub additional_files_for_assembler {
    my $assembler = shift;
    my $build = shift;
    my @addl_files;
    if( uc $assembler eq 'NEWBLER DE-NOVO-ASSEMBLE' ) {
        my @velvet_input_fastqs = velvet_input_fastq_files( $build );
        for( @velvet_input_fastqs ) {
            my $file_name = basename( $_ );
            push @addl_files, $file_name." edit_dir/$file_name";
        }
        push @addl_files, 'consed/phdball_dir/phd.ball.1 phdball_dir/phd.ball.1';
    }
    return @addl_files;
}

sub velvet_input_fastq_files {
    my $build = shift;
    my @input_fastqs = glob( $build->data_directory."/*input.fastq" );
    return @input_fastqs;
}

sub subdirs_for_assembler {
    my $assembler = shift;
    return (qw/ edit_dir chromat_dir phd_dir/) if uc $assembler eq 'VELVET ONE-BUTTON';
    return (qw/ edit_dir chromat_dir phd_dir phdball_dir /) if uc $assembler eq 'NEWBLER DE-NOVO-ASSEMBLE';
    return;
}

sub assemblers_edit_dir {
    my $assembler = shift;
    return 'edit_dir' if uc $assembler eq 'VELVET ONE-BUTTON';
    return 'consed/edit_dir' if uc $assembler eq 'NEWBLER DE-NOVO-ASSEMBLE';
    return;
}

sub map_inst_data_and_models_to_samples {
    my %h;
    # samples should always exist
    my $sample_parts = get_entity_class_parts( 'Sample' );
    die "Exiting .. no Genome::Sample project parts found in work-order\n" unless $sample_parts;
    for my $project_part( @$sample_parts ) {
        my $genome_obj = get_genome_class_obj( $project_part );
        # expect this to always exist
        die "Can't get genome sample for id, ".$project_part->entity_id."\n" unless $genome_obj;
        my $sample_name = $genome_obj->name;
        $h{ $sample_name }{ sample_name } = $sample_name;
    }
    # inst data may or may not exist
    my $inst_data = get_entity_class_parts( 'InstrumentData' );
    for my $project_part( @$inst_data ) {
        my $genome_obj = get_genome_class_obj( $project_part );
        # these can exist as wo parts but not actually
        # exist in GMS db because they haven't synced
        next unless $genome_obj;
        my $sample_name = $genome_obj->sample->name;
        die "Invalid instrument data in project .. found instrument data sample name, $sample_name, but sample does not exist as project part\n" unless
        exists $h{ $sample_name };
        push @{$h{ $sample_name }{ inst_data}}, $genome_obj->id;
    }
    # models may or may not exist
    my $models = get_entity_class_parts( 'Model::DeNovoAssembly' );
    for my $project_part( @$models ) {
        my $genome_obj = get_genome_class_obj( $project_part );
        next if model_is_newbler_assembly_with_multiple_inst_data( $genome_obj );
        die "Can't get genome model for id, ".$project_part->entity_id."\n" unless $genome_obj;
        # TODO - below could be taxon
        my $sample_name = $genome_obj->subject->name;
        next if $sample_name eq 'Unknown';
        die "Invalid model subject name, $sample_name was found but it does not exist as a project part\n" unless
        exists $h{ $sample_name };
        push @{$h{ $sample_name }{ models }}, $genome_obj->name;
    }
    #print Dumper \%h;

    print_project_parts( \%h );

    return \%h;
}

sub print_project_parts {
    my $parts = shift;

    my $spacing = "%-30s%-30s%-50s";
    my @headers = (qw/ SAMPLE_NAMES INSTRUMENT_DATA MODELS /);
    printf("$spacing\n", @headers );

    for my $sample_name( sort keys %$parts ) {
        my @output = ( $sample_name );
        if( $parts->{ $sample_name }->{inst_data} ) {
            push @output, join(',', @{$parts->{$sample_name}->{inst_data}});
        }
        else {
            push @output, 'NA';
        }
        if( $parts->{ $sample_name }->{models} ) {
            push @output, join( ',', @{$parts->{$sample_name}->{models}} );
        }
        else {
            push @output, 'NA';
        }
        printf( "$spacing\n", @output );
    }
}

sub get_entity_class_parts {
    my $entity_class = shift;
    my @parts = grep{ $_->entity_class_name eq 'Genome::'.$entity_class } @project_parts;
    return \@parts;
}

sub get_genome_class_obj {
    my $project_part = shift;
    my $genome_class = $project_part->entity_class_name;
    my $obj = $genome_class->get( $project_part->entity_id );
    return $obj if $obj;
    return;
}

