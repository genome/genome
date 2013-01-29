package Genome::Model::Build::MetagenomicComposition16s;

use strict;
use warnings;

use Genome;

require Carp;
require File::Basename;
require Mail::Sendmail;
use Switch;

class Genome::Model::Build::MetagenomicComposition16s {
    is => 'Genome::Model::Build',
    is_abstract => 1,
    subclassify_by => 'subclass_name',
    has => [
        subclass_name => { 
            is => 'String', len => 255, is_mutable => 0,
            calculate_from => ['model_id'],
            calculate => sub {
                my($model_id) = @_;
                return unless $model_id;
                my $model = Genome::Model->get($model_id);
                Carp::croak("Can't find Genome::Model with ID $model_id while resolving subclass for Build") unless $model;
                my $seq_platform = $model->sequencing_platform;
                Carp::croak("Can't subclass Build: Genome::Model id $model_id has no sequencing_platform") unless $seq_platform;
                return return __PACKAGE__ . '::' . Genome::Utility::Text::string_to_camel_case($seq_platform);
            },
        },
        map( { 
                $_ => { via => 'processing_profile' } 
            } Genome::ProcessingProfile::MetagenomicComposition16s->params_for_class 
        ),
        map ( {
                $_ => { is => 'Number', is_metric => 1, }
            } (qw/
                amplicons_attempted amplicons_processed amplicons_processed_success 
                amplicons_classified amplicons_classified_success amplicons_classification_error 
                amplicons_chimeric amplicons_chimeric_percent 
                reads_attempted reads_processed reads_processed_success 
                /)
        ),
    ],
};

sub length_of_16s_region {
    return 1542;
}

sub post_allocation_initialization {
    my $self = shift;
    return $self->create_subdirectories;
}

sub validate_instrument_data{ # overloaded
    my $self = shift;

    my @tags;
    my @instrument_data = $self->instrument_data;
    if ( not @instrument_data ) {
        push @tags, UR::Object::Tag->create(
            type => 'error',
            properties => ['instrument_data'],
            desc => 'No instrument data assigned to this build!',
        );
    }

    my $reads;
    for my $instrument_data ( @instrument_data ) {
        $reads += $instrument_data->read_count;
    }

    if ( not $reads ){
        push @tags, UR::Object::Tag->create(
            type => 'error',
            properties => [ 'instrument_data' ],
            desc => 'No reads for  instrument data assigned to build!',
        );
    }

    return @tags;
}

sub description {
    my $self = shift;

    return sprintf(
        'metagenomic composition 16s build (%s) for model (%s %s)',
        $self->id,
        $self->model->name,
        $self->model->id,
    );
}

sub amplicon_sets {
    my $self = shift;

    my %amplicon_set_names = map { $_->sequencing_platform => 1 } $self->instrument_data;
    my $amplicon_set_name = (keys %amplicon_set_names)[0];
    my %amplicon_set_names_and_primers = Genome::Model::Build::MetagenomicComposition16s::SetNamesAndPrimers->set_names_and_primers_for($amplicon_set_name);
    my @amplicon_sets;
    for my $set_name ( sort { $a cmp $b } keys %amplicon_set_names_and_primers ) {
        push @amplicon_sets, Genome::Model::Build::MetagenomicComposition16s::AmpliconSet->create(
            name => $set_name,
            primers => $amplicon_set_names_and_primers{$set_name},
            file_base_name => $self->file_base_name,
            directory => $self->data_directory,
            classifier => $self->classifier,
            chimera_detector => $self->chimera_detector,
        );
    }

    unless ( @amplicon_sets ) {
        $self->error_message("No amplicon sets found for ".$self->description);
        return;
    }

    return @amplicon_sets;
}

sub amplicon_sets_for_processing {
    my $self = shift;

    my @amplicon_sets = $self->amplicon_sets;
    return if not @amplicon_sets;

    if ( @amplicon_sets > 1 ) {
        push @amplicon_sets, Genome::Model::Build::MetagenomicComposition16s::AmpliconSet->create(
            name => 'none',
            directory => $self->data_directory,
            file_base_name => $self->file_base_name,
            classifier => $self->classifier,
        );
    }

    return @amplicon_sets;
}
#<>#

#< Dirs >#
sub create_subdirectories {
    my $self = shift;
    my @methods = (qw| classification_dir fasta_dir reports_dir |);
    push @methods, (qw/ chimera_dir /) if $self->processing_profile->chimera_detector;
    for my $method ( @methods ) {
        my $directory =  $self->$method;
        my $create_directory = eval{ Genome::Sys->create_directory($directory); };
        if ( not $create_directory or not -d $directory ) {
            $self->error_message($@) if $@;
            $self->error_message('Failed to create directory: '.$directory);
            return;
        }
    }
    return 1;
}

sub chimera_dir {
    return $_[0]->data_directory.'/chimera';
}

sub classification_dir {
    return $_[0]->data_directory.'/classification';
}

sub fasta_dir {
    return $_[0]->data_directory.'/fasta';
}

sub reports_dir {
    return $_[0]->data_directory.'/reports';
}
#<>#

#< Files >#
sub file_base_name {
    return Genome::Utility::Text::sanitize_string_for_filesystem( $_[0]->subject_name );
}

sub combined_input_fasta_file {
    return $_[0]->fasta_dir.'/'.$_[0]->file_base_name.'.'.'input.fasta';
}

sub processed_fasta_file { # returns them as a string (legacy)
    return join(' ', $_[0]->processed_fasta_files);
}

sub processed_fasta_files {
    my $self = shift;
    my @amplicon_sets = $self->amplicon_sets;
    return if not @amplicon_sets;
    return map { $_->processed_fasta_file } @amplicon_sets;
}

sub processed_qual_file { # returns them as a string (legacy)
    return join(' ', $_[0]->processed_qual_files);
}

sub processed_qual_files {
    my $self = shift;
    return map { $_->processed_qual_file } $self->amplicon_sets;
}

sub combined_original_fastq_file {
    my $self = shift;
    return sprintf(
        '%s/%s.%s.fastq',
        $self->fasta_dir,
        $self->file_base_name,
        'original',
    );
}

sub oriented_fasta_file { # returns them as a string
    return join(' ', $_[0]->oriented_fasta_files);
}

sub oriented_fasta_files {
    my $self = shift;
    return map { $_->oriented_fasta_file } $self->amplicon_sets;
}

sub oriented_qual_file { # returns them as a string (legacy)
    return join(' ', $_[0]->oriented_qual_files);
}

sub oriented_qual_files {
    my $self = shift;
    return map { $_->oriented_qual_file } $self->amplicon_sets;
}

sub classification_files_as_string {
    return join(' ', $_[0]->classification_files);
}

sub classification_files {
    my $self = shift;
    return map { $_->classification_file } $self->amplicon_sets;
}
#<>#

#< Process instrument data >#
sub sx_result_params_for_instrument_data {
    my ($self, @instrument_data) = @_;

    Carp::confess('No instrument data to get sx result params!') if not @instrument_data;

    my @amplicon_sets = $self->amplicon_sets_for_processing;
    return if not @amplicon_sets;

    my (@primers, @output_configs);
    for my $amplicon_set ( @amplicon_sets ) {
        my @set_primers = $amplicon_set->primers;
        for my $primer ( @set_primers ) {
            push @primers, $amplicon_set->name.'='.$primer;
        }
        push @output_configs, {
            basename => $amplicon_set->base_name_for('processed_fastq'),
            type => 'sanger',
        };
        my $writer_name = $amplicon_set->name;
        if ( $writer_name ) {
            if ( $writer_name eq 'none' ) {
                $writer_name = 'discard';
            }
            else {
                #strip off .F/.R for paired sets
                $writer_name =~ s/\.[FR]$//;
            }
        }
        $output_configs[$#output_configs]->{name} = $writer_name;
    }

    my @read_processor = ( 'rm-desc' );
    if ( @amplicon_sets > 1 ) {
        push @read_processor, 'bin by-primer --remove --primers '.join(',', @primers);
    }

    my @output_file_configs = map { Genome::Model::Tools::Sx::Functions->hash_to_config(%$_) } @output_configs;
    if ( @output_configs and ( not @output_file_configs or @output_file_configs != @output_configs ) ) {
        $self->error_message('Failed to get output file config as strings!');
        return;
    }

    if ( $self->processing_profile->amplicon_processor ) {
        push @read_processor, $self->processing_profile->amplicon_processor;
    }

    return (
        instrument_data_id => ( @instrument_data > 1 ? [ map { $_->id } @instrument_data ] : $instrument_data[0]->id ),
        read_processor => join(' | ', @read_processor),
        output_file_config => \@output_file_configs,
    );
}

sub process_instrument_data {
    my ($self, $instrument_data) = @_;
    $self->status_message('Process instrument data...');

    Carp::confess('No instrument data given to process!') if not $instrument_data;

    my %sx_result_params = $self->sx_result_params_for_instrument_data($instrument_data);
    return if not %sx_result_params;

    my $sx_result = Genome::InstrumentData::SxResult->get_or_create(%sx_result_params);
    if ( not $sx_result ) {
        $self->error_message('Failed to create SX result!');
        return;
    }

    $sx_result->add_user(user => $self, label => 'sx_result');

    $self->status_message('Process instrument data...OK');
    return 1;
}

sub merge_processed_instrument_data {
    my $self = shift;
    $self->status_message('Merge processed instrument data...');

    my @amplicon_sets = $self->amplicon_sets_for_processing;
    return if not @amplicon_sets;

    my @instrument_data = $self->instrument_data;
    my %sx_result_params = $self->sx_result_params_for_instrument_data(@instrument_data);
    return if not %sx_result_params;

    my @sx_results = Genome::InstrumentData::SxResult->get(%sx_result_params);
    if ( not @sx_results or @sx_results != @instrument_data ) {
        $self->error_message('Failed to find sx results for instrument data!');
        return;
    }

    $self->status_message('Merge SX results...');
    for my $amplicon_set ( @amplicon_sets ) {
        my @sx_processed_fastq_files;
        for my $sx_result ( @sx_results ) {
            my $sx_processed_fastq_file = $sx_result->output_dir.'/'.$amplicon_set->base_name_for('processed_fastq');
            next if not -s $sx_processed_fastq_file; # ok
            push @sx_processed_fastq_files, $sx_processed_fastq_file;
        }

        my $reader = Genome::Model::Tools::Sx::Reader->create(config => \@sx_processed_fastq_files);
        Carp::confess('Failed to open fastq reader for processed sx result fastqs! '.join(', ', @sx_processed_fastq_files)) if not $reader;

        my $processed_fasta_file = $amplicon_set->processed_fasta_file;
        my $writer = Genome::Model::Tools::Sx::PhredWriter->create(
            file => $processed_fasta_file,
            qual_file => $amplicon_set->processed_qual_file,
            mode => 'a',
        );
        Carp::confess("Failed to open phred writer for file! ".$amplicon_set->processed_fasta_file) if not $writer;
        while ( my $seqs = $reader->read ) { 
            for my $seq ( @$seqs ) { 
                $writer->write($seq); 
            }
        }
    }
    $self->status_message('Merge SX results...OK');

    $self->status_message('Instrument data processed metrics...');
    my $inst_data_metrics_fh = Genome::Sys->open_file_for_writing( $self->fasta_dir.'/metrics.processed' );
    my %metrics;
    for my $sx_result ( @sx_results ) {
        my $instrument_data = $sx_result->instrument_data;
        # link original data path
        my $original_data_path = 'NA';
        ORIGINAL_DATA_PATH: for my $attribute_label (qw/ bam_path sff_file archive_path /) {
            my $attribute = $instrument_data->attributes(attribute_label => $attribute_label);
            next if not $attribute;
            my $attribute_value = $attribute->attribute_value;
            next if not -e $attribute_value;
            last ORIGINAL_DATA_PATH;
            $original_data_path = $attribute_value;
        }

        # metrics: link files and get attempted and processed
        my %sx_metrics;
        for my $type (qw/ input output /) {
            my $file_method = 'read_processor_'.$type.'_metric_file';
            my $metrics_file_name = $sx_result->$file_method;
            my $metrics_file = $sx_result->output_dir.'/'.$metrics_file_name;
            my $metrics = Genome::Model::Tools::Sx::Metrics->from_file($metrics_file);
            Carp::confess("Failed to get $type metrics for SX result! ".$sx_result->id) if not $metrics;
            $sx_metrics{$type.'_bases'} += $metrics->bases;
            $sx_metrics{$type.'_count'} += $metrics->count;
        }
        $metrics{amplicons_attempted} += $sx_metrics{input_count};
        $metrics{amplicons_processed} += $sx_metrics{output_count};
        $inst_data_metrics_fh->print(
            join(' ', $metrics{input_count}, $sx_metrics{output_count}, $original_data_path)."\n"
        );
    }
    $inst_data_metrics_fh->close;
    $self->status_message('Instrument data processed metrics...OK');

    $self->amplicons_attempted( $metrics{amplicons_attempted} );
    $self->amplicons_processed( $metrics{amplicons_processed} );
    $self->amplicons_processed_success( 
        $metrics{amplicons_attempted} > 0 ?  sprintf('%.2f', $metrics{amplicons_processed} / $metrics{amplicons_attempted}) : 0
    );
    $self->status_message('Attempted: '.$self->amplicons_attempted);
    $self->status_message('Processed: '.$self->amplicons_processed);
    $self->status_message('Success:   '.($self->amplicons_processed_success * 100).'%');

    $self->status_message('Merge processed instrument data...');
    return 1;
}

sub prepare_instrument_data {
    my $self = shift;
    $self->status_message('Prepare instrument data...');

    my @instrument_data = $self->instrument_data;
    $self->status_message('Instrument data count: '.@instrument_data);
    unless ( @instrument_data ) {
        $self->error_message("No instrument data found for ".$self->description);
        return;
    }

    #writer original/unprocessed/combined fastq file
    if ( not $self->fastq_from_instrument_data ) {
        Carp::confess( "Failed to get fasta from instrument data" );
    }
    #read in orig fastq file
    my $orig_fastq = $self->combined_original_fastq_file;
    if ( not -s $orig_fastq ) {
        Carp::confess( "Original fasta file did not get created or is blank" );
    }

    my @cmd_parts = ( 'gmt sx rm-desc' );
    my @amplicon_sets = $self->amplicon_sets;
    my (@output, @primers);
    for my $amplicon_set ( @amplicon_sets ) {
        my @set_primers = $amplicon_set->primers;
        for my $primer ( @set_primers ) {
            push @primers, $amplicon_set->name.'='.$primer;
        }
        my $root_set_name = $amplicon_set->name;
        $root_set_name =~ s/\.[FR]$//; #strip off .F/.R for paired sets
        my $fasta_file = $amplicon_set->processed_fasta_file;
        my $qual_file = $amplicon_set->processed_qual_file;
        unlink $fasta_file, $fasta_file;
        my $output = 'file='.$fasta_file.':qual_file='.$qual_file.':type=phred';
        $output .= ':name='.$root_set_name if @set_primers;
        push @output, $output;
    }

    if ( @primers ) {
        my $none_amplicon_set = Genome::Model::Build::MetagenomicComposition16s::AmpliconSet->create(
            name => 'none',
            directory => $self->data_directory,
            file_base_name => $self->file_base_name,
            classifier => $self->classifier,
        );
        my $none_fasta_file = $none_amplicon_set->processed_fasta_file;
        my $none_qual_file = $none_amplicon_set->processed_qual_file;
        unlink $none_fasta_file, $none_qual_file;
        push @output, 'name=discard:file='.$none_fasta_file.':qual_file='.$none_qual_file.':type=phred';
        push @cmd_parts, 'gmt sx bin by-primer --remove --primers '.join(',', @primers);
    }

    #add amplicon processing sx commands
    if ( $self->processing_profile->amplicon_processor ) {
        push @cmd_parts, $self->processing_profile->amplicon_processor_commands;
    }

    # Add input and metrics to first cmd
    $cmd_parts[0] .= ' --input file='.$orig_fastq.':type=sanger';
    my $input_metrics_file = $self->fasta_dir.'/metrics.processed.in.txt';
    $cmd_parts[0] .= ' --input-metrics '.$input_metrics_file;

    # Add output and metrics to last cmd
    $cmd_parts[$#cmd_parts] .= ' --output '.join(',', @output);
    my $output_metrics_file = $self->fasta_dir.'/metrics.processed.out.txt';
    $cmd_parts[$#cmd_parts] .= ' --output-metrics '.$output_metrics_file;

    # Create, execute, check return
    my $cmd = join(' | ', @cmd_parts);
    $self->status_message('Run SX...');
    $self->status_message("SX comand: $cmd");
    my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
    if ( not $rv ) {
        $self->error_message('Failed to run sx command to create amplicon files.');
        return;
    }
    $self->status_message('Run SX...OK');

    # Rm empty output files
    $self->status_message('Remove empty otuput files...');
    for my $amplicon_set ( @amplicon_sets ) {
        for my $file_method (qw/ processed_fasta_file processed_qual_file /) {
            my $file = $amplicon_set->$file_method;
            my $sz = -s $file;
            unlink $file if not $sz or $sz == 0;
        }
    }
    $self->status_message('Remove empty otuput files...OK');

    # Set metrics
    $self->status_message('Get metrics...');
    my $input_metrics = Genome::Model::Tools::Sx::Metrics->from_file($input_metrics_file);
    if ( not $input_metrics ) {
        $self->error_message('Failed to get metrcis from file: '.$input_metrics_file);
        return;
    }
    my $attempted = $input_metrics->count;
    $attempted = 0 if not defined $attempted;
    my $reads_attempted = $attempted;

    my $output_metrics = Genome::Model::Tools::Sx::Metrics->from_file($output_metrics_file);
    if ( not $output_metrics ) {
        $self->error_message('Failed to get metrcis from file: '.$output_metrics_file);
        return;
    }
    my $processed = $output_metrics->count;
    $processed = 0 if not defined $processed;
    $self->status_message('Get metrics...OK');

    $self->amplicons_attempted($attempted);
    $self->amplicons_processed($processed);
    my $processed_success =  $attempted > 0 ?  sprintf('%.2f', $processed / $attempted) : 0;
    $self->amplicons_processed_success($processed_success);
    $self->reads_attempted($reads_attempted);
    $self->reads_processed($processed);
    $self->reads_processed_success( $reads_attempted > 0 ?  sprintf('%.2f', $processed / $reads_attempted) : 0 );

    $self->status_message('Attempted:  '.$self->amplicons_attempted);
    $self->status_message('Processed:  '.$self->amplicons_processed);
    $self->status_message('Success:    '.($self->amplicons_processed_success * 100).'%');

    $self->status_message('Prepare instrument data...OK');
    return 1;
}

sub original_fastq_writer {
   my $self = shift;

   return $self->{_fastq_writer} if $self->{_fastq_writer};
   
   my $fastq_file = $self->combined_original_fastq_file;
   unlink $fastq_file if -e $fastq_file;

   my $writer = Genome::Model::Tools::Sx::Writer->create( config => [$fastq_file.':type=sanger'] );
   Carp::confess("Failed to create fastq writer to write original fastq file") if not $writer;
   
   $self->{_fastq_writer} = $writer;
   
   return $self->{_fastq_writer};
}

#< Fasta/Qual Readers/Writers >#
sub fastq_from_instrument_data {
    my $self = shift;
    
    for my $inst_data ( $self->instrument_data ) {

        my $temp_dir = Genome::Sys->create_temp_directory;
        my @fastq_files;
        if ( @fastq_files = eval {$inst_data->dump_fastqs_from_bam( directory => $temp_dir );} ) {
            #solexa bam
            if ( not @fastq_files ) {
                Carp::confess( "Did not get any fastq files from instrument data bam path" );
            }
            $self->status_message( "Got fastqs from bam" );
            #write fastq to original.fastq.file
            if ( not $self->append_fastq_to_orig_fastq_file( @fastq_files ) ) {
                Carp::confess( "Attempt to get fastq from fastq via bam failed" );
            }
        }
        elsif (  my $rv = eval {Genome::Sys->shellcmd( cmd => "tar zxf " . $inst_data->archive_path ." -C $temp_dir" );} ) {
            #solexa archive
            my @fastq_files = glob $temp_dir.'/*';
            if ( not @fastq_files ) {
                Carp::confess( "Did not get any fastq files from instrument data archive path" );
            }
            $self->status_message( "Untarred fastqs from archive path" );
            if ( not $self->append_fastq_to_orig_fastq_file( @fastq_files ) ) {
                Carp::confess( "Attempt to get fastq from fastq via archive path failed" );
            }
        }
        elsif ( @fastq_files = eval{$inst_data->dump_sanger_fastq_files;} ) {
            #fastq from sff files .. 454
            if ( not @fastq_files ) {
                Carp::confess( "Did not get fastq files from inst data dump_sanger_fastq_file method" );
            }
            $self->status_message( "Got fastq from sff files" ); #better message?
            if ( not $self->append_fastq_to_orig_fastq_file( @fastq_files ) ) {
                Carp::confess( "Attempt to append to original fastq with with new fastq failed" );
            }
        }
        else {
            $self->status_message($inst_data->__display_name__.' does not have any sequences! Skipping!');
        }
    }
    return 1;
}

sub append_fastq_to_orig_fastq_file {
    my ( $self, @fastq_files ) = @_;

    my $writer = $self->original_fastq_writer;
    unless ( $writer ) {
        Carp::confess( "Failed to get fastq writer to write original fastq file" );
    }
    my $reader = Genome::Model::Tools::Sx::Reader->create(config => [ map { $_.':type=sanger' } @fastq_files ]);
    if ( not $reader ) {
        Carp::confess( "Did not get fastq reader for files: ".join(' ', @fastq_files) );
    }
    while ( my $fastqs = $reader->read ) {
        for my $fastq( @$fastqs ) {
            $writer->write( [$fastq] );
        }
    }
    $self->status_message( "Finished writing to original fastq file, files: ".join(' ', @fastq_files) );

    return 1;
}

#< Remove Chimeras >#
sub detect_and_remove_chimeras {
    my $self = shift;
    $self->status_message('Detect and remove chimeras...');

    if ( not $self->processing_profile->chimera_detector ) {
        $self->error_message('Tried to detect and remove chimeras, but there is not a chimera detector on the processing profile!');
        return;
    }

    my $amplicons_classified = $self->amplicons_classified;
    if ( not defined $amplicons_classified ) {
        $self->error_message('Cannot deteect and remove chimeras because "amplicons classified" metric is not set!');
        return;
    }

    if ( $amplicons_classified == 0 ) {
        $self->status_message("No amplicons were successfully classified, skipping detect and remove chimeras.");
        return 1;
    }

    my @amplicon_sets = $self->amplicon_sets;
    if ( not @amplicon_sets ) {
        $self->error_message('No amplicon sets for '.$self->description);
        return;
    }

    my $amplicons_processed = $self->amplicons_processed;
    if (! defined($amplicons_processed)) {
        $self->error_message("No value for amplicons processed.  Cannot remove chimeras!");
        return;
    } elsif ($amplicons_processed == 0) {
        $self->status_message('amplicons processed is 0.  Skipping chimera removal.');
        return 1;
    }

    $self->status_message('Chimera detector: '.$self->processing_profile->chimera_detector);
    $self->status_message('Chimera detector params: '.$self->processing_profile->chimera_detector_params);
    my $detector = $self->processing_profile->chimera_detector;
    my $chimera_reader_class = 'Genome::Model::Tools::'.Genome::Utility::Text::string_to_camel_case(join(' ', split('-', $detector))).'::Reader';

    my %metrics = ( input => 0, output => 0, );
    for my $amplicon_set ( @amplicon_sets ) {
        my $fasta_file = $amplicon_set->oriented_fasta_file;
        next if not -s $fasta_file; # ok, not error
        $self->status_message('Amplicon set'.$amplicon_set->name) if $amplicon_set->name;

        # DETECT
        $self->status_message('Detect chimeras...');
        my $fasta_base_name = File::Basename::basename($fasta_file);
        my $sequences = $amplicon_set->chimera_dir.'/'.$fasta_base_name;
        symlink($fasta_file, $sequences);
        if ( not -s $sequences ) {
            $self->error_message("Failed to link oriented fasta ($fasta_file) to chimera dir!");
            return;
        }
        my $chimera_file = $amplicon_set->chimera_file;
        my $cmd = "gmt $detector detect-chimeras --sequences $sequences --chimeras $chimera_file ".$self->processing_profile->chimera_detector_params;
        $self->status_message('Detect chimeras command: '.$cmd);
        my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
        if ( not $rv ) {
            $self->error_message($@) if $@;
            $self->error_message('Failed to detect chimeras!');
            return;
        }
        $self->status_message('Detect chimeras...OK');

        # REMOVE
        $self->status_message('Remove chimeras...');

        # chimera free seq writer
        my $writer = $amplicon_set->seq_writer_for('chimera_free');
        if ( not $writer ) {
            $self->error_message('Failed to get chimera free seq writer!');
            return;
        }

        # chimera free classification writer
        my $classification_file = $amplicon_set->chimera_free_classification_file;
        my $classification_fh = eval{ Genome::Sys->open_file_for_writing($classification_file); };
        if ( not $classification_fh ) {
            $self->error_message($@) if $@;
            $self->error_message("Failed to open chimera free classifcation file! $classification_file");
            return;
        }

        while ( my $amplicon = $amplicon_set->next_amplicon ) {
            $metrics{input}++;
            next if not $amplicon->{chimera_result} or not defined $amplicon->{chimera_result};
            next if $amplicon->{chimera_result}->{verdict} eq 'YES';
            $metrics{output}++;
            $writer->write($amplicon->{seq});
            $classification_fh->print( $amplicon->{classification_line}."\n" );
        }
        $classification_fh->close;

        $self->status_message('Remove chimeras...OK');
    }

    if ( $amplicons_classified != $metrics{input} ) {
        $self->error_message("Amplicons oriented ($amplicons_classified) and amplicons put through chimera detection ($metrics{input}) do not match!");
        return;
    }

    my $amplicons_chimeric = $metrics{input} - $metrics{output};
    $self->amplicons_chimeric($amplicons_chimeric);
    $self->amplicons_chimeric_percent( sprintf('%.2f', $amplicons_chimeric / $amplicons_classified) );

    $self->status_message('Classifed and oriented: '.$self->amplicons_classified);
    $self->status_message('Chimeric:               '.$self->amplicons_chimeric);
    $self->status_message('Chimeric percent:       '.($self->amplicons_chimeric_percent * 100).'%');

    $self->status_message('Detect and remove chimeras...OK');
    return 1;
}
#<>#

#< Orient >#
sub orient_amplicons {
    my $self = shift;
    $self->status_message('Orient amplicons...');

    my $amplicons_classified = $self->amplicons_classified;
    if ( not defined $amplicons_classified ) {
        $self->error_message('Cannot orient amplicons because "amplicons classified" metric is not set!');
        return;
    }

    if ( $amplicons_classified == 0 ) {
        $self->status_message("No amplicons were successfully classified, skipping orient.");
        return 1;
    }

    my @amplicon_sets = $self->amplicon_sets;
    if ( not @amplicon_sets ) {
        $self->error_message('No amplicon sets for '.$self->description);
        return;
    }

    my $no_classification = 0;
    for my $amplicon_set ( @amplicon_sets ) {
        next if not $amplicon_set->amplicon_iterator;
        my $writer = $amplicon_set->seq_writer_for('oriented');
        return if not $writer;

        while ( my $amplicon = $amplicon_set->next_amplicon ) {
            my $seq = $amplicon->{seq};
            next if not $seq; #OK - for now...

            if ( not defined $amplicon->{classification} ) {
                $no_classification++;
                next;
            }

            if ( $amplicon->{classification}->[1] eq '-' ) {
                $seq->{seq} = reverse $seq->{seq};
                $seq->{seq} =~ tr/ATGCatgc/TACGtacg/;
            }

            $writer->write($seq);
        }
    }

    my $classification_error = $self->amplicons_classification_error;
    if ( $no_classification != $classification_error ) {
        $self->error_message("Found $no_classification amplicons without classifications, but $classification_error amplicons failed to classify.");
    }

    $self->status_message('Orient amplicons...OK');
    return 1;
}

#< Classify >#
sub classify_amplicons {
    my $self = shift;
    $self->status_message('Classify amplicons...');

    my $attempted = $self->amplicons_attempted;
    if ( not defined $attempted ) {
        $self->error_message('No value for amplicons attempted set. Cannot classify.');
        return;
    }

    my @amplicon_sets = $self->amplicon_sets;
    if ( not @amplicon_sets ) {
        $self->error_message('No amplicon sets for '.$self->description);
        return;
    }

    $self->status_message('Classifier: '.$self->classifier);
    my $classifier_params = $self->processing_profile->classifier_params;
    $self->status_message('Classifier params: '.$classifier_params);

    my %metrics;
    @metrics{qw/ attempted success error total /} = (qw/ 0 0 0 0 /);
    for my $amplicon_set ( @amplicon_sets ) {
        my $fasta_file = $amplicon_set->processed_fasta_file;
        next if not $fasta_file or not -s $fasta_file;

        my $classification_file = $amplicon_set->classification_file;
        unlink $classification_file if -e $classification_file;

        # FIXME use $classifier
        my $cmd = "gmt metagenomic-classifier rdp --input-file $fasta_file --output-file $classification_file $classifier_params --metrics"; 
        my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
        if ( not $rv ) {
            $self->error_message('Failed to execute classifier command');
            return;
        }

        # metrics
        my $metrics_file = $classification_file.'.metrics';
        my $metrics_fh = eval{ Genome::Sys->open_file_for_reading($metrics_file); };
        if ( not $metrics_fh ) {
            $self->error_message("Failed to open metrics file ($metrics_file): $@");
            return;
        }
        while ( my $line = $metrics_fh->getline ) {
            chomp $line;
            my ($key, $val) = split('=', $line);
            $metrics{$key} += $val;
        }
    }

    $self->amplicons_processed($metrics{total});
    $self->amplicons_processed_success( 
        defined $attempted and $attempted > 0 ?  sprintf('%.2f', $metrics{total} / $attempted) : 0 
    );
    $self->amplicons_classified($metrics{success});
    $self->amplicons_classified_success( 
        $metrics{total} > 0 ?  sprintf('%.2f', $metrics{success} / $metrics{total}) : 0
    );
    $self->amplicons_classification_error($metrics{error});

    $self->status_message('Processed:  '.$self->amplicons_processed);
    $self->status_message('Classified: '.$self->amplicons_classified);
    $self->status_message('Error:      '.$self->amplicons_classification_error);
    $self->status_message('Success:    '.($self->amplicons_classified_success * 100).'%');

    $self->status_message('Classify amplicons...OK');

    return 1;
}

#< Reports >#
sub generate_reports {
    my $self = shift;
    
    my @report_names = (qw/ summary /);
    push @report_names, 'composition' if not $self->model->is_for_qc;

    for my $report_name ( @report_names ) {
        $self->_generate_and_save_report($report_name);
    }

    return 1;
}

sub _generate_and_save_report {
    my ($self, $name) = @_;

    $self->status_message(ucfirst($name).' report...');

    my $class = 'Genome::Model::MetagenomicComposition16s::Report::'.Genome::Utility::Text::string_to_camel_case($name);
    my $generator = $class->create(
        build_id => $self->id,
    );
    unless ( $generator ) {
        $self->error_message("Could not create $name report generator for ".$self->description);
        return;
    }
    my $report = $generator->generate_report;
    unless ( $report ) {
        $self->error_message("Could not generate $name report for ".$self->description);
        return;
    }

    unless ( $self->add_report($report) ) {
        $self->error_message("Could not save $name report for ".$self->description);
    }

    my @datasets = $report->get_datasets;
    unless ( @datasets ) { # not ok
        $self->error_message("No datasets in $name report for ".$self->description);
        return;
    }
    my $file_base = sprintf(
        '%s/%s',
        $self->reports_directory,
        $report->name_to_subdirectory($report->name),
    );

    for my $dataset ( @datasets ) {
        my $dataset_name = $dataset->name;
        my $file = sprintf(
            '%s/%s.%s.tsv',
            $file_base,
            $self->model->subject_name,
            $dataset_name,
        );
        unlink $file if -e $file;
        my $fh = Genome::Sys->open_file_for_writing($file)
            or return; # not ok
        my ($svs) = $dataset->to_separated_value_string(separator => "\t");
        unless ( $svs ) { # not ok
            $self->error_message("Could not get dataset ($dataset) for $name report for ".$self->description);
            return;
        }
        $fh->print($svs);
        $fh->close;
    }

    my $xsl_file = $generator->get_xsl_file_for_html;
    if ( -e $xsl_file ) {
        my $xslt = Genome::Report::XSLT->transform_report(
            report => $report,
            xslt_file => $xsl_file,
        );
        unless ( $xslt ) {
            $self->error_message("Can't transform report to html.");
            return;
        }
        # Save
        my $html_file = $report->directory.'/report.html';
        my $fh = Genome::Sys->open_file_for_writing($html_file);
        unless ( $fh ) {
        }
        $fh->print( $xslt->{content} );
    }

    $self->status_message(ucfirst($name).' report...OK');
    
    return $report;
}
#<>#

#< QC Email >#
sub perform_post_success_actions {
    my $self = shift;
    $self->status_message('Post success actions');

    $self->status_message('Check if model is for QC: '.$self->model_name);
    if ( not $self->model->is_for_qc ) {
        $self->status_message('Model is not for QC. Not sending confirmation email');
        return 1;
    }
    $self->status_message('This model is QC');

    # Make sure there is only one run and region represented
    my @instrument_data = $self->instrument_data;
    my %run_region_names = map { $_->run_name.' '.$_->region_number => 1 } @instrument_data;
    my @run_region_names = grep { defined } keys %run_region_names;
    if ( not @run_region_names ) {
        $self->error_message('No run names from instrument data for '.$self->description);
        return;
    }
    if ( @run_region_names > 1 ) {
        $self->error_message("Found multiple run region names (@run_region_names) for instrument data assigned to ".$self->description);
        return;
    }

    # Get the 454 instrument data for the run and region number
    my @instrument_data_for_run_region = Genome::InstrumentData::454->get(
        run_name => $instrument_data[0]->run_name,
        region_number => $instrument_data[0]->region_number,
    );
    if ( not @instrument_data_for_run_region ) {
        $self->error_message('No instrument data found for run region name: '.$run_region_names[0]);
        return;
    }

    # Filter out pooled and negative control samples
    my $expected_cnt = 0;
    for my $instrument_data_for_run_region ( @instrument_data_for_run_region ) {
        next if $instrument_data_for_run_region->sample->name =~ /^pool/i;
        next if $instrument_data_for_run_region->sample->name =~ /^n\-cn?trl$/i;
        $expected_cnt++;
    }
    if ( @instrument_data != $expected_cnt ) {
        $self->status_message('Not sending email for MC16s QC model. Have '.@instrument_data." instrument data, but expected $expected_cnt\.");
        return 1;
    }

    my $msg = "Hello,\n\nThis MC16s QC build is finished and has all instrument data included.\n\n";
    $msg .= 'Model name:      '.$self->model_name."\n";
    $msg .= 'Model id:        '.$self->model_id."\n";
    $msg .= 'Build id:        '.$self->id."\n";
    $msg .= 'Directory:       '.$self->data_directory."\n";
    $msg .= 'Run name:        '.$instrument_data[0]->run_name."\n";
    $msg .= 'Region number:   '.$instrument_data[0]->region_number."\n";
    $msg .= 'Inluded count:   '.@instrument_data."\n";
    $msg .= 'Expected count:  '.$expected_cnt."\n";
    $msg .= 'Attempted:       '.$self->amplicons_attempted."\n";
    $msg .= 'Processed:       '.$self->amplicons_processed."\n";
    $msg .= 'Success :        '.(100 * $self->amplicons_processed_success)."%\n";
    $msg .= "\n-APIPE";
    $self->status_message($msg);

    if ( not $ENV{UR_DBI_NO_COMMIT} ) { # do not send mail when in dev mode
        Mail::Sendmail::sendmail(
            To => 'esodergr@genome.wustl.edu, kmihindu@genome.wustl.edu', 
            #To => 'ebelter@genome.wustl.edu', 
            Cc => 'ebelter@genome.wustl.edu', 
            From => 'apipe@genome.wustl.edu', 
            Subject => 'MC16s QC Build is Done',
            Message => $msg,
        );
    }

    $self->status_message('Sent email to Erica (esodergren) and Kathie (kmihindu)');

    return 1;
}

#< calculate est kb usage >#
sub calculate_estimated_kb_usage {
    my $self = shift;

    my @instrument_data = $self->instrument_data;
    if ( not @instrument_data ) {
        Carp::confess( "No instrument data found for ".$self->description );
    }

    my $est_kb_usage = 0;
    for my $instrument_data ( @instrument_data ) {
        my $sequencing_platform = $instrument_data->sequencing_platform;
        switch ($sequencing_platform) {
               case '454'    { $est_kb_usage += $instrument_data->read_count * 5 * 1024 }
               case 'sanger' { $est_kb_usage += 30_000 }
               case 'solexa' { $est_kb_usage += 500_000 } # FIXME update once implemented
               else          { Carp::confess('Unknown sequencing platform! '.$sequencing_platform) }
           }
    }

    return ( $est_kb_usage >= 1024 ? $est_kb_usage : 1024 );
}

#< Diff >#
sub dirs_ignored_by_diff {
    return (qw{
        logs/
        reports/
        edit_dir/
        chromat_dir/
        classification/
    });
}

sub files_ignored_by_diff {
    return (qw/ build.xml /);
}

sub regex_for_custom_diff {
    return (
        gz => '\.gz$',
        rdp => '\.rdp1-[12]$',
        metrics => 'metrics\.processed\.(in|out)\.txt',
    );
}

sub diff_rdp {
    my ($self, $file1, $file2) = @_;

    my $reader1 = Genome::Model::Tools::MetagenomicClassifier::ClassificationReader->create(
        file => $file1,
    );
    return if not $reader1;

    my $reader2 = Genome::Model::Tools::MetagenomicClassifier::ClassificationReader->create(
        file => $file2,
    );
    return if not $reader2;

    my ($classification1_cnt, $classification2_cnt) = (qw/ 0 0 /);
    while ( my $classification1 = $reader1->read ) {
        $classification1_cnt++;
        my $classification2 = $reader2->read;
        last if not $classification2;
        $classification2_cnt++;
        if ( $classification1->{id} ne $classification2->{id} ) {
            $self->status_message("RDP differs at id: ".$classification1->{id}.' <=> '.$classification2->{id});
            return;
        }
        if ( $classification1->{complemented} ne $classification2->{complemented} ) {
            $self->status_message("RDP differs at complemented: ".$classification1->{complemented}.' <=> '.$classification2->{complemented});
            return;
        }
        for my $rank (qw/ domain phylum order class family genus /) {
            if ( $classification1->{$rank}->{id} ne $classification2->{$rank}->{id} ) {
                $self->status_message("RDP differs at $rank: ".$classification1->{$rank}->{id}.' <=> '.$classification2->{$rank}->{id});
                return;
            }
        }
    }

    if ( $classification1_cnt != $classification2_cnt ) {
        $self->error_message('Classification counts differ: '.$classification1_cnt.' <=> '.$classification2_cnt);
        return;
    }

    return 1;
}

sub diff_metrics {
    my ($self, $file1, $file2) = @_;

    my $metrics_from_file1 = Genome::Model::Tools::Sx::Metrics->from_file($file1);
    return if not $metrics_from_file1;

    my $metrics_from_file2 = Genome::Model::Tools::Sx::Metrics->from_file($file2);
    return if not $metrics_from_file2;

    for my $metric_name (qw/ bases count /) {
        if ( $metrics_from_file1->$metric_name ne $metrics_from_file2->$metric_name ) {
            $self->status_message("Metrics differ for $metric_name: ".$metrics_from_file1->$metric_name.' <=> '.$metrics_from_file2->$metric_name);
            return;
        }
    }

    return 1;
}
#<>#

1;

