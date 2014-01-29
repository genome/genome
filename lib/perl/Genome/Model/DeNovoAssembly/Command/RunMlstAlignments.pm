package Genome::Model::DeNovoAssembly::Command::RunMlstAlignments;

use strict;
use warnings;

use Genome;

class Genome::Model::DeNovoAssembly::Command::RunMlstAlignments {
    is => 'Genome::Command::Base',
    has => [
        sample_names => {
            is => 'Text',
            is_many => 1,
            doc => 'Instrument data sample name to run alignment on latest succeeded de-novo velvet build',
            is_optional => 1,
        },
        sample_list => {
            is => 'Text',
            doc => 'File of sample names to run alignments for',
            is_optional => 1,
        },
        query_build_ids => {
            is => 'Number',
            is_many => 1,
            doc => 'DeNovo assembly build id',
            is_optional => 1,
        },
        reference_build_id => {
            is => 'Number',
            doc => 'Imported reference sequence build id',
        },
        output_dir => {
            is => 'Text',
            doc => 'Directory to output alignments to',
        },
    ],
};

sub help_brief {
    'Command to run Nucmer MLST alignments on de-novo assemblies'
}

sub help_detail {
    return <<"EOS"
genome model de-novo-assembly run-mlst-alignments --sample-names MRSA-6461-N121,MRSA-6365-A50 --reference-build-id 121034788 --output-dir /gscmnt/100/MLSTs
genome model de-novo-assembly run-mlst-alignments --sample-list /gscmnt/100/sample.txt --reference-build-id 121034788 --output-dir /gscmnt/100/MLSTs
genome model de-novo-assembly run-mlst-alignments --query-build-ids 121482147,121483135 --reference-build-id 121034788 --output-dir /gscmnt/100/MLSTs
EOS
}

sub execute {
    my $self = shift;

    # validate output dir
    $self->debug_message('Invalid output dir name: '.$self->output_dir) and return if
        not -d $self->output_dir;

    # verify sample names and get builds, ids
    my ( $query_builds, $ref_seq_build ) = $self->_validate_samples_and_get_builds;

    my @output_files;

    for my $build ( @$query_builds ) {
        # run nucmer
        my $nucmer_out_file = $self->output_dir.'/'.$build->subject_name.'.delta';
        unlink $nucmer_out_file;
        my $nucmer = Genome::Model::Tools::Mummer::Nucmer->create(
            prefix    => $self->output_dir.'/'.$build->subject_name,
            query     => $build->contigs_bases_file,
            reference => $ref_seq_build->full_consensus_path('fa'),
        );
        if ( not $nucmer ) {
            $self->error_message('Failed to create nucmer tool for sample: '.$build->subject_name);
            return;
        }
        if ( not $nucmer->execute ) {
            $self->error_message('Failed to run nucmer for sample: '.$build->subject_name);
            return;
        }
        $self->debug_message('Successfully ran nucmer on sample: '.$build->subject_name);

        # run show-coords
        my $output_file = $self->output_dir.'/'.$build->subject_name.'.alignments.txt';
        push @output_files, $output_file;
        my $show_coords = Genome::Model::Tools::Mummer::ShowCoords->create(
            input_delta_file => $nucmer_out_file,
            output_file      => $output_file,
            r                => 1,
            c                => 1,
            l                => 1,
            T                => 1,
        );
        if ( not $show_coords ) {
            $self->error_message('Failed to create show-coords tool');
            return;
        }
        if ( not $show_coords->execute ) {
            $self->error_message('Failed to execute show-coords');
            return;
        }
        $self->debug_message('Successfully ran show-coords for sample: '.$build->subject_name);
    }

    # summarize output files
    if ( not $self->summarize_outputs(@output_files) ) {
        $self->debug_message('Failed to summarize alignment results');
        return;
    }

    return 1;
}

sub summarize_outputs {
    my ($self, @files) = @_;

    my $summary_file = $self->output_dir.'/alignments_summary.txt';
    unlink $summary_file if -e $summary_file;
    my $summary_fh = Genome::Sys->open_file_for_writing($summary_file);

    my @headers = qw/ Subject Locus Allel Qlength Rlength Identity Coverage /;
    $summary_fh->printf("%-35s%-12s%-12s%-12s%-12s%-12s%-12s\n", @headers );

    for my $file (@files ) {
        my $file_name = File::Basename::basename($file);
        my ($subject) = $file_name =~ /(\S+)\.alignments\.txt$/;
        my $align_fh = Genome::Sys->open_file_for_reading($file);
        while( my $line = $align_fh->getline ) {
            chomp $line;
            if ( $line =~ /^\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\.\d+/ ) {
                my @tmp = split( "\t", $line );
                # tmp[5] = query length
                # tmp[6] = identity
                # tmp[7] = ref length
                # tmp[9] = coverage
                # tmp[11] = locus_allel, eg, adk_1
                  # locus = adk
                  # allel = 1
                my($locus,$allel) = $tmp[11] =~ /^(\S+)_(\d+)$/;
                $locus = ( defined $locus ) ? $locus : 'undef';
                $allel = ( defined $allel ) ? $allel : 'undef';

                my @fields = ( $subject, $locus, $allel, $tmp[5], $tmp[7], $tmp[6], $tmp[9] );
                $summary_fh->printf("%-35s%-12s%-12s%-12s%-12s%-12s%-12s\n", @fields );
            }
        }
        $align_fh->close;
    }
    $summary_fh->close;
    return 1;
}

sub _validate_reference_seq_build {
    my $self = shift;

    my $build = Genome::Model::Build->get( $self->reference_build_id );
    $self->error_message('Failed to get build for reference build id: '.$self->reference_build_id) and return if
        not $build;
    $self->error_message('Expected ImportedReferenceSequence build but got: '.$build->subclass_name) and return if
        not $build->subclass_name =~ /ImportedReferenceSequence/;
    $self->error_message('Failed to get ref seq file from build or file is zero size') and return if
        not -s $build->full_consensus_path('fa');

    return $build;
}

sub _validate_samples_and_get_builds {
    my $self = shift;

    if ( ! $self->sample_names and ! $self->sample_list and ! $self->query_build_ids ) {
        $self->error_message('You must supply sample name(s) or query build id(s)');
        return;
    }

    $self->error_message('Failed to validate reference sequence build of id: '.$self->reference_build_id) and return if not
        my $ref_seq_build = $self->_validate_reference_seq_build;

    # get input sample names
    my @sample_names;
    if ( $self->sample_list ) {
        $self->error_message("Invalid fine name or file is zero size: ".$self->sample_list) and return if
            not -s $self->sample_list;
        my $fh = Genome::Sys->open_file_for_reading( $self->sample_list );
        while ( my $line = $fh->getline ) {
            chomp $line;
            push @sample_names, split( /\s+/, $line );
        }
        $fh->close;
    }
    push @sample_names, $self->sample_names if $self->sample_names;

    # verify samples and get builds
    my @valid_builds;
    my $expected_subclass = 'Genome::Model::Build::DeNovoAssembly::Velvet';
    if ( @sample_names ) {
        for my $name ( @sample_names ) {
            my $sample = Genome::Sample->get( name => $name );
            $self->debug_message("Skipping .. can not get genome sample for name: $name") and next if
                not $sample;
            $self->warning_message(
                    'Sample species name: '. $sample->species_name.' does not match ref seq species name: '.$ref_seq_build->species_name.', for sample '.$sample->name
                ) if not $sample->species_name eq $ref_seq_build->species_name;
            my @builds = Genome::Model::Build->get(
                subject_id    => $sample->id,
                subclass_name => $expected_subclass,
                status        => 'Succeeded',
            );
            $self->debug_message("Skipping .. NO Succeeded de-novo velvet build found for sample: $name") and next if
                not @builds;
            push @valid_builds, $builds[-1];#->id;
        }
    }

    # add query build ids if specified
    if ( $self->query_build_ids ) {
        for my $id ( $self->query_build_ids ) {
            my $build = Genome::Model::Build->get( $id );
            $self->debug_message("Skipping .. failed to get build for id: $id") and next if
                not $build;
            $self->debug_message("Skipping .. build $id is not of subclass: $expected_subclass") and next if
                not $build->subclass_name eq $expected_subclass;
            push @valid_builds, $build;
        }
    }

    return \@valid_builds, $ref_seq_build;
}

sub _validate_show_coords_params {
    my $self = shift;
    my %p;
    my @ps = ( split /-/, $self->show_coords_params );
    shift @ps if $ps[0] eq ''; # from split
    for my $param ( @ps ) {
        $param =~ s/\s+$//;
        $param =~ s/^\s+//;
        if ( $param =~ /\d+/ ) {
            # eg -I0.5, -I 0.5 or -I=0.5 or -M 500
            my ( $name, $value ) = $param =~ /^(\S)\s?\=?(\d+\.\d+|\d+)$/;
            if ( not $name and not $value ) {
                $self->debug_message("Failed to get valid name and value from param: $param");
                return;
            }
            $p{$name} = $value;
        } else {
            # boleans eg -rlcT
            my @ps = split( '', $param );
            for my $name ( @ps ) {
                $p{$name} = 1;
            }
        }
    }

    if ( not $self->validate_param_name(%p) ) {
        $self->error_message('Failed to validate show-coords param');
        return;
    }

    return %p;
}

sub validate_param_name {
    my ( $self, %params ) = @_;

    my $class = 'Genome::Model::Tools::Mummer::ShowCoords';
    my $class_meta = $class->__meta__; # dies if not
    if ( not $class_meta ) {
        $self->error_message("No genome class found for $class");
        return;
    }
 
    foreach my $name ( keys %params ) {
        my $property = $class_meta->property_meta_for_name($name);
        if ( not $property ) {
            $self->error_message("Invalid param, $name, specified for $class");
            return;
        }
        # TODO - validate values ?
    } 

    return 1;
}

1;
