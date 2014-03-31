package Genome::Model::Build::MetagenomicComposition16s::ProcessSangerInstrumentData; 

use strict;
use warnings;

use Genome;

class Genome::Model::Build::MetagenomicComposition16s::ProcessSangerInstrumentData {
    is => 'Command::V2',
    has_input => [
        input_build => {
            is => 'Genome::Model::Build::MetagenomicComposition16s',
        },
    ],
    has_output => [
        build => {
            is => 'Genome::Model::Build::MetagenomicComposition16s',
            via => '__self__',
            to => 'input_build',
        },
    ],
    has_constant => [
        amplicon_size => {
            calculate_from => [qw/ build /],
            calculate => q| 
                my $string = $build->processing_profile->amplicon_processor;
                return unless $string; #ok
                my ($amplicon_size) = $string =~ /filter\s+by-min-length\s+--?length\s+(\d+)/;
                return return $amplicon_size;
                |,
        },
        chromat_dir => {
            calculate => q| return $self->build->data_directory.'/chromat_dir'; |,
        },
        edit_dir => {
            calculate => q| return $self->build->data_directory.'/edit_dir'; |,
        },
    ],
    has_optional_transient => [
        _raw_reads_fasta_and_qual_writer => { },
        _processed_reads_fasta_and_qual_writer => { },
    ],
};

sub execute {
    my $self = shift;

    for my $subdir_name (qw/ chromat_dir edit_dir /) {
        Genome::Sys->create_directory($self->$subdir_name);
    }

    my @amplicon_sets = $self->build->amplicon_sets;
    return if not @amplicon_sets;
    if ( @amplicon_sets != 1 ) {
        $self->error_message('Expected only 1 amplicon set to process sanger inst data!');
        return;
    }

    my $link_ok = $self->_dump_and_link_instrument_data;
    return if not $link_ok;

    my $raw_reads_writer = $self->_create_raw_reads_writer;
    return if not $raw_reads_writer;

    my $processed_reads_writer = $self->_create_processed_reads_writer;
    return if not $processed_reads_writer;

    my $assembler_params_string = $self->build->processing_profile->assembler_params;
    my %assembler_params = Genome::Utility::Text::param_string_to_hash($assembler_params_string);
    unless ( %assembler_params ) {
        $self->error_message("Malformed assembler params: $assembler_params_string");
        return;
    }

    my $writer = $amplicon_sets[0]->seq_writer_for('processed');
    return if not $writer;

    my $amplicon_iterator = $self->_amplicon_iterator;
    return if not $amplicon_iterator;

    my ($attempted, $processed, $reads_attempted, $reads_processed) = (qw/ 0 0 0 /);
    while ( my $amplicon = $amplicon_iterator->() ) {
        $attempted++;
        $reads_attempted += @{$amplicon->{reads}};
        my $prepare_ok = $self->_prepare($amplicon);
        return if not $prepare_ok;

        my $trim_ok = $self->_trim($amplicon);
        return if not $trim_ok;

        my $assemble_ok = $self->_assemble($amplicon, %assembler_params);
        return if not $assemble_ok;

        $self->_clean_up($amplicon);

        $self->load_seq_for_amplicon($amplicon)
            or next; # ok
        $writer->write($amplicon->{seq});
        $processed++;
        $reads_processed += $amplicon->{reads_processed};
    }

    $self->_raw_reads_fasta_and_qual_writer->flush;
    $self->_processed_reads_fasta_and_qual_writer->flush;

    $self->build->amplicons_attempted($attempted);
    $self->build->amplicons_processed($processed);
    $self->build->amplicons_processed_success( $attempted > 0 ?  sprintf('%.2f', $processed / $attempted) : 0 );
    $self->build->reads_attempted($reads_attempted);
    $self->build->reads_processed($reads_processed);
    $self->build->reads_processed_success( $reads_attempted > 0 ?  sprintf('%.2f', $reads_processed / $reads_attempted) : 0 );

    return 1;
}

sub _amplicon_iterator {
    my $self = shift;

    # open chromat_dir
    my $dh = Genome::Sys->open_directory( $self->chromat_dir );
    unless ( $dh ) {
        $self->error_message("Can't open chromat dir to get reads. See above error.");
        return;
    }

    my @all_read_names;
    while ( my $read_name = $dh->read ) {
        next if ($read_name eq '.' || $read_name eq '..');
        $read_name =~ s/\.gz$//;
        push @all_read_names, $read_name;
    }
    # make sure we got some 
    unless ( @all_read_names ) {
        $self->error_message( "No reads found in chromat dir! ".$self->chromat_dir);
        return;
    }
    #sort
    @all_read_names = sort { $a cmp $b } @all_read_names;

    my $amplicon_name_reqgexp = qr/^(.+)\.[bg]\d+$/;
    my $pos = 0;
    return sub{
        AMPLICON: while ( $pos < $#all_read_names ) {
            # Get amplicon name
            $all_read_names[$pos] =~ $amplicon_name_reqgexp;
            my $amplicon_name = $1;
            unless ( $amplicon_name ) {
                Carp::confess('Could not determine amplicon name for read: '.$all_read_names[$pos]);
            }
            # Start reads list
            my @read_names = ( $all_read_names[$pos] );
            READS: while ( $pos < $#all_read_names ) {
                # incremnent
                $pos++;
                # Get amplicon name
                $all_read_names[$pos] =~ $amplicon_name_reqgexp;
                my $read_amplicon_name = $1;
                unless ( $read_amplicon_name ) {
                    Carp::confess sprintf(
                        'Could not determine amplicon name from read name (%s) for build (%s)',
                        $all_read_names[$pos],
                        $self->build->id,
                    );
                }
                unless ( $read_amplicon_name eq $amplicon_name ) { 
                    # go on to filtering
                    last READS; 
                }
                push @read_names, $all_read_names[$pos]; # add read
            }

            # Create amplicon object
            my $amplicon = {
                name => $amplicon_name,
                reads => \@read_names,
            };

            # Processed oseq
            $self->load_seq_for_amplicon($amplicon); # dies on error

            return $amplicon;
        }
    };
}

sub _create_raw_reads_writer {
    my $self = shift;

    my $fasta_file = $self->build->fasta_dir.'/'.$self->build->file_base_name.'.reads.raw.fasta';
    my $qual_file = $fasta_file.'.qual';
    my $writer = Genome::Model::Tools::Sx::PhredWriter->create(
        file => $fasta_file,
        qual_file => $qual_file,
    );
    if ( not $writer ) {
        $self->error_message('Failed to create raw reads writer');
        return;
    my $amplicon_name_reqgexp = qr/^(.+)\.[bg]\d+$/;
    }

    return $self->_raw_reads_fasta_and_qual_writer($writer);
}

sub _create_processed_reads_writer {
    my $self = shift;

    my $fasta_file = $self->build->fasta_dir.'/'.$self->build->file_base_name.'.reads.processed.fasta';
    my $qual_file = $fasta_file.'.qual';
    my $writer = Genome::Model::Tools::Sx::PhredWriter->create(
        file => $fasta_file,
        qual_file => $qual_file,
    );
    if ( not $writer ) {
    my $amplicon_name_reqgexp = qr/^(.+)\.[bg]\d+$/;
        $self->error_message('Failed to create processed reads writer');
        return;
    }

    return $self->_processed_reads_fasta_and_qual_writer($writer);
}

sub _dump_and_link_instrument_data {
    my $self = shift;

    my @instrument_data = $self->build->instrument_data;
    unless ( @instrument_data ) { # should not happen
        $self->error_message('No instrument data found for '.$self->build->description);
        return;
    }

    my $chromat_dir = $self->chromat_dir;
    for my $instrument_data ( @instrument_data ) {
        my $full_path = $instrument_data->full_path;
        if ( not $full_path or not -s $full_path ) {
            $self->error_message(
                'Instrument data full path does not exist! Pleasee contact support for assistance. '.$instrument_data->id,
            );
            return;
        }

        # link
        my $instrument_data_dir = $instrument_data->full_path;
        my $dh = Genome::Sys->open_directory($instrument_data_dir);
        return if not $dh;

        my $cnt = 0;
        while ( my $trace = $dh->read ) {
            next if ($trace eq '.' || $trace eq '..');
            $cnt++;
            my $target = sprintf('%s/%s', $instrument_data_dir, $trace);
            my $link = sprintf('%s/%s', $chromat_dir, $trace);
            next if -e $link; # link points to a target that exists
            unlink $link if -l $link; # remove - link exists, but points to something that does not exist
            Genome::Sys->create_symlink($target, $link)
                or return;
        }

        unless ( $cnt ) {
            $self->error_message("No traces found in instrument data directory ($instrument_data_dir)");
            return;
        }
    }

    return 1;
}

sub _prepare {
    #< Scf to Fasta via Phred >#
    my ($self, $amplicon) = @_;

    # scfs file
    my $scfs_file = $self->edit_dir.'/'.$amplicon->{name}.'.scfs';
    unlink $scfs_file;
    my $scfs_fh = Genome::Sys->open_file_for_writing($scfs_file);
    return if not $scfs_fh;
    for my $scf ( @{$amplicon->{reads}} ) { 
        $scfs_fh->print($self->chromat_dir."/$scf.gz\n");
    }
    $scfs_fh->close;

    # scf 2 fasta
    my $fasta_file = $self->edit_dir.'/'.$amplicon->{name}.'.fasta';
    unlink $fasta_file;
    my $qual_file =  $fasta_file.'.qual';
    unlink $qual_file;
    my $command = sprintf(
        'phred -if %s -sa %s -qa %s -nocall -zt /tmp',
        $scfs_file,
        $fasta_file,
        $qual_file,
    );

    my $rv = eval{ Genome::Sys->shellcmd(cmd => $command); };
    if ( not $rv ) {
        $self->error_message('Failed to convert '.$amplicon->{name}.' SCFs to FASTA: '.$@);
        return;
    }

    # write the 'raw' read fastas
    my $reader = Genome::Model::Tools::Sx::PhredReader->create(
        file => $fasta_file, 
        qual_file => $qual_file,
    );
    return if not $reader;
    my $writer = $self->_raw_reads_fasta_and_qual_writer;
    while ( my $seq = $reader->read ) {
        $writer->write($seq) or return;
    }
    
    return 1;
}

sub _trim {
    my ($self, $amplicon) = @_;

    my $fasta_file = $self->edit_dir.'/'.$amplicon->{name}.'.fasta';
    return unless -s $fasta_file; # ok

    my $trim3 = Genome::Model::Tools::Fasta::Trim::Trim3->create(
        fasta_file => $fasta_file,
        min_trim_quality => 10,
        min_trim_length => 100,
    );
    unless ( $trim3 ) { # not ok
        $self->error_message("Can't create trim3 command for amplicon: ".$amplicon->name);
        return;
    }
    $trim3->execute; # ok

    next unless -s $fasta_file; # ok

    my $screen = Genome::Model::Tools::Fasta::ScreenVector->create(
        fasta_file => $fasta_file,
    );
    unless ( $screen ) { # not ok
        $self->error_message("Can't create screen vector command for amplicon: ".$amplicon->name);
        return;
    }
    $screen->execute; # ok

    next unless -s $fasta_file; # ok

    my $reader = Genome::Model::Tools::Sx::PhredReader->create(
        file => $fasta_file,
        qual_file => $fasta_file.'.qual',
    );
    return if not $reader;
    while ( my $seq = $reader->read ) {
        $self->_processed_reads_fasta_and_qual_writer->write($seq) or return;
    }
 
    return 1;
}

sub _assemble {
    my ($self, $amplicon, %assembler_params) = @_;

    my $fasta_file = $self->edit_dir.'/'.$amplicon->{name}.'.fasta';
    next unless -s $fasta_file; # ok

    my $phrap = Genome::Model::Tools::PhredPhrap::Fasta->create(
        fasta_file => $fasta_file,
        %assembler_params,
    );
    unless ( $phrap ) { # bad
        $self->error_message(
            "Can't create phred phrap command for build's (".$self->id.") amplicon (".$amplicon->{name}.")"
        );
        return;
    }
    $phrap->dump_status_messages(1);
    $phrap->execute; # no check

    return 1;
}

sub _clean_up {
    my ($self, $amplicon) = @_;

    for my $ext (qw/
        fasta.contigs fasta.contigs.qual 
        fasta.log fasta.singlets
        fasta.phrap.out fasta.memlog
        fasta.problems fasta.problems.qual
        fasta.preclip fasta.qual.preclip 
        fasta.prescreen fasta.qual.prescreen
        scfs
        /) {
        my $file = sprintf('%s/%s.%s', $self->edit_dir, $amplicon->{name}, $ext);
        unlink $file if -e $file;
    }

    return 1;
}

sub load_seq_for_amplicon {
    my ($self, $amplicon) = @_;

    die "No amplicon to load seq." unless $amplicon;

    # get contig from acefile
    my $acefile = $self->edit_dir.'/'.$amplicon->{name}.'.fasta.ace';
    return unless -s $acefile; # ok

    my $ace = Genome::Model::Tools::Consed::AceReader->create(file => $acefile);
    my $contig;
    my $amplicon_size = ( $self->amplicon_size ) ? $self->amplicon_size : 0;
    while ( my $obj = $ace->next ) {
        next if $obj->{type} ne 'contig';
        next unless $obj->{read_count} > 1;
        $obj->{unpadded_consensus} = $obj->{consensus};
        $obj->{unpadded_consensus} =~ s/\*//g;
        next if length $obj->{unpadded_consensus} < $amplicon_size;
        $contig = $obj;
        last;
    }
    return unless $contig; # ok

    my $seq = {
        id => $amplicon->{name},
        seq => $contig->{unpadded_consensus},
        qual => join('', map { chr($_ + 33) } @{$contig->{base_qualities}}),
    };
    if ( $seq->{seq} !~ /^[ATGCNX]+$/i ) {
        Carp::confess('Illegal characters in sequence for amplicon: '.$amplicon->{name}."\n".$seq->{seq}); 
    }

    if ( length $seq->{seq} !=  length $seq->{qual} ) {
        Carp::confess('Unequal lengths of sequence and quality for amplicon: '.$amplicon->{name}."\n".$seq->{seq}."\n".$seq->{qual});
    }

    $amplicon->{seq} = $seq;
    $amplicon->{reads_processed} = $contig->{read_count};

    return $seq;
}

1;

