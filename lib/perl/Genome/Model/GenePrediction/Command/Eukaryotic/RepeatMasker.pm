package Genome::Model::GenePrediction::Command::Eukaryotic::RepeatMasker;

# There seems to be a Genome::Model::Tools::RepeatMasker command tree that
# could replace this module, with some minor modifications.

use strict;
use warnings;

use Genome;
use Carp 'confess';
use File::Temp 'tempdir';
use File::chdir;
use File::Basename 'fileparse';

class Genome::Model::GenePrediction::Command::Eukaryotic::RepeatMasker {
    is => 'Command',
    has => [
        fasta_file => { 
            is  => 'FilePath',
            is_input => 1,
            doc => 'Fasta file to be masked',
        },
    ],
    has_optional => [
        raw_output_directory => {
            is => 'DirectoryPath',
            doc => 'Directory in which raw output from this tool should be placed',
        },
        version => {
            is => 'Text',
            doc => 'Version of repeat masker to use',
            default => '3.2.9',
            valid_values => [
                '3.0.5',
                '3.0.8',
                '3.1.0',
                '3.1.5',
                '3.1.8',
                '3.2.7',
                '3.2.9',
                '12122000',
                '20020505',
                '20030921',
                '2004March6',
            ],
        },
        masked_fasta => { 
            is => 'FilePath',
            is_input => 1,
            is_output => 1,
            doc => 'Masked sequence is placed in this file (fasta format)' 
        },
        make_ace => {
            is => 'Boolean',
            is_input => 1,
            default => 1,
            doc => 'If set, repeat masker will create an ace file in addition to normal output',
        },
        ace_file_location => {
            is => 'FilePath',
            is_input => 1,
            is_output => 1,
            doc => 'If ace file is generated, it goes here',
        },
        make_gff => {
            is => 'Boolean',
            is_input => 1,
            default => 1,
            doc => 'If set, repeat masker will create a gff file in addition to normal output',
        },
        gff_file_location => {
            is => 'FilePath',
            is_input => 1,
            is_output => 1,
            doc => 'If gff file is generated, it goes here',
        },
        repeat_library => {
            is => 'FilePath',
            is_input => 1,
            doc => 'Repeat library to pass to RepeatMasker', 
        },
	    species	=> {
            is  => 'Text',
            is_input => 1,
            doc => 'Species name',
		},
	    xsmall	=> {
            is  => 'Boolean',
            is_input => 1,
            default => 0,
            doc => 'If set, masked sequence is marked by lowercasing bases instead of using N',
		},
        temp_working_directory => {
            is => 'DirectoryPath',
            is_input => 1,
            doc => 'Temporary working files are written here',
        },
        skip_masking => {
            is => 'Boolean',
            is_input => 1,
            default => 0,
            doc => 'If set, masking is skipped',
        },
        exclude_overly_masked => {
            is => 'Boolean',
            is_input => 1,
            default => 0,
            doc => 'If set, sequences with N content greater than maximum_percent_masked are excluded from output fasta',
        },
        maximum_percent_masked => {
            is => 'Number',
            is_input => 1,
            doc => 'Sequences with an N content greater than this are omitted. Value is a percent',
        },
        overly_masked_sequence_fasta => {
            is => 'FilePath',
            is_input => 1,
            is_output => 1,
            doc => 'If overly masked sequences are excluded from output fasta, those sequences are placed into this file',
        },
    ], 
};

sub help_brief { return "RepeatMask the contents of the input file and write the result to the output file" };
sub help_synopsis { return help_brief() };
sub help_detail {
    return "Masks out repetitive regions of sequences from input fasta and places results in output fasta. " .
        "Can optionally create an ace file detailing masked regions";
}

sub version_to_path_mapping {
    return (
        '3.0.5' => '/gsc/scripts/pkg/bio/repeatmasker/RepeatMasker-open-3-0-5/RepeatMasker',
        '3.0.8' => '/gsc/scripts/pkg/bio/repeatmasker/RepeatMasker-open-3-0-8/RepeatMasker',
        '3.1.0' => '/gsc/scripts/pkg/bio/repeatmasker/RepeatMasker-open-3-1-0/RepeatMasker',
        '3.1.5' => '/gsc/scripts/pkg/bio/repeatmasker/RepeatMasker-open-3-1-5/RepeatMasker',
        '3.1.8' => '/gsc/scripts/pkg/bio/repeatmasker/RepeatMasker-open-3-1-8/RepeatMasker',
        '3.2.7' => '/gsc/scripts/pkg/bio/repeatmasker/RepeatMasker-3.2.7/RepeatMasker',
        '3.2.9' => '/gsc/scripts/pkg/bio/repeatmasker/RepeatMasker-open-3.2.9/RepeatMasker',
        '12122000' => '/gsc/scripts/pkg/bio/repeatmasker/RepeatMasker12122000/RepeatMasker',
        '20020505' => '/gsc/scripts/pkg/bio/repeatmasker/RepeatMasker20020505/RepeatMasker',
        '20030921' => '/gsc/scripts/pkg/bio/repeatmasker/RepeatMasker20030921/RepeatMasker',
        '2004March6' => '/gsc/scripts/pkg/bio/repeatmasker/RepeatMasker2004March6/RepeatMasker',
    );
}

sub path_for_version {
    my ($self, $version) = @_;
    my %mapping = $self->version_to_path_mapping;
    return $mapping{$version};
}

# TODO For now, just allowing use of the params/flags I care about. This will
# need to be expanded if other groups are to use this wrapper.
sub flags {
    return qw/
        xsmall
        ace
        gff
    /;
}

sub is_flag {
    my ($self, $option) = @_;
    for my $flag ($self->flags) {
        return 1 if $flag eq $option;
    }
    return 0;
}

sub _generate_params {
    my $self = shift;
    my %params;

    $params{'dir'} = $self->temp_working_directory;

    if ($self->repeat_library and $self->repeat_library ne '') {
        $params{'lib'} = $self->repeat_library;
    }
    if ($self->species) {
        $params{'species'} = $self->species;
    }
    if ($self->make_ace) {
        $params{'ace'} = 1;
    }
    if ($self->xsmall) {
        $params{'xsmall'} = 1;
    }
    if ($self->make_gff) {
        $params{'gff'} = 1;
    }

    return %params;
}

sub execute {
    my $self = shift;

    $self->debug_message("Preparing to execute repeat masker.");

    $self->_check_input_fasta;
    $self->_set_masked_fasta unless defined $self->masked_fasta;
    $self->_set_ace_file_location if $self->make_ace and not defined $self->ace_file_location;
    $self->_set_gff_file_location if $self->make_gff and not defined $self->gff_file_location;
    $self->_set_excluded_sequence_file_location if $self->exclude_overly_masked and not defined $self->overly_masked_sequence_fasta;

    # Prepare for skip copies the input fasta to the output location, which is why the input fasta
    # needs to be checked first and the masked fasta location needs to be figured out. Also, the
    # ace file and gff file locations need to be set so they can be passed on to later steps. Workflow 
    # doesn't allow links to not have values, so even though the files don't exist in the case that
    # execution is skipped, they still gotta be defined.
    if ($self->skip_masking) {
        $self->_prepare_for_skip;
        return 1;
    }

    $self->_validate_species_and_library;
    $self->_set_temp_working_directory unless defined $self->temp_working_directory;

    my $executable = $self->path_for_version($self->version);
    $self->_validate_executable($executable);

    my %params = $self->_generate_params;
    my $cmd_string = join(' ', $executable, $self->_params_to_string(%params));

    $self->_execute_repeat_masker($cmd_string);
    $self->debug_message("Successfully executed repeat masker!");
    return 1;
}

sub _execute_repeat_masker {
    my ($self, $cmd_string) = @_;
    
    # Repeat masker dumps files next to the input fasta file, which may not be
    # writable. Instead use the temp working directory and copy files to where
    # you want them.
    my $rv;
    {
        local $CWD = $self->temp_working_directory;
        $self->debug_message("Executing repeat masker command $cmd_string from " . 
            "working directory " . $self->temp_working_directory);
        $rv = system($cmd_string);
    }

    unless ($rv == 0) {
        confess "Could not execute repeat masker command: $cmd_string";
    }

    $self->_move_files_from_working_directory;
    return 1;
}

sub _move_files_from_working_directory {
    my $self = shift;
    my $working_dir = $self->temp_working_directory;

    # All of the output files use the base name of the input fasta file
    # and then append a suffix.
    my ($fasta_name) = fileparse($self->fasta_file);

    my %suffix_to_property_mapping = (
        'out' => 'raw_output_directory',
        'masked' => 'masked_fasta',
        'ace' => 'ace_file_location',
        'gff' => 'gff_file_location',
    );

    my @output_files = glob("$working_dir/*");
    for my $output_file (@output_files) {
        my (undef, undef, $suffix) = fileparse($output_file, keys %suffix_to_property_mapping);
        next unless $suffix;
        next unless exists $suffix_to_property_mapping{$suffix};
        
        my $property = delete $suffix_to_property_mapping{$suffix};
        my $target = $self->$property;
        next unless defined $target;

        Genome::Sys->copy_file($output_file, $self->$property);

        # Filter out sequences from fasta file that are more than X% masked
        if ($property eq 'masked_fasta' and -e $self->$property) {
            $self->_filter_masked_sequences_from_fasta_file;
        };
    }

    if (%suffix_to_property_mapping) {
        my @missing = sort keys %suffix_to_property_mapping;
        $self->warning_message("Found no output files in working directory $working_dir " .
            "with suffixes: " . join(" ", @missing));
    }

    return 1;
}

# Filter sequences that are greater than X% masked from output fasta, where X
# is whatever maximum_percent_masked was set to. 
sub _filter_masked_sequences_from_fasta_file {
    my $self = shift;
    return 1 unless $self->exclude_overly_masked;
    my $max_percent = $self->maximum_percent_masked;
    return 1 unless defined $max_percent;

    my $temp_file_fh = File::Temp->new();
    my $temp_file = $temp_file_fh->filename;
    $temp_file_fh->close;

    my $temp_excluded_fh = File::Temp->new();
    my $temp_excluded = $temp_excluded_fh->filename;
    $temp_excluded_fh->close;

    my $input_fasta_io = Bio::SeqIO->new(
        -format => 'fasta',
        -file => $self->masked_fasta,
    );
    my $output_fasta_io = Bio::SeqIO->new(
        -format => 'fasta',
        -file => ">$temp_file",
    );
    my $excluded_fasta_io = Bio::SeqIO->new(
        -format => 'fasta',
        -file => ">$temp_excluded",
    );

    while (my $seq_obj = $input_fasta_io->next_seq()) {
        my $seq = $seq_obj->seq;
        my $n_count = $seq =~ tr/n|N//;
        my $total_length = length($seq);
        my $n_percent = ($n_count / $total_length) * 100;
        # Need some record of what was excluded
        if ($n_percent > $self->maximum_percent_masked) {
            $self->debug_message("Sequence " . $seq_obj->display_name . " $n_percent % masked, excluding from masked fasta file");
            $excluded_fasta_io->write_seq($seq_obj);
            next;
        }
        $output_fasta_io->write_seq($seq_obj);
    }
    $input_fasta_io->close;

    # Replace original fasta file with filtered fasta
    unlink $self->masked_fasta or die "Could not remove " . $self->masked_fasta;
    Genome::Sys->copy_file($temp_file, $self->masked_fasta);
    Genome::Sys->copy_file($temp_excluded, $self->overly_masked_sequence_fasta);
    $output_fasta_io->close;

    return 1;
}

sub _params_to_string {
    my ($self, %params) = @_;
    
    my @words;
    for my $param (sort keys %params) {
        if ($self->is_flag($param)) {
            push @words, '-' . $param;
        }
        else {
            push @words, '-' . $param, $params{$param};
        }
    }

    push @words, $self->fasta_file;
    
    if ($self->raw_output_directory) {
        my $file = $self->raw_output_directory . '/repeat_masker.output';
        push @words, ">$file 2>&1";
    }
    else {
        push @words, ">/dev/null 2>&1";
    }

    return join(' ', @words);
}

sub _check_input_fasta {
    my $self = shift;
    my $rv = eval { Genome::Sys->validate_file_for_reading($self->fasta_file) };
    if ($@ or not $rv) {
        my $msg = "Cannot read from file " . $self->fasta_file;
        $msg .= ", reason: $@" if $@;
        confess $msg;
    }
    return 1;
}

sub _validate_species_and_library {
    my $self = shift;
    if ((not defined $self->repeat_library or $self->repeat_library eq '') and not defined $self->species) {
        confess "Either repeat library or species must be defined!";
    }
    elsif ((defined $self->repeat_library and $self->repeat_library ne '') and defined $self->species) {
        $self->warning_message("Both repeat library and species are specified, choosing repeat library!");
    }

    if (defined $self->repeat_library and $self->repeat_library ne '') {
        my $rv = eval { Genome::Sys->validate_file_for_reading($self->repeat_library) };
        if ($@ or not $rv) {
            confess "Could not validate repeat library " . $self->repeat_library . " for reading!";
        }
    }

    return 1;
}

sub _set_masked_fasta {
    my $self = shift;
    my $masked_fasta = $self->fasta_file . ".repeat_masked";
    $self->masked_fasta($masked_fasta);
    $self->debug_message("Masked fasta file path not given, defaulting to " . $self->masked_fasta);
    return 1;
}

sub _set_ace_file_location {
    my $self = shift;
    my $default_ace_file = $self->fasta_file. ".repeat_masker.ace";
    $self->ace_file_location($default_ace_file);
    $self->debug_message("Ace file is being generated and location not given, defaulting to $default_ace_file");
    return 1;
}

sub _set_gff_file_location {
    my $self = shift;
    my $default_gff_file = $self->fasta_file . ".repeat_masker.gff";
    $self->gff_file_location($default_gff_file);
    $self->debug_message("Gff is file is being generated and location not given, defaulting to $default_gff_file");
    return 1;
}

sub _set_excluded_sequence_file_location {
    my $self = shift;
    my $default_file = $self->fasta_file . ".overly_masked";
    $self->overly_masked_sequence_fasta($default_file);
    $self->debug_message("Overly masked fasta file not given, defaulting to " . $self->overly_masked_sequence_fasta);
    return 1;
}

sub _prepare_for_skip {
    my $self = shift;
    $self->debug_message("skip_masking flag is set, copying input fasta to masked fasta location");
    my $rv = Genome::Sys->copy_file($self->fasta_file, $self->masked_fasta);
    confess "Trouble executing copy of " . $self->fasta_file . " to " . $self->masked_fasta unless defined $rv and $rv;
    $self->debug_message("Copy of input fasta at " . $self->fasta_file . " to masked fasta path at " .
        $self->masked_fasta . " successful, exiting!");
    return 1;
}

sub _set_temp_working_directory {
    my $self = shift;
    my $temp_dir = tempdir(
        'RepeatMasker-XXXXXX',
        DIR => '/tmp/',
        CLEANUP => 0,
        UNLINK => 0,
    );
    chmod(0775, $temp_dir);
    $self->temp_working_directory($temp_dir);
    $self->debug_message("Not given temp working directory, using $temp_dir");
    return 1;
}

sub _validate_executable {
    my ($self, $exe) = @_;
    my $rv = eval { Genome::Sys->validate_file_for_execution($exe) };
    if ($@ or not $rv) {
        my $msg = "Cannot execute file $exe";
        $msg .= ", reason: $@" if $@;
        confess $msg;
    }
    return 1;
}

sub _generate_ace_file {
    my $self = shift;
    my @ace_files = glob($self->temp_working_directory."/*ace");
    $self->debug_message("Ace file $ace_files[0]");

    $self->debug_message("Fasta " . $self->masked_fasta);
	my ( $masked_file_name, $dir ) = fileparse($self->masked_fasta);
	my @masked_file_name = split(/repeat_masker/, $masked_file_name);

    #local $CWD = $self->temp_working_directory;

	my $new_ace_file_name = $self->ace_file_location;
	$self->debug_message("Ace file name: ". $new_ace_file_name);
	my $ace_file_fh = IO::File->new($new_ace_file_name, "a");
	$self->error_message("Could not get handle for ace file($ace_file_fh): $!") and return unless $ace_file_fh;

    foreach my $ace_file (@ace_files) {
		my @file_name = split(/\./, $ace_file);
		my $masked_fh = IO::File->new($file_name[0], "r");
		$self->error_message("Could not get handle for masked file($file_name[0]): $!") and return 0 unless $masked_fh;
		undef my $contig_name;
		while (my $line = $masked_fh->getline) {
			chomp $line;
			if ($line =~ m/^>(.*)/) {
		  		$contig_name = $1;
				$self->debug_message("contig_name: ". $contig_name);
				last;
			}
		}
		$masked_fh->close;

		my $ace_fh = IO::File->new($ace_file, "r");
		$self->error_message("Could not get handle for ace file($ace_file): $!") and return unless $ace_fh;

		undef my $line_count;
		$line_count = 0;
		while (my $line = $ace_fh->getline) {
			$ace_file_fh->print("Sequence $contig_name\n") if ($line_count == 0);
            $line = ~s/\s+(\+|\-)\s+/ /;
			$ace_file_fh->print($line);
			$line_count++;
		}
		$ace_file_fh->print("\n");
	}
    $ace_file_fh->close;
    return 1;
}

1;
