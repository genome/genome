package Genome::Model::Tools::Sx::Trim::Flexbar;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Sx::Trim::Flexbar {
    is => 'Genome::Model::Tools::Sx::ExternalCmdBase',
    has => [
        __PACKAGE__->_cmd_properties,
        adapter => {
            is => 'Text',
            is_optional => 1,
            doc => 'Adaptor sequence to be removed.',
        },
        remove_revcomp => {
            is => 'Boolean',
            is_optional => 1,
            default => 0,
            doc => 'Remove the reverse compl;ement of the single adapter.',
        },
        version => {
            is => 'Text',
            valid_values => [ __PACKAGE__->cmd_versions ],
            doc => 'Verion of flexbar to use.',
        },
        _tmp_flexbar_inputs => { is => 'Array', is_optional => 1, },
    ],
};

sub help_brief {
    return 'Trim with the [flex]ible [bar]code detection and adapter removal';
}

sub cmd_display_name {
    return 'flexbar';
}

sub _cmd_versions {
    return (
        '229' => '/usr/bin/flexbar229',
	'230' => '/gsc/bin/flexbar',
    );
}

sub cmd_versions {
    my $self = shift;
    my %versions = $self->_cmd_versions;
    return sort keys %versions;
}

sub _cmd_properties {
    # Removed algorithm, cut_off
    return (
        adapters => {
            is => 'Text',
            is_optional => 1,
            doc => 'Fasta file of adapter sequences to be removed.',
        },
        adapter_min_overlap => {
            is => 'Text',
            is_optional => 1,
            doc => ' Minimum overlap of adapter and read in base pairs.',
        },
         adapter_threshold  => {
             is => 'Number',
             is_optional => 1,
             doc => 'Allowed mismatches and indels per 10 bases for adapter',
         },
         adapter_trim_end => {
             is => 'Text',
             valid_values => [qw/ ANY RIGHT LEFT RIGHT_TAIL LEFT_TAIL  /],
             doc => 'Decides on which end adapter removal is performed.',
         },
         adapter_tail_length => {
             is => 'Number',
             is_optional => 1,
             doc => ' Number of bases for tail trim-end types, default: adapter length',
         },
         adapter_no_adapt => {
             is => 'Boolean',
             is_optional => 1,
             doc => 'Do not adapt min-overlap to overlap length at ends.',
         },
         adapter_match => {
             is => 'Number',
             is_optional => 1,
             doc => 'Match score.',
         },
         adapter_mismatch => {
             is => 'Number',
             is_optional => 1,
             doc => 'Mismatch score.',
         },
         adapter_gap_cost => {
             is => 'Number',
             is_optional => 1,
             doc => 'Gap score.',
         },
         min_readlength => {
             is => 'Text',
             is_optional => 1,
             doc => 'Minimum readlength in basepairs after adapter removal or read will be discarded.',
         },
        min_read_length => {
             is => 'Text',
             is_optional => 1,
             doc => 'Minimum readlength in basepairs after adapter removal or read will be discarded, use in version 230 only',
         },
         max_uncalled => {
             is => 'Text',
             is_optional => 1,
             doc => 'Number of allowed uncalled bases in a read.',
         },
         no_length_dist => {
             is => 'Boolean',
             is_optional => 1,
#            default_value => '1',
             doc => 'Prevent writing length distributions for read output files.',
         },
	post_trim_length =>  {
	    is => 'Number',
	    is_optional => 1,
	    doc => 'Trim to specified read length from 3\' end after removal, v230 only',
	},
	pre_trim_left => {
	    is => 'Number',
	    is_optional => 1,
	    doc => 'Trim given number of bases on 5\' read end before detection, v230 only.',
	},
	pre_trim_right => {
	    is => 'Number',
	    is_optional => 1,
	    doc => 'Trim specified number of bases on 3\' end prior to detection, v230 only.',
	},
	pre_trim_phred => {
	    is => 'Number',
	    is_optional => 1,
	    doc => 'Trim 3\' end until specified or higher quality reached, v230 only',
	},
         removal_tag => {
             is => 'Boolean',
             is_optional => 1,
             doc => 'Tag reads for which adapter or barcode is removed.',
         },
         threads => {
             is => 'Text',
             is_optional => 1,
             default_value => 1,
             doc => 'Number of threads to use.',
         },
     );
 }

 sub execute {
     my $self = shift;

     my $resolve_adapters = $self->_resolve_adapters;
     return if not $resolve_adapters;

     my @input_params = $self->_resolve_input_params;
     return if not @input_params;

     my $output = $self->_init_ouptut;
     return if not $output;

     $self->debug_message('Run flexbar...');

     my $cmd = $self->build_command;

     # set input param name .. either reads or source
     my $input_param_name = $self->_input_data_param_name();
     $cmd .= ' --'.$input_param_name.' '.$input_params[0]->{file};
     $cmd .= ' --'.$input_param_name.'2 '.$input_params[1]->{file} if $input_params[1];

     # set input format ..sanger or sanger-fastq
     my $input_format = $self->_input_data_format();
     $cmd .= ' --format '.$input_format;

     # set output file name
     $cmd .= ' --target '.$self->_tmpdir.'/output.fastq',

     # set versions env variable
     my $env_var = $self->_env_var_for_version();
     $cmd = $env_var.' '.$cmd if $env_var;

    my $rv = $self->_run_command($cmd);
    return if not $rv;
    $self->debug_message('Run flexbar...OK');

    my @fastq_files = glob $self->_tmpdir.'/output*.fastq';
    my @outputs = grep { $_ !~ /single/ } @fastq_files;
    if ( not @outputs ) {
        $self->error_message('Failed to find output files! Files in output directory: ');
        return;
    }
    $self->debug_message('Output: '.join(' ', @outputs));
    
    my $output_reader = Genome::Model::Tools::Sx::Reader->create(
        config => [ map { $_.':type=sanger' } @outputs ],
    );
    if ( not $output_reader ) {
        $self->error_message('Failed to open reader for flexbar output!');
        return;
    }

    $self->debug_message('Processing flexbar output...');
    while ( my $seqs = $output_reader->read ) {
        $output->write($seqs);
    }

    $self->_rm_tmpdir;

    return 1;
}

sub _input_data_param_name {
    my $self = shift;

    my $name = 'source';
    if( $self->version eq '230' ) {
	$name = 'reads';
    }   
    return $name;
}

sub _input_data_format {
    my $self = shift;

    my $format = 'fastq-sanger';
    if( $self->version eq '230' ) {
	$format = 'sanger';
    }
    return $format;
}

sub _env_var_for_version {
    my $self = shift;

    my $env_var;
    if( $self->version eq '230' ) {
	$env_var = 'LD_LIBRARY_PATH=/gscmnt/gc3001/assembly/Downloads/Flexbar_2.4/flexbar_v2.4_linux64';
	# this is needed to use multiple nodes for v230
    }
    return $env_var;
}

sub _resolve_adapters {
    my $self = shift;

    if ( $self->adapters and $self->adapter ) {
        $self->error_message('Cannot specify both adpaters file (adapters) and single adapter (adapter) params!');
        return;
    }

    if ( $self->adapters ) { # using adapters file
        $self->error_message('Cannot specify both adpaters file (adapters) and single adapter (adapter) params!') if $self->adapter;
        $self->error_message('Cannot specify to remove reverse complement adapter (remove_revcomp) and adpaters file (adapters) params!') if $self->remove_revcomp;
        return 1 if -s $self->adapters;
        $self->error_message('Adapters file does not exist! '.$self->adapters);
        return;
    }

    my $adapters_file = $self->_tmpdir.'/adapters.fasta';
    $self->adapters($adapters_file);
    my $writer = Genome::Model::Tools::Sx::PhredWriter->create(file => $adapters_file);
    if ( not $writer ) {
        $self->error_message('Failed to create fasta writer!');
        return;
    }

    my $adapter = uc($self->adapter);
    $writer->write({ id => 'Adapter', seq => $adapter, });
    if ( $self->remove_revcomp ) {
        my $revcomp_adapter = reverse $adapter;
        $revcomp_adapter =~ tr/TCGA/AGCT/;
        $writer->write({ id => 'Adapter-Rev-Complement', seq => uc( $revcomp_adapter ), });
    }

    return 1;
}

sub _required_type_and_counts_for_inputs {
    return ( 'sanger', [qw/ 1 2 /], );
}

1;

