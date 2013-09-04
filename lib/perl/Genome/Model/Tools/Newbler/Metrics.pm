package Genome::Model::Tools::Newbler::Metrics;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Newbler::Metrics {
    is => 'Genome::Model::Tools::Newbler',
    has => [
        assembly_directory => {
            is => 'Text',
            doc => 'Path to soap assembly',
        },
    ],
    has_optional => [
        first_tier => {
            is => 'Number',
            doc => 'First tier value',
        },
        second_tier => {
            is => 'Number',
            doc => 'Second tier value',
        },
        major_contig_length => {
            is => 'Number',
            default_value => 500,
            doc => 'Cutoff value for major contig length',
        },
        output_file => {
            is => 'Text',
            doc => 'Stats output file',
        },
        _metrics => { is_transient => 1, },
    ],
};

sub help_brief {
    return 'Produce metrics for newbler assemblies'
}

sub __errors__ {
    my $self = shift;

    my @errors = $self->SUPER::__errors__(@_);
    return @errors if @errors;

    if ( $self->assembly_directory ) {
        if ( not -d $self->assembly_directory ) {
            push @errors, UR::Object::Tag->create(
                type => 'invalid',
                properties => [qw/ assembly_directory /],
                desc => 'The assembly_directory is not a directory!',
            );
            return @errors;
        }
        if ( not defined $self->output_file ) {
            my $create_dir = $self->create_consed_edit_dir;
            return if not $create_dir;
            $self->output_file( $self->stats_file );
        }
    }
    elsif ( not $self->output_file ) { 
        push @errors, UR::Object::Tag->create(
            type => 'invalid',
            properties => [qw/ output_file /],
            desc => 'No output file given and no assembly_directory given to determine the output file!',
        );
    }

    my @reads_files = grep { -s } $self->input_fastq_files;
    if ( not @reads_files ) {
        push @errors, UR::Object::Tag->create(
            type => 'invalid',
            properties => [qw/ assembly_directory /],
            desc => 'No input reads files found in assembly_directory: '.$self->assembly_directory,
        );
    }

    for my $file_method ( qw/ all_contigs_fasta_file newb_read_status_file / ) {
        my $file = $self->$file_method;
        if ( not -s $file ) {
            my $file_name = File::Basename::basename( $file );
            push @errors, UR::Object::Tag->create(
                type => 'invalid',
                properties => [qw/ assembly_directory /],
                desc => "No newbler $file_name file found in assembly_directory: ".$self->assembly_directory,
           ); 
        }
    }

    # 454Scaffolds.txt 
    # assembly not scaffolded if missing

    return @errors;
}

sub execute {
    my $self = shift;
    $self->status_message('Newbler metrics...');

    # tier values
    my ($t1, $t2);
    if ($self->first_tier and $self->second_tier) {
        $t1 = $self->first_tier;
        $t2 = $self->second_tier;
    }
    else {
        my $est_genome_size = -s $self->all_contigs_fasta_file;
        $t1 = int ($est_genome_size * 0.2);
        $t2 = int ($est_genome_size * 0.2);
    }
    $self->status_message('Tier one: 1 to '.$t1);
    $self->status_message('Tier two: '.$t1.' to '.($t1 + $t2));

    # metrics
    my $metrics = Genome::Model::Tools::Sx::Metrics::Assembly->create(
        major_contig_threshold => $self->major_contig_length,
        tier_one => $t1,
        tier_two => $t2,
    );
    $self->_metrics($metrics);

    # add contigs
    $self->status_message('Add metrics from 454AllContigs.fna file');
    my $add_contigs_ok = $self->_metrics_from_newbler_contigs_file( $metrics );
    return if not $add_contigs_ok;

    # input reads
    for my $reads_file ( $self->input_fastq_files ) {
        $self->status_message('Add reads file: '.$reads_file);
        my $add_reads = $metrics->add_reads_file_with_q20($reads_file);
        return if not $add_reads;
    }

    # assembled reads
    $self->status_message('Add read status metrics: ');
    my $add_read_depth = $self->_metrics_from_newb_read_status_file($metrics);
    return if not $add_read_depth;

    # transform metrics
    my $text = $metrics->transform_xml_to('txt');
    if ( not $text ) {
        $self->error_message('Failed to transform metrics to text!');
        return;
    }

    # write file
    my $output_file = $self->output_file;
    unlink $output_file if -e $output_file;
    $self->status_message('Write output file: '.$output_file);
    my $fh = eval{ Genome::Sys->open_file_for_writing($output_file); };
    if ( not $fh ) {
        $self->error_message('Failed to open metrics output file!');
        return;
    }
    $fh->print($text);
    $fh->close;

    $self->status_message('Velvet metrics...DONE');
    return 1;
}

sub _scaffolding_info_from_newbler {
    my $self = shift;

    my %scaffolds;
    if ( not -s $self->scaffolds_agp_file ) { 
        $self->status_message('No newbler scaffolds file found, assuming assembly is not scaffolded');
        return \%scaffolds;
    }

    my $fh = Genome::Sys->open_file_for_reading( $self->scaffolds_agp_file );
    while ( my $line = $fh->getline ) {
        next if $line =~ /\s+fragment\s+/;
        chomp $line;
        my @tmp = split( /\s+/, $line );
        #tmp[0] = scaffold name
        #tmp[5] = contig name
        #tmp[7] = contig length
        my $scaffold_number = $tmp[0];
        $scaffold_number =~ s/scaffold//;
        $scaffold_number =~ s/^0+//;
        $scaffolds{$tmp[5]} = $scaffold_number;
    }
    $fh->close;

    $self->{_SCAFFOLDS} = \%scaffolds;

    return \%scaffolds;
}

sub _metrics_from_newbler_contigs_file {
    my ( $self, $metrics ) = @_;

    #empty hashref if not scaffolds in assembly
    my $scaffolds = $self->_scaffolding_info_from_newbler;

    my ( $content_at, $content_gc, $content_nx, $contigs_5k, $assembly_length ) = ( 0,0,0,0,0 );

    my $reader = Genome::Model::Tools::Sx::Reader->create(config => [$self->all_contigs_fasta_file.':type=fasta']);
    while (my $seqs = $reader->read ) {
        for my $seq ( @$seqs ) {
            # contigs lengths
            $self->_metrics->contigs->{$seq->{id}} = length $seq->{seq};
            $assembly_length += length $seq->{seq};

            # supercontigs lengths
            my $supercontig_id = ( exists $scaffolds->{$seq->{id}} ) ?
                $scaffolds->{$seq->{id}} : #part of scaffold
                $seq->{id};                #it's own scaffold
            $self->_metrics->supercontigs->{$supercontig_id} += length $seq->{seq};

            # gt 5k contigs
            $contigs_5k += length $seq->{seq} if length $seq->{seq} >= 5000;

            # genome content metrics
            my %base_counts = ( a => 0, t => 0, g => 0, C => 0, n => 0, x => 0 );
            foreach my $base ( split('', $seq->{seq}) ) {
                $base_counts{lc $base}++;
            }
            $content_at += $base_counts{a} + $base_counts{t};
            $content_gc += $base_counts{g} + $base_counts{c};
            $content_nx += $base_counts{n} + $base_counts{x};
        }
    }

    #set major/minor lengths
    for my $type ( 'contigs', 'supercontigs' ) {
        my ( $major, $minor, $count ) = ( 0,0,0 );
        for my $type_name ( keys %{$self->_metrics->{$type}} ) {
            $count++;
            my $length = $self->_metrics->{$type}->{$type_name};
            if ( $length >= $self->major_contig_length ) {
                $major += $length;
            } else {
                $minor += $length;
            }
        }
        $metrics->set_metric($type.'_count', $count);
        $metrics->set_metric($type.'_length', $major + $minor);
        $metrics->set_metric($type.'_major_length', $major);
        $metrics->set_metric($type.'_minor_length', $minor);
    }

    # set at/gc/nx contents
    $metrics->set_metric('content_at', $content_at);
    $metrics->set_metric('content_gc', $content_gc);
    $metrics->set_metric('content_nx', $content_nx);

    # set 5k contigs
    $metrics->set_metric('contigs_length_5k', $contigs_5k);

    # assembly length
    $metrics->set_metric('assembly_length', $assembly_length);

    return 1;
}

sub _metrics_from_newb_read_status_file {
    my ( $self, $metrics ) = @_;

    my %covered;
    my %reads_in_contigs;
    my $reads_assembled = 0;
    
    my $fh = Genome::Sys->open_file_for_reading( $self->newb_read_status_file );
    # exclude header
    my $header = $fh->getline;
    while ( my $line = $fh->getline ) {
        chomp $line;
        my @tmp = split( /\s+/, $line );
        #$tmp[0] = read name
        #$tmp[1] = status, assembled or not
        #$tmp[2] = contig name
        #$tmp[3] = 3' position
        #$tmp[4] = 3' complimented or uncomp
        #$tmp[5] = contig name
        #$tmp[6] = 5' position
        #$tmp[7] = 5' complimented or uncomp
        if ( $tmp[1] eq 'Assembled' or $tmp[1] eq 'PartiallyAssembled' ) {
            # no split reads .. shouldn't happend anymore
            next if not $tmp[2] eq $tmp[5];
            # mal-formed line die?
            next if not $tmp[3] =~ /^\d+$/ and not $tmp[6] =~ /^\d+$/;
            # ignore contigs not part of main assembly
            next if not exists $self->_metrics->contigs->{$tmp[2]};

            my $start = ( $tmp[3] >= $tmp[6] ) ? $tmp[6] : $tmp[3];
            my $stop = ( $tmp[3] > $tmp[6] ) ? $tmp[3] : $tmp[6];
            for ( $start .. $stop ) {
                @{$covered{$tmp[2]}}[$_]++;
            }
            
            #track # of reads in contigs
            $reads_in_contigs{$tmp[2]}++;
            $reads_assembled++;
        }
    }
    $fh->close;

    # determine coverage
    my ( $zero_x, $one_x, $two_x, $three_x, $four_x, $five_x ) = ( 0,0,0,0,0,0 );
    for my $contig ( keys %covered ) {
        # ignore contigs not in assembly
        #next if not exists $self->_metrics->contigs->{$contig};
        shift @{$covered{$contig}}; #[0] is always undef
        for ( @{$covered{$contig}} ) {
            if ( not defined $_ ) {
                $zero_x++;
                next;
            }
            $one_x++   if $_ > 0;
            $two_x++   if $_ > 1;
            $three_x++ if $_ > 2;
            $four_x++  if $_ > 3;
            $five_x++  if $_ > 4;
        }
    }
    $metrics->set_metric('coverage_5x', $five_x);
    $metrics->set_metric('coverage_4x', $four_x);
    $metrics->set_metric('coverage_3x', $three_x);
    $metrics->set_metric('coverage_2x', $two_x);
    $metrics->set_metric('coverage_1x', $one_x);
    $metrics->set_metric('coverage_0x', $zero_x);

    # major/minor contigs/supercontigs read counts
    my $scaffolds = $self->{_SCAFFOLDS};

    my $major_contigs_read_count = 0;
    my $minor_contigs_read_count = 0;
    my $major_supercontigs_read_count = 0;
    my $minor_supercontigs_read_count = 0;

    for my $contig_id ( keys %reads_in_contigs ) {
        # ignore contigs not in final assembly
        if ( exists $self->_metrics->contigs->{$contig_id} ) {
            my $read_count = $reads_in_contigs{$contig_id};

            # contigs read count
            my $contig_length = $self->_metrics->contigs->{$contig_id};
            if ( $contig_length >= $self->major_contig_length ) {
                $major_contigs_read_count += $read_count;
            } else {
                $minor_contigs_read_count += $read_count
            }
            
            # supercontig read count
            my $supercontig_id = ( exists $scaffolds->{$contig_id} ) ?
                $scaffolds->{$contig_id} :
                $contig_id;

            my $supercontig_length = $self->_metrics->supercontigs->{$supercontig_id};
            if ( $supercontig_length >= $self->major_contig_length ) {
                $major_supercontigs_read_count += $read_count;
            } else {
                $minor_supercontigs_read_count += $read_count;
            }
        }
    }

    $metrics->set_metric('contigs_major_read_count', $major_contigs_read_count);
    $metrics->set_metric('supercontigs_major_read_count', $major_supercontigs_read_count);
    my $major_contigs_read_percent = sprintf( "%.1f", $major_contigs_read_count/ $reads_assembled * 100);
    $metrics->set_metric('contigs_major_read_percent', $major_contigs_read_percent);
    my $major_supercontigs_read_percent = sprintf( "%.1f", $major_supercontigs_read_count/ $reads_assembled * 100);
    $metrics->set_metric('supercontigs_major_read_percent', $major_supercontigs_read_percent);

    
    $metrics->set_metric('contigs_minor_read_count', $minor_contigs_read_count);
    $metrics->set_metric('supercontigs_minor_read_count', $minor_supercontigs_read_count);
    my $minor_contigs_read_percent = sprintf( "%.1f", $minor_contigs_read_count/ $reads_assembled * 100);
    $metrics->set_metric('contigs_minor_read_percent', $minor_contigs_read_percent);
    my $minor_supercontigs_read_percent = sprintf( "%.1f", $minor_supercontigs_read_count/ $reads_assembled * 100);
    $metrics->set_metric('supercontigs_minor_read_percent', $minor_supercontigs_read_percent);



    # all contigs read counts
    $metrics->set_metric('reads_assembled', $reads_assembled);
    $metrics->set_metric('reads_assembled_unique', $reads_assembled);
    my $reads_processed = $metrics->get_metric('reads_processed');
    $metrics->set_metric('reads_assembled_success', sprintf('%.3f', $reads_assembled / $reads_processed));
    my $reads_not_assembled = $reads_processed - $reads_assembled;
    $metrics->set_metric('reads_not_assembled', $reads_not_assembled);

    return 1;
}

1;

