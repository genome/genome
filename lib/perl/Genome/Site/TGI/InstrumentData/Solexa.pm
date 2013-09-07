
# review jlolofie:
# 0. "project" property in class definition- cant find research_project_name():  $self->research_project_name
# 1. calculate_alignment_estimated_kb_usage() could estimate usage instead of harcoded value
# 2. resolve_full_path - remove return that comments say is not needed?
# 3. add comments to resolve_full_path- why it is trying to find a collection of paths and die if multiple- do the other
#    file types still get used?
# 4. dump_to_filesystem() ?
# 5. resolve_fastq_filenames- funny Temp/ hack not needed anymore
# 6. resolve_adapter_file - make properties out of paths



package Genome::Site::TGI::InstrumentData::Solexa;

use strict;
use warnings;

use Genome;

use File::Basename;

class Genome::Site::TGI::InstrumentData::Solexa {
    is => ['Genome::Site::TGI::InstrumentData', 'Genome::Sys'],
    table_name => <<EOS
        (
            select
                --to_char(s_rev.seq_id) id,
                to_char(i.analysis_id) id,

                'solexa' sequencing_platform,

                i.research_project project_name,

                i.target_region_set_name,

                --s_rev.sample_id,
                lib.sample_id,

                i.library_id,

                --s_rev.run_name,
                fc.run_name,

                fc.flow_cell_id,
                i.lane,

                r2.read_length,
                r2.filt_error_rate_avg,

                --(case when s_rev.run_type = 'Paired End Read 2' then s_rev.filt_error_rate_avg else null end) rev_filt_error_rate_avg,
                (case when r1.seq_id is not null then r2.filt_error_rate_avg else null end) rev_filt_error_rate_avg,

                --(case when s_fwd.run_type = 'Paired End Read 1' then s_fwd.filt_error_rate_avg else null end) fwd_filt_error_rate_avg,
                r1.filt_error_rate_avg fwd_filt_error_rate_avg,

                --(case when s_rev.run_type = 'Paired End Read 2' then s_rev.filt_aligned_clusters_pct else null end) rev_filt_aligned_clusters_pct,
                (case when r1.seq_id is not null then r2.filt_aligned_clusters_pct else null end) rev_filt_aligned_clusters_pct,

                --(case when s_fwd.run_type = 'Paired End Read 1' then s_fwd.filt_aligned_clusters_pct else null end) fwd_filt_aligned_clusters_pct,
                r1.filt_aligned_clusters_pct fwd_filt_aligned_clusters_pct,

                --(case when s_rev.run_type = 'Paired End Read 2' then s_rev.seq_id else null end) rev_seq_id,
                (case when r1.seq_id is not null then r2.sls_seq_id else null end) rev_seq_id,

                --(case when s_fwd.run_type = 'Paired End Read 1' then s_fwd.seq_id else null end) fwd_seq_id,
                r1.sls_seq_id fwd_seq_id,

                --(case when s_rev.run_type = 'Paired End Read 2' then s_rev.read_length else null end) rev_read_length,
                --(case when r1.seq_id is not null then r2.read_length else null end) rev_read_length,
                (case when r1.seq_id is not null then r2.read_length else -1 end) rev_read_length,

                (case when r1.seq_id is not null then r2.kilobases_read else -1 end) fwd_kilobases_read,
                (case when r2.seq_id is not null then r2.kilobases_read else -1 end) rev_kilobases_read,

                --(case when s_fwd.run_type = 'Paired End Read 1' then s_fwd.read_length else null end) fwd_read_length,
                --r1.read_length fwd_read_length,
                nvl(r1.read_length,-1) fwd_read_length,

                --(case when s_rev.run_type = 'Paired End Read 2' then s_rev.run_type else null end) rev_run_type,
                (case when r1.seq_id is not null then 'Paired End Read 2' else null end) rev_run_type,

                --(case when s_fwd.run_type = 'Paired End Read 1' then s_fwd.run_type else null end) fwd_run_type,
                (case when r1.seq_id is not null then 'Paired End Read 1' else null end) fwd_run_type,

                --(case when s_rev.run_type = 'Paired End Read 2' then 'Paired' else 'Standard' end) run_type,
                (case when r1.seq_id is not null then 'Paired' else 'Standard' end) run_type,

                --s_rev.gerald_directory,
                i.gerald_directory,

                --s_rev.median_insert_size,
                i.median_insert_size,

                --s_rev.sd_above_insert_size,
                i.sd_above_insert_size,

                --s_rev.sd_below_insert_size,
                i.sd_below_insert_size,

                --s_rev.is_external,
                0 is_external,

                --archive.path archive_path,
                archive2.path archive_path,
                gerald_bam.path gerald_bam_path,

                collect_gc_bias.path collect_gc_bias_path,
                fastqc.path fastqc_path,

                --adaptor.path adaptor_path,
                --adaptor2.path adaptor_path,
                '/gscmnt/sata114/info/medseq/adaptor_sequences/solexa_adaptor_pcr_primer'
                    || (case when sample_type = 'rna' then '_SMART' else '' end) adaptor_path,

                --(case when s_fwd.run_type = 'Paired End Read 1' then s_fwd.FILT_CLUSTERS else null end) fwd_filt_clusters,
                (case when r1.seq_id is not null then i.filt_clusters else null end) fwd_filt_clusters,

                --(case when s_rev.run_type = 'Paired End Read 2' then s_rev.FILT_CLUSTERS else null end) rev_filt_clusters,
                (case when r1.seq_id is not null then i.filt_clusters else null end) rev_filt_clusters,

                --(nvl(s_fwd.FILT_CLUSTERS,0) + s_rev.FILT_CLUSTERS) filt_clusters, 	-- s_rev.FILT_CLUSTERS is still the expected value for fragment reads
                i.filt_clusters,

                --s_rev.analysis_software_version,
                i.analysis_software_version,
                --Analysis Project
                swo.analysis_project_id,

                i.index_sequence

                --from solexa_lane_summary s_rev
                --join read_illumina r2 on r2.sls_seq_id = s_rev.seq_id --and r1.read_number = 1
                from index_illumina i
                    join flow_cell_illumina fc on fc.flow_cell_id = i.flow_cell_id
                    join read_illumina r2
                        on i.seq_id = r2.ii_seq_id
                        and (
                            (fc.run_type = 'Paired End' and r2.read_number = 2)
                            or
                            (fc.run_type = 'Fragment' and r2.read_number = 1)
                        )
                    left join seq_fs_path archive2 on archive2.seq_id = i.seq_id
                        and archive2.data_type = 'illumina fastq tgz'
                    left join seq_fs_path gerald_bam on gerald_bam.seq_id = i.seq_id
                        and gerald_bam.data_type = 'gerald bam'
                    left join seq_fs_path collect_gc_bias on collect_gc_bias.seq_id = i.seq_id
                        and collect_gc_bias.data_type = 'collect gc bias'
                    left join seq_fs_path fastqc on fastqc.seq_id = i.seq_id
                        and fastqc.data_type = 'fastqc'
                    left join read_illumina r1
                        on run_type = 'Paired End'
                        and r1.ii_seq_id = i.seq_id
                        and r1.read_number = 1
                  join GSC.library_summary lib on lib.library_id = i.library_id
                  join GSC.organism_sample sam on sam.organism_sample_id = lib.sample_id
                  join GSC.woi_sequence_product wsp on wsp.seq_id = i.seq_id
                  join work_order_item\@oltp woi on wsp.woi_id = woi.woi_id
                  join setup_work_order\@oltp swo on swo.setup_wo_id = woi.setup_wo_id
            /*
                    left join solexa_lane_summary s_fwd on s_fwd.sral_id = s_rev.sral_id and s_fwd.run_type = 'Paired End Read 1'
                    left join seq_fs_path archive on archive.seq_id = s_rev.seq_id
                        and archive.data_type = 'illumina fastq tgz'
                    left join seq_fs_path adaptor on adaptor.seq_id = s_rev.seq_id
                        and adaptor.data_type = 'adaptor sequence file'
                    where s_rev.run_type in ('Standard','Paired End Read 2')
                        and s_rev.flow_cell_id = '617ER'
            */
        )
        solexa_detail
EOS
    ,
    has_constant => [
        sequencing_platform => { value => 'solexa' },
    ],
    has_optional => [
        flow_cell_id                    => { }, # = short name
        flow_cell                       => { is => 'Genome::InstrumentData::FlowCell', id_by => 'flow_cell_id' },
        lane                            => { },
        gc_bias_path                    => { column_name => 'collect_gc_bias_path' },
        fastqc_path                     => { },
        index_sequence                  => { },
        read_length                     => { },
        fwd_read_length                 => { },
        rev_read_length                 => { },
        fwd_kilobases_read              => { },
        rev_kilobases_read              => { },
        #TODO These three columns will point to "read_length" or whatever name is decided
        #(see also https://gscweb.gsc.wustl.edu/wiki/Software_Development:Illumina_Indexed_Runs_Warehouse_Schema)
        _sls_read_length                => { calculate => q| return $self->read_length + 1| },
        _sls_fwd_read_length            => { calculate => q| return $self->fwd_read_length + 1| },
        _sls_rev_read_length            => { calculate => q| return $self->rev_read_length + 1| },
        cycles                          => { calculate => q| return $self->read_length + 1| }, #TODO point to an actual "cycles" column
        run_type                        => { },
        fwd_run_type                    => { },
        rev_run_type                    => { },
        gerald_directory                => { },
        median_insert_size              => { },
        sd_above_insert_size            => { },
        sd_below_insert_size            => { },
        is_external                     => { },
        adaptor_path                    => { },
        archive_path                    => { },
        bam_path                        => { column_name => 'gerald_bam_path'},
        analysis_software_version       => { },
        clusters                        => { column_name => 'filt_clusters' },
        fwd_clusters                    => { column_name => 'fwd_filt_clusters' },
        rev_clusters                    => { column_name => 'rev_filt_clusters' },
        fwd_seq_id                      => { },
        rev_seq_id                      => { },
        filt_error_rate_avg             => { },
        fwd_filt_error_rate_avg         => { },
        rev_filt_error_rate_avg         => { },
        fwd_filt_aligned_clusters_pct   => { },
        rev_filt_aligned_clusters_pct   => { },
        target_region_set_name          => { },
        analysis_project_id             => { },

        short_name => {
            doc => 'The essential portion of the run name which identifies the run.  The rest is redundent information about the instrument, date, etc.',
            is => 'Text',
            calculate_from => ['run_name'],
            calculate => q|($run_name =~ /_([^_]+)$/)[0]|
        },

        is_paired_end                   => {
                                            calculate_from => ['run_type'],
                                            calculate => q| if (defined($run_type) and $run_type =~ m/^Paired$/) {
                                                                return 1;
                                                             }
                                                             else {
                                                                 return 0;
                                                             } |
                                        },
        project_name => { },
        project => {
            is => "Genome::Site::TGI::Project",
            calculate => q|Genome::Site::TGI::Project->get(name => $self->research_project_name)|
        },
        _run_lane_solexa => {
            doc => 'Solexa Lane Summary from LIMS.',
            is => 'GSC::RunLaneSolexa',
            calculate => q| GSC::RunLaneSolexa->get($id); |,
            calculate_from => ['id']
        },
        # Index Illumina
        index_illumina => {
            doc => 'Index Illumina from LIMS.',
            is => 'GSC::IndexIllumina',
            calculate => q| GSC::IndexIllumina->get(analysis_id=>$id); |,
            calculate_from => [ 'id' ]
        },
        # basic relationship to the "source" of the lane

        sample_source       => { via => 'sample', to => 'source' },
        sample_source_name  => { via => 'sample_source', to => 'name' },

        # indirect via the sample source, but we let the sample manage that
        # since we sometimes don't know the source, it also tracks taxon directly
        taxon               => { via => 'sample', to => 'taxon', is => 'Genome::Site::TGI::Taxon' },
        species_name        => { via => 'taxon' },
    ],
};

sub __display_name__ {
    my $self = $_[0];
    return $self->flow_cell_id . '/' . $self->subset_name;
}

sub _calculate_paired_end_kb_usage {
    return Genome::InstrumentData::Solexa::_calculate_paired_end_kb_usage(@_);
}

sub _calculate_non_paired_end_kb_usage {
    return Genome::InstrumentData::Solexa::_calculate_non_paired_end_kb_usage(@_);
}

sub calculate_alignment_estimated_kb_usage {
    return Genome::InstrumentData::Solexa::calculate_alignment_estimated_kb_usage(@_);
}

sub resolve_full_path {
    my $self = shift;

    my @fs_path = GSC::SeqFPath->get(
        seq_id => $self->genome_model_run_id,
        data_type => [qw/ duplicate fastq path unique fastq path /],
    )
        or return; # no longer required, we make this ourselves at alignment time as needed

    my %dirs = map { File::Basename::dirname($_->path) => 1 } @fs_path;

    if ( keys %dirs > 1) {
        $self->error_message(
            sprintf(
                'Multiple directories for run %s %s (%s) not supported!',
                $self->run_name,
                $self->lane,
                $self->genome_model_run_id,
            )
        );
        return;
    }
    elsif ( keys %dirs == 0 ) {
        $self->error_message(
            sprintf(
                'No directories for run %s %s (%s)',
                $self->run_name,
                $self->lane,
                $self->id,
            )
        );
        return;
    }

    my ($full_path) = keys %dirs;
    $full_path .= '/' unless $full_path =~ m|\/$|;

    return $full_path;
}

#< Dump to File System >#
sub dump_to_file_system {
    #$self->warning_message("Method 'dump_data_to_file_system' not implemented");
    return 1;
}

sub dump_illumina_fastq_files {
    return Genome::InstrumentData::Solexa::dump_illumina_fastq_files(@_);
}

sub dump_solexa_fastq_files {
    return Genome::InstrumentData::Solexa::dump_solexa_fastq_files(@_);
}

sub dump_sanger_fastq_files {
    return Genome::InstrumentData::Solexa::dump_sanger_fastq_files(@_);
}

sub _unprocessed_fastq_filenames {
    my $self = shift;
    my @fastqs;
    if ($self->is_external) {
        @fastqs = $self->resolve_external_fastq_filenames(@_);
    } else {
        @fastqs = @{$self->resolve_fastq_filenames(@_)};
    }
    return @fastqs;
}

sub desc {
    return Genome::InstrumentData::Solexa::desc(@_);
}

sub read1_fastq_name {
    return Genome::InstrumentData::Solexa::read1_fastq_name(@_);
}

sub read2_fastq_name {
    return Genome::InstrumentData::Solexa::read2_fastq_name(@_);
}

sub fragment_fastq_name {
    return Genome::InstrumentData::Solexa::fragment_fastq_name(@_);
}

sub resolve_fastq_filenames {
    my $self = shift;
    my $lane = $self->lane;
    my $desc = $self->desc;

    my %params = @_;
    my $paired_end_as_fragment = delete $params{'paired_end_as_fragment'};
    my $requested_directory = delete $params{'directory'} || Genome::Sys->base_temp_directory;

    my @illumina_output_paths;
    my @errors;

    # First check the archive directory and second get the gerald directory
    for my $dir_type (qw(archive_path gerald_directory)) {
        $self->status_message("Now trying to get fastq from $dir_type for $desc");

        my $directory = $self->$dir_type;
        $directory = $self->validate_fastq_directory($directory, $dir_type);
        next unless $directory;

        if ($dir_type eq 'archive_path') {
            $directory = $self->dump_illumina_fastq_archive($requested_directory); #need eval{} here ?
            $directory = $self->validate_fastq_directory($directory, 'dump_fastq_dir');
            next unless $directory;
        }

        eval {
            #handle fragment or paired-end data
            if ($self->is_paired_end) {
                if (!$paired_end_as_fragment || $paired_end_as_fragment == 1) {
                    if (-e "$directory/" . $self->read1_fastq_name) {
                        push @illumina_output_paths, "$directory/" . $self->read1_fastq_name;
                    }
                    elsif (-e "$directory/Temp/" . $self->read1_fastq_name) {
                        push @illumina_output_paths, "$directory/Temp/" . $self->read1_fastq_name;
                    }
                    else {
                        die "No illumina forward data in directory for lane $lane! $directory";
                    }
                }
                if (!$paired_end_as_fragment || $paired_end_as_fragment == 2) {
                    if (-e "$directory/" . $self->read2_fastq_name) {
                        push @illumina_output_paths, "$directory/" . $self->read2_fastq_name;
                    }
                    elsif (-e "$directory/Temp/" . $self->read2_fastq_name) {
                        push @illumina_output_paths, "$directory/Temp/" . $self->read2_fastq_name;
                    }
                    else {
                        die "No illumina reverse data in directory for lane $lane! $directory";
                    }
                }
            }
            else {
                if (-e "$directory/" . $self->fragment_fastq_name) {
                    push @illumina_output_paths, "$directory/" . $self->fragment_fastq_name;
                }
                elsif (-e "$directory/Temp/" . $self->fragment_fastq_name) {
                    push @illumina_output_paths, "$directory/Temp/" . $self->fragment_fastq_name;
                }
                else {
                    die "No fragment illumina data in directory for lane $lane! $directory";
                }
            }
        };

        push @errors, $@ if $@;
        last if @illumina_output_paths;
    }
    unless (@illumina_output_paths) {
        $self->error_message("No fastq files were found for $desc");
        $self->error_message(join("\n",@errors)) if @errors;
        die $self->error_message;
    }
    return \@illumina_output_paths;
}


sub dump_illumina_fastq_archive {
    my ($self, $dir) = @_;

    my $archive = $self->archive_path;
    $dir = Genome::Sys->base_temp_directory unless $dir;

    #Prevent unarchiving multiple times during execution
    #Hopefully nobody passes in a $dir expecting to overwrite another set of FASTQs coincidentally from the same lane number
    my $already_dumped = 0;

    if($self->is_paired_end) {
        if (-s $dir . '/' . $self->read1_fastq_name and -s $dir . '/' . $self->read2_fastq_name) {
            $already_dumped = 1;
        }
    }
    else {
        if (-s $dir . '/' . $self->fragment_fastq_name) {
            $already_dumped = 1;
        }
    }

    unless($already_dumped) {
        my $cmd = "tar -xzf $archive --directory=$dir";
        unless (Genome::Sys->shellcmd(
            cmd => $cmd,
            input_files => [$archive],
        )) {
            $self->error_message('Failed to run tar command '. $cmd);
            return;
            #die($self->error_message); Should try to get fastq from gerald_directory instead of dying
        }
    }
    return $dir;
}

sub validate_fastq_directory {
    my ($self, $dir, $dir_type) = @_;

    my $msg_base = "$dir_type : $dir";

    unless ($dir) {
        $self->error_message("$msg_base is null");
        return;
    }

    unless (-e $dir) {
        $self->error_message("$msg_base not existing in file system");
        return;
    }

    unless ($dir_type eq 'archive_path') {
        my @files = glob("$dir/*"); #In scalar context, a glob functions as an iterator--we instead want to check the number of files
        unless (scalar @files) {
            $self->error_message("$msg_base is empty");
            return;
        }
    }

    return $dir;
}


sub resolve_external_fastq_filenames {
    my $self = shift;

    my @fastq_pathnames;
    my $fastq_pathname = $self->create_temp_file_path('fastq');
    unless ($fastq_pathname) {
        die "Failed to create temp file for fastq!";
    }
    return ($fastq_pathname);
}

sub _calculate_total_read_count {
    my $self = shift;

    if($self->is_external) {
        my $data_path_object = Genome::MiscAttribute->get(entity_id => $self->id, property_name=>'full_path');
        my $data_path = $data_path_object->value;
        my $lines = `wc -l $data_path`;
        return $lines/4;
    }
    if ($self->clusters <= 0) {
        die('Impossible value '. $self->clusters .' for clusters field for solexa lane '. $self->id);
    }

    return $self->clusters;
}

sub resolve_quality_converter {

    # old stuff needed sol2sanger, new stuff all uses sol2phred, but
    # we dont care what the version is anymore

    my $self = shift;

    my %analysis_software_versions = (
                                     'GAPipeline-0.3.0'       => 'sol2sanger',
                                     'GAPipeline-0.3.0b1'     => 'sol2sanger',
                                     'GAPipeline-0.3.0b2'     => 'sol2sanger',
                                     'GAPipeline-0.3.0b3'     => 'sol2sanger',
                                     'GAPipeline-1.0'         => 'sol2sanger',
                                     'GAPipeline-1.0-64'      => 'sol2sanger',
                                     'GAPipeline-1.0rc4'      => 'sol2sanger',
                                     'GAPipeline-1.1rc1p4'    => 'sol2sanger',
                                     'SolexaPipeline-0.2.2.5' => 'sol2sanger',
                                     'SolexaPipeline-0.2.2.6' => 'sol2sanger',
                                 );

    my $analysis_software_version = $self->analysis_software_version;
    unless ($analysis_software_version) {
        die('No analysis_software_version found for instrument data '. $self->id);
    }

    return $analysis_software_versions{$analysis_software_version} || 'sol2phred';
}

sub resolve_adaptor_file {
    my $self = shift;

    #these are constants and should probably be defined in class properties...TODO
    my $dna_primer_file = '/gscmnt/sata114/info/medseq/adaptor_sequences/solexa_adaptor_pcr_primer';
    my $rna_primer_file = '/gscmnt/sata114/info/medseq/adaptor_sequences/solexa_adaptor_pcr_primer_SMART';

    my $adaptor_file;
    if ( $self->sample_type eq 'rna' ) {
        $adaptor_file = $rna_primer_file;
    }
    else {
        $adaptor_file = $dna_primer_file;
    }
    unless (-f $adaptor_file) {
        $self->error_message('Specified adaptor file'. $adaptor_file .' does not exist.');
        die($self->error_message);
    }
    return $adaptor_file;
}

sub create_mock {
    my $class = shift;
    my $self = $class->SUPER::create_mock(@_);
    return unless $self;

    for my $method (qw/
        dump_sanger_fastq_files
        resolve_fastq_filenames
        _calculate_total_read_count
        resolve_adaptor_file
        run_identifier
    /) {
        my $ref = $class->can($method);
        die "Unknown method $method on " . $class . ".  Cannot make a pass-through for mock object!" unless $ref;
        $self->mock($method,$ref);
    }

    return $self;
}

sub run_start_date_formatted {
    my $self = shift;

    my ($y, $m, $d) = $self->run_name =~ m/^(\d{2})(\d{2})(\d{2})_.*$/;

    my $dt_format = UR::Time->config('datetime');
    #UR::Time->config(datetime => '%a %b %d %T %Z %Y');
    UR::Time->config(datetime => '%Y-%m-%d');
    my $dt = UR::Time->numbers_to_datetime(0, 0, 0, $d, $m, "20$y");
    UR::Time->config(datetime => $dt_format);

    return $dt;
}

sub total_bases_read {
    my $self = shift;
    my $filter = shift; # optional?
    if(!defined($filter))
    {
        $filter = 'both';
    }
    my $total_bases; # unused?


    my $count;
    if ($self->is_paired_end) {
        # this changed in case we only want the fwd or rev counts...
        $count += ($self->fwd_read_length * $self->fwd_clusters)  unless $filter eq 'reverse-only';
        $count += ($self->rev_read_length * $self->rev_clusters) unless $filter eq 'forward-only';
    } else {
        $count += ($self->read_length * $self->clusters);
    }

    return $count;
}

sub summary_xml_content {
    my $self = shift;
    my $rls = $self->_run_lane_solexa;
    unless ($rls) { return; }
    return $rls->summary_xml_content;
}

sub run_identifier {
    my $self = shift;
    return $self->flow_cell_id;
}

1;

#$HeaderURL$
#$Id: Solexa.pm 61055 2010-07-16 19:30:48Z boberkfe $
