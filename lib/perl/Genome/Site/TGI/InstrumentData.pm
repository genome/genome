package Genome::Site::TGI::InstrumentData;

use strict;
use warnings;

use Genome;

use Data::Dumper;
require Genome::Sys;

class Genome::Site::TGI::InstrumentData {
    id_by => ['id'],
    table_name => <<EOS
    (
         SELECT to_char(solexa.analysis_id) id,
               fc.run_name,
               'Genome::Site::TGI::InstrumentData::Solexa' subclass_name,
               'solexa' sequencing_platform,
               (
                    case
                        when solexa.index_sequence is null then to_char(solexa.lane)
                        else to_char(solexa.lane) || '-' || solexa.index_sequence
                    end
               ) subset_name,
               library_id
          FROM index_illumina solexa
          JOIN flow_cell_illumina fc on fc.flow_cell_id = solexa.flow_cell_id
     UNION ALL
            SELECT
               to_char(case when ri.index_sequence is null then ri.region_id else ri.seq_id end) id,
               r.run_name,
               'Genome::Site::TGI::InstrumentData::454' subclass_name,
               '454' sequencing_platform,
               (
                case
                    when ri.index_sequence is null then to_char(r.region_number)
                    else to_char(r.region_number) || '-' || ri.index_sequence
                end
               ) subset_name,
               (
                case
                    when ri.index_sequence is null then r.library_id
                    else ri.library_id
                end
               ) library_id
           FROM run_region_454 r
           JOIN region_index_454 ri on ri.region_id = r.region_id
    ) idata
EOS
    ,
    is_abstract => 1,
    subclassify_by => 'subclass_name',
    has => [
        subclass_name       => { is => 'Text', len => 255 },
        sequencing_platform => { is => 'Text', len => 255 },
        run_name            => { is => 'Text', len => 500, is_optional => 1 },
        subset_name         => { is => 'Text', len => 32, is_optional => 1, },
        
        # TODO: see if this stuff is used and if not delete it -ss
        full_name => { calculate_from => ['run_name','subset_name'], calculate => q|"$run_name/$subset_name"| },        
    ],
    has_optional => [        
        library_id          =>  { is => 'VARCHAR2', len => 15 },
        library             =>  { is => 'Genome::Library', id_by => 'library_id' },
        library_name        =>  { via => 'library', to => 'name' },

        # Library id is overridden in some subclasses, having this alias allows the actual db value to be accessed in that case.
        library_summary_id  =>  { calculate_from => 'library_id', calculate => q{ return $library_id } },

        sample_id           =>  { is => 'Number', via => 'library' },
        sample              =>  { is => 'Genome::Site::TGI::Sample', id_by => 'sample_id' },
        sample_name         =>  { via => 'sample', to => 'name' },
	
        # TODO: see if this stuff is used and if not delete it -ss
        events => { is => 'Genome::Model::Event', is_many => 1, reverse_id_by => "instrument_data" },
        # TODO: see if this stuff is used and if not delete it -ss
        full_path => {
            via => 'attributes',
            to => 'value', 
            where => [ property_name => 'full_path' ],
            is_mutable => 1,
        },
    ],
    has_many_optional => [
        attributes => {
            is => 'Genome::MiscAttribute',
            reverse_as => '_instrument_data',
        },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::Dwrac',
};

#< UR Methods >#
sub delete {
    my $self = shift;

    for my $attr ( $self->attributes ) {
        $attr->delete;
    }
    $self->SUPER::delete;

    return $self;
}

sub _resolve_subclass_name {
	my $class = shift;

	if (ref($_[0]) and $_[0]->isa(__PACKAGE__)) {
		my $sequencing_platform = $_[0]->subclass_name;
		return $class->_resolve_subclass_name_for_sequencing_platform($sequencing_platform);
	}
    elsif (my $sequencing_platform = $class->define_boolexpr(@_)->value_for('subclass_name')) {
        return $class->_resolve_subclass_name_for_sequencing_platform($sequencing_platform);
    }
	else {
		return;
	}
}

sub seq_id { shift->id }

sub _resolve_subclass_name_for_sequencing_platform {
    my ($class,$sequencing_platform) = @_;
    my @type_parts = split(' ',$sequencing_platform);

    my @sub_parts = map { ucfirst } @type_parts;
    my $subclass = join('',@sub_parts);

    my $class_name = join('::', 'Genome::Site::TGI::InstrumentData' , $subclass);
    return $class_name;
}

sub _resolve_sequencing_platform_for_subclass_name {
    my ($class,$subclass_name) = @_;
    my ($ext) = ($subclass_name =~ /Genome::Site::TGI::InstrumentData::(.*)/);
    return unless ($ext);
    my @words = $ext =~ /[a-z]+|[A-Z](?:[A-Z]+|[a-z]*)(?=$|[A-Z])/g;
    my $sequencing_platform = lc(join(" ", @words));
    return $sequencing_platform;
}

#< Paths >#
sub create_data_directory_and_link {
    my $self = shift;

    my $data_path = $self->resolve_full_path;
    Genome::Sys->create_directory($data_path)
        or return;
    
    my $link = $self->data_link;
    unlink $link if -l $link;
    Genome::Sys->create_symlink($data_path, $link)
        or return;

    return $data_path;
}

sub _links_base_path {
    return '/gscmnt/839/info/medseq/instrument_data_links/';
}

sub data_link {
    return sprintf('%s/%s', _links_base_path(), $_[0]->id);
}

sub _data_base_path {
    return '/gscmnt/sata835/info/medseq/instrument_data/';
}

sub resolve_full_path{
    my $self = shift;

    return $self->full_path if $self->full_path;

    return $self->full_path( $self->_default_full_path );
}

sub _default_full_path {
    my $self = shift;
    sprintf('%s/%s', $self->_data_base_path, $self->id)
}

#< Dump to File System >#
sub dump_to_file_system {
    my $self = shift;
    $self->warning_message("Method 'dump_data_to_file_system' not implemented");
    return 1;
}

sub allocations {
    my $self = shift;

    my @allocations = Genome::Disk::Allocation->get(
                                                    owner_class_name => $self->class,
                                                    owner_id => $self->id,
                                                );
    return @allocations;
}

sub calculate_alignment_estimated_kb_usage {
    my $self = shift;
    return;
}

sub sample_type {
    my $self = shift;
    my $dna_type;
    my @dna = GSC::DNA->get(dna_name => $self->sample_name);
    if (@dna == 1) {
        
        if ($dna[0]->dna_type eq 'genomic dna') {
            return 'dna';
        } elsif ($dna[0]->dna_type eq 'pooled dna') {
            return 'dna';
        } elsif ($dna[0]->dna_type eq 'rna') {
            return 'rna';
        }
    }
    return;
}


######
# Needed??
# WHY NOT USE RUN_NAME FROM THE DB????
sub old_name {
    my $self = shift;

    my $path = $self->full_path;

    my($name) = ($path =~ m/.*\/(.*EAS.*?)\/?$/);
    if (!$name) {
	   $name = "run_" . $self->id;
    }
    return $name;
}

sub create_mock {
    my $class = shift;
    return $class->SUPER::create_mock(subclass_name => 'Genome::Site::TGI::InstrumentData', @_);
}

sub run_identifier  {
    die "run_identifier not defined in instrument data subclass.  please define this. this method should " . 
         "provide a unique identifier for the experiment/run (eg flow_cell_id, ptp barcode, etc).";
}

sub dump_fastqs_from_bam {
    my $self = shift;
    my %p = @_;

    die "cannot call bam path" if (!$self->can('bam_path'));
    
    unless (-e $self->bam_path) {
	$self->error_message("Attempted to dump a bam but the path does not exist:" . $self->bam_path);
	die $self->error_message;
    }
    
    my $temp_dir = Genome::Sys->create_temp_directory('unpacked_bam');

    my $subset = (defined $self->subset_name ? $self->subset_name : 0);

    my %read_group_params;

    if (defined $p{read_group_id}) {
        $read_group_params{read_group_id} = delete $p{read_group_id};
        $self->status_message("Using read group id " . $read_group_params{read_group_id});
    } 

    my $fwd_file = sprintf("%s/s_%s_1_sequence.txt", $temp_dir, $subset);
    my $rev_file = sprintf("%s/s_%s_2_sequence.txt", $temp_dir, $subset);
    my $fragment_file = sprintf("%s/s_%s_sequence.txt", $temp_dir, $subset);
    my $cmd = Genome::Model::Tools::Picard::SamToFastq->create(input=>$self->bam_path, fastq=>$fwd_file, fastq2=>$rev_file, fragment_fastq=>$fragment_file, no_orphans=>1, %read_group_params);
    unless ($cmd->execute()) {
        die $cmd->error_message;
    }

    if ((-s $fwd_file && !-s $rev_file) ||
        (!-s $fwd_file && -s $rev_file)) {
        $self->error_message("Fwd & Rev files are lopsided; one has content and the other doesn't. Can't proceed"); 
        die $self->error_message;
    }

    my @files;
    if (-s $fwd_file && -s $rev_file) { 
        push @files, ($fwd_file, $rev_file);
    }
    if (-s $fragment_file) {
        push @files, $fragment_file;
    }
   
    return @files; 
}

1;

#$HeadURL$
#$Id$
