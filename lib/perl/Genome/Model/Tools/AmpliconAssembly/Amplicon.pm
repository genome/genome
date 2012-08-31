package Genome::Model::Tools::AmpliconAssembly::Amplicon;

use strict;
use warnings;

use Bio::SeqIO;
use Bio::Seq::Quality;
use Carp 'confess';
use Data::Dumper 'Dumper';
use File::Grep 'fgrep';
require Genome::Utility::MetagenomicClassifier::SequenceClassification;
use Storable;

class Genome::Model::Tools::AmpliconAssembly::Amplicon {
    is => 'UR::Object',
    has => [
    name => {
        is => 'Text',
        doc => 'Name of amplicon',
    },
    directory => {
        is => 'Text',
        doc => 'Directory where the amplicon fasta and qual files are located.',
    },
    ],
    has_many => [
    reads => {
        is => 'Text',
        doc => 'Reads for the amplicon',
    },
    ],
};

#:jpeck The code looks clean.  

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_)
        or return;

    for my $attr (qw/ name reads /) {
        next if $self->$attr;
        $self->error_message("Attribute ($attr) is required to create");
        $self->delete;
        return;
    }

    unless ( Genome::Sys->validate_existing_directory( $self->directory ) ) {
        $self->delete;
        return;
    }

    return $self;
}

sub DESTROY {
    my $self = shift;

    $self->{_ace_factory}->disconnect if $self->{_ace_factory};
    
    return $self->SUPER::DESTROY;
}

sub _fatal_msg {
    my ($self, $msg) = @_;

    confess ref($self)." - ERROR: $msg\n";
}

#< Successfully Assembled Reqs >#
sub successful_assembly_length { successfully_assembled_length(@_); }
sub successfully_assembled_length {
    return 1150;
}

sub successfully_assembled_read_count {
    return 2;
}

sub successfully_assembled_requirements_as_string {
    return 'length >= 1150, reads >= 2';
}

#< Basic Accessors >#
sub get_name {
    return $_[0]->name;
}

sub get_directory {
    return $_[0]->directory;
}

sub get_reads {
    return [ $_[0]->reads ];
}

sub get_read_count {
    return scalar(@{$_[0]->get_reads});
}

#< Files >#
# name.fasta
# name.fasta.preclip
# name.fasta.view
# name.fasta.ace  
# name.fasta.prescreen
# name.phds
# name.fasta.contigs 
# name.fasta.problems
# name.reads.fasta
# name.fasta.contigs.qual
# name.fasta.problems.qual
# name.reads.fasta.qual
# name.fasta.log  
# name.fasta.qual  
# name.scfs
# name.fasta.memlog 
# name.fasta.qual.preclip
# name.fasta.phrap.out
# name.fasta.singlets
sub _file {
    return $_[0]->get_directory.'/'.$_[0]->get_name.'.'.$_[1];
}

#< Cleanup >#
sub remove_unneeded_files {
    my $self = shift;

    my @unneeded_file_exts = (qw/
        fasta.view fasta.log fasta.singlets
        fasta.problems fasta.problems.qual
        fasta.phrap.out fasta.memlog
        fasta.preclip fasta.qual.preclip 
        fasta.prescreen fasta.qual.prescreen
        scfs phds
        /);
    for my $ext ( @unneeded_file_exts ) {
        my $file = $self->_file($ext);
        #print "$file\n";
        unlink $file if -e $file;
    }

    return 1;
}

#< Aux Files #>
sub scfs_file { # .scfs
    return _file($_[0], 'scfs');
}

sub create_scfs_file {
    my $self = shift;

    my $scfs_file = $self->scfs_file;
    unlink $scfs_file if -e $scfs_file;
    my $scfs_fh = Genome::Sys->open_file_for_writing($scfs_file)
        or return;
    for my $scf ( @{$self->get_reads} ) { 
        $scfs_fh->print("$scf\n");
    }
    $scfs_fh->close;

    if ( -s $scfs_file ) {
        return $scfs_file;
    }
    else {
        unlink $scfs_file;
        return;
    }
}

sub phds_file { # .phds
    return $_[0]->_file('phds');
}

sub classification_file {
    return $_[0]->_file('classification.stor');
}

#< Sequence/Qual Files >#
sub fasta_file_for_type {
    my ($self, $type) = @_;
    my $method = $type.'_fasta_file';
    return $self->$method;
}

sub qual_file_for_type {
    my ($self, $type) = @_;
    my $method = $type.'_qual_file';
    return $self->$method;
}

sub fasta_file { # main fasta file 
    return $_[0]->_file('fasta');
}

sub qual_file { # main qual file
    return $_[0]->_file('fasta.qual');
}

sub reads_fasta_file { # .fasta - after phred
    return $_[0]->_file('reads.fasta');
}

sub reads_qual_file { # .fasta.qual - after phred
    return $_[0]->_file('reads.fasta.qual');
}

sub processed_fasta_file { # processed.fasta - after screen and trim
    return $_[0]->_file('fasta');
    #return $_[0]->_file('processed.fasta');
}

sub processed_qual_file { # processed.fasta.qual - after screen and trim
    return $_[0]->_file('fasta.qual');
    #return $_[0]->_file('processed.fasta.qual');
}

sub assembly_fasta_file { return contigs_fasta_file(@_); }
sub contigs_fasta_file { # .fasta.contigs - fasta file of assembly 
    return $_[0]->_file('fasta.contigs');
}

sub assembly_qual_file { return contigs_qual_file(@_); }
sub contigs_qual_file { # .fasta.contigs.qual - qual file of assembly
    return $_[0]->_file('fasta.contigs.qual');
}

sub ace_file { # .fasta.ace - ce file from phrap
    return $_[0]->_file('fasta.ace');
}

sub oriented_fasta_file { # .oriented.fasta - post assembly, then oriented
    return $_[0]->_file('oriented.fasta');
}

sub oriented_qual_file { # .oriented.fasta.qual - post assembly, then oriented
    return $_[0]->_file('oriented.fasta.qual');
}

#< Assembled Sequence/Qual and Info >#
sub get_bioseq { # gets the oriented or assembly bioseq and sets the other info about it
    return $_[0]->_get_bioseq_info->{bioseq};
}

sub get_bioseq_source {
    return $_[0]->_get_bioseq_info->{source};
}

sub get_assembled_reads {
    return $_[0]->_get_bioseq_info->{assembled_reads};
}

sub get_assembled_read_count {
    return scalar(@{$_[0]->get_assembled_reads});
}

sub was_assembled_successfully {
    return ( @{$_[0]->get_assembled_reads} ) ? 1 : 0;
}

sub is_bioseq_oriented {
    return $_[0]->_get_bioseq_info->{oriented};
}

sub get_assembly_bioseq { # this just gets the assembly bioseq, and doesn't set the bioseq info
    my $self = shift;

    my %assembly_bioseq_info = $self->_get_bioseq_info_from_assembly
        or return;

    return $assembly_bioseq_info{bioseq};
}

sub get_oriented_bioseq { # this just gets the oriented bioseq, and doesn't set the bioseq info
    my $self = shift;

    my %oriented_bioseq_info = $self->_get_bioseq_info_from_oriented_fasta
        or return;

    return $oriented_bioseq_info{bioseq};
}

sub _get_bioseq_info {
    my $self = shift;

    unless ( $self->{_bioseq_info} ) {
        my %info;
        # 2009apr30 No longer getting bioseq if not assembled and has contig length > 1150
        for my $method (qw/ _get_bioseq_info_from_oriented_fasta _get_bioseq_info_from_assembly _get_bioseq_info_from_nothing /) {
            %info = $self->$method;
            last if %info;
        }

        $self->{_bioseq_info} = \%info;
    }

    return $self->{_bioseq_info};
}

sub confirm_orientation {
    my ($self, $complement) = @_;

    my $bioseq = $self->get_bioseq;
    return unless $bioseq;

    if ( $complement ) {
        my $revcom_bioseq;
        eval { $revcom_bioseq = $bioseq->revcom; };
        unless ( $revcom_bioseq ) {
            print "$!\n";
            $self->_fatal_msg("Can't reverse complement bioseq for amplicon: ".$self->get_name);
        }
        $bioseq = $revcom_bioseq;
    }

    my $fasta_file = $self->oriented_fasta_file;
    unlink $fasta_file if -e $fasta_file;
    my $fasta_o = Bio::SeqIO->new(
        '-format' => 'fasta',
        '-file' => "> $fasta_file",
    )
        or $self->_fatal_msg("Can't open fasta file ($fasta_file)");
    $fasta_o->write_seq($bioseq);

    my $qual_file = $self->oriented_qual_file;
    unlink $qual_file if -e $qual_file;
    my $qual_o = Bio::SeqIO->new(
        '-format' => 'qual',
        '-file' => "> $qual_file",
    )
        or $self->_fatal_msg("Can't open qual file ($qual_file)");
    $qual_o->write_seq($bioseq);

    $self->{_bioseq_info}->{bioseq} = $bioseq;
    $self->{_bioseq_info}->{oriented} = 1;
    
    return 1;
}
    
sub _get_bioseq_info_from_oriented_fasta {
    my $self = shift;

    # Fasta
    my $fasta_file = $self->oriented_fasta_file;
    return unless -s $fasta_file;
    my $fasta_reader = Bio::SeqIO->new(
        '-file' => "< $fasta_file",
        '-format' => 'fasta',
    )
        or return;
    my $fasta = $fasta_reader->next_seq
        or return;

    # 2009apr30 no longer supported, all are from assembly, this will remove ones from previous dirs
    #  if source is read - remove oriented file?
    #  if length is too short, skip
    # Check if source is read
    my ($source_string) = $fasta->desc =~ /source=(\S+)/;
    if ( not defined $source_string or $source_string eq 'read' ) {
        #print $self->get_name." is from read\n";
        #unlink $fasta_file;
        #unlink $self->oriented_qual_file;
        return;
    }
    # Check length
    unless ( $fasta->length >= $self->successfully_assembled_length ) {
        #print $self->get_name." is too short\n";
        #unlink $fasta_file;
        #unlink $self->oriented_qual_file;
        return;
    }

    my ($read_string) = $fasta->desc =~ /reads=(\S+)/;
    #my ($read_string) = $fasta->desc =~ /reads=([\w\d\.\-\_\,]+)/;
    unless ( $read_string ) {
        $self->_fatal_msg('No reads found in description of oriented fasta for amplicon: '.$self->get_name);
    }
    my @reads = split(',', $read_string);
    #print Dumper([$self->get_name, $fasta->desc, \@reads]);
    
    # Qual
    my $qual_file = $self->oriented_qual_file;
    my $qual_reader = Bio::SeqIO->new(
        '-file' => "< $qual_file",
        '-format' => 'qual',
    )
        or return;
    my $qual = $qual_reader->next_seq
        or return;

    $self->_validate_fasta_and_qual_bioseq($fasta, $qual)
        or return;
    
    return (
        bioseq => Bio::Seq::Quality->new(
            '-id' => $self->get_name,
            '-desc' => "source=oriented reads=$read_string",
            '-alphabet' => 'dna',
            '-force_flush' => 1,
            '-seq' => $fasta->seq,
            '-qual' => $qual->qual,
        ), 
        assembled_reads => \@reads,
        source => 'oriented',
        oriented => 1,
    );
}

sub get_successfully_assembled_contig_from_assembly {
    my $self = shift;

    my $acefile = $self->ace_file;
    return unless -s $acefile; # ok
    my $ace_reader = Genome::Model::Tools::Consed::AceReader->create(
        file => $acefile,
    );
    return if not $ace_reader;

    my $contig;
    while ( my $obj = $ace_reader->next ) {
        next if $obj->{type} ne 'contig';
        next unless $obj->{read_count} > 1;
        $obj->{unpadded_consensus} = $obj->{consensus};
        $obj->{unpadded_consensus} =~ s/\*//g;
        next if length $obj->{unpadded_consensus} < $self->successfully_assembled_length;
        $contig = $obj;
        last;
    }

    return if not $contig;

    while ( my $obj = $ace_reader->next ) {
        last if $obj->{type} eq 'contig';
        next if $obj->{type} ne 'read';
        push @{$contig->{read_names}}, $obj->{name};
    }

    return $contig;
}

sub get_reads_from_successfully_assembled_contig {
    my $self = shift;

    my $contig = $self->get_successfully_assembled_contig_from_assembly;
    return if not $contig;

    return $contig->{read_count};
}

sub _get_bioseq_info_from_assembly { # was _get_bioseq_info_from_longest_contig
    my $self = shift;

    my $contig = $self->get_successfully_assembled_contig_from_assembly
        or return;

    my @read_names = sort { $a cmp $b } @{$contig->{read_names}};

    # Bioseq
    my $bioseq = Bio::Seq::Quality->new(
        '-id' => $self->get_name,
        '-desc' => 'source=assembly reads='.join(',', @read_names),
        '-alphabet' => 'dna',
        '-force_flush' => 1,
        '-seq' => $contig->{unpadded_consensus},
        '-qual' => join(' ', @{$contig->{base_qualities}}),
    );

    $self->_validate_fasta_and_qual_bioseq($bioseq, $bioseq)
        or return;

    return (
        bioseq => $bioseq,
        assembled_reads => \@read_names,
        source => 'assembly',
        oriented => 0,
    );
}

sub _get_bioseq_info_from_nothing {
    my $self = shift;

    return (
        bioseq => undef,
        assembled_reads => [],
        oriented => 0,
        source =>undef, 
    );
}

#< Read Bioseq >#
# all
sub get_bioseqs_for_raw_reads {
    my ($self, $read_name) = @_;

    return $self->_get_bioseqs_for_reads('reads');
}

sub get_bioseqs_for_processed_reads {
    my ($self, $read_name) = @_;

    return $self->_get_bioseqs_for_reads('processed');
}

sub _get_bioseqs_for_reads {
    my ($self, $type) = @_;

    my $fasta_file_method = $type.'_fasta_file';
    my $fasta_file = $self->$fasta_file_method;
    unless ( -e $fasta_file ) { # ok
        print "No fasta file ($fasta_file) for type ($type)\n";
        return;
    }
    
    my $qual_file_method = $type.'_qual_file';
    my $qual_file = $self->$qual_file_method;
    unless ( -e $qual_file ) { # not ok
        $self->_fatal_msg("Found fasta file, but not quality file for type ($type)");
        return;
    }

    my @bioseqs;
    my $fasta_reader = Bio::SeqIO->new(
        '-file' => "< $fasta_file",
        '-format' => 'fasta',
    );
    while ( my $fasta = $fasta_reader->next_seq ) {
        push @bioseqs, Bio::Seq::Quality->new(
            '-id' => $fasta->id,
            '-desc' => 'source='.$type,
            '-seq' => $fasta->seq,
        );
    }
    unless ( @bioseqs ) { # not ok
        $self->error_message("No reads found in fasta file ($fasta_file) for type ($type)\n");
        return;
    }

    my $qual_reader = Bio::SeqIO->new(
        '-file' => "< $qual_file",
        '-format' => 'qual',
    );
    while ( my $qual = $qual_reader->next_seq ) {
        # should be in order, but just in case
        my ($bioseq) = grep { $_->id eq $qual->id } @bioseqs;
        $self->_validate_fasta_and_qual_bioseq($bioseq, $qual)
            or return;
        $bioseq->qual( $qual->qual );
    }

    return @bioseqs;
}

# by read
sub get_bioseq_for_raw_read {
    my ($self, $read_name) = @_;

    $self->_validate_read_name($read_name);
    
    return $self->_get_bioseq_for_read_and_type($read_name, 'reads');
}

sub get_bioseq_for_processed_read {
    my ($self, $read_name) = @_;

    $self->_validate_read_name($read_name);
    
    return $self->_get_bioseq_for_read_and_type($read_name, 'processed');
}

sub _validate_read_name {
    my ($self, $read_name) = @_;

    unless ( defined $read_name ) {
        $self->_fatal_msg("No read name given");
    }

    unless ( grep { $_ eq $read_name } @{$self->get_reads} ) {
        $self->_fatal_msg("Read ($read_name) is not in the list of attempted reads: ".join(', ', @{$self->get_reads}));
    }
    
    return 1;
}

sub _get_bioseq_for_read_and_type {
    my ($self, $read_name, $type) = @_;

    my $fasta_file_method = $type.'_fasta_file';
    my $fasta_file = $self->$fasta_file_method;
    unless ( -e $fasta_file ) { # ok
        print "No fasta file ($fasta_file) for type ($type)\n";
        return;
    }
    
    my $qual_file_method = $type.'_qual_file';
    my $qual_file = $self->$qual_file_method;
    unless ( -e $qual_file ) { # not ok
        $self->_fatal_msg("Found fasta file, but not quality file for type ($type)");
        return;
    }

    my ($fasta, $qual);
    my $fasta_reader = Bio::SeqIO->new(
        '-file' => "< $fasta_file",
        '-format' => 'fasta',
    );
    while ( $fasta = $fasta_reader->next_seq ) {
        last if $fasta->id eq $read_name;
    }
    unless ( $fasta ) { # ok
        print "Fasta not found for read ($read_name) and type ($type)\n";
        return;
    }
    
    my $qual_reader = Bio::SeqIO->new(
        '-file' => "< $qual_file",
        '-format' => 'qual',
    );
    while ( $qual = $qual_reader->next_seq ) {
        last if $qual->id eq $read_name;
    }
    unless ( $qual ) { # not ok
        $self->_fatal_msg("Found fasta, but not quality for read ($read_name) and type ($type)");
    }

    $self->_validate_fasta_and_qual_bioseq($fasta, $qual)
        or return;

    return Bio::Seq::Quality->new(
        '-id' => $read_name,
        '-desc' => 'source='.$type,
        '-seq' => $fasta->seq,
        '-qual' => $qual->qual,
    );
}

#< Classification >#
sub get_classification {
    my $self = shift;

    my $classification_file = $self->classification_file;
    return unless -s $classification_file;

    return retrieve($classification_file);
}

sub save_classification {
    my ($self, $classification) = @_;

    unless ( $classification ) {
        $self->_fatal_msg('No classification to save for amplicon: '.$self->get_name);
    }
    
    my $classification_file = $self->classification_file;
    unlink $classification_file if -e $classification_file;
    Storable::nstore($classification, $classification_file);
    
    return 1;
}

#< Bio::Seq Helpers >#
sub _validate_fasta_and_qual_bioseq {
    my ($self, $fasta, $qual) = @_;

    #print "Validating ".$fasta->id.' and '.$qual->id."\n";
    
    unless ( $fasta->seq =~ /^[ATGCNX]+$/i ) {
        $self->error_message(
            sprintf(
                "Illegal characters found in fasta (%s) seq:\n%s",
                $fasta->id,
                $fasta->seq,
            )
        );
        return;
    }

    unless ( $fasta->length == $qual->length ) {
        $self->error_message(
            sprintf(
                'Unequal length for fasta (%s) and quality (%s)',
                $fasta->id,
                $qual->id,
            )
        );
        return;
    }
    
    return 1;
}

1;

=pod

=head1 Name

Genome::Model::Tools::AmpliconAssembly::Amplicon

=head1 Synopsis

A package for interacting with amplicons in amplicon assembly models 

=head1 Usage

Although an amplicon can be created directly, it is best to use the 'get_amplicons' method on a amplicon assembly build:

Create directly:

 my $amplicon = Genome::Model::Tools::AmpliconAssembly::Amplicon->new(
    name => 'HMPB-aaa01a01', # base name of the amplicon
    directory => $build->edit_dir, # dir w/ fasta and quals for the amplicon
    reads => [qw/ HMPB-aaa01a01.b1 HMPB-aaa01a01.b2 HMPB-aaa01a01.g1 /], # the trace names for the amplicon that were used in attempting to assemble
 );

 ...


=head1 Attribute Getters

These are getters for the attributes that are used to create the amplicon.

=over

=item B<get_name>

=item B<get_directory>

=item B<get_reads> 

=item B<get_read_count> - simple read count

=back

=head1 Files

The methods get the file names

=over

=item B<scfs_file>

=item B<create_scfs_file>

=item B<phds_file>

=item B<classification_file>

=item B<fasta_file_for_type>

=item B<qual_file_for_type>

=item B<reads_fasta_file>

=item B<reads_qual_file>

=item B<processed_fasta_file>

=item B<processed_qual_file>

=item B<assembly_fasta_file> - same as contigs_fasta_file

=item B<assembly_qual_file> - same as contigs_qual_file

=item B<contigs_fasta_file>

=item B<contigs_qual_file>

=item B<ace_file>

=item B<oriented_fasta_file>

=item B<oriented_qual_file>

=back

=head1 Main Bioseq Methods

The methods refer to the main bioseq methods, which when called we set the bioseq and it's properties.  These methods try to first get the oriented fasta/qual, but if that does not exist, the ace assembly will be used.  If the amplicon reads failed to assemble, or do not meet the successfully_assembled_length, no bioseq info will be returned.

=head2 get_bioseq

=over

=item I<Synopsis>   Trys to get the assembled bioseq from oriented fasta, then the assembly fasta.  If not assembled no bioseq is returned.

=item I<Arguments>  none

=item I<Returns>    bioseq object (Bio::Seq::Quality)

=back

=head2 get_bioseq_source

=over

=item I<Synopsis>   The source of the main bioseq.  Currently oriented or assembly

=item I<Arguments>  none

=item I<Returns>    source (string)

=back

=head2 get_assembled_reads

=over

=item I<Synopsis>   The reads used in successful assembly of the amplicon

=item I<Arguments>  none

=item I<Returns>    read names (array ref of strings)

=back


=head2 get_assembled_read_count

=over

=item I<Synopsis>   The number of reads used in successful assembly of the amplicon

=item I<Arguments>  none

=item I<Returns>    read count (integer)

=back

=head2 was_assembled_successfully

 if ( $amplicon->was_assembled_successfully ) {
    ...
 }

=over

=item I<Synopsis>   Tells wheter or not the amplicon was assembled successfully

=item I<Arguments>  none

=item I<Returns>    boolean (true for success)

=back

=head2 is_bioseq_oriented

=over

=item I<Synopsis>   Tells if the bioseq was oriented

=item I<Arguments>  none

=item I<Returns>    boolean (true for success)

=back

=head1 Bioseq Methods for Speicifc Types

=head2 get_assembly_bioseq

 my $bioseq = $amplicon->get_assembly_bioseq;

=over

=item I<Synopsis>   Gets the assembly bioseq, straight from the acefile, if assembled and if the assembly meets the successful assembly length.

=item I<Arguments>  none

=item I<Returns>    bioseq object (Bio::Seq::Quality)

=back

=head2 get_oriented_bioseq

 my $bioseq = $amplicon->get_oriented_bioseq;

=over

=item I<Synopsis>   Gets the oriented bioseq, if it exists

=item I<Arguments>  none

=item I<Returns>    bioseq object (Bio::Seq::Quality)

=back

=head2 get_bioseq_for_raw_read

 my $bioseq = $amplicon->get_bioseq_for_raw_readi($read_name);

=over

=item I<Synopsis>   Gets the bioseq for a raw, unprocessed read

=item I<Arguments>  read name (string)

=item I<Returns>    bioseq object (Bio::Seq::Quality)

=back

=head2 get_bioseqs_for_raw_reads

 my @bioseqs = $amplicon->get_bioseqs_for_raw_reads;

=over

=item I<Synopsis>   Gets the oriented bioseq, if it exists

=item I<Arguments>  none

=item I<Returns>    array of bioseq objects (Bio::Seq::Quality)

=back

=head2 get_bioseq_for_processed_read

 my $bioseq = $amplicon->get_bioseq_for_processed_read($read_name);

=over

=item I<Synopsis>   Gets the bioseq for a processed read

=item I<Arguments>  read name (string)

=item I<Returns>    bioseq object (Bio::Seq::Quality)

=back

=head2 get_bioseqs_for_processed_reads

 my @bioseqs = $amplicon->get_bioseqs_for_processed_reads;

=over

=item I<Synopsis>   Gets the bioseqs for all of the processed reads

=item I<Arguments>  none

=item I<Returns>    array of bioseq objects (Bio::Seq::Quality)

=back

=head1 Classification

=head2 get_classification

=over

=item I<Synopsis>   Gets the classification

=item I<Arguments>  none

=item I<Returns>    classification object (Genome::Utility::MetagenomicClassifier::SequenceClassification)

=back

=head2 save_classification

=over

=item I<Synopsis>   Saves the classification

=item I<Arguments>  classification object (Genome::Utility::MetagenomicClassifier::SequenceClassification)

=item I<Returns>    boolean (true for success)

=back

=head1 Misc

=head2 successfully_assembled_length

=over

=item I<Synopsis>   Gets the length that an assembly must be greater than or eqaul to to be considered a success

=item I<Arguments>  none

=item I<Returns>    1150 (integer)

=back

=head2 confirm_orientation

 if ( $amplicon->confirm_orientation ) {
    print "Yay! Saved oriented fasta for amplicon\n";
 }

=over

=item I<Synopsis>   Confirms the orientation of the amplicon.  This will then create 'oriented' fasta and quality files.

=item I<Arguments>  complement (boolean) - whether or not the assembly bioseq must be complemented before saving.

=item I<Returns>    boolean (true for success)

=back

=head2 remove_unneeded_files
 
 $amplicon->remove_unneeded_files;

=over

=item I<Synopsis>   Removes the redundant, uneeded files for an amplicon from the directory

=item I<Arguments>  none

=item I<Returns>    boolean (true for success)

=back

=head1 See Also

=over

=item B<Genome::Model::Tools::AmpliconAssembly>

=back

=head1 Disclaimer

Copyright (C) 2005 - 2009 Genome Center at Washington University in St. Louis

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied
 warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

=head1 Author(s)

B<Eddie Belter> I<ebelter@genome.wustl.edu>

=cut

#$HeadURL$
#$Id$
