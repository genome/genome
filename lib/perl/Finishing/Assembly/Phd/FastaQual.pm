package Finishing::Assembly::Phd::Schema;

use strict;
use warnings;
no warnings 'reserved';

use Finfo::Std;

use Finishing::Assembly::Cache;
use Finishing::Assembly::Phd::Reader;
#use Finishing::Assembly::Phd::Writer;
use IO::Dir;
use IO::File;
use IO::String;
use Compress::Zlib;
use IO::String;
use Storable;

#(_reader _writer _input _input_directory _input_file _index conserve_memory);
#
my %input :name(input:r);
my %in_type :name(input_type:r) :isa([ 'in_list', __PACKAGE__->input_types ]);

my %cache :name(_cache:p);

sub connect
{
    my ($class, $input, $type) = @_;

    return $class->new
    (
        input => $input,
        input_type => $type,
    );
}

sub input_types
{
    return (qw/ directory ball fasta_qual /);
}

sub reader
{
    return Finishing::Assembly::Phd::Reader->instance;
}

sub START
{
    my $self = shift;

    if ( $self->input_type eq 'directory' )
    {
        $self->_build_cache_for_directory;
    }
    elsif ( $self->input_type eq 'fasta_qual' )
    {
        $self->_connect_to_index;
    }
    else
    {
        $self->_build_cache_for_ball;
    }
    
    return $self;
}

sub _split_phd_name : PRIVATE
{
    my ($self, $name) = @_;

    return split(/\.phd\./, $name);
}

sub _build_cache_for_directory : PRIVATE
{
    my $self = shift;

    my $dh = IO::Dir->new($self->input)
        or $self->fatal_msg( sprintf( 'Can\'t open directory (%s): %s', $self->input, $!) );
    my $cache = Finishing::Assembly::Cache->new();
    while ( my $file = $dh->read )
    {
        next unless $file =~ /\.phd\./;
        my $name = File::Basename::basename($file);
        #$name =~ s/\.gz//;
        $cache->set($name, 'offset', 0);
    }

    $self->_cache($cache);

    return 1;
}

sub _build_cache_for_ball : PRIVATE
{
    my $self = shift;

    my $fh = IO::File->new('<' . $self->input)
        or $self->fatal_msg( sprintf( 'Can\'t open file (%s): %s', $self->input, $!) );
    my $cache = Finishing::Assembly::Cache->new();

	my $temp=[qw(0)];
	while( my $line = $fh->getline )
	{
        next unless $line =~ /BEGIN_SEQUENCE/;
        chomp $line;
        my @tokens = split(/ /, $line);
        my $offset = $fh->tell - length($line);
        $cache->set($tokens[1], 'offset', $offset);
    }

    $self->_cache($cache);

    return 1;
}

#- INFO -#
sub phd_names
{
    my $self = shift;

    return $self->_cache->ids;
}

sub phd_exists
{
	my ($self, $name) = @_;

	return $self->_cache->id_exists($name);
}

#- GET PHD -#
sub get_phd
{
	my ($self, $name) = @_;
	
    my $file;
    if ( $self->input_type eq 'directory' )
    {
        $file = sprintf('%s/%s', $self->input, $name);
        # handle gzip??
        #my $fh = IO::File->new($self->_input_directory.$name.".gz");
        #my $sPhdContent = join( '', <$fh> );
        #my $PHDFILEDATA = Compress::Zlib::memGunzip($sPhdContent);
        #my $PHDFILE = new IO::String($PHDFILEDATA);
        #my $phd = $reader->read($PHDFILE);
        #close $fh;
        #close $PHDFILE;
    }
    else
    {
        $file = $self->input;
    }

    my $fh = IO::File->new("< $file")
        or $self->fatal_msg( sprintf( 'Can\'t open file (%s): %s', $file, $!) );
    my $position = $self->_cache->get($name, 'offset');
    $fh->seek($position, 0);
    my $phd = $self->reader->execute($fh);
    $fh->close;

    #return Finihing::Assembly::Phd->new(%phd);
    return $phd;
}

sub get_latest_phd
{
	my ($self, $name) = @_;
	
    $self->fatal_msg("Need name to get next phd name") unless $name;

    my $phd_name = $self->get_latest_phd_name($name);
    
    return unless defined $phd_name and -s $self->_input_directory . $phd_name;
    
	return $self->get_phd($phd_name);
}

sub get_latest_phd_name
{
	my ($self, $name) = @_;
	
    $self->fatal_msg("Need read name to get latest phd name") unless $name;

    my $iteration = $self->_get_lastest_iteration($name);

    return unless $iteration;
    
	return "$name.phd.$iteration";
}

sub get_next_phd_name
{
	my ($self, $name) = @_;
    
    $self->fatal_msg("Need read name to get next phd name") unless $name;

    my $iteration = $self->_get_lastest_iteration($name);

	return sprintf('%s.phd.%s', $name, ++$iteration);
}

sub get_phd_version_names
{
    my ($self, $name) = @_;

    $self->fatal_msg("Need read name to get next phd name") unless $name;

    my $iteration = $self->_get_lastest_iteration($name);

    return unless $iteration;
    
    return map { "$name.phd.$_" } (1..$iteration);
}

sub _get_lastest_iteration
{
    my ($self, $name) = @_;

    $self->fatal_msg("Need read name to get lastest iteration") unless $name;
    
    my $iteration = 0;
    while ( 1 )
    {
        $iteration++;
        my $phd_name = sprintf('%s.phd.%d', $name, $iteration);
        last unless defined $self->_cache->get($phd_name, 'offset');
	}

    $iteration--; # zero means not found
    
	return $iteration;
}

sub write_phd
{
	my ($self, $phd, $phd_name) = @_;

	my $phd_dir = $self->_input_directory;
	my $writer = $self->_writer;
	
    die "Need a phd_dir to write new phd files to.\n" unless -d $phd_dir;
	
    $phd_name = $self->get_next_phd_name( $phd->name ) unless defined $phd_name;
    
    die "Need a phd anem to write phd info to" unless defined $phd_name;
    
	$writer->write( IO::File->new('> ' . $phd_dir . $phd_name), $phd );
}

sub write_phd_ball
{
	my ($self, %params) = @_;
	
	my $writer = $self->_writer;
	
	if(exists $params{output})
	{
		$self->_output ( $params{output});				
	}
	elsif(exists $params{output_file})
	{
		$self->_output_file ($params{output_file});
		$self->_output(IO::File->new(">".$self->_output_file));
	}
	elsif(defined $self->_output)
	{
		$self->_output->seek(0,0) or die "Could not seek to beginning of write file\n";		
	}
	else
	{
		die "Could not find find file to write to.\n";
	}
	
	#my ($dir_name, $file_name) = $output_file =~ /(.+\?)(.+)/;
	my $fh = IO::File->new(">/tmp/phd.ball");
	my $index = $self->_index;
	my $input = $self->_input;
	foreach my $value (values %$index)
	{
		if($value->[0] == -1)
		{
			$writer->write($fh, $value->[1]);		
		}
		else
		{
			#write out phd
		    my $phd_string;
    		$input->seek($value->[0],0);
    		$input->read($phd_string, $value->[1]);
    		print $fh $phd_string;
		}	
	}
	$fh->seek(0,0);
	$input->close;
	my $output = $self->_output;
	$output->seek(0,0);
	while(<$fh>)
	{
		print $output $_;
	
	}
	
}

sub add_phd
{
	my ($self, $phd) = @_;

	my $index = $self->_index;
	my ($read_name, $version) = $phd->name =~ /(.+)\.phd\.\d+$/;
	
	if((defined $index && $version == 1) )
	{
		$$index{$read_name} = [-1, $phd];	
	}
	elsif(-e $self->_input_directory)
	{
		$self->_writer->write(">IO::File->new($self->_input_directory$phd->name)",$phd);
		
	}
	else
	{
		die "Could not find appropriate place to write Phd file.\n";
	}	
	
}


sub rename_phds_and_traces_to_read_names
{
# TODO put in phd schema
    my ($self, %p) = @_;
    
    $self->info_msg("rename_phds_and_traces_to_read_names not fully tested, but should work.");
    
    my $base_dir = delete $p{base_dir};

    my $phd_dir = $base_dir . '/phd_dir';
    return unless Finfo::Validate->validate
    (
        attr => 'phd_dir',
        value => $phd_dir,
        isa => 'dir_w',
        msg => 'fatal',
    );

    my $phd_obj = Finishing::Assembly::Phd->new(input_directory => $phd_dir);

    my $chromat_dir = $base_dir . '/chromat_dir';
    return unless Finfo::Validate->validate
    (
        attr => 'chromat_dir',
        value => $chromat_dir,
        isa => 'dir_w',
        msg => 'fatal',
    );

    my $current_name_cb = delete $p{current_name_cb};
    return unless Finfo::Validate->validate
    (
        attr => 'current name callback',
        value => $current_name_cb,
        isa => 'code',
        msg => 'fatal',
    );

    $self->fatal_msg("Can't create Phd object") unless $phd_obj;

    foreach my $contig_name ( $self->contig_names )
    {
        my $contig = $self->get_contig($contig_name);

        my $reads = $contig->reads;
        foreach my $new_name ( keys %$reads )
        {
            my $current_name = $current_name_cb->($new_name);

            next unless defined $current_name;
            
            my $current_trace_name = sprintf('%s/%s.gz', $chromat_dir, $current_name);
            my $new_trace_name = sprintf('%s/%s.gz', $chromat_dir, $new_name);
            unless ( -e $new_trace_name )
            {
                #$self->info_msg("move($current_trace_name, $new_trace_name)");
                unless ( move($current_trace_name, $new_trace_name) )
                {
                    $self->error_msg("Can't move $current_trace_name to $new_trace_name\: $!");
                    next;
                }
            }
            
            my $ext = $phd_obj->get_lastest_iteration_for_read_name($current_name);
            $self->info_msg("Could not get phd name for read name ($current_name)")
                and next unless $ext;

            my $current_phd_name = sprintf('%s/%s.phd.%d', $phd_dir, $current_name, $ext);
            my $new_phd_name = sprintf('%s/%s.phd.%d', $phd_dir, $new_name, $ext);

            next if -e $new_phd_name;

            #$self->info_msg("move($current_phd_name, $new_phd_name)");
            unless ( move($current_phd_name, $new_phd_name) )
            {
                $self->error_msg("Can't move $current_phd_name to $new_phd_name\: $!");
                next;
            }
        }
    }

    return 1;
}

1;		

=pod

=head1 NAME

 Finishing::Assembly::Phd
 
  > Object oriented phd/phd.ball file reader/writer

=head1 SYNOPSIS

 my $phd_object = Finishing::Assembly::Phd->new
 (
    input_directory => "inputdirname",
 );

 my @phd_names = $phd_object->get_phd_names();
 my $phd = $phd_object->get_phd("vef07");

    
=head1 DESCRIPTION

Finishing::Assembly::Phd takes either a Phd file, and allows the user to get Contig objects from the ace file, edit them, and write the file back to the hard disk when finished.

=head1 METHODS

=head2 get_phd_names

=head2 get_phd

=head2 get_latest_phd

=head2 get_latest_phd_name

=head2 get_next_phd_name

=head2 write_phd

=head2 write_phd_ball

=head2 add_phd

=head1 Disclaimer

=head1 Author(s)

=cut

#$HeadURL: svn+ssh://svn/srv/svn/gscpan/perl_modules/branches/jschindl/Assembly/Phd.pm $
#$Id: Phd.pm 20568 2007-01-09 04:47:08Z jschindl $
