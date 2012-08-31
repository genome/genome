package Genome::Model::Tools::Pcap::Phd;

our $VERSION = 0.01;
my $pkg = 'Genome::Model::Tools::Pcap::Phd';

use strict;
use warnings;

use base qw(Class::Accessor::Fast);

use Carp;
use Genome::Model::Tools::Pcap::Phd::Reader;
use Genome::Model::Tools::Pcap::Phd::Writer;
use IO::File;
use IO::String;
use Compress::Zlib;
use IO::String;
use Storable;

Genome::Model::Tools::Pcap::Phd->mk_accessors(qw(_reader _writer _input _input_directory _input_file _index using_db));

sub new
{
    croak("$pkg:new:no class given, quitting") if @_ < 1;
    my ($caller, %params) = @_;
	my $caller_is_obj = ref($caller);
    my $class = $caller_is_obj || $caller;
    my $self = {};
    bless ($self, $class);
	
	if(exists $params{input})
	{
		$self->_input ( $params{input});		
	}
	elsif(exists $params{input_file}&&defined $params{input_file}&&-e $params{input_file})
	{
		$self->_input_file ($params{input_file});
		$self->_input(IO::File->new($self->_input_file));
	}
	
	if(exists $params{input_directory})
	{
        my $dir = $params{input_directory};
        $dir .= '/' unless $dir =~ /\/$/;
		$self->_input_directory($dir);	
	}
		
	if(!defined($self->_input) && 
	   !defined($self->_input_directory))
	{
		die "Need either a phd.ball, or a phd_dir\n";
	}
	
	if(exists $params{using_db})
	{
		$self->using_db($params{using_db});	
	}
	else
	{
		$self->using_db(0);
	}
	
	$self->_reader( Genome::Model::Tools::Pcap::Phd::Reader->new );
	$self->_writer( Genome::Model::Tools::Pcap::Phd::Writer->new );
	if(exists $params{index_file}&& defined $params{index_file} && -e $params{index_file})
	{
		my $fh = IO::File->new($params{index_file});
		$self->_load_index_from_file($fh);
	}
	elsif(exists $params{index} && defined $params{index})
	{
		my $fh = $params{index};
		$self->_load_index_from_file($fh);
    }
    elsif(defined $self->_input_file && -e $self->_input_file.".index" && $self->using_db)
    {
        $self->_index(Storable::retrieve($self->_input_file.".index"));
    }
	elsif(defined $self->_input)
	{
		$self->_build_index();
        if($self->using_db && !-e $self->_input_file.".index")
        {
            Storable::nstore $self->_index, $self->_input_file.".index";
        }
	}

    return $self;
}

sub _build_index
{
	my ($self) = @_;	

	my $input = $self->_input;
	my %index;
	my $temp=[qw(0)];
	while(<$input>)
	{
       	if(/BEGIN_SEQUENCE/)
        {
            my @tokens = split / /;
            my $offset = tell ($input) - length ($_);
            chomp $tokens[1];
			$$temp[1] = $offset-$$temp[0];#calculate length
			$temp = [$offset];
			$index{$tokens[1]} = $temp;                
        }
	}
	$$temp[1] = ((stat ($input))[7])-$$temp[0];#calculate length
	$self->_index(\%index);

}

sub _load_index_from_file
{
	my ($self, $fh) = @_;
		
	my %index;
	my $input = $self->_input;
	my $temp=[qw(0)];
	while(<$fh>)
	{
		chomp;
		my @tokens = split / /;
		$$temp[1] = $tokens[1]-$$temp[0];#calculate length
		$temp = [$tokens[1]];
		$index{$tokens[0]} = $temp;
	}
	$$temp[1] = ((stat ($input))[7])-$$temp[0];#calculate length
	$self->_index(\%index);		
}

sub get_phd_size
{
	my ($self, $name) = @_;
	return $self->_index->{$name}[1];
}

sub phd_exists
{
	my ($self, $name) = @_;
	return exists $self->_index->{$name};
}

sub get_phd_names
{
	my ($self) = @_;

	return [keys %{$self->_index}];
}

sub get_phd
{
	my ($self, $name) = @_;
	
	my $phd_ball = $self->_input;
	my $index = $self->_index;
	my $reader = $self->_reader;
	my ($read_name, $version) = ($name =~ /(.+)\.phd\.(\d+)$/);

	if(exists $$index{$read_name}&&$version==1)
	{
		$phd_ball->seek($$index{$read_name}[0],0);
		#my $phd_string;
		#read $phd_ball, $phd_string, $$index{$read_name}[1];	
		#my $fh = IO::String->new($phd_string);
		return $reader->read($phd_ball);
	}
	elsif(defined $self->_input_directory && -e $self->_input_directory.$name)
	{
		my $fh = IO::File->new($self->_input_directory.$name);
		my $phd = $reader->read($fh);
		close $fh;
		return $phd;
	}
	elsif(defined $self->_input_directory && -e $self->_input_directory.$name.".gz")
	{
		my $fh = IO::File->new($self->_input_directory.$name.".gz");
		my $sPhdContent = join( '', <$fh> );
		my $PHDFILEDATA = Compress::Zlib::memGunzip($sPhdContent);
		my $PHDFILE = new IO::String($PHDFILEDATA);
		my $phd = $reader->read($PHDFILE);
		close $fh;
		close $PHDFILE;
		return $phd;
	}	
	else
	{
		die "Could not locate $name\n";
	}	
}

sub get_latest_phd
{
	my ($self, $read_name) = @_;
	
    die "Need read name to get next phd name\n" unless defined $read_name;

    my $phd_name = $self->get_latest_phd_name($read_name);
    
    return unless defined $phd_name and -s $self->_input_directory . $phd_name;
    
	return $self->get_phd($phd_name);
}

sub get_latest_phd_name
{
	my ($self, $read_name) = @_;
	
    die "Need read name to get next phd name\n" unless defined $read_name;

    my $int = $self->_get_lastest_iteration_for_read_name($read_name);

    return unless defined $int;
    
	return "$read_name.phd.$int";
}

sub get_next_phd_name
{
	my ($self, $read_name) = @_;
    
    die "Need read name to get next phd name\n" unless defined $read_name;

    my $int = $self->_get_lastest_iteration_for_read_name($read_name);

    $int++; # should increment to 1 if not defined
    
	return "$read_name.phd.$int";
}

sub get_phd_version_names
{
    my ($self, $read_name) = @_;

    die "Need read name to get next phd name\n" unless defined $read_name;

    my @version_names;

    my $int = $self->_get_lastest_iteration_for_read_name($read_name);

    return unless defined $int;
    
    push @version_names, "$read_name.phd.$_" for (1..$int);
    return @version_names;    
}

sub _get_lastest_iteration_for_read_name
{
    my ($self, $read_name) = @_;

    die "Need read name to get lastest iteration\n" unless defined $read_name;
    
	opendir my $DIR, $self->_input_directory;
	my @all_files = readdir $DIR;
    closedir $DIR;

	my $int = 0;
	for (my $i = 0; $i < @all_files; $i++)
	{
		if ($all_files[$i] =~ /$read_name/)
		{
			my ($ext) = $all_files[$i] =~ /phd\.(\d+)$/;
            next unless defined $ext;
            $int = $ext if $ext > $int;
        }
	}
    
    return if $int == 0; # This means the file was not found

	return $int;
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

1;		

=pod

=head1 NAME

 Genome::Model::Tools::Pcap::Phd
 
  > Object oriented phd/phd.ball file reader/writer

=head1 SYNOPSIS

 my $phd_object = Genome::Model::Tools::Pcap::Phd->new
 (
    input_directory => "inputdirname",
 );

 my @phd_names = $phd_object->get_phd_names();
 my $phd = $phd_object->get_phd("vef07");

    
=head1 DESCRIPTION

Genome::Model::Tools::Pcap::Phd takes either a Phd file, and allows the user to get Contig objects from the ace file, edit them, and write the file back to the hard disk when finished.

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

#$HeadURL$
#$Id$
