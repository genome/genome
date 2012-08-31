package Genome::Model::Tools::Pcap::Ace;

our $VERSION = 0.01;

use strict;
use warnings;
use Carp;
use base qw(Class::Accessor::Fast);
use Genome::Model::Tools::Pcap::Ace::Reader;
use Genome::Model::Tools::Pcap::Ace::Writer;
use Genome::Model::Tools::Pcap::Item;
use Genome::Model::Tools::Pcap::SequenceItem;
use Genome::Model::Tools::Pcap::Contig;
use Genome::Model::Tools::Pcap::Read;
use Genome::Model::Tools::Pcap::Tag;
use Genome::Model::Tools::Pcap::TagParser;
use Genome::Model::Tools::Pcap::Sources::Ace::Contig;
use Genome::Model::Tools::Pcap::Config;
use Genome::Model::Tools::Pcap::FHManager;
use File::Temp;
use IO::File;
use IO::String;
use Storable;

use File::Basename;
use Cwd 'abs_path';

no warnings 'once';
local $Storable::Deparse = 1;
local $Storable::Eval = 1;
#local $Storable::forgive_me = 1;

Genome::Model::Tools::Pcap::Ace->mk_accessors(qw(_reader _writer _input_files _show_progress _cache_dir _output _output_file)); #_input removed, ability to parse stream to a index/db will be a future feature, as needed



sub new {
    croak("__PACKAGE__:new:no class given, quitting") if @_ < 1;
    my ($caller, %params) = @_;
    my $caller_is_obj = ref($caller);
    my $class = $caller_is_obj || $caller;
    my $self = {};
    bless ($self, $class);
    if( exists $params{input_file})
    {
        push @{$params{input_files}}, $params{input_file};
    }
    if (exists $params{input_files})
    {
        @{$params{input_files}} = map { abs_path($_); } @{$params{input_files}};
	}    
    $self->_reader ( Genome::Model::Tools::Pcap::Ace::Reader->new);

#TODO:make it check for dir before building index
    $self->_init_cache_dir(%params);
    
    return $self;
}

sub _init_cache_dir
{
    my ($self,%params) = @_;
    if(defined $params{cache_dir})
    {
        $self->_cache_dir($params{cache_dir});
        return 1 if(-d $self->_cache_dir);
        #since the cache dir does not exist, we'll need to create it and initialize it below
        `mkdir -p $params{cache_dir}`;
    }
    else
    {        
        $self->{temp_cache} = File::Temp->newdir('/tmp/temp_ace_cacheXXXXX');
        $self->_cache_dir($self->{temp_cache}->dirname);
    }
    my $cache_dir = $self->_cache_dir;   
    `touch $cache_dir/VERSION_1_0`;
    if(defined $params{input_files})
    {
        foreach my $input_file (@{$params{input_files}})
        {
            if (! -d $input_file.'.idx')
            {
                my $result;
                #TODO: make the index either go to a temp location, or derive itself from the acefile
                $result = $self->_build_ace_index(file_name => $input_file, dir => $input_file.'.idx') ;
                print "There was an error indexing the ace file\n" and return unless $result;
            }
            #now add the ace file's index to the cache
            my @files = `ls $input_file.idx/*.ci`;
            chomp @files;
            foreach my $contig_index (@files)
            {
                `ln -s $contig_index $cache_dir/.`;
            }
            `touch $input_file.idx/assembly_tags` if (! -e "$input_file.idx/assembly_tags");
            `cat $input_file.idx/assembly_tags >> $cache_dir/assembly_tags`;            
        }
    }       
}

sub _build_index
{
    _build_ace_index(@_);
}

sub _build_ace_index
{
	my ($self, %params) = @_;
    #TODO: get this working with file handles too
    my $file_name = $params{file_name};
    my $fh = Genome::Model::Tools::Pcap::FHManager->get_fh($file_name);
	my %contigs;
	my @assembly_tag_indexes;
    my $old_hash = undef;
    my $contig = undef;
	my $co_count = 0;
	my $index_dir;
    my $get_assembly_tags = $params{get_assembly_tags};
    $index_dir = $params{dir};
    if(! -d $index_dir) { `mkdir -p $index_dir`;}
    my $reader = $self->_reader;
    $reader->{input} =$fh;
    `touch $index_dir/VERSION_1_0`;
    
    my %tag_hash;
	while(my $line = <$fh>)
	{
		my $first_three = substr($line, 0,3);        
		if($first_three eq "CO ")
		{
			$contig = $self->_build_contig_index($file_name,$line);
            Storable::nstore $contig, "$index_dir/$contig->{name}.ci";
		}
		elsif($first_three eq "WA ")
		{
			my $offset = (tell $fh) - length $line;
            $old_hash = { offset => $offset, length => length($line) };
            next unless $get_assembly_tags;            
			push (@assembly_tag_indexes, $old_hash);
		}
		elsif($first_three eq "CT{")
		{
			my $offset = (tell $fh) - length $line;
			my $first_length = length $line;
			$line = <$fh>;
			$line =~ s/^\s*// if $line =~ /\w/;
			my @tokens = split(/[ {]/,$line);	            
			$old_hash = {offset => $offset, length => (length($line) + $first_length)};
			my $contig;
			$contig = $contigs{$tokens[0]};

			if(!defined $tag_hash{$tokens[0]}) {$tag_hash{$tokens[0]} = [];}
            push (@{$tag_hash{$tokens[0]}}, $old_hash); 
		}				        
        else
        {            
            $old_hash->{length} += length($line) if(defined($old_hash));
        }
	}
    
    #TODO: add contig tags to contigs
    foreach my $contig_name (keys %tag_hash)
    {
        my $ci = retrieve "$index_dir/$contig_name.ci";
        push @{$ci->{tags}},@{$tag_hash{$contig_name}};
        delete $tag_hash{$contig_name};
        Storable::nstore $ci, "$index_dir/$contig_name.ci";
    }
    
    #write assembly tags to file
    my $fh_assembly_tags = IO::File->new(">>$index_dir");
	foreach my $assembly_tag_index (@assembly_tag_indexes)
	{		
		$fh->seek($assembly_tag_index->{offset},0);
		my $tag_string;
        $fh->read($tag_string, $assembly_tag_index->{length});
        print $fh_assembly_tags $tag_string;
	}
    
    return 1;    
}

sub _build_contig_index
{
	my ($self, $file_name, $line, $file_offset) = @_;
	my %contigs;
	my @assembly_tag_indexes;
    my $old_hash = undef;
    my $contig = undef;
    my $first_bs = 1;
    my $found_bq = 0;
    my $found_qa = 0;
    my $found_af = 0;
    my $found_rd = 0;
	my $co_count = 0;
    my $fh = Genome::Model::Tools::Pcap::FHManager->get_fh($file_name);
    $fh->seek($file_offset,0) if defined $file_offset;
    
	
	#build contig data structures
	my @tokens = split(/[ {]/,$line);
	my $offset = (tell $fh) - length $line;            

	$old_hash = { offset => $offset,
                      name => $tokens[1],
                      read_count => $tokens[3],
                      base_segments => {line_count => $tokens[4]},#\@base_segments,
                      contig_tags => [] ,
                      contig_loaded => 0 ,
					  contig_indexed => 1,
                      length => length($line),
					  sequence_length => $tokens[2],
                      contig_length => length($line),
                      file_name => $file_name
                      };
    $contig = $old_hash;    
    $found_af = 0;
    $found_rd = 0;
	$fh->seek($tokens[2],1);            							  
		
	while(my $line = <$fh>)
	{
		my $first_three = substr($line, 0,3);        
		if($first_three eq "BS ")
		{           
            $offset = (tell $fh) - length $line;
            $old_hash = undef;            
            $contig->{base_segments}{offset} = $offset;
			for(my $i = 0;$i<$contig->{base_segments}{line_count};$i++)
			{
				my $line = <$fh>;			
			}            
			
			my $offset2 = tell $fh;
	        $contig->{base_segments}{length} = $offset2 - $contig->{base_segments}{offset};
            
		}
		elsif($first_three eq "AF ")
		{
			if($found_bq == 1)
            {
                $contig->{base_qualities}{length} = (tell $fh) - length($line) - $contig->{base_qualities}{offset};
                $found_bq = 0;
            }
			my @tokens = split(/[ {]/,$line);
            my $end = (tell $fh);			
			$offset = $end - length $line;            
            $old_hash = { name => $tokens[1], offset => $offset };
            $contig->{reads}{$tokens[1]}{read_position}= $old_hash;
            $contig->{af_start} = $offset;
			for(my $i=1;$i<$contig->{read_count};$i++)
			{
				my $line = <$fh>;
				my @tokens = split(/[ {]/,$line);
            	my $end = (tell $fh);			
				my $offset = $end - length $line;            
            	$old_hash = { name => $tokens[1], offset => $offset };
            	$contig->{reads}{$tokens[1]}{read_position}= $old_hash;				
            }		
			
            $contig->{af_end} = (tell $fh);
            
		}
		elsif($first_three eq "RD ")
		{
			my @tokens = split(/[ {]/,$line);
			$offset = (tell $fh) - length $line;
			if(!$found_rd)
			{
				$contig->{rd_start} = $offset ;			
            	$found_rd = 1;
			}
			$old_hash = { offset => $offset,
                          name => $tokens[1],
                          read_tags => [],
                          length => length($line)};
                        
                        
            $old_hash->{sequence}{offset} = (tell $fh)+1;
			$fh->seek($tokens[2],1);
			$contig->{reads}{$tokens[1]}{read} = $old_hash;			
			
		}        
		elsif($first_three eq "CO ")
		{
            $offset = (tell $fh) - length $line;
			$contig->{contig_length} = $offset - $contig->{offset};
			$contig->{rd_end} = $offset;
			$fh->seek(- length($line),1);
			last;
		
		}		
		elsif(substr($first_three,0,2) eq "BQ")
        {
            $offset = (tell $fh) - length $line;
            $contig->{base_qualities}{offset} = $offset; 
            $contig->{base_sequence}{length} = ($offset-1)-$contig->{offset};
			$found_bq = 1;
            #It looks like this was intended to get though the base qual lines
            #quickly by skipping the seq_length * 2 seek positions but this
            #causes problems when ratio of pads to bases is more then 1/3 which
            #causes seek pos to end up in the middle AF lines thus preventing 
            #proper parsing of AF lines
			#$fh->seek($contig->{sequence_length}*2,1);
			
        }
		elsif($first_three eq "WA ")
		{
			$offset = (tell $fh) - length $line;
            $contig->{contig_length} = $offset - $contig->{offset};
			$contig->{rd_end} = $offset;
            $fh->seek(-length($line),1);
			last;
		}
		elsif($first_three eq "CT{"  )
		{
			$offset = (tell $fh) - length $line;
            $contig->{contig_length} = $offset - $contig->{offset};
			$contig->{rd_end} = $offset;
			$fh->seek(-length ($line),1);
			last;
			
		}
		elsif($first_three eq "DS ")
		{
			$offset = (tell $fh) - length $line;
        	$old_hash->{ds}{offset} = $offset;  
        	$old_hash->{length} = $offset - $old_hash->{offset};							
		}
		elsif($first_three eq "QA ")
		{
        	$offset = (tell $fh) - length $line;
        	$old_hash->{qa}{offset} = $offset;
        	$old_hash->{sequence}{length} = ($offset - 1) - $old_hash->{sequence}{offset};
		    $old_hash->{length} = $offset - $old_hash->{offset};
        }
		elsif($first_three eq "RT{")
		{
			$offset = (tell $fh) - length $line;
			my $first_length = length $line;
			$line = <$fh>;
			my @tokens = split(/[ {]/,$line);            
			$old_hash = {offset => $offset, length => (length($line)+ $first_length)};
			push (@{$contig->{reads}{$tokens[0]}{read}{read_tags}} ,$old_hash );	
			
		}		        
        else
        {            
            $old_hash->{length} += length($line) if(defined($old_hash));
        }
	}
    

    $offset =  (stat($fh))[7] if eof $fh;        
    $contig->{contig_length} = $offset - $contig->{offset};
	$contig->{rd_end} = $offset;
	return $contig;

}

sub get_contig_names
{
    my $self = shift;
    
    my $cache_dir = $self->_cache_dir;
    my @contig_names = `ls $cache_dir/*.ci`;
    chomp @contig_names;
    @contig_names = map { /.*(Contig.*)\.ci/ } @contig_names;
       
    return [sort { _cmptemp($a, $b) } @contig_names ] ;
}

sub get_contig
{
    my ($self, $contig_name, $load) = @_;
    
    confess "No contig name given.\n" unless defined $contig_name;    
        
    my $contig_index;
    my $cache_dir = $self->_cache_dir;
    my $contig_file = "$cache_dir/$contig_name.ci";
        
    if(-e $contig_file)
    {
        $contig_index = retrieve $contig_file;        
    }                
    else
    {
        confess "Could not get contig for name: $contig_name\n";
    }
    
    if($contig_index->{contig_loaded})
    {          
        #my $contig = Storable::dclone($contig_index->{contig_object});
		my $contig = $contig_index->{contig_object};
        #TODO: make sure that the contig is able to get a handle to it's input file from the contig index, add the path to the ace file to the contig index when indexing the ace file
        $contig->thaw($self);
        return $contig;       
    }
    else 
    {
		if($load)
		{
			return $self->get_contig_old($contig_index);		
		}
		else
		{    
        	my $reader = $self->_reader;
        	my $contig_callback = Genome::Model::Tools::Pcap::Sources::Ace::Contig->new(name => $contig_index->{name},
        	index => $contig_index, reader => $self->_reader); 
            $contig_callback->thaw;#TODO: Figure out a less hacky way to do this   
        	return Genome::Model::Tools::Pcap::Contig->new(callbacks => $contig_callback);   
    	}
    }
}

sub get_contig_old
{
	my ($self, $contig_index) = @_;    		
    
    my $input = Genome::Model::Tools::Pcap::FHManager->get_fh($contig_index->{file_name});
    my $reader = $self->_reader;
    $reader->{input} = $input;
    my $ace_contig;
    my %reads; #contins read_positions and read_tags
    my @base_segments;
    my @contig_tags;
    my $result = $input->seek($contig_index->{offset},0);
    $ace_contig = $reader->next_object;
    #grab reads
    foreach my $read_index (values %{$contig_index->{reads}})
    {
        $input->seek($read_index->{read}{offset},0);
        my $ace_read = $reader->next_object;
        $reads{$ace_read->{name}} = Genome::Model::Tools::Pcap::Read->new(ace_read => $ace_read);
		#grab read_tags
		foreach my $read_tag_index (@{$read_index->{read}{read_tags}})
		{
			$input->seek($read_tag_index->{offset},0);
			my $read_tag = $self->_build_read_tag($reader->next_object);
			$reads{$read_tag->parent}->add_tag($read_tag);	
		}	
    }
	#grab read_positions
	foreach my $read_position_index (values %{$contig_index->{reads}})
	{
		$input->seek($read_position_index->{read_position}{offset},0);
		my $ace_read_position = $reader->next_object;
		$reads{$ace_read_position->{read_name}}->ace_read_position ($ace_read_position);	
	}	
		
	#grab contig_tags
	foreach my $contig_tag_index (@{$contig_index->{contig_tags}})
	{
		$input->seek($contig_tag_index->{offset},0);
		push @contig_tags, Genome::Model::Tools::Pcap::TagParser->new()->parse($input);
	}
	#grab base_segments
	$input->seek($contig_index->{base_segments}{offset},0);
	while(my $obj = $reader->next_object)
	{
		last if ($obj->{type} ne "base_segment");
		push @base_segments, $obj;	
	}	
	#glue everything together
    my $contig = Genome::Model::Tools::Pcap::Contig->new(ace_contig => $ace_contig,
                                                reads => \%reads,
                                                contig_tags => \@contig_tags,
                                                base_segments => \@base_segments);
	
	return $contig;	
}

sub add_contig
{
	my ($self, $contig) = @_;
	
	my $contig_index;
    $contig->freeze;
    my $cache_dir = $self->_cache_dir;
    my $contig_name = $contig->name;
    my $contig_file = "$cache_dir/$contig_name.ci";
#TODO: check to see if copy is still necessary, probably not
#TODO: add option to flush all data to a new ace file, and re-index it
#TODO: the cache should store frozen contig objects, not their indices alone
    $contig_index = { offset => -1, contig_loaded => 1, name => $contig->name, contig_object => $contig->copy($contig) };
    unlink $contig_file if(-e $contig_file);
    Storable::nstore $contig_index, $contig_file;

    $contig->thaw($self);
}

sub remove_contig
{
    my ($self, $contig_name) = @_;
#TODO: delete contig ace file and delete ci file
    my $cache_dir = $self->_cache_dir;
    my $contig_index_file = "$cache_dir/$contig_name.ci";
    if(-e $contig_index_file)
    {
        unlink $contig_index_file;    
    }
    my $contig_data_file = "$cache_dir/$contig_name.ace";
    if(-e $contig_data_file)
    {
        unlink $contig_data_file;
    }
}

sub get_assembly_tags
{
    my ($self) = @_;
    my $reader = $self->_reader;
#TODO:  parse tags file
    my $cache_directory = $self->_cache_dir;

    my @assembly_tags;
    my $assembly_tags_file = $cache_directory."/assembly_tags";
    
    return \@assembly_tags unless -s $assembly_tags_file; #No tags present

    $reader->{input} = IO::File->new("$assembly_tags_file");

    while(my $obj=$reader->next_object)
    {       
        my $assembly_tag = $self->_build_assembly_tag($obj);
        push @assembly_tags, $assembly_tag;
    }
    $reader->{input} = undef;
    return \@assembly_tags; 
    
}

sub set_assembly_tags
{
    my ($self, $assembly_tags) = @_;
    
    my $cache_directory = $self->_cache_dir;
    my $assembly_tags_file = $cache_directory."/assembly_tags";
    my $writer = Genome::Model::Tools::Pcap::Ace::Writer->new(IO::File->new(">$assembly_tags_file"));
    foreach my $obj (@{$assembly_tags})
    {
        $writer->write_obj($obj);
    }
}

#TODO: need to create different version of write_file
#1.  writes cache to ace file(s)
#2.  writes contig to ace file
#3.  writes

sub write_directory
{
    my ($self, %params) = @_;
    my $output_directory = $params{output_directory};
    Carp::confess("output_directory not defined\n") and return unless defined $output_directory;
    `mkdir -p $output_directory` unless -d $output_directory;
    my $prefix = $params{prefix}||'';
    my $number = $params{number}||1;
    my $contig_names =$self->get_contig_names;
    if($number > 1)
    {
        for(my $i=0;$i<$number;$i++)
        {
        #TODO: come up with better names for the confusing index/number params below
            $self->write_file(output_file => "$output_directory/$prefix$i.ace", index => $i, number => $number);
        }
    }
    else
    {
        $self->write_file(output_file => "$output_directory/$prefix".'0.ace');
    }    
}
sub write_file
{
    my ($self, %params) = @_;
    #reopen output file
    
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
    my $index = $params{index};
    my $number = $params{number};
    $self->_writer (Genome::Model::Tools::Pcap::Ace::Writer->new($self->_output));
    #first, come up with a list of read and contig counts
    my $read_count=0;
    my $contig_count=0;
    my @contig_names;
    
    my @contig_index_files;
    my $cache_dir = $self->_cache_dir;
    @contig_index_files = `ls $cache_dir/*.ci`;
    chomp @contig_index_files;
    @contig_index_files = map { /(.*)\.ci/ } @contig_index_files;#chop off ci extension for sorting       
    @contig_index_files = sort { _cmptemp($a, $b) } @contig_index_files ;#sort by scaffold_num.contig_num
    @contig_index_files = map { $_.'.ci' } @contig_index_files ;#re-add extension
    
    if(defined $index && defined $number)
    {
        my @temp_ci;
        foreach my $ci_file (@contig_index_files)
        {
            my ($contig_num) = $ci_file =~ /.*Contig(\d+).*/;
            #print "writing $contig_num\n";
            next unless (($contig_num%$number) == $index);
            push @temp_ci,$ci_file;
        }
        @contig_index_files = @temp_ci;
    }
    $contig_count = @contig_index_files;
#TODO: continue to work on ace file writing code
    foreach my $ci_file (@contig_index_files)
    {
        my $contig = retrieve $ci_file;
        
        if($contig->{offset} != -1)
        {
            $read_count += scalar keys %{$contig->{reads}}; 
        }
        else
        {
            $contig->{contig_object}->thaw;
            $read_count += $contig->{contig_object}->read_count;
        }
    }
      
    my $ace_assembly = { type => 'assembly', 
                          contig_count => $contig_count,
                          read_count => $read_count };
    $self->_writer->write_object($ace_assembly);
    #initialize contig tags array
    $self->{contig_tags} = [];
    #write out contigs
    #my $ci_files =     
    foreach my $ci_file (@contig_index_files)
    {
        my $contig_index = retrieve $ci_file;
        if($contig_index->{offset} == -1)
        {
            $self->_write_contig_from_object($contig_index->{contig_object});           
        }
        else
        {
            $self->_write_contig_from_file($contig_index);
        }   
    }
    
    #write out contig tags
    if(defined $self->{contig_tags})
    {
        foreach my $tag (@{$self->{contig_tags}})
        {
            $self->_writer->write_object($tag);     
        }
        $self->{contig_tags} = undef;   
    }
    

    
    #write out assembly tags
    #TODO: copy assembly tags from cache to output file
    #$self->_output->close; 
    my $assembly_tags = $self->get_assembly_tags;
    foreach my $assembly_tag (@{$assembly_tags})
    {
        $self->_writer->write_object($assembly_tag);
    }
    $self->_output->autoflush(1);
}

sub _write_contig_from_object
{
	my ($self, $contig) = @_;
	my $writer = $self->_writer;
    my $output = $self->_output;
	#my %reads = %{$contig->reads}; #contins read_positions and read_tags
	my @contig_tags = @{$contig->{tags}||[]} ;
    my $contig_index = $contig->{callbacks}{index};
	
    my $reader = $self->_reader;
	$contig->thaw();
    my $input = $contig->{callbacks}{fh};
    $reader->{input} = $input;
	#first write contig	hash
    if($contig->loaded||$contig->check_data_changed("contig")||$contig->check_data_changed("padded_base_string"))
    {
        my @tokens = ("CO",$contig->name,$contig->base_count,$contig->read_count,$contig->base_segment_count,
                      $contig->complemented?"C":"U");
        my $line = join " ",@tokens;
        $line .= "\n";
        print $output $line;
    }
    else
    {
        $input->seek($contig_index->{offset},0);
        my $string = <$input>;
        print $output $string;    
    }
    #next, write contig sequence
    if($contig->loaded||$contig->check_data_changed("padded_base_string"))
    {
        my $consensus = $contig->padded_base_string;
		$writer->_write_sequence($output, $consensus);
    }
    else
    {
        my $consensus;
        $input->seek($contig_index->{offset},0); 
        my $line = <$input>;
        my $length = $contig_index->{base_sequence}{length} - length($line);
        $input->read($consensus, $length);
        print $output $consensus;
        
    }
    #write contig base qualities
	if($contig->loaded||$contig->check_data_changed("padded_base_quality"))
    {
        print $output "\n\nBQ";

        my $width = 50;#$self->width();
        my @bq = @{$contig->unpadded_base_quality};
        for (my $i = 0; $i < @bq; $i += 1) {
            if ($i % $width == 0) {
                print $output "\n";
            }
            print $output " $bq[$i]";
        }
        print $output "\n\n";
    }
    else
    {
		print $output "\n\n";
        my $string;
        $input->seek($contig_index->{base_qualities}{offset},0);
        $input->read($string, $contig_index->{base_qualities}{length});
        print $output $string;      
    }
    #write out the contig's read positions
	my @reads;
    if($contig->loaded||$contig->check_data_loaded("children"))
    {
        #if it's been loaded, we need to check each read to see if it's changed
        #yeah, I know, this is kind of slow, I'll add a better callback mechanism in 
        #the future
        my $children = $contig->children;
		#@reads = sort { $a->position <=> $b->position } values %{$children};
        foreach my $read (values %{$children})
        {
            if($contig->loaded||$read->check_data_changed("self"))
            {
                if($contig->loaded||$read->check_data_changed("read_position"))
                {
                    print $output "AF ".$read->name;
                    if($read->complemented)
                    {
                        print $output " C ".$read->position."\n";
                    }
                    else
                    {
                        print $output " U ".$read->position."\n";
                    }                
                }
                else
                {
                   $input->seek($contig_index->{reads}{$read->name}{read_position}{offset},0);
                   my $line = <$input>;
                   print $output $line;                
                }            
            }
            else
            {
                $input->seek($contig_index->{reads}{$read->name}{read_position}{offset},0);
                my $line = <$input>;
                print $output $line;            
            }        
        }    
    }
    else
    {
        my $string;
        $input->seek($contig_index->{af_start},0);
        $input->read($string,$contig_index->{af_end}-$contig_index->{af_start});
        print $output $string;
    }
    #write out the contigs base segments
    if($contig->loaded||$contig->check_data_changed("base_segments"))
    {	
		my @base_segments = @{$contig->base_segments};
        foreach my $base_segment (@base_segments)
        {
            $writer->write_object($base_segment);
        } 
		print $output "\n";       
    }
    else
    {
        my $string;
        $input->seek($contig_index->{base_segments}{offset},0);
        $input->read($string, $contig_index->{base_segments}{length});
        print $output $string;    
    }
    #write out the contig's read lines
    if($contig->loaded||$contig->check_data_loaded("children"))
    {
        #if it's been loaded, we need to check each read to see if it's changed
        #yeah, I know, this is kind of slow, I'll add a better callback mechanism in 
        #the future
        my $children = $contig->children;
        foreach my $read (values %{$children})
        {
            if($contig->loaded||$read->check_data_changed("self"))
            {
                if($contig->loaded||$read->check_data_changed("read"))
                {
                    my @tokens = ("RD",$read->name,$read->length,$read->info_count,scalar @{$read->tags});
                    my $string = join " ",@tokens,"\n";
                    print $output $string;
                
                }
                else
                {
                    my $string;
                    $input->seek($contig_index->{reads}{$read->name}{read}{offset},0);
                    $string = <$input>;;
                    print $output $string;                
                }
                
                if($contig->loaded||$read->check_data_changed("padded_base_string"))
                {
                    my $sequence = $read->padded_base_string;
                    my $seq_len = length($sequence);
                    my $width = 50;
                    for (my $i = 0;$i < $seq_len; $i += $width ) {
                        print $output substr($sequence, $i, $i + ($width-1) < $seq_len ? $width : $seq_len - $i) . "\n";
                    }
                    print $output "\n";                
                }
                else
                {
                    my $string;
                    $input->seek($contig_index->{reads}{$read->name}{read}{sequence}{offset},0);
                    $input->read($string, $contig_index->{reads}{$read->name}{read}{sequence}{length});
                    print $output $string; 
					print $output "\n";               
                }
                if($contig->loaded||$read->check_data_changed("qa"))
                {
                    print $output ("QA ", $read->qual_clip_start, " ",$read->qual_clip_end, " ", $read->align_clip_start, " ", $read->align_clip_end,"\n");                
                }
                else
                {
                    $input->seek($contig_index->{reads}{$read->name}{read}{qa}{offset},0);
                    my $string = <$input>;
                    print $output $string;                
                }
                if($contig->loaded||$read->check_data_changed("ds"))
                {
                    print $output ("DS CHROMAT_FILE: ",$read->chromat_file," PHD_FILE: ",$read->phd_file,
                    " CHEM: ", $read->chemistry, " TIME: ", $read->time); 
					print $output " DYE: ", $read->dye if ($read->dye);
					print $output "\n\n";               
                }
                else
                {
                    $input->seek($contig_index->{reads}{$read->name}{read}{ds}{offset},0);
                    my $string = <$input>;
                    print $output $string;
					print $output "\n";                
                }            
            }
			else #check if data has changed
			{
				my $string;
				$input->seek($contig_index->{reads}{$read->name}{read}{offset},0);
				$input->read($string, $contig_index->{reads}{$read->name}{read}{length});
				print $output $string;

			}        
        }    
    }
    else
    {
        my $string;
        $input->seek($contig_index->{rd_start},0);
        $input->read($string,$contig_index->{rd_end}-$contig_index->{rd_start});
        print $output $string;    
    }
    if($contig->loaded||$contig->check_data_changed("tags"))
    {
        #store contig tags for writing, also convert
        #tag to low level format for writing
        my @contig_tags = @{$contig->tags};
        foreach my $tag (@contig_tags)
        {
            #$self->_write_contig_tag($contig_tag);
            my $contig_tag =  {
                type => 'contig_tag',
                tag_type => $tag->type,
                date => $tag->date,
                program => $tag->source,
                contig_name => $tag->parent,
                scope => 'ACE',
                start_pos => $tag->start,
                end_pos => $tag->stop,
                data => $tag->text,
                no_trans => $tag->no_trans,
            };
            push @{$self->{contig_tags}}, $contig_tag;
        }    
    }
    else
    {
        my @contig_tags_index = @{$contig_index->{contig_tags}};         
          
        #write out contig tags
        foreach my $contig_tag_index (@contig_tags_index)
        {
            $input->seek($contig_tag_index->{offset},0);
            push @{$self->{contig_tags}},
            map
            {
                {
                    type => 'contig_tag',
                    tag_type => $_->type,
                    date => $_->date,
                    program => $_->source,
                    contig_name => $_->parent,
                    scope => 'ACE',
                    start_pos => $_->start,
                    end_pos => $_->stop,
                    data => $_->text,
                    no_trans => $_->no_trans,
                }
            }       
            Genome::Model::Tools::Pcap::TagParser->new()->parse($input);
        }    
    }
}

sub _write_contig_from_file
{
	my ($self, $contig_index) = @_;
	
    my $output = $self->_output;
	my $writer = $self->_writer;
	my @contig_tags_index = @{$contig_index->{contig_tags}};
	my $input_file = $contig_index->{file_name};
    my $input = Genome::Model::Tools::Pcap::FHManager->get_fh($input_file);

    #write out contig
    my $contig_string;
    $input->seek($contig_index->{offset},0);
    $input->read($contig_string, $contig_index->{contig_length});
    print $output $contig_string;	
	#write out contig tags
	foreach my $contig_tag_index (@contig_tags_index)
	{
	    $input->seek($contig_tag_index->{offset},0);
	    push @{$self->{contig_tags}},
	    map
	    {
		    {
		        type => 'contig_tag',
		        tag_type => $_->type,
		        date => $_->date,
		        program => $_->source,
		        contig_name => $_->parent,
		        scope => 'ACE',
		        start_pos => $_->start,
		        end_pos => $_->stop,
		        data => $_->text,
		        no_trans => $_->no_trans,
		    }
        }		
	    Genome::Model::Tools::Pcap::TagParser->new()->parse($input);
	}
}

sub _build_assembly_tag {
    my ($self, $obj) = @_;

    my $tag = new Genome::Model::Tools::Pcap::Tag(
        type => $obj->{tag_type},
        date => $obj->{date},
        source => $obj->{program},
        text => $obj->{data},
    );
    return $tag;
}

sub _build_read_tag {
    my ($self, $obj) = @_;
    my $tag = new Genome::Model::Tools::Pcap::Tag(
        type => $obj->{tag_type},
        date => $obj->{date},
        source => $obj->{program},
        parent => $obj->{read_name},
        scope => 'ACE',
        start => $obj->{start_pos},
        stop => $obj->{end_pos},
    );
    return $tag;
}

sub _copy_assembly_tag
{
	my ($self, $assembly_tag) = @_;
	return Storable::dclone($assembly_tag);
}

sub _write_assembly_tag
{
    my ($self, $tag) = @_;

	
	my $ace_tag = { type => 'assembly_tag',
					tag_type => $tag->type,
					program => $tag->source,
					date => $tag->date,
					data => $tag->text,
				  };
    
	$self->_writer->write_object($ace_tag);    
        
    return;
}

sub _write_read_tag
{
    my ($self, $tag) = @_;
		
	my $read_tag =  {
		type => 'read_tag',
        tag_type => $tag->type,
        date => $tag->date,
        program => $tag->source,
        read_name => $tag->parent,
        scope => 'ACE',
        start_pos => $tag->start,
        end_pos => $tag->stop,
    };    
    
    $self->_writer->write_object($read_tag);
    return;
}

sub _write_contig_tag
{
    my ($self, $tag) = @_;

    my $contig_tag =  {
		type => 'contig_tag',
        tag_type => $tag->type,
        date => $tag->date,
        program => $tag->source,
        contig_name => $tag->parent,
        scope => 'ACE',
        start_pos => $tag->start,
        end_pos => $tag->stop,
		data => $tag->text,
		no_trans => $tag->no_trans,
    };
	
	$self->_writer->write_object($contig_tag);   	

    return;
}

sub get_num
{
    my ($name) = @_;
    my ($ctg_num) = $name =~ /Contig(\d+)\.\d+/;
	($ctg_num) = $name =~ /Contig(\d+)/ if(!defined $ctg_num);
	($ctg_num) = $name =~ /.?(\d+)/ if(!defined $ctg_num);
	
    return $ctg_num;
}

sub get_ext
{
    my ($name) = @_;
    my ($ctg_ext) = $name =~ /Contig\d+\.(\d+)/;
    return $ctg_ext;
}

sub _cmptemp
{
    my ($a, $b) = @_;
    my $num1 = get_num($a);
    my $num2 = get_num($b);
    if($num1 > $num2)
    {
        return 1;
    }elsif($num2 > $num1)
    {
        return -1;
    }
    my $ext1 = get_ext($a);$ext1 = 0 if(!defined $ext1);
    my $ext2 = get_ext($b);$ext2 = 0 if(!defined $ext2);
    if($ext1 > $ext2)
    {
        return 1;
    }elsif($ext2 > $ext1)
    {
        return -1;
    }
    return 0;
    
}





1;

=pod

=head1 NAME

Ace - Object oriented ace file reader/writer

=head1 SYNOPSIS

my $ace_object = Genome::Model::Tools::Pcap::Ace->new(input_file => "inputfilename", output_file => "outputfilename", using_db => 1, input_file_index => "inputfileindex");

 my @contig_names = $ace_object->get_contig_names();
 my $contig = $ace_object->get_contig("Contig0.1");
 $ace_object->remove_contig("Contig0.1");

 $ace_object->write_file;
    
=head1 DESCRIPTION

Genome::Model::Tools::Pcap::Ace indexes an ace file, and allows the user to get Contig objects from the ace file, edit them, and write the file back to the hard disk when finished.

=head1 METHODS

=head1 new 

my $ace_object = new Genome::Model::Tools::Pcap::Ace(input_file => $input_file, output_file => $output_file);

input_file - required, the name of the input ace file.

output_file - option, the name of output ace file.  You can give ace_object the file handle when you create it, or later when you write it.  If you are reading, then you don't need to specify the file handle.

=head1  get_contig_names 

my @contig_names = $ace_object->get_contig_names();

returns a list of contig names in the ace file.

=head1 get_contig 

my $contig = $ace_object->get_contig("Contig0.1");
    
returns a Contig object to the user.

=head1 add_contig 

 my $contig = $ace_object->get_contig("Contig0.1");
 ...
 $ace_object->add_contig($contig);
    
inserts a contig into the ace file.  If a contig with that name already exists, then it is overwritten by the data in the newly added contig.

=head1 remove_contig 

$ace_object->remove_contig("Contig0.1");
    
returns a Contig from the ace file.

=head1 get_assembly_tags 

my @assembly_tags = $ace_object->get_assembly_tags;
    
returns an array off assembly tags to the user.

=head1 set_assembly_tags 

$ace_object->set_assembly_tags(\@assembly_tags);
    
replaces the current array of assembly tags in the ace file with a new list of assembly tags.

=head1 write_file

$ace_object->write_file;
    
This function will write the ace object in it's current state to the output ace file specified during object construction.

=head1 Author(s)

 Jon Schindler <jschindl@watson.wustl.edu>

=cut

#$HeadURL: svn+ssh://svn/srv/svn/gscpan/perl_modules/trunk/Genome/Assembly/Pcap/Ace.pm $
#$Id: Ace.pm 58706 2010-05-13 17:55:53Z jschindl $
