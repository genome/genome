package Genome::Model::Tools::Pcap::Proxy;

#use UNIVERSAL qw(can);
*can = \&UNIVERSAL::can;

sub new
{
	croak("new:no class given, quitting") if @_ < 1;
    my ($caller, %args) = @_;
    my $caller_is_obj = ref($caller);
    my $class = $caller_is_obj || $caller;
    my $self = {};
    bless ($self, $class);		
	     
    return $self;
}

sub name
{
    my ($self, $value) = @_;
    if(@_ > 1)
    {
        $self->{callbacks}{name} = $value;
        $self->been_changed("name",1);
    }
    return $self->{callbacks}{name};
}

sub callbacks
{
    my ($self, $callbacks) = @_;
    if(@_ > 1)
    {
        $self->{callbacks} = $callbacks;
    }
    return $self->{callbacks};
}

sub check_and_load_data
{
	my ($self,$name,$value,$just_load) = @_;
	
	if (@_ > 2)
	{
		$self->{$name} = $value;
		$self->already_loaded($name,1);
		$self->been_changed($name,1)if(!$self->{just_load});
        #$self->been_changed($name,0)if($self->{just_load});
	}
	
	if(!$self->already_loaded($name))
	{
		$self->load_data($name);
	}
	return $self->{$name};

}

sub load_data
{
	my ($self,$name) = @_;
    my @callbacks = ($self->{callbacks});
	return if(!defined $self->{callbacks});
	foreach(@callbacks)
	{
		if(!defined $_)
		{
			my $num = 0;
		}
		if(can($_,$name) && $_->$name)
		{
			$_->$name($self);
            last;			
		}
	}
	return;
	
}

sub already_loaded
{
	my ($self,$name,$been_loaded) = @_;
	if(@_ > 2)
	{
		$self->{data}{$name}[0] = $been_loaded;
        $self->{data}{self}[0] = 1;
	}

	return $self->{data}{$name}[0];

}

#probably remove this one
sub check_map {
	my ($self,$name) = @_;
	my @callbacks = ($self->{callbacks});
	my $map = $callbacks[0]->get_map();
	my @props_list = @{$map->{$name}};

	foreach my $prop_name (@props_list)
	{
		return 1 if($self->been_changed($prop_name));
	}
	return 0;
}

sub check_data_loaded{
    my ($self,$name) = @_;
    my @callbacks = ($self->{callbacks});
    my $map = $callbacks[0]->get_map();
    my @props_list = @{$map->{$name}};

    foreach my $prop_name (@props_list)
    {
        return 1 if($self->already_loaded($prop_name));
    }
    return 0;
}

sub check_data_changed{
    my ($self,$name) = @_;
    my @callbacks = ($self->{callbacks});
	if(!defined $callbacks[0])
	{
		return 1;
	}
    my $map = $callbacks[0]->get_map();
    my @props_list = @{$map->{$name}};

    foreach my $prop_name (@props_list)
    {
        return 1 if($self->been_changed($prop_name));
    }
    return 0;
}

sub been_changed {
    my ($self,$name,$changed) = @_;
    if(@_ > 2)
    {
        $self->{data}{$name}[1] = $changed;
        $self->{data}{self}[1] = 1;        
    }
    return $self->{data}{$name}[1];
}

sub freeze {
	my ($self) = @_;
	my @callbacks = ($self->{callbacks});
	foreach(@callbacks)
	{
		can($_,'freeze')&&$_->freeze;
	}
}

sub thaw {
	my ($self) = @_;	
	my @callbacks = ($self->{callbacks});
	foreach(@callbacks)
	{
		can($_,'thaw')&&$_->thaw;
	}
}

1;
