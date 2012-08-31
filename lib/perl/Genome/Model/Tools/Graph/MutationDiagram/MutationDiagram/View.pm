#----------------------------------
# $Author: bshore $ 
# $Date: 2008-07-01 09:11:56 -0500 (Tue, 01 Jul 2008) $
# $Revision: 36100 $
# $URL: svn+ssh://svn/srv/svn/gscpan/perl_modules/trunk/MG/MutationDiagram/View.pm $
#----------------------------------
package Genome::Model::Tools::Graph::MutationDiagram::MutationDiagram::View;
#------------------------------------------------
our $VERSION = '1.0';
#------------------------------------------------
use strict;
use warnings;
use SVG;

# If you know about and want to use inheritance:
#use Base::Class;
#our @ISA = qw( Base::Class );
our $subview_number = 0;

sub new {
 	my ($class, %arg) = @_;
	my $self = {};
    unless(exists($arg{width}) && exists($arg{height})) {
        return undef;
    }
    else {
        my %create_parameters;
        $self->{_id} = $create_parameters{id} = $arg{id};
        $self->{_size} =  { width => $arg{width},
                          height => $arg{height},
                      };
        $create_parameters{width} = $arg{width};
        $create_parameters{height} = $arg{height};
        
        $self->{_parent} = $arg{parent};       
        $self->{_style} = $arg{style};
        if(exists($arg{right_margin})) {
            $self->{_margins}{right_margin} = $arg{right_margin};
        }
        if(exists($arg{left_margin})) {
            $self->{_margins}{left_margin} = $arg{left_margin};
        }
        if(exists($arg{top_margin})) {
            $self->{_margins}{top_margin} = $arg{top_margin};
        }
        if(exists($arg{bottom_margin})) {
            $self->{_margins}{bottom_margin} = $arg{bottom_margin};
        }
				$self->{_margins}{right_margin} ||= 0;
				$self->{_margins}{left_margin} ||= 0;
				$self->{_margins}{top_margin} ||= 0;
				$self->{_margins}{bottom_margin} ||= 0;
        if(exists($arg{x}) && exists($arg{y})) { 
            $self->{_origin}{x} = $arg{x};
            $self->{_origin}{y} = $arg{y};
            $create_parameters{x} = $arg{x};
            $create_parameters{y} = $arg{y};
        }
				$self->{_origin}{x} ||= 0;
				$self->{_origin}{y} ||= 0;
        if(exists($arg{viewport})) {
            $self->{_viewport} = $arg{viewport};
            $create_parameters{viewBox} = 
                sprintf "%d %d %d %d", @{$arg{viewport}}{('x','y','width','height')}; 
        }
        if(defined($self->{_parent})) {
            $subview_number++;
            $self->{_svg} = $self->{_parent}->svg;
            $self->{_svg} = $self->{_svg}->svg(%create_parameters); 
            $self->{_svg}->group(id => "view$subview_number");
        }
        else {
            #root svg
            $self->{_svg} = SVG->new(%create_parameters);
        }
    }
    bless($self, ref($class) || $class);
    return $self;
}
#-------------------------------------------------
sub size {
  my ($self,) = @_;
  return $self->{_size};
}
#----------------------------------
sub margins {
    my ($self,) = @_;
    return $self->{_margins};
}
#----------------------------------
sub svg {
    my ($self,) = @_;
    return $self->{_svg};
}
#----------------------------------
sub content_view {
    my ($self,) = @_;
    my %dimensions;
    #expectation is that margins are in viewport coordinates
    $dimensions{x} = $self->{_origin}{x} + $self->{_margins}{left_margin};
    $dimensions{y} = $self->{_origin}{y} + $self->{_margins}{top_margin};
    $dimensions{width} = $self->{_size}{width} - $self->{_margins}{right_margin} - $self->{_margins}{left_margin};
    $dimensions{height} = $self->{_size}{height} - $self->{_margins}{top_margin} - $self->{_margins}{bottom_margin};
    return %dimensions;
}
#----------------------------------
sub draw {
    my ($self,) = @_;
    my $svg = $self->{_svg};
	$svg->rectangle(
														 id => $self->{id},
														 x => $self->{_viewport}{x} , 
                                                         y => $self->{_viewport}{y},
                                                         width => $self->{_viewport}{width},
                                                         height => $self->{_viewport}{height},
														 style => $self->{_style} 
                                                     );
    
}
=head1 AUTHOR

David Larson, E<lt>dlarson@watson.wustl.eduE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2008 by David Larson

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8 or,
at your option, any later version of Perl 5 you may have available.


=cut
1;

