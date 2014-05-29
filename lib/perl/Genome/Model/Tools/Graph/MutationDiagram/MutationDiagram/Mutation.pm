#----------------------------------
# $Author: dlarson $ 
# $Date: 2008-09-16 16:33:54 -0500 (Tue, 16 Sep 2008) $
# $Revision: 38655 $
# $URL: svn+ssh://svn/srv/svn/gscpan/perl_modules/trunk/MG/MutationDiagram/Mutation.pm $
#----------------------------------
package Genome::Model::Tools::Graph::MutationDiagram::MutationDiagram::Mutation; 
#------------------------------------------------
our $VERSION = '1.0';
#------------------------------------------------
use strict;
use warnings;
use SVG;
use Carp;

# If you know about and want to use inheritance:
#use Base::Class;
#our @ISA = qw( Base::Class );

sub new {
 	my ($class, %arg) = @_;
	my $self = {
							_backbone => $arg{backbone},
							_frequency => $arg{frequency},
							_color => $arg{color},
							_style => $arg{style},
							_amino_acid => $arg{start_aa},
							_text => $arg{text},
							_id => $arg{id},
                            _max_frequency => $arg{max_freq},
                            _shape => $arg{shape},
                            _only_label_max => $arg{only_label_above_max_freq},
						 };
    $self->{_feature_length} = 1;

    #set up default layout attributes
    my $backbone = $self->{_backbone};
    #add to subgroup if passed in

    my $mutation_top = 43;
    if($self->{_frequency} > 1) {
        $mutation_top += ($self->{_frequency} - 1) * 13; 
    }

    if($self->{_max_frequency}) {
        if($self->{_frequency} > $self->{_max_frequency}) {
            my $gutter_till_label = 14;
            my $first_mutation_location = 30;
            $mutation_top = $first_mutation_location + ($self->{_max_frequency} + 1) * 13 + $gutter_till_label;

            #add in the allele count to the label
            $self->{_text} = $self->{_text} . " (" . $self->{_frequency} . ")";
        }
        elsif($self->{_only_label_max}) {
            $self->{_text} = q{};
        }
    }
    
    my @line_path_y =
    ($backbone->dimensions->{y} - $mutation_top,$backbone->dimensions->{y}-20,$backbone->dimensions->{y}-10,$backbone->dimensions()->{y},$backbone->dimensions()->{y}+$backbone->dimensions()->{height});
    my @line_path_x = ($backbone->x_coordinate_for_amino_acid_pos($self->{_amino_acid})) x scalar(@line_path_y);

    #store attributes
    $self->{_line_path_y} = \@line_path_y;
    $self->{_line_path_x} = \@line_path_x;
    
    #store lollipop center
    $self->{_lollipop}{x} = $backbone->x_coordinate_for_amino_acid_pos($self->{_amino_acid});
    $self->{_lollipop}{y} = $backbone->dimensions->{y}-30;

    #store label offsets
    my $text_height = 7;
    $self->{_text_label}{x} = $backbone->x_coordinate_for_amino_acid_pos($self->{_amino_acid}) + $text_height/2;
    $self->{_text_label}{y} = $backbone->dimensions->{y}-$mutation_top;

    bless($self, ref($class) || $class);

    if($arg{shape}) {
        $self->{_shape} = "_" . lc($self->{_shape});
        unless($self->can($self->{_shape})) {
            croak "Shape $arg{shape} is currently unsupported";
        }
    }
    
    return $self;
}
#----------------------------------
sub get_label_position {
    my ($self) = @_;
    #assume that everything is being programmatically accounted for
    #return the x coordinate of the lollipop
    return $self->{_lollipop}{x};
}
#----------------------------------
sub set_label_position {
    my ($self,$new_position) = @_;
    @{$self->{_line_path_x}}[(0,1)] = $new_position;
    $self->{_lollipop}{x} = $new_position;
    $self->{_text_label}{x} = $new_position;
}
#----------------------------------
sub translate_label_position {
    my ($self,$translation) = @_;
    @{$self->{_line_path_x}}[0] += $translation;
    @{$self->{_line_path_x}}[1] += $translation; 
    $self->{_lollipop}{x} += $translation;
    $self->{_text_label}{x} += $translation;
}
#----------------------------------
sub length {
    my ($self,) = @_;
    return $self->{_feature_length};
}
#----------------------------------
sub amino_acid {
    my ($self) = @_;
    return $self->{_amino_acid};
}
#----------------------------------
sub vertically_align_to {
    my ($self, $mutation) = @_;
    my @new_line_path_y = @{$mutation->{_line_path_y}};
    $self->{_line_path_y} = \@new_line_path_y;
    $self->{_text_label}{y} = $mutation->{_text_label}{y};
    return 1;
}
#----------------------------------
sub draw {
    my ($self,$group) =  @_;
    #add to subgroup if passed in
    my $svg = defined($group) ? $group : $self->{_backbone}->svg;
    $svg = $svg->group(id => $self->{_id} . "_group");
    #This appears to evaluate the closed parameter using an if. Make sure to
    #pass along a value that perl evaluate as false
    #My solution is to not pass it at all :-)
    my $path = $svg->get_path(
        x=>$self->{_line_path_x}, y=>$self->{_line_path_y},
        -type=>'path',
    );
    $svg->path(
        %$path,   	
        id => $self->{_id}, 
        style => $self->{_style}, 
    );

    my $shape_func = $self->{_shape};
    my $drawable_frequency = defined $self->{_max_frequency} && $self->{_frequency} > $self->{_max_frequency} ? $self->{_max_frequency} : $self->{_frequency};                                        
    for (my $i = 0; $i < $drawable_frequency; $i++) {
        #draw lollipop                                            
        my $shape = $self->$shape_func($svg, $self->{_id} . "_lollipop_" . $i, $self->{_lollipop}{x}, $self->{_lollipop}{y} - ($i * 13), 5, {fill => $self->{_color}, stroke => 'black'});  
    }

    if(defined $self->{_max_frequency} && $self->{_frequency} > $self->{_max_frequency}) {
        #frequency is truncated
        $self->_broken_count_indicator($svg, $self->{_lollipop}{x}, $self->{_lollipop}{y} - ($drawable_frequency * 13), 5);
        my $shape = $self->$shape_func($svg, $self->{_id} . "_lollipop_" . "gutter", $self->{_lollipop}{x}, $self->{_lollipop}{y} - (($drawable_frequency + 1) * 13), 5, {fill => $self->{_color}, stroke => 'black'});  
    }

    #add text
    my $transform = sprintf("matrix(0 -1 1 0 %f %f)",$self->{_text_label}{x},$self->{_text_label}{y});
    my $text = $svg->text(id => $self->{_id}."_label",
        transform => $transform,

    )->cdata($self->{_text});
                                                
}

sub _broken_count_indicator {
    my ($self, $svg, $x, $y, $radius) = @_;
    my @x = ($x-$radius,$x+$radius,$x+$radius,$x-$radius);
    my @y = ($y+$radius/2,$y,$y-$radius/2,$y);

    my $group = $svg->group(id => $self->{_id} . "_count_break_group");

    my $path = $group->get_path( x => \@x, y => \@y, -type => 'polygon' );
    $group->polygon( %$path, id => $self->{_id} . "_count_break",style => {fill => 'white', stroke => 'white'});

    $group->line( id => $self->{_id} . "_count_break" . "_bottom_line", x1 => $x[0], x2 => $x[1], y1 => $y[0], y2 => $y[1], style => {stroke => 'black'});
    $group->line( id => $self->{_id} . "_count_break" . "_top_line", x1 => $x[2], x2 => $x[3], y1 => $y[2], y2 => $y[3],style=>{stroke => 'black'});
    return $group;
}

sub _circle {
    my ($self, $svg, $id, $x, $y, $radius, $style) = @_;
    my $circle = $svg->circle( cx => $x, 
        cy => $y, 
        r => $radius,
        id => $id,
        style => $style);
    return $circle;
}

sub _diamond {
    my ($self, $svg, $id, $x, $y, $radius, $style) = @_;
    my @x = ($x-$radius,$x,$x+$radius,$x);
    my @y = ($y,$y+$radius,$y,$y-$radius);
    my $path = $svg->get_path( x => \@x, y => \@y, -type => 'polygon' );
    my $diamond = $svg->polygon( %$path,
        id => $id,
        style => $style);
    return $diamond;
}

sub _square {
    my ($self, $svg, $id, $x, $y, $radius, $style) = @_;
    my $rect = $svg->rect( x => $x - $radius, 
        y => $y - $radius, 
        width => $radius * 2, 
        height => $radius * 2,
        id => $id,
        style => $style);
    return $rect;
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

