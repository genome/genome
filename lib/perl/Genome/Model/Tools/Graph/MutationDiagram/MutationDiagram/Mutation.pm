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
						 };
    $self->{_feature_length} = 1;

    #set up default layout attributes
    my $backbone = $self->{_backbone};
    #add to subgroup if passed in

    my $mutation_top = 70;
    if($self->{_frequency} > 3) {
        $mutation_top += ($self->{_frequency} - 3) * 13; 
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
		for (my $i = 0; $i < $self->{_frequency};$i++) {
			#draw lollipop                                            
			my $circle = $svg->circle( cx => $self->{_lollipop}{x}, 
																 cy => $self->{_lollipop}{y} - ($i * 13), 
																 r => 5,
																 id=>$self->{_id}."_lollipop_" . $i,
																 style => {fill => $self->{_color}, stroke => 'black'});  
		}
    #add text
    my $transform = sprintf("matrix(0 -1 1 0 %f %f)",$self->{_text_label}{x},$self->{_text_label}{y});
    my $text = $svg->text(id => $self->{_id}."_label",
                          transform => $transform,

                      )->cdata($self->{_text});
                                                
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

