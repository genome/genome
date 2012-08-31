#----------------------------------
# $Author: bshore $ 
# $Date: 2008-06-26 17:39:39 -0500 (Thu, 26 Jun 2008) $
# $Revision: 35986 $
# $URL: svn+ssh://svn/srv/svn/gscpan/perl_modules/trunk/MG/MutationDiagram/Legend.pm $
#----------------------------------
package Genome::Model::Tools::Graph::MutationDiagram::MutationDiagram::Legend; 
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
							_values => $arg{values},
							_style => $arg{style},
							_object => $arg{object},
							_x => $arg{x},
							_id => $arg{id},
						 };
    #set up default layout attributes
    my $backbone = $self->{_backbone};

    $self->{_x} = $backbone->x_coordinate_for_amino_acid_pos($self->{_x});
    #add to subgroup if passed in
    my @line_path_y =
    ($backbone->dimensions->{y}-70,$backbone->dimensions->{y}-20,$backbone->dimensions->{y}-10,$backbone->dimensions()->{y},$backbone->dimensions()->{y}+$backbone->dimensions()->{height});
    my @line_path_x = ($self->{_amino_acid}) x scalar(@line_path_y);

    bless($self, ref($class) || $class);
    return $self;
}
#----------------------------------
sub draw {
    my ($self,$group) =  @_;
    #add to subgroup if passed in
    my $svg = defined($group) ? $group : $self->{_backbone}->svg;
    $svg = $svg->group(id => $self->{_id} . "_group");
		my $backbone = $self->{_backbone};
		my $y = $backbone->dimensions->{y} + $backbone->dimensions->{height} + 80;
		foreach my $label (sort (keys %{$self->{_values}})) {
			$y += 15;
			my $color = $self->{_values}{$label};
			if ($self->{_object} eq 'circle') {
				#draw lollipop
				my $circle = $svg->circle( cx => $self->{_x} + 5, 
																	 cy => $y - 4,
																	 r => 5,
																	 id=>$self->{_id}."_" . $label,
																	 style => {fill => $color, stroke => 'black'});  
			} elsif ($self->{_object} eq 'rectangle') {
				$svg->rectangle(
												id => $self->{_id} . '_' . $label, 
												x => $self->{_x} + 5,
												y => $y - 9,
												width => 10,
												height => 10,
												style => {fill => $color, stroke => 'black'}
											 );
			} else {
				$svg->rectangle(
												id => $self->{_id} . '_' . $label, 
												x => $self->{_x} + 5,
												y => $y - 9,
												width => 10,
												height => 10,
												style => {fill => $color, stroke => 'black'}
											 );
			}
			#add text
			my $text = $svg->text(id => $self->{_id}."_label" . $label,
														x => $self->{_x} + 25,
														y => $y,
														style => {
																			'text-align' => 'center',
																		 }
																			
													 )->cdata($label);
		}
                                                
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

