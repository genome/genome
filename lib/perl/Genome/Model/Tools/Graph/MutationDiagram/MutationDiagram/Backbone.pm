#----------------------------------
# $Author: dlarson $ 
# $Date: 2008-09-16 16:33:54 -0500 (Tue, 16 Sep 2008) $
# $Revision: 38655 $
# $URL: svn+ssh://svn/srv/svn/gscpan/perl_modules/trunk/MG/MutationDiagram/Backbone.pm $
#----------------------------------
package Genome::Model::Tools::Graph::MutationDiagram::MutationDiagram::Backbone; 
#------------------------------------------------
our $VERSION = '1.0';
#------------------------------------------------
use strict;
use warnings;
use SVG;
use Genome::Model::Tools::Graph::MutationDiagram::MutationDiagram::View;

# If you know about and want to use inheritance:
#use Base::Class;
our @ISA = qw( Genome::Model::Tools::Graph::MutationDiagram::MutationDiagram::View);

sub new {
 	my ($class, %arg) = @_;
    #set up viewport
    #$arg{viewport} = { x => 0,
    #                 y => 0,
    #                 height => $arg{height}, #default to available space
    #                 width => $arg{protein_length}, #expected to be in amino acids
    #             };
                 
    my $self = Genome::Model::Tools::Graph::MutationDiagram::MutationDiagram::View->new(%arg);
    $self->{_backbone}{height} = $arg{backbone_height};
	$self->{_gene} = $arg{gene};
    $self->{_protein_length} = $arg{protein_length}; 
    #we want to center the backbone within the available space
    $self->{_backbone}{x} = 0;
    $self->{_backbone}{y} = ($self->{_size}{height} - $self->{_backbone}{height}) / 2;

    $self->{_aa_scaling} = ($self->{_size}{width} - $self->{_backbone}{x}) / $self->{_protein_length};
                         
    bless($self, ref($class) || $class);
    $self->{_backbone}{length} = $self->x_coordinate_for_amino_acid_pos($self->{_protein_length});
    return $self;
}
#----------------------------------
sub dimensions {
    my ($self,) = @_;
    return $self->{_backbone};
}
#----------------------------------
sub x_coordinate_for_amino_acid_pos {
    my ($self, $amino_acid) = @_;
    return $self->{_backbone}{x} + $amino_acid * $self->{_aa_scaling};
}
#-------------------------------------------------
sub draw {
    my ($self,) = @_;
    my $svg = $self->{_svg};
	$svg->rectangle(
														 id => $self->{id},
														 x => $self->{_backbone}{x} + $self->{_backbone}{x} , 
                                                         y => $self->{_backbone}{y},
                                                         width => $self->{_backbone}{length},
                                                         height => $self->{_backbone}{height},
														 style => $self->{_style},
                                                         id => "protein_backbone_outline",
                                                         
                       
                              );

		
		my $coord_increment = int($self->{_protein_length}/5/100) * 100;
		$coord_increment ||= 50;
		my $coord_i;
		for ($coord_i = 0;
				 $coord_i < $self->{_protein_length};
				 $coord_i += $coord_increment) {
			$svg->line(
								 x1 => $self->x_coordinate_for_amino_acid_pos($coord_i),
								 y1 => $self->{_backbone}{y} + $self->{_backbone}{height},
								 x2 => $self->x_coordinate_for_amino_acid_pos($coord_i),
								 y2 => $self->{_backbone}{y} + $self->{_backbone}{height} + 5,
								 style => {
													 fill => 'none',
													 stroke => 'black'
													}
								);
			
#			if ($coord_i % (2 * $coord_increment)) {
				my $coord_label = $coord_i;
				$svg->text(
									 id => 'coordinate_' . $coord_label,
									 x => $self->x_coordinate_for_amino_acid_pos($coord_i),
									 y => $self->{_backbone}{y} + $self->{_backbone}{height} + 20,
									 style => {
														 'text-align' => 'center',
														 'text-anchor' => 'middle'
														}
									)->cdata($coord_label);
#			}
		}
		if ($coord_i - int($coord_increment/2) < $self->{_protein_length}) {
			$coord_i = $self->{_protein_length};
			$svg->line(
								 x1 => $self->x_coordinate_for_amino_acid_pos($coord_i),
								 y1 => $self->{_backbone}{y} + $self->{_backbone}{height},
								 x2 => $self->x_coordinate_for_amino_acid_pos($coord_i),
								 y2 => $self->{_backbone}{y} + $self->{_backbone}{height} + 5,
								 style => {
													 fill => 'none',
													 stroke => 'black'
													}
								);
			my $coord_label = $coord_i;
			$svg->text(
								 id => 'coordinate_' . $coord_label,
								 x => $self->x_coordinate_for_amino_acid_pos($coord_i),
								 y => $self->{_backbone}{y} + $self->{_backbone}{height} + 20,
								 style => {
													 'text-align' => 'center',
													 'text-anchor' => 'middle'
													}
								)->cdata($coord_label);
		}

		my $coord_scale =
			$svg->text(
								 id => 'coordinate_scale',
								 x => $self->x_coordinate_for_amino_acid_pos($self->{_protein_length})/2,
								 y => $self->{_backbone}{y} + $self->{_backbone}{height} + 40,
								 style => {
													 'text-align' => 'center',
													 'text-anchor' => 'middle',
													}
								)->cdata('Scale (AA)');

		my $gene_name =
			$svg->text(
								 id => 'gene_name',
								 x => $self->x_coordinate_for_amino_acid_pos($self->{_protein_length})/2,
								 y => $self->{_backbone}{y} + $self->{_backbone}{height} + 70,
								 style => {
													 'text-align' => 'center',
													 'text-anchor' => 'middle',
													 'font-size' => 24,
													 'font-style' => 'italic',
													 'font' => 'Myriad-Roman',
													}
							 
							)->cdata($self->{_gene});
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

