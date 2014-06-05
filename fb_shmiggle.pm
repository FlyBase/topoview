package Bio::Graphics::Glyph::fb_shmiggle;

# [ to see proper formatting set tab==2 ]
#
# 2009-2010 Victor Strelets, FlyBase.org 

use constant DEBUG => 0;

#use strict;
use GD;
use vars '@ISA';
use Bio::Graphics::Glyph::generic;
BEGIN{ @ISA = 'Bio::Graphics::Glyph::generic';
  local($colors_selected,$black,$white,$yellow,$red,$lightgrey,$grey,$darkgrey);
	@colors= ();
	@colcycle= (
		# red/orange
		'A62A2A', # + 166,42,42
		'CC3299', # + 204,50,153
		'FF7F00', # + 255,127,0
		# yellow/brown
		'B87333', # + 184,115,51
		# greenish
		'6B8E23', # + 107,142,35
		'238E68', # + 35,142,104
		# blue/violet
		'6982CA', # + 105,130,202
		'2222AD', # + 34,34,173
		);
	%Indices= ();
	local(*DATF);
	}
	
#--------------------------
sub draw {
  my $self = shift;
  my $gd = shift;
  my ($left,$top,$partno,$total_parts) = @_;
	my $ft= $self->feature;
	my $ftid= $ft->{id};
	my $pn= $self->panel;
  my $fgc = $self->fgcolor;
  my $bgc = $self->bgcolor;
	my($r,$g,$b);

	unless( $colors_selected ) {
		$black= $gd->colorClosest(0,0,0);
		$white= $gd->colorClosest(255,255,255);
		$yellow= $gd->colorClosest(255,255,0);
		$red= $gd->colorClosest(255,0,0);
		$lightgrey= $gd->colorClosest(225,225,225);
		$grey= $gd->colorClosest(200,200,200);
		$darkgrey= $gd->colorClosest(125,125,125);
		foreach my $ccc ( @colcycle ) {
			if( $ccc=~/^[#]?(\S\S)(\S\S)(\S\S)$/ ) {
				($r,$g,$b)= ($1,$2,$3);
				eval('$r=hex("0x'.$r.'");'); 
				eval('$g=hex("0x'.$g.'");'); 
				eval('$b=hex("0x'.$b.'");');
				}
			my $c= $gd->colorClosest($r,$g,$b);
			push(@colors,$c);
			}
		$colors_selected= 1;
		}

	my($pnstart,$pnstop)= ($pn->start,$pn->end); # in seq coordinates
	my $nseqpoints= $pnstop - $pnstart + 1;
  my($xf1,$yf1,$xf2,$yf2)= $self->calculate_boundaries($left,$top);
	my $leftpad= $pn->pad_left;
	my $serverdir= $ENV{SERVER_PATH} || $ENV{DOCUMENT_ROOT};	
	$serverdir.= '/' unless $serverdir=~m|/$|;
	my $datadir= $serverdir . $ft->{datadir};
	
	my($start,$stop)= $pn->location2pixel(($ft->{start},$ft->{end}));
	my $ftscrstop= $stop + $leftpad;
	my $ftscrstart= $start + $leftpad;

	my $chromosome= $ft->{ref};
	#warn("pn start stop leftpad nseq dir = $pnstart $pnstop $leftpad $datadir $chromosome\n");
	my $flipped= $self->{flip} ? 1 : 0;
	my($subsets,$subsetsnames,$signals)= $self->getData($ft,$datadir,$chromosome,$pnstart,$pnstop,$xf1,$xf2,$flipped);

  my $poly_pkg = $self->polygon_package;

	my @orderedsubsets= @{$subsets};
	my $nsets= $#orderedsubsets+1;
	my($xstep,$ystep)= (2,int(100.0/$nsets));
	$ystep= 7 unless $ystep>=7; # empiricaly found - to read lines of tiny fonts
	$ystep= 12 if $ystep>12; # empirically found - to preserve topo feel when number of subsets is small
	$ystep= 7 if $ystep>7; # tmp unification
	my($xw,$yw)= ( $nsets*$xstep, ($nsets-1)*$ystep );

  my $polybg= $poly_pkg->new();
  $polybg->addPt($xf1,$yf2-$yw);
  $polybg->addPt($xf2,$yf2-$yw);
  $polybg->addPt($xf2-$xw, $yf2); 
  $polybg->addPt($xf1-$xw, $yf2); 
  $gd->filledPolygon($polybg,$lightgrey); # background
	for( my $xx= $xf1+2; $xx<$xf2; $xx+=6 ) { $gd->line($xx,$yf2-$yw,$xx-$xw,$yf2,$grey); } # grid-helper
	
	my $xshift= 0;
	my $yshift= $nsets * $ystep;
	($r,$g,$b)= (10,150,80);
	my $colcycler= -1; 
	my @screencoords= @{$signals->{screencoords}};
	my $max_signal= 30;
	my $koeff= 4;
	if( exists $signals->{max_signal} ) {
		$max_signal= $signals->{max_signal};
		$koeff= 80.0/$max_signal;
		}
	my $predictor_cutoff= int($max_signal*0.95); # empirically found
	my @prevx= ();
	my @prevy= ();
	my @prevvals= ();
	my $profilen= 0;
	my %SPEEDUP= ();
	foreach my $subset ( @orderedsubsets ) {
		my $edgecolor= ($xshift==0) ? $black : $yellow;
		#my $color= ($xshift==0) ? $darkgrey : $gd->colorClosest(oct($r),oct($g),oct($b));
		my $color= ($profilen==0) ? $darkgrey : $colors[$colcycler];
		$xshift -= $xstep;
		$yshift -= $ystep;
		my @values= @{$signals->{$subset}};
		unless( exists $signals->{$subset} ) { 
			$colcycler= 0; $color= $black; # to use same colors to mark same subsets on both strands
			$gd->setThickness($ystep+1);
			$gd->line($xf1+$xshift,$yf2-$yshift+1,$xf2+$xshift+1,$yf2-$yshift+1,$white); 
			$gd->setThickness(1);
			next;
			} 
		my($xold,$yold)= ($xf1+$xshift,$yf2-$yshift+1);
		my $xpos= 0;
  	my $poly= $poly_pkg->new();
    $poly->addPt($xold,$yold+1);
		my @allx= ($xold);
		my @ally= ($yold);
		my @allvals= (0);
		my $runx= $xf1 + $xshift;
		foreach my $val ( @values ) {
			$scrx += $leftpad;
			my $x=  $screencoords[$xpos] + $xshift;
			my $visval;
			if( exists $SPEEDUP{$val} ) { $visval= $SPEEDUP{$val}; }
			else { $visval= int($val*$koeff); $SPEEDUP{$val}= $visval; }
			my $y= $yf2 - $yshift - $visval;
			push(@allx,$x);
			push(@ally,$y);
			push(@allvals,$visval);
			if( $xpos>0 ) {
    		$poly->addPt($x,$y+1);
				}
			($xold,$yold)= ($x,$y);
			$xpos++;
			}
    $poly->addPt($xf2+$xshift, $yf2-$yshift+1); 
  	$gd->filledPolygon($poly,$color) unless $profilen==0; # not on MAX predictor
		($xold,$yold)= ($allx[0],$ally[0]);
		for( my $en=1; $en<=$#allx; $en++ ) {
			my $x= $allx[$en];
			my $y= $ally[$en];
			$gd->line($xold,$yold,$x,$y,$edgecolor);
			($xold,$yold)= ($x,$y);
			} 
		if( $profilen==0 ) { # drawing mRNA (cutoff-based) predictor on MAX subset
			my($xxx,$yyy)= ($allx[1]-1,$yf2-$yw);
			$gd->line($xxx-4,$yyy,$xxx-2,$yyy,$black);
			$gd->string(GD::Font->Tiny,$xxx-12, $yyy-3,'0',$black);
			$gd->line($xxx-2,$yyy,$xxx-2,$yyy-50,$black);
			$gd->line($xxx-4,$yyy-47,$xxx-2,$yyy-50,$black);
			$gd->line($xxx,$yyy-47,$xxx-2,$yyy-50,$black);
			$gd->line($xxx-4,$yyy-44,$xxx-2,$yyy-44,$black);
			$gd->string(GD::Font->Tiny,$xxx-18, $yyy-47,$max_signal,$black);
			my($inexon,$exstart,$exend,$ymax)= (0,0,0,999);
			for( my $en=1; $en<=$#allx; $en++ ) {
				my $y= $ally[$en];
				$ymax= $y if $y < $ymax;
				if( $allvals[$en]>=$predictor_cutoff ) {
					my $x= $allx[$en];
					unless( $inexon ) { $inexon= 1; $exstart= $x; } # start exon
					$exend= $x; 
					}
				elsif( $inexon ) { # end exon and draw it
					$inexon= 0;
					$ymax -= 6;
					my $allowedymax= $yf2-$yshift-45;
					$ymax= $allowedymax if $ymax < $allowedymax; # set limit for huge peaks
					$gd->line($exstart,$ymax,$exstart,$ymax+2,$red);
					$gd->line($exstart,$ymax,$exend,$ymax,$red);
					$gd->line($exend,$ymax,$exend,$ymax+2,$red);
					($inexon,$exstart,$exend,$ymax)= (0,0,0,999);
					}
				}
			if( $inexon ) { # exon which ends beyond this screen
					$ymax -= 6;
					my $allowedymax= $yf2-$yshift-45;
					$ymax= $allowedymax if $ymax < $allowedymax; # set limit for huge peaks
					$gd->line($exstart,$ymax,$exstart,$ymax+2,$red);
					$gd->line($exstart,$ymax,$exend,$ymax,$red);
					}
			}
		if( 0 && $profilen>1 ) { # blocked - drawing roof lines doesn't work well with this view..
			for( my $en=3; $en<=$#allx; $en+=6 ) {
				next if $allvals[$en]==0 || $prevvals[$en]==0;
				my $y= $ally[$en];
				$yold= $prevy[$en];
				if( $yold<$y ) {
					my $x= $allx[$en];
					$xold= $prevx[$en];
					$gd->line($xold,$yold,$x,$y,$darkgrey);
					}
				}
			}
		$gd->string(GD::Font->Tiny,$xf2+$xshift+3, $yf2-$yshift-5,$subsetsnames->{$subset},$color);
		$colcycler++;
		$colcycler= 0 if $colcycler>$#colors;
		unless( $profilen==0 ) { @prevx= @allx; @prevy= @ally; @prevvals= @allvals; }
		$profilen++;
		}
	 
	return;
}

#--------------------------
sub getData {
  my $self = shift;
  my($ft,$datadir,$chromosome,$start,$stop,$scrstart,$scrstop,$flipped) = @_;
	my %Signals= ();
	$self->openDataFiles($datadir);
	my @subsets= (exists $ft->{'subsetsorder'}) ? @{$ft->{'subsetsorder'}} : sort split(/\t+/,$Indices{'subsets'});
	shift(@subsets) if $subsets[0] eq 'MAX';
	warn("subsets: @subsets\n") if DEBUG;
	my %SubsetsNames= map { $_, $_ } @subsets;
	%SubsetsNames=  %{$ft->{'subsetsnames'}} if(exists $ft->{'subsetsnames'});
	$SubsetsNames{MAX}= 'MAX'; 
	my $screenstep= ($scrstop-$scrstart+1) * 1.0 / ($stop-$start+1);
	my $donecoords= 0;
	foreach my $subset ( @subsets ) {
		my $nstrings= 0;
		# scan seq ranges offsets to see where to start reading
		my $key= $subset.':'.$chromosome;
		my $poskey= $key.':offsets';
		my $ranges_pos= (exists $Indices{$poskey}) ? int($Indices{$poskey}) : -1;
		if( $ranges_pos == -1 ) { next; } # no such signal..
		warn("  positioning for $poskey starts at $ranges_pos\n") if DEBUG;
		if( $start>=1000000 ) {  
			my $bigstep= int($start/1000000.0);
			if( exists $Indices{$key.':offsets:'.$bigstep} ) {
				my $jumpval= $Indices{$key.':offsets:'.$bigstep}; 
				warn("  jump in offset search to $jumpval\n") if DEBUG;
				$ranges_pos= int($jumpval); }
			}
		seek(DATF,$ranges_pos,0);
		my($offset,$offset1)= (0,0);
		my $lastseqloc= -999999999;
		my $useoffset= 0;
		while( (my $strs=<DATF>) ) {
			$nstrings++ if DEBUG;
			if( DEBUG ) {
				chop($strs); warn("  	positioning read for coord $start ($strs)\n"); }
			last unless $strs=~m/^(-?\d+)[ \t]+(\d+)/;
			my($seqloc,$fileoffset)= ($1,$2);
			if( DEBUG ) {
				warn("  positioning read for $poskey => $seqloc, $fileoffset ($strs)\n"); }
			$offset1= $offset;
			$offset= $fileoffset;
			$lastseqloc= $seqloc;
			if( $seqloc > $start ) { $useoffset= int($offset1); last; } 
			}
		warn("  will use offset $useoffset\n") if DEBUG;
		warn("  	(scanned $nstrings offset strings)\n") if DEBUG;
		if( $useoffset==0 ) { # data offset cannot be 0 - means didn't find where to read required data..
			next;
			my @emptyvals= ();
			for( my $ii= $scrstart; $ii++ <= $scrstop; ) { push(@emptyvals,0); }
			$Signals{$subset}= \@emptyvals;
			}
		$nstrings= 0;
		# read signal profile 
		seek(DATF,$useoffset,0);
		$lastseqloc= -999999999;
		my $lastsignal= 0;
		my($scrx,$scrxold)= ($scrstart,$scrstart-1);
		my $runmax= 0;
		my @values= ();
		my @xscreencoords= ();
		while( (my $str=<DATF>) ) {
			$nstrings++ if DEBUG;
			unless( $str=~m/^(-?\d+)[ \t]+(\d+)/ ) {
				warn("  header read: $str") if DEBUG;
				last; # because no headers were indexed at the beginning of data packs
				}
			my($seqloc,$signal)= ($1,$2);
			warn("  signal read: $seqloc, $signal 		line: $str") if DEBUG;
			last if $lastseqloc > $seqloc; # just in case, as all sits merged in one file..
			if( $seqloc>=$start ) { # current is the next one after the one we need to start from..
				unless( $lastseqloc== -999999999 ) { # expand previous
					$lastseqloc= $start-2 if $lastseqloc<$start; # limit empty steps (they may start from -200000)
					while( $lastseqloc < $seqloc ) { # until another (one we just retrieved) wiggle reading
						last if $lastseqloc > $stop; # end of subset data 
						next if $lastseqloc++ < $start; 
						# we have actual new seq position in our required range
						my $scrpos= int($scrx);
						$runmax= $lastsignal if $runmax < $lastsignal;
						if( $scrpos != $scrxold ) { # we have actual new seq _and_ screen position
							push(@values,$runmax);
							push(@xscreencoords,$scrpos) unless $donecoords;
							$scrxold= $scrpos;
							$runmax= 0;
							}
						$scrx += $screenstep; # remember - it is not integer
						}
					}
				}
			($lastseqloc,$lastsignal)= ($seqloc,$signal);
			last if $seqloc > $stop; # end of subset data
			}
		if( $lastseqloc < $stop ) { # if on the end of signal profile, but still in screen range
			# just assume that we are getting one more reading with signal == 0
			my $signal= 0;
			while( $lastseqloc++ < $stop ) {
				my $scrpos= int($scrx);
				if( $scrpos != $scrxold ) { # we have actual new seq _and_ screen position
							push(@values,$signal);
							push(@xscreencoords,$scrpos) unless $donecoords;
							$scrxold= $scrpos;
							}
				$scrx += $screenstep;
				}
			}
		warn("  	(scanned $nstrings signal strings)\n") if DEBUG;
		$nstrings= 0;
		if( $flipped ) {
			my @ch= reverse @values; @values= @ch;
			}
		warn("  ".$subset."=> ".@values." values @values\n") if DEBUG && $#values<1000;
		$Signals{$subset}= \@values;
		$Signals{screencoords}= \@xscreencoords unless $donecoords;
		$donecoords= 1;
		} # foreach my $subset ( @subsets ) {
	if( exists $Indices{max_signal} ) {
		$Signals{max_signal}= $Indices{max_signal};
		warn("  max_signal=> ".$Indices{max_signal}." \n") if DEBUG;
		}
	# prepare MAX profile - will be used as a base for exon/UTR prediction
	my @maxprofile= ();
	my @ruler= @{$Signals{screencoords}};
	for( my $npos= 0; $npos<=$#ruler; $npos++ ) {
		my $maxval= 0;
		foreach my $subset ( @subsets ) {
			my $p= $Signals{$subset};
			my $val= $p->[$npos];
			$maxval= $val if $maxval < $val;
			}
		push(@maxprofile,$maxval);
		}
	$Signals{MAX}= \@maxprofile;
	warn("  MAX=> ".@maxprofile." values @maxprofile\n") if DEBUG && $#maxprofile<1000;
	unshift(@subsets,'MAX');
	
	return(\@subsets,\%SubsetsNames, \%Signals);
}

#--------------------------
sub openDataFiles {
  my $self = shift;
  my $datadir= shift;
	$datadir.= '/' unless $datadir=~m|/$|;
	my $datafile= $datadir.'data.cat';
	open(DATF,$datafile) || warn("cannot open $datafile\n");
	use BerkeleyDB; # caller should already used proper 'use lib' command with path
  my $bdbfile= $datadir . 'index.bdbhash';
	tie %Indices, "BerkeleyDB::Hash", -Filename => $bdbfile, -Flags => DB_RDONLY || warn("can't read BDBHash $bdbfile\n"); 
	if( DEBUG ) { foreach my $kk ( sort keys %Indices ) { warn("	$kk => ".$Indices{$kk}."\n"); } }
	return;
}

1;

