# FlyBase TopoView glyph for GBrowse 2.x


*Warning:* This software is still in the developmental stage and is distributed
"as is", without packaging and in the same exact condition as currently used by FlyBase.
You are free to use it, modify and develop further, but proper reference
to the original author and FlyBase is required.

[TopoView](fb_shmiggle.pm) glyph (AKA shmiggle) was developed for fast
3D-like demonstration of RNA-seq data consisting of multiple
individual subsets. Main purposes were to compact presentation as
much as possible (in one reasonably sized track) and
to allow easy visual detection of coordinated behavior
of the expression profiles of different subsets.


It was found that log2 conversion dramatically changes
perception of expression profiles and kind of illuminates
coordinated behavior of different subsets. Glyph and data
indexer/formatter were in fact modified with the assumption 
that final data produced by indexer/formatter will always
be a log2 conversion of the original coverage, therefore
represented by short integer with values in range of 0-200 
or so.


Comparing performance (retrieval of several Kbp of data profiles
for several subsets of some RNA-seq experiment) of wiggle binary
method and of several possible alternatives, it was discovered that
one of the approaches remarkably outperforms wiggle bin method
(although it requires several times more space for formatted data 
storage). Optimal storage/retrieval method stores all experiment
data (all subsets of the experiment) in one text file, where
structure of the file in fact is one of the most simple wiggle
(coverage files) formats with the addition of some positioning
data (two-column format, without runlength specification, without
omission of zero values). This is the only format which glyph is able
to handle (there are many reasons for that) so any modification
of indexer/formatter **must** produce exact equivalent of that
format. In my experience, 90% of the debugging with new incoming
data was related to 
the problems of that exact format conversion. Example of the formatted
data:


```
# subset=BS107_all_unique chromosome=2LHet
-200000 0
0       0
19955   1
19959   0
19967   2
19972   0
19977   2
20027   0
20031   2
20035   0
20043   1
20045   0
20049   1
20055   0
20062   2
20069   0
20073   2
20082   0
20097   3
20115   0
20125   3
20127   0
20134   3
20139   0
20140   3
20144   0
20145   3
20150   0
20157   3
20162   0
20172   3
20183   0
```

Glyph is supplied with a [data indexer/formatter](index_cov_files.pl)
which is converting original coverage (wiggle) files into data structure which will
be used for fast retrieval. You should run this script in some separate directory,
containing original coverage files (gzipped form works too). After it finishes,
directory will contain two new files: data.cat and index.bdbhash. Both files required
for data retrieval by glyph. Files can be moved freely between different directories 
or even operational systems (Mac and PC included, I think). Content of the dat file
is subject of accurate check - this is if you want to avoid long debugging sessions
on the level of running GBrowse. Size of files is quite big, but in my experience it
is like twice less than gzipped size of all initial coverage files - which is quite 
acceptable.


Example of GBrowse conf file insert (shows actual FlyBase config sections for
Baylor and modENCODE RNA-seq tracks):

```
[baylor_wiggle]
feature       = RNAseq_profile:Baylor
glyph         = fb_shmiggle
height        = 124
bgcolor       = sub { my $f= shift;
        $f->{datadir}= '/.data/genomes/dmel/current/rnaseq-gff/baylor/'; # trick it this way..
        my @subsetsorder= qw(
                E2-4hr
                E2-16hr
                E2-16hr100
                E14-16hr
                L
                L3i
                L3i100
                P
                P3d
                MA3d
                FA3d
                A17d
                );
        $f->{subsetsorder}= \@subsetsorder;
        return 'lightgrey';
        }
key           = Baylor group RNA-seq coverage by subsets (devel.stages) [log2 converted]
category      = RNA-seq data
label         = ""
title         = ""
link = sub { my $f= shift;
  my $id= $f->{'id'};
  my $lnk="javascript:void(0);";
  "$lnk\" id=\"$id\" onmouseover=\"showdata_description('Baylor');return false;\" onmouseout=\"delsumm_overlib();";
  }

[celniker_wiggle]
feature       = RNAseq_profile:Celniker 
glyph 				= fb_shmiggle
height      	= 250
bgcolor       = sub { my $f= shift;
	$f->{datadir}= '/.data/genomes/dmel/current/rnaseq-gff/celniker/'; # trick it this way..
	my @subsetsorder= qw(
		BS40_all_unique
		BS43_all_unique
		BS46_all_unique
		BS49_all_unique
		BS54_all_unique
		BS55_all_unique
		BS58_all_unique
		BS62_all_unique
		BS66_all_unique
		BS67_all_unique
		BS71_all_unique
		BS73_all_unique
		BS107_all_unique
		BS111_all_unique 
		BS113_all_unique 
		BS196_all_unique 
		BS200_all_unique 
		BS203_all_unique 
		BS129_all_unique 
		BS133_all_unique 
		BS136_all_unique 
		BS137_all_unique 
		BS140_all_unique 
		BS143_all_unique 
		BS150_all_unique 
		BS156_all_unique 
		BS162_all_unique 
		BS153_all_unique 
		BS159_all_unique 
		BS165_all_unique 
		);
	$f->{subsetsorder}= \@subsetsorder;
	my %subsetsnames= qw(
		BS40_all_unique em0-2hr
		BS43_all_unique em2-4hr
		BS46_all_unique em4-6hr
		BS49_all_unique em6-8hr
		BS54_all_unique em8-10hr
		BS55_all_unique em10-12hr
		BS58_all_unique em12-14hr
		BS62_all_unique em14-16hr
		BS66_all_unique em16-18hr
		BS67_all_unique em18-20hr
		BS71_all_unique em20-22hr
		BS73_all_unique em22-24hr
		BS107_all_unique L1
		BS111_all_unique L2
		BS113_all_unique L3_12hr
		BS196_all_unique L3_PS1-2
		BS200_all_unique L3_PS3-6
		BS203_all_unique L3_PS7-9
		BS129_all_unique WPP
		BS133_all_unique WPP_12hr
		BS136_all_unique WPP_24hr
		BS137_all_unique WPP_2days
		BS140_all_unique WPP_3days
		BS143_all_unique WPP_4days
		BS150_all_unique AdM_Ecl_1days
		BS156_all_unique AdM_Ecl_5days
		BS162_all_unique AdM_Ecl_30days
		BS153_all_unique AdF_Ecl_1days
		BS159_all_unique AdF_Ecl_5days
		BS165_all_unique AdF_Ecl_30days
		);
	$f->{subsetsnames}= \%subsetsnames;
	return 'lightgrey';
	}
key           = modENCODE Transcription Group RNA-seq coverage (unique reads only) by subsets (devel. stages) [log2 converted]
category      = RNA-seq data
label         = "" 
title         = ""
link = sub { my $f= shift;
  my $id= $f->{'id'};
	my $lnk="javascript:void(0);";
	"$lnk\" id=\"$id\" onmouseover=\"showdata_description('Celniker');return false;\" onmouseout=\"delsumm_overlib();";
	}
```

In configuration, it is very important to set 'datadir' variable (relative
to server DOCUMENT_ROOT) so that glyph will know where to take data and index.

Setting 'subsetsorder' allows you to display expression profiles of subsets in
some predefined order. If setting omitted, glyph will display sets in alphabetical 
order of the initial subsets names.


Setting 'subsetsnames' allows to rename subsets (very important as in most cases
workflow names of subsets are unsutable for intelligent data display to end users).
If setting omitted, initial subsets names will be used for display.

For the glyph to be properly activated, you need to insert in all of your GFF files
(ones for which you have RNA-seq data) virtual contig-long features which will activate
expression data display. To cover whole range of the contig (chromosome arm), it is
better to use coordinates presented in 'sequence-region' definition at the top of GFF file.
Example of such feature lines for FlyBase data is shown below:

```
2LHet   Baylor  RNAseq_profile  1       368874  .       +       .       Comment=This is a reference feature for RNAseq wiggle tracks
2L      Baylor  RNAseq_profile  1       23011544        .       +       .       Comment=This is a reference feature for RNAseq wiggle tracks
2RHet   Baylor  RNAseq_profile  1       3288763 .       +       .       Comment=This is a reference feature for RNAseq wiggle tracks
2R      Baylor  RNAseq_profile  1       21146708        .       +       .       Comment=This is a reference feature for RNAseq wiggle tracks
3LHet   Baylor  RNAseq_profile  1       2555493 .       +       .       Comment=This is a reference feature for RNAseq wiggle tracks
3L      Baylor  RNAseq_profile  1       24543557        .       +       .       Comment=This is a reference feature for RNAseq wiggle tracks
3RHet   Baylor  RNAseq_profile  1       2517509 .       +       .       Comment=This is a reference feature for RNAseq wiggle tracks
3R      Baylor  RNAseq_profile  1       27905053        .       +       .       Comment=This is a reference feature for RNAseq wiggle tracks
4       Baylor  RNAseq_profile  1       1351857 .       +       .       Comment=This is a reference feature for RNAseq wiggle tracks
XHet    Baylor  RNAseq_profile  1       204113  .       +       .       Comment=This is a reference feature for RNAseq wiggle tracks
X       Baylor  RNAseq_profile  1       22422827        .       +       .       Comment=This is a reference feature for RNAseq wiggle tracks
YHet    Baylor  RNAseq_profile  1       347040  .       +       .       Comment=This is a reference feature for RNAseq wiggle tracks
```

Questions about TopoView glyph should be directed to Victor Strelets (strelets@bio.indiana.edu).
