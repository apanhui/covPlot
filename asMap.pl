#!/Bio/bin/perl
#-----------------------------------------------+
#    [APM] This script was created by amp.pl    |
#    [APM] Created time: 2018-11-18 21:24:36    |
#-----------------------------------------------+
# name: asMap.pl
# func: 
# version: 1.0

use strict;
use warnings;

use SVG;
use Config::General;
use Getopt::Long;
use File::Basename qw/basename dirname/;
use FindBin qw($Bin);
use List::Util qw/min max sum/;

use lib "$Bin/lib";
use Font;
use CONF qw/load_conf/;
use General qw/timeLOG WARN ERROR/;

use lib "/home/aipeng/work/develepment/SBV/lib";
use SBV::STAT qw(dividing bezier3_xy);

my $samtools = "/Bio/bin/samtools-1.9";

#-------------------------------------------------------------------------------
#  set the options
#-------------------------------------------------------------------------------
my %OPTS = (conf=>"$Bin/asMap.conf");
GetOptions(\%OPTS,
    'conf:s','loci:s','gene:s','bams:s','beds:s',
    'intron_fix:i','intron_scale:s','intron_background:s',
    'width:i','from_top:s','help');

&usage if ($OPTS{'help'});

unless ($OPTS{conf}){
    WARN("[FATA ERROR] you must defined the conf file!\n");
    &usage;
}

# load conf file 
my $conf = load_conf($OPTS{conf});
check_conf($conf);

# re defined some important options as var
my $legend = $conf->{legend};
my $gene   = $conf->{gene};

#-------------------------------------------------------------------------------
# load font 
my $font_path = "/Bio/User/aipeng/bin/OSGO/fonts/afm/";
load_font($font_path);

#-------------------------------------------------------------------------------
#  fetch the data
#-------------------------------------------------------------------------------
# fetch the coordinate
my %loci = fetch_gene_loci($conf->{gtf},$conf->{gene});
my ($chr,$start,$end,$strand) = ($loci{chr},$loci{start},$loci{end},$loci{strand});
timeLOG "fetch the gene coordinate done ... ";

# fetch the depth
my %depths = fetch_depths($conf,$chr,$start,$end);

my @samples = $conf->{samples} ? split /,/ , $conf->{samples} : sort keys %depths;
my %samples = map { $_ => 1 } @samples;

timeLOG "fetch the depths done ... ";

# fetch the junction
my %juncs = fetch_junctions($conf,$chr,$start,$end,$conf->{gene});
my $max_junc = fetch_max_junc(%juncs);
timeLOG "fetch the junctions done ... ";

# fetch marker info
my @marker = fetch_marker($conf,$chr,$start,$end);

#-------------------------------------------------------------------------------
# prepare to draw
#-------------------------------------------------------------------------------
my @flags   = ("exon");
my %colors  = (five_prime_utr=>"#C01C30","CDS"=>"#29B473","three_prime_utr"=>"#2A3890",exon=>"#000000");
my %heights = (five_prime_utr=>14,CDS=>18,three_prime_utr=>14,exon=>18);

@flags = reverse @flags if ($strand eq "-");

my %flags = exon_or_intron(\%loci,\@flags,$conf);
my %coord = rebuild_gene_coord(\%loci,\@flags,$conf);
my $count = $coord{$end};
timeLOG "rebuild the region coordinate system done ... ";

#-------------------------------------------------------------------------------
#  start to draw AS plot 
#-------------------------------------------------------------------------------
# define some default options
my $font = Font->new("font-size:12");
my $font_height = $font->fetch_text_height;
my $font_style = $font->toStyle;

my $sr_font = Font->new("font-size:12");
my $sr_font_h = $sr_font->fetch_text_height;
my $sr_font_style = $sr_font->toStyle;

my $margin = 20;
my $module_margin = 20;
my $ctrl_height = 20;
my $samples_spacing = 10;
my $width  = $conf->{width};
my $height = 1000;
my $spacing = 6;
my $lw = 40; # legend item width

my $max_flag_width = $font->fetch_max_text_width(\@flags);
# my $max_sample_width = $font->fetch_max_text_width(\@samples);
my $legend_width = $legend ?  $max_flag_width + $spacing*4 + $lw : 0;
my $sum_width = $width + $legend_width;

my $svg = SVG->new(width=>$sum_width,height=>$height,id=>"as_map");
my $ox = $margin;
my $oy = $margin;

## create the defs.
my $defs   = $svg->defs(id=>"defs_1");

## create backgroud color for text with a filter
my $textbg = $defs->filter(x=>0,y=>0,width=>1,height=>1,id=>"textbg");
$textbg->fe(-type=>"flood","flood-color"=>"white");
$textbg->fe(-type=>"composite",in=>"SourceGraphic");

## groups for intron backgroud
my $intron_bg = $svg->g(id=>"intron_bg");

## draw ref genes
my @tsids = keys %{$loci{loci}};
my $gene_width = $font->fetch_max_text_width( [$gene,@tsids] );
my $chr_width = $font->fetch_text_width($chr);
my $end_width = $font->fetch_text_width($end);

my $gox = $ox + $gene_width + $spacing;
my $gex = $width - $margin;

# draw chr
$svg->text(x=>$gox-$chr_width-$spacing,y=>$oy+$font_height,style=>$font_style)->cdata($chr);
$svg->text(x=>$gox,y=>$oy+$font_height,style=>$font_style)->cdata($start);
$svg->text(x=>$gex-$end_width,y=>$oy+$font_height,style=>$font_style)->cdata($end);

# draw gene 
my $geneid_width = $font->fetch_text_width($gene);
$svg->text(x=>$gox-$spacing-$geneid_width,y=>$oy+$font_height+$spacing+$font_height,style=>$font_style)->cdata($gene);
my $gene_arrowid = "gene_arrow";
end_arrow($defs,$gene_arrowid,10,6,-zh=>0,-style=>"stroke-width:1;stroke:#000;fill:#000");

my $goy = $oy + $font_height + $spacing + $font_height/2;

my ($gx1,$gx2) = ($gox,$gex);
($gx1,$gx2) = ($gx2,$gx1) if ($strand eq "-");
$svg->line(x1=>$gx1,x2=>$gx2,y1=>$goy,y2=>$goy,style=>"stroke-width:2;stroke:#000000",
    'marker-end'=>"url(#$gene_arrowid)");

# draw transcript
my $tsh = 24;
my $toy = $goy + $font_height/2 + $spacing*3;
foreach my $tsid (sort keys %{$loci{loci}}){
    my $ts_width = $font->fetch_text_width($tsid);
    $svg->text(x=>$gox-$spacing-$ts_width,y=>$toy+$tsh/2+$font_height/2,style=>$font_style)->cdata($tsid);
    $svg->line(x1=>$gox,x2=>$gex,y1=>$toy+$tsh/2,y2=>$toy+$tsh/2,style=>"stroke-width:1;stroke:#000000;stroke-dasharray:4 4 2 4");
    
    my $lastx = 0;
    foreach my $flag (@flags)
    {
        next unless $loci{loci}{$tsid}{$flag};
        my @recoreds = @{$loci{loci}{$tsid}{$flag}};
        for my $record (@recoreds)
        {
            my ($rs,$re) = @$record;
            $lastx = $rs if ($rs > $lastx);
        }
    }

    foreach my $flag (@flags)
    {
        next unless $loci{loci}{$tsid}{$flag};
        my @recoreds = @{$loci{loci}{$tsid}{$flag}};
        my $color = $colors{$flag};
        my $sw = $heights{$flag}/2;
        
        my @order_recoreds = sort {$a->[0] <=> $b->[0]} @recoreds;
        foreach my $record (@order_recoreds)
        {
            my ($rs,$re) = @$record;
            my $startx = $gox + fetch_axis_dis($coord{$rs},1,$count,$gex-$gox);
            my $endx = $gox + fetch_axis_dis($coord{$re},1,$count,$gex-$gox);
            $svg->rect(x=>$startx,y=>$toy+$tsh/2-$heights{$flag}/2,
                width=>$endx-$startx,height=>$heights{$flag},style=>"stroke-width:0;fill:$color;");
        }
    }

    $toy += $tsh + $spacing;
}

timeLOG("the gene and transcript info was drew ...");

# draw marker 
if (@marker){
    my $color = $conf->{marker_color};
    my $height = $conf->{marker_height};
    
    foreach my $mark (@marker) {
        my ($msta,$mend,$mid) = @$mark;
        my $startx = $gox + fetch_axis_dis($coord{$msta},1,$count,$gex-$gox);
        my $endx = $gox + fetch_axis_dis($coord{$mend},1,$count,$gex-$gox);
        $svg->rect(x=>$startx,y=>$toy+$tsh/2-$height/2,
            width=>$endx-$startx,height=>$height,style=>"stroke-width:0;fill:$color;");
    }
    
    push @flags , "TE";
    $heights{"TE"} = $height;
    $colors{"TE"} = $color;
    $toy += $tsh + $spacing;

    timeLOG("the marker info was drew ...");
}

# draw legend for TS
if ($legend)
{
    my $ly = $goy;
    my $lx = $gex + $spacing * 2;
    my $lh = $tsh;
    foreach my $item (@flags)
    {
        my $item_h = $heights{$item};
        $svg->rect(x=>$lx,y=>$ly+$lh/2-$item_h/2,width=>$lw,height=>$item_h,style=>"stroke-width:0;fill:$colors{$item};");
        $svg->text(x=>$lx+$lw+$spacing,y=>$ly+$lh/2+$font_height/2)->cdata($item);
        $ly += $lh + $spacing;
    }
    
    timeLOG("the Legend was drew ...");
}


# re set the height
$height = $toy - $tsh/2 - $spacing + $module_margin;

#-------------------------------------------------------------------------------
#  draw reads depth distribution and the splice reads info
#-------------------------------------------------------------------------------
my @regions = fetch_regions(\%flags,$start,$end);
my %styles = parse_styles($conf);

# draw connect from top
my $from_top = $conf->{from_top};

# unit height of depth 
my $uhd = $conf->{unit_depth_height};

# unit height of depth spacing
my $uhds = 20;

# ratio of depth height (with $uhd)
my $rdh = 0.8;

# calc depths for each samples with window
my %rdepths;
my $all_max_depth = 0;
for my $sample (@samples) {
    ERROR("${sample}'s info is not defined!") unless $styles{$sample}{window};
    my @temp = calc_depths($depths{$sample},$styles{$sample}{window},\@regions,\%coord);
    $rdepths{$sample}{depths} = \@temp;
    my @nums = map { $_->[1] } @temp;
    $rdepths{$sample}{max_depth} = max(@nums);
    $all_max_depth = $rdepths{$sample}{max_depth} if ($all_max_depth < $rdepths{$sample}{max_depth});
}

foreach my $sample (@samples){
    $height += $uhd;

    my @depths = @{$rdepths{$sample}{depths}};
    my $max_num  = $conf->{fix_axis} && $conf->{max_depth} ? $conf->{max_depth}  :
                                        $conf->{fix_axis}  ? $all_max_depth/$rdh : $rdepths{$sample}{max_depth}/$rdh;

    my $dividing = dividing(0,$max_num,-ntrue=>1,-xtrue=>1);
    my ($min,$max,$step) = split /\s/,$dividing;

    # draw depths distribution
    my @xv = ($gox);
    my @yv = ($height);
    foreach my $depth (@depths){
        my ($pos,$num) = @$depth;
        push @xv , $gox + fetch_axis_dis($pos,1,$count,$gex-$gox);
        push @yv , $height - fetch_axis_dis($num,$min,$max,$uhd);
    }
    push @xv , $gex;
    push @yv , $height;
    my $points = $svg->get_path(x=>\@xv,y=>\@yv,-type=>'polygon');
    $svg->polygon(%$points,style=>{fill=>$styles{$sample}{color},'stroke-width'=>0},id=>$sample);

    # draw y axis 
    $svg->line(x1=>$gox,x2=>$gox,y1=>$height,y2=>$height-$uhd,style=>"stroke-width:1;stroke:#000");
    my $i;
    for ($i=$min;$i<=$max;$i+=$step){
        my $ticky = $height - fetch_axis_dis($i,$min,$max,$uhd);
        $svg->line(x1=>$gox,x2=>$gox+4,y1=>$ticky,y2=>$ticky,style=>"stroke-width:1;stroke:#000");
        my $tick_width = $font->fetch_text_width($i);
        $svg->text(x=>$gox-$spacing-$tick_width,y=>$ticky+$font_height/2,style=>$font_style)->cdata($i)
    }
    
    # show samples 
    my $sample_font = $font;
    $sample_font->setAttr("fill:$styles{$sample}{sr_color};weight:bold");
    my $sample_label_width = $sample_font->fetch_text_width($sample);
    $svg->text(x=>$gex-$sample_label_width,y=>$height-$uhd+$font_height,style=>$sample_font->toStyle)->cdata($sample);
    
    # draw junctions 
    if ($juncs{$sample}){
        my @juncs = @{$juncs{$sample}};
        my %sites_count;
        my $bottom_flag = 0;
        my $pre_sta = 0;
        my $pre_end = 0;
        foreach my $junc ( sort { $a->[0] <=> $b->[0] || $a->[1] <=> $a->[1] } @juncs){
            my ($jsta,$jend,$jnum,$jid) = @$junc;
            
            next if ($jnum < $conf->{min_junction_depth});
            my $jx1 = $gox + fetch_axis_dis($coord{$jsta},1,$count,$gex-$gox);
            my $jx2 = $gox + fetch_axis_dis($coord{$jend},1,$count,$gex-$gox);
            my $cx1 = $jx1;
            my $cx2 = $jx2;
            my ($jy1,$jy2,$cy);
            $from_top = isoverlap($jsta,$jend,$pre_sta,$pre_end,$from_top,$conf->{from_top});
            if ( $from_top ){
                my $jsta_depth = fetch_site_depth($coord{$jsta},@depths);
                my $jend_depth = fetch_site_depth($coord{$jend},@depths);
                $jy1 = $height - fetch_axis_dis($jsta_depth,$min,$max,$uhd);
                $jy2 = $height - fetch_axis_dis($jend_depth,$min,$max,$uhd);
                $cy  = min($jy1,$jy2) - $ctrl_height;
            } else {
                $jy1 = $height;
                $jy2 = $height;
                $cy  = $height + $ctrl_height;
                $bottom_flag = 1;
            }
            
            # add junction connect line 
            my $pathd = "M $jx1 $jy1 C $cx1 $cy, $cx2 $cy, $jx2 $jy2";
            my $sr_size = sprintf "%.2f" , $styles{$sample}{sr_size} * $jnum / $max_junc;
            $sr_size = 1 if ($sr_size < 1);
            $svg->path(d=>$pathd,style=>"stroke-width:$sr_size;fill:none;stroke:$styles{$sample}{sr_color}",class=>"bezier");
            
            # add the splice reads num [text]
            my $jlabel_w = $sr_font->fetch_text_width($jnum);
            my $jlabel_x = ($jx1 + $jx2)/2 - $jlabel_w/2;
            my $jlabel_y = bezier3_xy([$jx1,$jy1],[$cx1,$cy],[$cx2,$cy],[$jx2,$jy2],0.5) + $sr_font_h/2;
            
            if ($conf->{text_filter}){
                $svg->text(x=>$jlabel_x,y=>$jlabel_y,style=>$sr_font_style,filter=>"url(#textbg)")->cdata($jnum);
            } else {
                my $sr_font_w = $sr_font->fetch_text_width($jnum);
                $svg->rect(x=>$jlabel_x,y=>$jlabel_y-$sr_font_h,width=>$sr_font_w,height=>$sr_font_h,style=>"stroke-width:0;fill:#fff");
                $svg->text(x=>$jlabel_x,y=>$jlabel_y,style=>$sr_font_style)->cdata($jnum);
            }

            $pre_sta = $jsta;
            $pre_end = $jend;
        }

        $height += $ctrl_height if ($bottom_flag);
    }
    
    $height +=  $samples_spacing;
}

timeLOG("the depth and junction reads number was drew ...");

#-------------------------------------------------------------------------------
#  draw the junction info, optional but is the main part
#-------------------------------------------------------------------------------
if ($conf->{show_junctions}){
    # &plot_junctions($svg);
}

# add the intron_background
if ($conf->{intron_background}){
    for(my$i=1; $i<= $#regions-2; $i+=2){
        my $intron_sta = $regions[$i]+1;
        my $intron_end = $regions[$i+1]-1;
        my $intron_x1 = $gox + fetch_axis_dis($coord{$intron_sta},1,$count,$gex-$gox);
        my $intron_x2 = $gox + fetch_axis_dis($coord{$intron_end},1,$count,$gex-$gox);
        
        $intron_bg->rect(x=>$intron_x1,y=>$goy,width=>$intron_x2-$intron_x1,
            height=>$height-$goy,style=>"stroke-width:0;fill:#ddd;fill-opacity:0.2");
    }
}

##  save figure
# reset the svg heigt
my $root = $svg->getElementByID("as_map");
$root->setAttribute("height",$height+$margin);

open OUT,">$gene.ASplot.svg" or die $!;
print OUT $svg->xmlify;
close OUT;
timeLOG("the figure [$gene.ASplot.svg] was finished :)");

#-------------------------------------------------------------------------------
#  sub functions
#-------------------------------------------------------------------------------
sub fetch_gene_loci {
    my ($gtf,$gene) = @_;

    my %flags = (exon=>1);
    my @loci;
    my %loci;

    open GTF,$gtf or die "can't open gtf file, $gtf $!";
    while(<GTF>){
        chomp;
        next unless $_;
        next if /^#/;
        
        my ($chr,$flag,$start,$end,$strand,$attrs) = (split /\t/)[0,2,3,4,6,8];
        next unless $flags{$flag};

        my $tsid     = $attrs =~ /transcript_id \"([\w\-\.]+)\";/ ? $1 : die;
        my $geneid   = $attrs =~ /gene_id \"([\w\-\.]+)\";/       ? $1 : die;
        my $genename = $attrs =~ /gene_name \"([\w\-\.]+)\";/     ? $1 : "";
        next unless ($geneid eq $gene || ($genename && $genename eq $gene) );

        push @{$loci{loci}{$tsid}{$flag}} , [$start,$end];
        $loci{chr} = $chr unless $loci{chr};
        $loci{strand} = $strand unless $loci{strand};

        push @loci , $start;
        push @loci , $end;
    }
    close GTF;
    
    ERROR("gene_not_in_gtf","$gene is not exists in gtf file!") if ($#loci == -1);

    $loci{start} = min(@loci);
    $loci{end} = max(@loci);
    
    return %loci;
}

# check the conf and set some default values 
sub check_conf {
    my $conf = shift;
    
    default_OPTS($conf,"gtf");
    default_OPTS($conf,"loci");
    default_OPTS($conf,"gene");
    default_OPTS($conf,"intron_fix");
    default_OPTS($conf,"intron_scale");
    default_OPTS($conf,"width");
    default_OPTS($conf,"unit_depth_height");
    default_OPTS($conf,"intron_background");
    $conf->{legend} = 1 if ($OPTS{legend});

    ERROR("[no_gtf_ERROR] gtf file must be defined") unless $conf->{gtf};
    ERROR("[no_loci_ERROR] the target loci must be defined") unless $conf->{loci} || $conf->{gene};
    
    default_set($conf,"width",1000);
    default_set($conf,"unit_depth_height",100);
    default_set($conf,"intron_background",1);
    default_set($conf,"intron_scale","auto");
    default_set($conf,"min_junction_depth",3);
    default_set($conf,"marker_color","#FF0000");
    default_set($conf,"marker_height",18);
    default_set($conf,"format","tophat2");

    return $conf;
}

# the OPTS defined will cover the conf file defined
sub default_OPTS {
    my ($conf,$attr) = @_;
    $conf->{$attr} = $OPTS{$attr} if ($OPTS{$attr});
}

# set the default value which not defined 
sub default_set {
    my ($conf,$attr,$val) = @_;
    $conf->{$attr} = $val unless (defined $conf->{$attr} && $conf->{$attr} ne "");
}

sub inherit_opts {
    my $parent = shift;
    my $child = shift;

    foreach my $attr (keys %$parent){
        $child->{$attr} = $parent->{$attr} unless exists $child->{$attr};
    }
}

# parse styles 
sub parse_styles {
    my $conf = shift ;
    my %styles = ();
    
    my @ases = ref $conf->{ases}->{as} eq "HASH"  ? ( $conf->{ases}->{as}  ) : 
               ref $conf->{ases}->{as} eq "ARRAY" ? @{ $conf->{ases}->{as} } : die;
   
    foreach my $as (@ases){
        inherit_opts($conf->{ases},$as);
        default_set($as,"sr_color",$as->{color});
        default_set($as,"sr_size",3);
        $styles{$as->{label}} = $as;
    }

    return %styles;
}

sub exon_or_intron {
    my ($loci,$flags,$conf) = @_;
    my ($start,$end) = ($loci->{start},$loci->{end});

    my %coord = map { $_ => 0 } $start .. $end;
    foreach my $tsid (keys %{$loci->{loci}}){
        foreach my $flag (@$flags){
            next unless $loci->{loci}{$tsid}{$flag};
            my @recoreds = @{$loci->{loci}{$tsid}{$flag}};
            foreach my $record (sort {$a->[0] <=> $b->[0]} @recoreds){
                map { $coord{$_} = 1 } $record->[0] .. $record->[1];
            }
        }
    }
    
    return %coord;
}

#-------------------------------------------------------------------------------
# rebuild the gene coordinate
# rebuild the CDS and UTR coord
# intron_fix, set all introns as length [intron_fix]
# intron_scale, zoom the intron size as [intron_scale] ratio to raw size
# intron_fix is prior to intron_scale
#-------------------------------------------------------------------------------
sub rebuild_gene_coord {
    my ($loci,$flags,$conf) = @_;
    my ($start,$end) = ($loci->{start},$loci->{end});

    my %coord = map { $_ => 0 } $start .. $end;
    foreach my $tsid (keys %{$loci->{loci}}){
        foreach my $flag (@$flags){
            next unless $loci->{loci}{$tsid}{$flag};
            my @recoreds = @{$loci->{loci}{$tsid}{$flag}};
            foreach my $record (sort {$a->[0] <=> $b->[0]} @recoreds){
                map { $coord{$_} = 1 } $record->[0] .. $record->[1];
            }
        }
    }
   

    my $count = 0;
    if ($conf->{intron_fix}){
        my $flag = 1;
        my $intron_scale = 1;
        
        foreach my $i ( $start .. $end ){
            if (1 == $coord{$i}){
                $count ++;
                $flag = 1;
            } elsif (0 == $coord{$i} && 1 == $flag){
                $intron_scale = fetch_intron_scale(\%coord,$i,$conf->{intron_fix});
                $count += $intron_scale;
                $flag = 0;
            } elsif (0 == $coord{$i} && 0 == $flag){
                $count += $intron_scale;
            }
            $coord{$i} = $count;
        }

    } elsif ($conf->{intron_scale} ne "auto") {
        map {
                $count = $coord{$_} == 1 ? $count + 1 : $count + $conf->{intron_scale};
                $coord{$_} = $count;
            } $start .. $end;
    } elsif ( $conf->{intron_scale} eq "auto" ) {
        my $exon_len = 0;
        map { $exon_len ++ if $coord{$_} == 1 } $start .. $end;
        $conf->{intron_scale} = $end-$start+1-$exon_len == 0 ? 1 : $exon_len/($end-$start+1-$exon_len);
        timeLOG("  the 'intron_scale' will be set as [$conf->{intron_scale}] auto ...");
        map {
                $count = $coord{$_} == 1 ? $count + 1 : $count + $conf->{intron_scale};
                $coord{$_} = $count;
            } $start .. $end;
        
    }

    return %coord;
}

sub fetch_intron_scale {
    my ($coord,$i,$intron_fix) = @_;
    my $intron_len = 0;
    while(exists $coord->{$i}){
        if (0 == $coord->{$i}){
            $intron_len ++;
            $i++
        } else {
            last;
        }
    }

    return  $intron_fix/$intron_len;
}

sub end_arrow
{
    my ($defs,$markid,$w,$h,%par) = @_;
    my $zh = defined $par{'-zh'} ? $par{'-zh'} : 1;
    $zh = 1 if ($zh < 0 || $zh > 1);
    my $zhx = $zh * $w;
    my $halfh = $h/2;
    my $marker = $defs->marker(
            refX => $w, refY => $h/2, # the basic point coord, is very important
            markerWidth=>$w,markerHeight=>$h,orient => 'auto',
            markerUnits => 'strokeWidth',id=>$markid);
    $marker->setAttribute('style',$par{-style}) if ($par{-style});
    $marker->setAttribute('class',$par{-class}) if ($par{-class});

    $marker->path(d => "M0,0 L$w,$halfh L0,$h L$zhx,$halfh z");
}

sub fetch_depths {
    my ($conf,$chr,$sta,$end) = @_;
    my %depths;
    
    my @ases = ref $conf->{ases}->{as} eq "HASH"  ? ( $conf->{ases}->{as}  ) : 
               ref $conf->{ases}->{as} eq "ARRAY" ? @{ $conf->{ases}->{as} } : die;

    foreach my $as (@ases){
        my $bam_file = $as->{bam};
        $depths{$as->{label}} = fetch_depth($bam_file,$chr,$sta,$end);
    }

    return %depths;
}

sub fetch_depth {
    my ($bam,$chr,$sta,$end) = @_;
    my %depth;

    ERROR("bam file is not exists, [$bam]") unless -e $bam;
    ERROR("cannot found the index for [$bam]") unless -e "$bam.bai";

    my $text = `$samtools depth -d 0 -r "$chr:$sta-$end" $bam`;
    my @lines = split /\n/ , $text;

    foreach my $line (@lines){
        my ($scf,$site,$num) = split /\t/ , $line;
        $depth{$site} = $num;
    }
    
    my @depths = map { $depth{$_} || 0 } $sta .. $end;
    return \@depths;
}

sub calc_depths {
    my ($depths,$window,$regions,$coord) = @_;
    my @region_depths = ();

    my $start = $regions->[0];

    # the first point
    my $index = 0;
    my $pos   = $index;
    my $depth = 0;
    
    for (my$i=0; $i<= $#regions-1; $i+=2){
        my $exon_sta = $regions->[$i] - $start;
        my $exon_end = $regions->[$i+1] - $start;
        
        # the exon start site
        push @region_depths , [ $exon_sta + $start , $depths->[$exon_sta] ];
        
        # the exon region with window
        my $j;
        for ($j=$exon_sta;$j<$exon_end-$window;$j+=$window){
            $depth = sum( @{$depths}[$j+1 .. $j+$window] )/$window;
            push @region_depths , [ $j + $window + $start , $depth ];
        }
        
        # besides exon region (the exon stop site)
        if ($j+$window > $exon_end){
            my @nums = @{$depths}[ $j+1 .. $exon_end ];
            $depth = sum (@nums)/($#nums+1);
            push @region_depths , [ $exon_end + $start , $depth ];
        }
        
        last if ($i == $#regions - 1); # the last exon

        ##  intron region 
        my $intron_sta = $exon_end + 1;
        my $intron_end = $regions->[$i+2] - 1 - $start;
        
        # the intron start site 
        push @region_depths , [ $intron_sta + $start , $depths->[$intron_sta] ];
        
        # the intron region with window
        for ($j=$intron_sta;$j<$intron_end-$window;$j+=$window){
            $depth = sum( @{$depths}[$j+1 .. $j+$window] )/$window;
            push @region_depths , [ $j + $window + $start , $depth ];
        }           

        # besides intron region (the intron stop site)
        if ($j+$window > $intron_end){
            my @nums = @{$depths}[ $j+1 .. $intron_end ];
            $depth = sum (@nums)/($#nums+1);
            push @region_depths , [ $intron_end + $start , $depth ];
        }
    }
    
    @region_depths = map { [ $coord->{$_->[0]} , $_->[1] ] } @region_depths;
    return @region_depths;
}

sub fetch_junctions {
    my ($conf,$chr,$sta,$end,$aimgene) = @_;
    my %junctions;
    
    my @ases = ref $conf->{ases}->{as} eq "HASH"  ? ( $conf->{ases}->{as}  ) : 
               ref $conf->{ases}->{as} eq "ARRAY" ? @{ $conf->{ases}->{as} } : die;

    foreach my $as (@ases){
        next unless $samples{$as->{label}};
        my $junc_file = $as->{junction} || next;
        $junctions{$as->{label}} = $conf->{format} eq "tophat2" ? fetch_junction($junc_file,$chr,$sta,$end) : fetch_bed($junc_file,$aimgene);
    }

    return %junctions;
}

sub fetch_junction {
    my ($junc_file,$aim_chr,$aim_sta,$aim_end) = @_;
    my @juncs;
    
    open my $fh_junc , $junc_file or die $!;
    <$fh_junc>;
    while(<$fh_junc>){
        chomp;
        my @arr = split /\t/;
        my ($chr,$sta,$end) = @arr[0,1,2];
        if ( $chr eq $aim_chr && $sta >= $aim_sta && $end <= $aim_end )
        {
            my ($junc_id,$num,$junc_len,$junc_sta) = @arr[3,4,10,11];
            my @junc_lens = split /,/ , $junc_len;
            my @junc_stas = split /,/ , $junc_sta;
            
            for my $i (1 .. $#junc_lens){
                my $jsta = $sta + $junc_stas[$i-1] + $junc_lens[$i-1];
                my $jend = $sta + $junc_stas[$i] + 1;
                push @juncs , [$jsta,$jend,$num,$junc_id];
            }
        }
    }
    close $fh_junc;

    return \@juncs;
}

sub fetch_bed {
    my ($file,$gene) = @_;
    my @juncs;

    open my $fh , $file or die $!;
    while(<$fh>){
        chomp;
        my ($chr,$sta,$end,$id,$geneid,$num) = split /\t/;
        next unless $geneid eq $gene;

        push @juncs , [$sta,$end,$num,$id];
    }
    close $fh;

    return \@juncs;
}

# fetch marker region
sub fetch_marker {
    my ($conf,$chr,$sta,$end) = @_;
    my @marker;
    
    return @marker unless $conf->{marker} ;

    open my $fh_marker , $conf->{marker} or die $!;
    while(<$fh_marker>){
        chomp;
        my ($mchr,$msta,$mend,$id) = split /\t/;
        if ($mchr eq $chr && $msta >= $sta && $mend <= $end){
            push @marker , [$msta,$mend,$id];
        }
    }
    close $fh_marker;

    return @marker;
}

# fetch the regions with exon_or_intron flags hash
sub fetch_regions {
    my ($hash,$sta,$end) = @_;
    my @exons;
    
    my @regions = ($sta);
    for my $i ( $sta+1 .. $end ){
        if ($hash->{$i} != $hash->{$i-1}){
            push @regions , $i;
        }
    }

    push @regions , $end;
    return @regions;
}

sub fetch_axis_dis{
    my ($val,$min,$max,$axis_len) = @_;
    my $dis = sprintf ( "%.4f" , (($val-$min)/($max-$min))*$axis_len );
    return $dis;
}

sub fetch_max_junc {
    my %hash = @_;
    my @maxes = map {
        my @nums = map { $_->[2] } @$_;
        @nums ? max(@nums) : 0;
        } values %hash;
    return max(@maxes);
}

sub fetch_site_depth {
    my ($site,@depths) = @_;
    
    my $res;
    foreach (@depths) {
        my ($pos,$num) = @$_;
        if ($site == $pos){
            $res = $num;
            last;
        } elsif ($site > $pos) {
            $res = $num;
            next;
        } else {
            last;
        }
    }
    return $res;
}

sub isoverlap {
    my ($s1,$e1,$s2,$e2,$flag,$rawflag) = @_;

    my $isoverlap = ($s1 >= $s2 && $s1 <= $e2) || ($e1 >= $s2 && $e1 <= $e2) ? 1 : 0;
    if ($isoverlap) { $flag = not $flag; return $flag } else { return $rawflag }
}

sub usage {
    print <<HELP;
Usage:   perl $0 [options] --conf your.conf [ --gene target_gene ]

Import Options:
         --help                      print the simple usage info
         --conf   <FILE>             the configuration file, ["$Bin/asMap.conf"]
         --gene   <STR>              set the gene name which you want to draw (must be in the gtf file)
         --gtf    <FILE>             set the gtf file

Styles Options (it's better to be set in conf file):
         --intron_fix         <INT>      fix the intron sizse as [200] nt
         --intron_scale       <FLOAT>    zoom the intron size as [auto] ratio to raw size, [scale the introns total size == exons]
         --intron_background             display the intron backgroud or not
         --width              <INT>      the width size of the figure, [1000]
         --from_top                      the connect junction lines draw from top first, otherwise will from bottom
         --unit_depth_height  <INT>      the unit depth height of each sample, [100]

HELP
    exit 1;
}
