package CONF;
#-------------------------------------------------+
#    [APM] This moudle was generated by amp.pl    |
#    [APM] Created time: 2015-11-11 17:52:30      |
#-------------------------------------------------+
=pod

=head2 v1.0

Date: 2015-11-11 17:52:30

=head1 Name

CONF

=head1 Synopsis

This module is not meant to be used directly

=head1 Feedback

Author: Peng Ai
Email:  aipeng0520@163.com

=head1 Version

Version history

=cut


use strict;
use warnings;
require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(load_conf fetch_rawdata fetch_insert_size fetch_path groups2file fetch_groups fetch_group_names groups2file_qiime fetch_group_samples);

use Config::General;
use File::Basename qw(dirname basename);

use FindBin;
use lib "$FindBin::RealBin";
use lib "$FindBin::RealBin/lib";
use lib "$FindBin::RealBin/../";
use lib "$FindBin::RealBin/../lib";

use DEBUG;

sub load_conf 
{
	my $confF = shift;
	my %opts = @_;
	$confF = check_path($confF,"$confF.conf");
	
	if ($confF eq "")
	{
		ERROR('no_conf_path');
	}
	
	timeLOG("found conf file: $confF");

	my @confPath = (
		dirname($confF),
		dirname($confF)."etc",
		"$FindBin::RealBin",
		"$FindBin::RealBin/etc",
		"$FindBin::RealBin/../etc",
		"$FindBin::RealBin/.."
	);

	my $conf = Config::General->new(
		-SplitPolicy       => 'equalsign',
		-ConfigFile        => $confF,
		-AllowMultiOptions => 1,
		-LowerCaseNames    => 1,
		-IncludeAgain      => 1,
		-ConfigPath        => \@confPath,
		-AutoTrue => 1
	);	
	
	my $conf_root = { $conf->getall };
	
	&check_conf($conf_root) if ($opts{'-check'});
	&init_conf($conf_root) if ($opts{'-init'});

	return $conf_root;
}

sub check_conf
{
	my $conf = shift;
	
	# check some needed options
	my @needed = ("project_id","samples","rawdata","soft","db");
	foreach (@needed)
	{
		ERROR("[$_] is not defined in conf file!") unless $conf->{$_};
	}
}

sub init_conf 
{
	my $conf = shift;

	default($conf,"out_dir","pipe");
	default($conf,"sequence_target","16S rDNA");
	default($conf,"sequence_region","V3 + V4");
	default($conf,"length_range","300-490");

	default($conf,"phred",33);
	
	default($conf,"reads_filter","yes");
	default($conf,"filter_options","2 20 0.4 0.1");
	default($conf,"overlap_options","-m 10 -x 0.2");
	default($conf,"tags_filter","yes");
	default($conf,"tags_filter_options","-lqual 3 -conl 3 -hqual 3 -minp 75");
	default($conf,"remove_chimera","yes");
	
	default($conf,"fix_tag_num","-1");
	default($conf,"cut_tag_num","-1");
	
	default($conf,"cluster_soft","uparse");
	default($conf,"queue","all.q");
	
	my @attrs = qw/two_groups_diff multi_groups_diff otus_groups_venn otus_samples_venn adonis_groups_diff metastats_groups_diff lefse_groups_diff/;
	for (@attrs) { default($conf,$_,"") }

	init_diffs($conf);
}

sub init_diffs
{
	my $conf = shift;
	
	my %groups = fetch_group_samples($conf);
	
	# check the repeat compares
	$conf->{two_groups_diff}   = check_diffs($conf->{two_groups_diff});
	$conf->{multi_groups_diff} = check_diffs($conf->{multi_groups_diff});
	my $all_diffs = join "," , ($conf->{two_groups_diff},$conf->{multi_groups_diff});
	chop $all_diffs if ($all_diffs =~ /,$/);
	$conf->{all_diffs} = $all_diffs;

	timeLOG("all two groups diff: $conf->{two_groups_diff}");
	timeLOG("all multi groups diff: $conf->{multi_groups_diff}");

	# check the number in groups and fetch the compares whose group contains >3 samples all
	my $two_pass_diffs = fetch_replicate_diffs($conf->{two_groups_diff},\%groups);
	my $multi_pass_diffs = fetch_replicate_diffs($conf->{multi_groups_diff},\%groups);
	my $all_pass_diffs = "$two_pass_diffs,$multi_pass_diffs";
	$all_pass_diffs =~ s/,$//;
	$conf->{all_pass_diffs} = $all_pass_diffs;

	$conf->{two_pass_diffs}   = $two_pass_diffs;
	$conf->{multi_pass_diffs} = $multi_pass_diffs;

	timeLOG(">=3 replicates two groups diff: $two_pass_diffs");
	timeLOG(">=3 replicates multi groups diff: $multi_pass_diffs");

	# init otus_venn_diff 
	$conf->{otus_groups_venn} ||= $all_diffs;
	timeLOG("otus_samples_venn: $conf->{otus_samples_venn}");
	timeLOG("otus_groups_venn: $conf->{otus_groups_venn}");

	# init Anosim and Adonis 
	$conf->{adonis_groups_diff} ||= $all_pass_diffs;
	timeLOG("diffs for Anosim and Adonis: $conf->{adonis_groups_diff}");

	# init metastat diff
	$conf->{metastats_groups_diff} ||= "$conf->{two_groups_diff}";
	timeLOG("diffs for Metastats: $conf->{metastats_groups_diff}");

	# init Lefse diff 
	$conf->{lefse_groups_diff} ||= $all_pass_diffs;
	timeLOG("diffs for Lefse: $conf->{lefse_groups_diff}");
}

sub fetch_replicate_diffs
{
	my ($diffs,$groups_hash) = @_;
	my @diffs = split /,/ , $diffs;
	
	my @newdiffs;
	foreach my $diff (@diffs)
	{
		my @groups = split /&/ , $diff;
		my %groups = map { $_ => 1 } @groups;
		my $group_num = scalar keys %groups;
		ERROR("your diffs' name contains the same group",$diff) if ($group_num < $#groups+1);
		
		my $flag = 1;
		for (@groups)
		{
			my $snum = scalar @{$groups_hash->{$_}};
			if ($snum < 3)
			{
				$flag = 0;
				last;
			}
		}

		push @newdiffs , $diff if ($flag);
	}

	return ${[ join ",",@newdiffs ]}[0];
}

sub check_diffs
{
	my $diffs = shift;

	my @diffs = split /,/ , $diffs;
	my %diffs;
	for (@diffs) { $diffs{$_} ++ }
	
	my @multis = grep { $diffs{$_} > 1 } keys %diffs;

	if ($#multis >= 0)
	{
		WARN("there are some diffs are duplicated","[@multis]");
		
		my %count;
		my @newdiffs;
		for (@diffs)
		{
			push @newdiffs , $_ unless $count{$_};
			$count{$_} = 1;
		}

		return ${[ join ",",@newdiffs ]}[0];
	}
	else 
	{
		return $diffs;
	}
}

# set the default value
sub default
{
	my $hash = shift;
	my $tag = shift;
	my $val = shift;

	$hash->{$tag} = $val if (! defined $hash->{$tag});
}

sub fetch_rawdata
{
	my $conf = shift;
	my $sample = shift;
	
	my ($read1,$read2);
	if ($conf->{rawdata})
	{
		my $read_dir = check_path($conf->{rawdata});
		
		if (-e "$read_dir/$sample\_1.fq.gz" and -e "$read_dir/$sample\_2.fq.gz")
		{
			$read1 = "$read_dir/$sample\_1.fq.gz";
			$read2 = "$read_dir/$sample\_2.fq.gz";
		}
		elsif (-e "$read_dir/$sample\_1.fq" and -e "$read_dir/$sample\_2.fq")
		{
			$read1 = "$read_dir/$sample\_1.fq";
			$read2 = "$read_dir/$sample\_2.fq";
		}
		else 
		{
			ERROR("The reads is not defined!",$sample);
		}
	}
	elsif ($conf->{reads})
	{
		ERROR("The reads is not defined!",$sample) unless $conf->{reads}->{$sample};
		my @tmp = split /\s+/,$conf->{reads}->{$sample};
		$read1 = check_path($tmp[0]);
		$read2 = check_path($tmp[1]);
	}
	else 
	{
		ERROR("The reads is not defined!",$sample);
	}

	return ($read1,$read2);
}

sub fetch_insert_size
{
	my $conf = shift;
	my @samples = @_;
	my %hash;

	my @lib_size = split /[;,]/,$conf->{insert_size};

	if ($#lib_size == $#samples)
	{
		%hash = map { $samples[$_] => $lib_size[$_] } 0 .. $#samples;
	}
	elsif ($#lib_size == 0)
	{
		%hash = map { $samples[$_] => $lib_size[0] } 0 .. $#samples;
	}

	return %hash;
}

sub fetch_path
{
	my $conf = shift;
	my $name = shift;
	my %opts = @_;

	my $lcname = lc $name;
	
	if (! $conf->{soft}->{$lcname})
	{
		ERROR("the [$name] is not defined in file [software.conf]");
	}
	
	if (defined $opts{'-check'} && $opts{'-check'} == 0)
	{
		return $conf->{soft}->{$lcname};
	}

	my $path = check_path($conf->{soft}->{$lcname});

	if ($path eq "")
	{
		ERROR("the path of [$name] defined in [software.conf] is not exists!");
	}
	
	return $path;
}

sub fetch_taxon_db
{
	my $conf = shift;
	
	my %dbs = ("greengene"=>1,"silva"=>1,"rdp"=>1,"unite"=>1);

	my $type = $conf->{sequence_target};
	
	my $db = $conf->{taxon_db};
	my $dbname;

	if ($type =~ /16S/i)
	{
		if (! $db)
		{
			WARN("taxon_db is not set for 16S, greengene will be used");
		}
		if ( $db && $db eq "unite" )
		{
			WARN("unite is for ITS, so greengene will be used");
			$dbname = "greengene";
		}
		elsif ($dbs{$db})
		{
			$dbname = $db;
		}
		else 
		{
			ERROR("the 'taxon_db' must be in [greengene, silva, rdp]","[$db]");
		}
	}
	elsif ($type =~ /its/i)
	{
		$dbname = "unite";
	}
	
	return $dbname;
}

sub groups2file
{
	my $conf = shift;
	my $group_file = shift;

	open my $ofh_group , ">" , $group_file or die $!;
	
	foreach my $index (sort {$a<=>$b} keys %{$conf->{groups}})
	{
		my ($group,$samples) = split /:/ , $conf->{groups}->{$index};
		print $ofh_group "$group = $samples\n";
	}

	close $ofh_group;
}

sub groups2file_qiime
{
	my $conf = shift;
	my $group_file = shift;

	open my $ofh_group , ">" , $group_file or die $!;
	print $ofh_group "#SampleID\tGroup\n";
	foreach my $index (sort {$a<=>$b} keys %{$conf->{groups}})
	{
		my ($group,$samples) = split /:/ , $conf->{groups}->{$index};
		my @samples = split /,/ , $samples;
		print $ofh_group "$_\t$group\n" for (@samples);
	}
	close $ofh_group;
}

sub fetch_group_names
{
	my $conf = shift;
	my @group_names;
	foreach my $index (sort {$a<=>$b} keys %{$conf->{groups}})
	{
		my ($group,$samples) = split /:/ , $conf->{groups}->{$index};
		push @group_names , $group;
	}

	return @group_names;
}

sub fetch_groups
{
	my $conf = shift;
	my %groups;
	foreach my $index (sort {$a<=>$b} keys %{$conf->{groups}})
	{
		my ($group,$samples) = split /:/ , $conf->{groups}->{$index};
		my @samples = split /,/ , $samples;
		for (@samples) { $groups{$_} = $group }
	}
	return %groups;
}

sub fetch_group_samples
{
	my $conf = shift;
	my %groups;

	foreach my $index (sort {$a<=>$b} keys %{$conf->{groups}})
	{
		my ($group,$samples) = split /:/ , $conf->{groups}->{$index};
		my @samples = split /,/ , $samples;
		$groups{$group} = \@samples;
	}
	return %groups;
}