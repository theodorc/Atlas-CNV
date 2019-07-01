#!/hgsc_software/perl/latest/bin/perl
# --------------------------------------------------------------------------------
# main atlas_cnv.pl v0, Sep 11, 2018. Ted Chiang
# Copyright 2016-2018, Baylor College of Medicine Human Genome Sequencing Center.
# All rights reserved.
# --------------------------------------------------------------------------------
use strict;
#use warnings;
#use diagnostics;
use Getopt::Long;
use Pod::Usage;
use Config::General;

pod2usage(-exitstatus => 0) if (@ARGV == 0);
my %options = ();
GetOptions (
    \%options,
    'config=s',
    'panel=s',
    'sample=s',
    'threshold_del=f',
    'threshold_dup=f',
    'help|?',
    'man'
) or pod2usage(1);
pod2usage(-exitstatus => 0, -verbose => 1) if defined($options{help});
pod2usage(-exitstatus => 0, -verbose => 2) if defined($options{man});

# -------
# logfile
# -------
open (LOG, ">>atlas_cnv.log") or die "Can't open fcd_parser.log file.\n";
my $date = `date`;
print LOG "\n$date#---------------------------\n";
chomp( my $pwd = `pwd`);


# --------------------
# Validate the inputs.
# --------------------
my $path_where_this_script_lives = $0;
$path_where_this_script_lives =~ s/\/atlas_cnv.pl//;

if (!-e $options{'config'}) {
	print LOG "No config provided: $options{'config'}. Looking for a config file in CWD or in atlas_cnv dir.\n";
	if (-e "config") {
		$options{'config'} = 'config';
		print "Found local config: $pwd/$options{'config'}\n";
		print LOG "Found config: $pwd/$options{'config'}\n";
	}
	elsif (-e "$path_where_this_script_lives/config") {
		$options{'config'} = $path_where_this_script_lives.'/config';
		print "Found config in atlas_cnv code dir: $options{'config'}\n";
		print LOG "Found config: $options{'config'}\n";
	}
	else { 
		print "No config found, please define one.. Exit.\n";
		print LOG "No config found, exiting.\n";
		exit(0);
	}
}
else { print LOG "config: $options{'config'}\n"; }

if (!-e $options{'panel'}) {
	print LOG "No panel design file provided: $options{'panel'}. Looking for panel design file in atlas_cnv dir.\n";
	if (-e "$path_where_this_script_lives/panel_design") {
		$options{'panel'} = $path_where_this_script_lives.'/panel_design';
		print "Found a panel_design file in atlas_cnv code dir: $options{'panel'}\n";
		print LOG "Found panel: $options{'panel'}\n";
	}
	else {
		print "No panel design file found, please define one.. Exit.\n";
		print LOG "No panel design file, exiting.\n";
		exit(0);
	}
}
else { print LOG "panel: $options{'panel'}\n"; }

if (!-e $options{'sample'}) {
	print "No sample file. $options{'sample'}.\n";
	exit(0);
}
else { print LOG "sample: $options{'sample'}\n"; }

if (!exists $options{'threshold_del'}) {
	print "No user defined cnv threshold_del provided, so defaults will be used.\n";
	print LOG "CNV threshold_del: defaults used, as none provided by user.\n";
}
else { print LOG "CNV threshold_del: $options{'threshold_del'}\n"; }

if (!exists $options{'threshold_dup'}) {
	print "No user defined cnv threshold_dup provided, so defaults will be used.\n";
	print LOG "CNV threshold_dup: defaults used, as none provided by user.\n";
}
else { print LOG "CNV threshold_dup: $options{'threshold_dup'}\n"; }

# -------------------------------------------------
# Load configuration file of paths and executables.
# -------------------------------------------------
# keys: 
# 	ATLASCNV
# 	LIMS
# 	JAVA
# 	RPATH
# 	RSCRIPT
my %config = new Config::General($options{'config'})->getall;


# -------------------------------------------------------------------------------
# (1a). Process sample file.
# (1b). Process the panel design file.
# (2).  Create midpool dirs.
# (3).  Copy GATK *.DATA.sample_interval_summary files to which midpool dirs.
# (4).  convert_GATK_DoC.R (GATK_DoC --> RPKM).
# (5).  Create a R rpkm matrix for each midpool.
# (6).  Call CNVS using atlas_cnv.R on each R rpkm matrix.
# (7).  Copy .cnv and .pdf files to the sample dir after making CNV dir. (SPECIFIC TO IPIPE)
# -------------------------------------------------------------------------------
# (1a).
my %sample;
my %midpool;
open (IN, "$options{'sample'}") or die "Can't open the CNV sample file: <$options{'sample'}>\n";
while (<IN>) {
	chomp;
	(my $sample, my $gender, my $midpool) = split /\t/, $_;
	(my $id, my $fclbc) = split /_/, $sample;
	$sample{ $sample } = $gender."\t".$midpool;
	$midpool{ $midpool }++;
}
close(IN);

# (1b).
my @panel_targets;
my %AUTO_GENES; my @auto_genes;
my %SEX_GENES;  my @sex_genes;
open (PANEL, "$options{'panel'}");  # matrix header from panel file $options{'panel'}
while (<PANEL>) { 
	chomp; 
        if ($_ !~ /\d/) { next; }  # skip header
	(my $exon_coor, my $gene_exon, my $callcnv, my $refseq) = split /\t/, $_; 
	if ($callcnv eq 'Y') { push @panel_targets, $exon_coor; }
        if ($exon_coor =~ /X:/ || $exon_coor =~ /Y:/) { $SEX_GENES { $gene_exon } = $exon_coor; push @sex_genes,  $gene_exon; }
        else                                          { $AUTO_GENES{ $gene_exon } = $exon_coor; push @auto_genes, $gene_exon; }
}

# (2).
for (keys %midpool) { system("mkdir $_") == 0 or die "Failed to create midpool dir: $_ ($options{'sample'}). Maybe already exists? Exit and die.\n"; }
for my $sample (keys %sample)  {
	(my $gender, my $midpool) = split /\t/, $sample{$sample};
	(my $id, my $fclbc) = split /_/, $sample;
	my $this_sample = $config{GATKDIR};           # Sample_[FCLBC]/GATK_DoC/[SAMPLE_FCLBC].DATA.sample_interval_summary
	$this_sample =~ s/\[FCLBC\]/$fclbc/;          # HMWCYBCXX-1-IDMB1
	$this_sample =~ s/\[SAMPLE_FCLBC\]/$sample/;  # 1000015113_HMWCYBCXX-1-IDMB1
	# (3)
	system("cp $this_sample $midpool/") == 0 or print LOG "Failed to copy $sample GATK-DoC file to $midpool.\n";
	# (4)
	system("$config{RSCRIPT} $config{ATLASCNV}/convert_GATK_DoC.R  $midpool/$sample.DATA.sample_interval_summary  $options{'panel'}  $gender  $midpool") ==0 or print LOG "Failed to convert_GATK_DoC.R on $sample.\n";
}

# (5)
my $matrix_header = join "\t", 'RPKM_data', @panel_targets;
for my $midpool (keys %midpool) {
	open (OUT, ">$midpool/RPKM_matrix.$midpool") or die "Can't open an outputfile RPKM_matrix.$midpool ($options{'sample'}).\n";
	print OUT "$matrix_header\n";
	opendir (MIDPOOL, "$midpool/") or die "Can't read directory: $midpool ($options{'sample'})\n";
	my @files = readdir (MIDPOOL);
	closedir (MIDPOOL);
	for (@files) {
		if ($_ =~ /(.*).rpkm.txt$/) {
			my $sample = $1;
			open (IN, "$midpool/$sample.rpkm.txt") or die "Can't open rpkm.txt file for sample: $sample.\n";
			$sample =~ s/.DATA.sample_interval_summary//; print OUT "$sample";
			while (<IN>) {
				chomp;
				if ($_ =~ /Gene_Exon/) { next; }  # skip header
				my @line = split /\t/, $_;
				my $rpkm = sprintf ( "%.0f", $line[2] );
				print OUT "\t$rpkm";
			}
			print OUT "\n";
		}
	}
	close (OUT);
	# (6)
	if (!exists $options{'threshold_del'} && !exists $options{'threshold_dup'}) {
	  system("$config{RSCRIPT} $config{'ATLASCNV'}/atlas_cnv.R  --rpkm_file $midpool/RPKM_matrix.$midpool  --panel_file $options{'panel'}  >> atlas_cnv.log") ==0 or print LOG "Failed to run atlas_cnv.R on $midpool/RPKM_matrix.$midpool.\n";
	}
	elsif (exists $options{'threshold_del'} && exists $options{'threshold_dup'}) {
	  system("$config{RSCRIPT} $config{'ATLASCNV'}/atlas_cnv.R  --rpkm_file $midpool/RPKM_matrix.$midpool  --panel_file $options{'panel'}  --threshold_del $options{'threshold_del'}  --threshold_dup $options{'threshold_dup'}  >> atlas_cnv.log") ==0 or print LOG "Failed to run atlas_cnv.R on $midpool/RPKM_matrix.$midpool.\n";
	}
	elsif (exists $options{'threshold_del'} && !exists $options{'threshold_dup'}) {
	  system("$config{RSCRIPT} $config{'ATLASCNV'}/atlas_cnv.R  --rpkm_file $midpool/RPKM_matrix.$midpool  --panel_file $options{'panel'}  --threshold_del $options{'threshold_del'}  >> atlas_cnv.log") ==0 or print LOG "Failed to run atlas_cnv.R on $midpool/RPKM_matrix.$midpool.\n";
	}
	elsif (!exists $options{'threshold_del'} && exists $options{'threshold_dup'}) {
	  system("$config{RSCRIPT} $config{'ATLASCNV'}/atlas_cnv.R  --rpkm_file $midpool/RPKM_matrix.$midpool  --panel_file $options{'panel'}  --threshold_dup $options{'threshold_dup'}  >> atlas_cnv.log") ==0 or print LOG "Failed to run atlas_cnv.R on $midpool/RPKM_matrix.$midpool.\n";
	}
}
# (7)
for my $sample (keys %sample)  {
	(my $gender, my $midpool) = split /\t/, $sample{$sample};
# 	my $sample_dir = $sample;   # $sample = 1000015132_HMWCYBCXX-1-IDMB20
# 	$sample_dir =~ s/.*_/Sample_/;
	system("if [ ! -e $midpool/CNV ]; then mkdir -p $midpool/CNV; fi") == 0 or print LOG "Failed to create CNV subdir in CNV dir. Maybe already exists?\n";
	system("cp $midpool/$sample.cnv* $midpool/CNV/") == 0 or print LOG "Failed to copy $midpool/$sample.cnv to $midpool/CNV/.\n";
	system("cp $midpool/$sample.pdf $midpool/CNV/") == 0 or print LOG "Failed to copy $midpool/$sample.pdf to $midpool/CNV/.\n";
}






# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------
=head1 NAME

atlas_cnv.pl - Call panel cnvs.

=head1 SYNOPSIS

atlas_cnv.pl [options] [file ...]

Options:

 --config             path to the configuration file [default = config]
 --panel              path to the panel design file (4 tab-delimited columns)
                      Exon_Target	Gene_Exon	Call_CNV	RefSeq
		      1:1220087-1220186	SNP_1		N		rs2144440

 --sample             path to the CNV sample tab-delimited file containing sample_id, gender, midpool
 --threshold_del      log2 ratio threshold to call a deletion [uses defaults in atlas_cnv.R]
 --threshold_dup      log2 ratio threshold to call a duplication [uses default in atlas_cnv.R]
 --help|?             prints a brief help message
 --man                prints a man page

=head1 OPTIONS

=over 8

=item B<--config-path>

The path to the file containing configuration items.

=item B<--panel>

The path to the panel design tab-delimited file. Exon_Target, Gene_Exon, Call_CNV, RefSeq

=item B<--sample>

The path to the tab-delimited file containing sample_id, gender, and midpool

=item B<--threshold_del>

The log2 ratio threshold to call a deletion log2(sample/median) [uses defaults in atlas_cnv.R]

=item B<--threshold_dup>

The log2 ratio threshold to call a duplication log2(sample/median) [uses defaults in atlas_cnv.R]

=item B<--help>

Print a brief help message and exits.

=item B<--man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<atlas_cnv.pl> will call cnvs based on a midpool set of samples using a log2 ratio of sample/median approach.

=cut
# --------------------------------------------------------------------------------------------------

