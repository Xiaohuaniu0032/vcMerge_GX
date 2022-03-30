use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use FindBin qw/$Bin/;

my ($report_dir,$outdir) = @ARGV;

my @Snvindel_tsv = glob "$report_dir/*/outputs/VcMetricActor-00/Snvindel.tsv"; # same as alleles.xls
my @cov_vcf; # SARS vcf
for my $file (@Snvindel_tsv){
	# get sample name
	my $sample_name = (split /\//, $file)[-4]; # /ChipLane2/training-test2_LibPrep60/outputs/VcMetricActor-00/Snvindel.tsv
	#my $lane = (split /\//, $file)[-5]; # ChipLane2
	#my $sample_name_final = "$lane\_$sample_name"; # ChipLane2_training-test2_LibPrep60
	my $if_cov2_vcf = &check_if_cov2_vcf($file);
	if ($if_cov2_vcf eq "YES"){
		print "[SARS-CoV-2 VCF]: $file\n";
		push @cov_vcf, $file;
	}else{
		next;
	}
}

my %sample_vars; # 记录特定样本的变异
my %all_vars;    # 记录这个run下所有样本的变异
my %sample;      # 样本名称(barcode编号)

for my $vcf (@cov_vcf){
	my $sample_name = (split /\//, $file)[-4];
	$sample{$sample_name} = 1;

	open TSV, "$vcf" or die;
	while (<TSV>){
		chomp;
		next if /^$/;
		next if /User Classification/;
		my @arr = split /\t/;
		next if (/2019-nCoV/); # keep only 2019-nCoV chrom
		next if (/ABSENT/);
		my $freq = $arr[6];
		my $var = "$arr[1]\:$arr[7]\:$arr[8]"; # 2019-nCoV:210:G:T
		$sample_vars{$sample_name}{$var} = $freq;
		my $pos = (split /\:/, $arr[1])[1];
		push @{$all_vars{$pos}}, $var; # 1) one pos may has diff variants; 2) may contain dup vars
	}
	close TSV;

	my $outfile = "$outdir/$dir\.TSVC_variants.merged.vcf.xls";
	print "[Merge VCF is]: $outfile\n";
	# Chrom/Position/Ref/Variant/IonXpress_001.Freq/IonXpress_002.Freq/.../
	# 2019-nCoV/210/G/T/1/0.98/.../
	
	open O, ">$outfile" or die;
	print O "Chrom\tPosition\tRef\tVariant";

	my @barcode = keys %sample;
	my @barcode_sort = sort {$a cmp $b} @barcode;

	for my $s (@barcode_sort){
		print O "\t$s";
	}
	print O "\n";

	# sort variant by pos
	foreach my $pos (sort { $a <=> $b } keys %all_vars){
		my $var_aref = $all_vars{$pos}; # one pos may have multi var type
		my %uniq_var;
		for my $var (@{$var_aref}){
			$uniq_var{$var} = 1;
		}

		for my $var (keys %uniq_var){ # 2019-nCoV:210:G:T
			my @var = split /\:/, $var;
			print O "$var[0]\t$var[1]\t$var[2]\t$var[3]";
			for my $s (@barcode_sort){
				my $freq;
				if (exists $sample_vars{$s}{$var}){
					$freq = $sample_vars{$s}{$var};
				}else{
					$freq = "NA";
				}
				print O "\t$freq";
			}
			print O "\n";
		}
	}
	close O;
}


sub check_if_cov2_vcf{
	my $vcf = $_[0];
	my $res;
	open IN, "$vcf" or die;
	<IN>;
	<IN>;
	my $first_line = <IN>;
	close IN;

	if ($first_line =~ /2019-nCoV/){
		$res = "YES";
	}else{
		$res = "NO";
	}

	return($res);
}




