#!/usr/bin/perl -w
#adapted the original perl script from Adrian Tan
 
use warnings;
use strict;
use POSIX;
use Getopt::Long;
use File::Path;
use File::Basename;
use Pod::Usage;
 
=head1 NAME
 
generate_simple_stuff_makefile
 
=head1 SYNOPSIS
 
 generate_simple_stuff_makefile [options]
 
  -o     output directory : location of all output files
  -l     local or slurm mode
  -m     output make file
  -j     jobname
  -i     sample id
 
 example: ./generate_simple_stuff_makefile.pl
 
=head1 DESCRIPTION
 
=cut
 
#option variables
my $help;
my $verbose;
my $debug;
my $outputDir = "/net/snowwhite/home/khlin/topmed";
my $makeFile = "makefile";
my $launchMethod = "local";
my $id = "";
my $jobName = "runmake.time.log";

#initialize options
Getopt::Long::Configure ('bundling');
 
if(!GetOptions ('h'=>\$help, 'v'=>\$verbose, 'd'=>\$debug,
                'o:s'=>\$outputDir,
                'l:s'=>\$launchMethod,
                'm:s'=>\$makeFile,
                'j:s'=>\$jobName,
                'i:s'=>\$id)
  || !defined($outputDir)
  || scalar(@ARGV)!=0)
{
    if ($help)
    {
        pod2usage(-verbose => 2);
    }
    else
    {
        pod2usage(1);
    }
}
 
if ($launchMethod ne "local" && $launchMethod ne "slurm")
{
    print STDERR "Launch method has to be local or slurm\n";
    exit(1);
}
 
##############
#print options
##############
printf("Options\n");
printf("\n");
printf("sample ID        : %s\n", $id);
printf("output directory : %s\n", $outputDir);
printf("launch method    : %s\n", $launchMethod);
printf("output jobName   : %s\n", $jobName);
printf("makefile         : %s\n", $makeFile);
printf("\n");

 
#arrays for storing targets, dependencies and commands
my @tgts = ();
my @deps = ();
my @cmds = ();
 
#temporary variables
my $tgt;
my $dep;
my @cmd;
 
mkpath($outputDir);
mkpath("output");
mkpath("1000g/$id");
my $slurmScriptsDir = "$outputDir/slurm_scripts";
mkpath($slurmScriptsDir);
my $slurmScriptNo = 0;
my $toolsDir = '/net/snowwhite/home/khlin/tools';

my $inputFiles = "";
my $inputFilesOK = "";
my $inputFile = "";
my $outputFile = "";


#############################################
#RFMix local ancestry inference using all 7 ancestral populations
#############################################


######################
#1.0. extract chromosomes of sample
######################
for my $chr (1..22)
{
    $tgt = "$outputDir/1000g/${id}/${id}_filtered_phased_chr${chr}.vcf.gz.OK";
    $dep = "";
    @cmd = ("bcftools view $outputDir/1000g/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \\
        --types snps --max-alleles 2 --exclude-uncalled -f PASS --phased --regions ${chr} --force-samples -s ${id} \\
        --output-type z --output-file $outputDir/1000g/${id}/${id}_filtered_phased_chr${chr}.vcf.gz && \\
    bcftools index -t -f $outputDir/1000g/${id}/${id}_filtered_phased_chr${chr}.vcf.gz");
    makeJob($launchMethod, $tgt, $dep, @cmd);
}

######################
#1.1. ind. common snp between sample and HGDP
######################
for my $chr (1..22)
{
    $tgt = "$outputDir/common_site/1000g_chr${chr}_HGDP_common.txt.OK";
    $dep = "$outputDir/1000g/${id}/${id}_filtered_phased_chr${chr}.vcf.gz.OK";
    @cmd = ("/net/snowwhite/home/khlin/bin/vcftools --gzvcf $outputDir/1000g/${id}/${id}_filtered_phased_chr${chr}.vcf.gz --gzdiff $outputDir/HGDP_938/HGDP_938_chr${chr}_phased.vcf.gz \\
        --diff-site --stdout | awk '{if(\$4 == \"B\")  print \$1 \"\\t\" \$2}' > $outputDir/common_site/1000g_chr${chr}_HGDP_common.txt");
    makeJob($launchMethod, $tgt, $dep, @cmd);
}

######################
#1.2. log the start time
######################

$tgt = "$outputDir/log/start.runmake.${id}_HGDP_RFMix.OK";
$dep = "";
@cmd = ("date | awk '{print \"Local ancestry pipeline\\n\\nstart: \"\$\$0}' > $outputDir/log/runmake_${id}_HGDP_RFMix_time.log");
makeJob("local", $tgt, $dep, @cmd);

######################
#1.3. merge sample and HGDP
######################
for my $chr (1..22)
{
    $tgt = "$outputDir/1000g/${id}/${id}_HGDP_chr${chr}_filtered_phased.vcf.gz.OK";
    $dep = "$outputDir/common_site/1000g_chr${chr}_HGDP_common.txt.OK";
    @cmd = ("bcftools merge -O z -o $outputDir/1000g/${id}/${id}_HGDP_chr${chr}_filtered_phased.vcf.gz -R $outputDir/common_site/1000g_chr${chr}_HGDP_common.txt \\
        $outputDir/1000g/${id}/${id}_filtered_phased_chr${chr}.vcf.gz $outputDir/HGDP_938/HGDP_938_chr${chr}_phased.vcf.gz && \\
        bcftools index -t -f $outputDir/1000g/${id}/${id}_HGDP_chr${chr}_filtered_phased.vcf.gz");
    makeJob($launchMethod, $tgt, $dep, @cmd);
}

######################
#1.4. make files needed for RFMix allele, location, classes
######################
for my $chr (1..22)
{
    $tgt = "$outputDir/1000g/${id}/${id}_HGDP_chr${chr}_RFMix_files_prep.OK";
    $dep = "$outputDir/1000g/${id}/${id}_HGDP_chr${chr}_filtered_phased.vcf.gz.OK";
    @cmd = ("bcftools view $outputDir/1000g/${id}/${id}_HGDP_chr${chr}_filtered_phased.vcf.gz \\
      | bcftools query -f '[%GT]\\n' - | sed 's/|//g' > $outputDir/1000g/${id}/${id}_HGDP_chr${chr}.alleles && \\
    bcftools view $outputDir/1000g/${id}/${id}_HGDP_chr${chr}_filtered_phased.vcf.gz \\
        | bcftools query -f '%POS\\n' - | Rscript $outputDir/utilities/generate_markerLocations_file.R stdin $outputDir/genetic_map_GRCh37/genetic_map_chr${chr}_combined_b37.txt $outputDir/1000g/${id}/${id}_HGDP_chr${chr}.locations && \\
    bcftools query -l $outputDir/1000g/${id}/${id}_HGDP_chr${chr}_filtered_phased.vcf.gz | Rscript $outputDir/utilities/make_classes_file.R \\
        stdin $outputDir/1000g/${id}/${id}_HGDP_chr${chr}.classes");
    makeJob($launchMethod, $tgt, $dep, @cmd);
}

######################
#1.5. run RFMix
######################
for my $chr (1..22)
{
    $inputFiles .= " ${outputDir}/1000g/${id}/${id}_HGDP_chr${chr}.0.Viterbi.txt";
    $inputFiles .= " $outputDir/common_site/1000g_chr${chr}_HGDP_common.txt";
    $inputFilesOK .= " $outputDir/1000g/${id}/${id}_HGDP_chr${chr}_RFMix_run.OK";
    $tgt = "$outputDir/1000g/${id}/${id}_HGDP_chr${chr}_RFMix_run.OK";
    $dep = "$outputDir/1000g/${id}/${id}_HGDP_chr${chr}_RFMix_files_prep.OK";
    @cmd = ("cd /net/snowwhite/home/khlin/tools/RFMix_v1.5.4 && \\
    python ./RunRFMix.py PopPhased $outputDir/1000g/${id}/${id}_HGDP_chr${chr}.alleles \\
        ${outputDir}/1000g/${id}/${id}_HGDP_chr${chr}.classes \\
        ${outputDir}/1000g/${id}/${id}_HGDP_chr${chr}.locations \\
        -o ${outputDir}/1000g/${id}/${id}_HGDP_chr${chr} --forward-backward --num-threads 1");
    makeJob($launchMethod, $tgt, $dep, @cmd);
}

######################
#1.6. make local ancestry plot
######################
$outputFile = "${outputDir}/output/${id}_HGDP_local_ancestry_RFMix.png";
$tgt = "${outputDir}/output/${id}_HGDP_local_ancestry_RFMix.png.OK";
$dep = "$inputFilesOK";
@cmd = ("Rscript $outputDir/utilities/plot_local_ancestry_pipeline.R ${id}_HGDP_RFMix $outputFile $inputFiles 40000");
makeJob($launchMethod, $tgt, $dep, @cmd);

######################
#1.7. log end time
######################
$tgt = "$outputDir/log/end.runmake.${id}_HGDP_RFMix.OK";
$dep = "${outputDir}/output/${id}_HGDP_local_ancestry_RFMix.png.OK";
@cmd = ("date | awk '{print \"\\nend: \"\$\$0}' >> $outputDir/log/runmake_${id}_HGDP_RFMix_time.log");
makeJob("local", $tgt, $dep, @cmd);


#############################################
#RFMix local ancestry inference using onlu europe, africa, native american
#############################################


######################
#2.0. Europe, Africa, Native America subset of HGDP
######################
for my $chr (1..22)
{
    $tgt = "$outputDir/HGDP_938/HGDP_938_chr${chr}_subset_phased.vcf.gz.OK";
    $dep = "";
    @cmd = ("bcftools view $outputDir/HGDP_938/HGDP_938_chr${chr}_phased.vcf.gz \\
        --force-samples -S $outputDir/HGDP_938/HGDP_europe_africa_native_america.txt \\
        --output-type z --output-file $outputDir/HGDP_938/HGDP_938_chr${chr}_subset_phased.vcf.gz; \\
    bcftools index -t -f $outputDir/HGDP_938/HGDP_938_chr${chr}_subset_phased.vcf.gz");
    makeJob($launchMethod, $tgt, $dep, @cmd);
}

######################
#2.1. log the start time
######################
$tgt = "$outputDir/log/start.runmake.${id}_HGDP_subset_RFMix.OK";
$dep = "";
@cmd = ("date | awk '{print \"Local ancestry pipeline using HGDP subset \\n\\nstart: \"\$\$0}' > $outputDir/log/runmake_${id}_HGDP_subset_RFMix_time.log");
makeJob("local", $tgt, $dep, @cmd);

######################
#2.2. merge sample and HGDP subset
######################
for my $chr (1..22)
{
    $tgt = "$outputDir/1000g/${id}/${id}_HGDP_subset_chr${chr}_filtered_phased.vcf.gz.OK";
    $dep = "$outputDir/common_site/1000g_chr${chr}_HGDP_common.txt.OK";
    @cmd = ("bcftools merge -O z -o $outputDir/1000g/${id}/${id}_HGDP_subset_chr${chr}_filtered_phased.vcf.gz -R $outputDir/common_site/1000g_chr${chr}_HGDP_common.txt \\
        $outputDir/1000g/${id}/${id}_filtered_phased_chr${chr}.vcf.gz $outputDir/HGDP_938/HGDP_938_chr${chr}_subset_phased.vcf.gz && \\
        bcftools index -t -f $outputDir/1000g/${id}/${id}_HGDP_subset_chr${chr}_filtered_phased.vcf.gz");
    makeJob($launchMethod, $tgt, $dep, @cmd);
}

######################
#2.3. make files needed for RFMix allele, location, classes
######################
for my $chr (1..22)
{
    $tgt = "$outputDir/1000g/${id}/${id}_HGDP_subset_chr${chr}_RFMix_files_prep.OK";
    $dep = "$outputDir/1000g/${id}/${id}_HGDP_subset_chr${chr}_filtered_phased.vcf.gz.OK";
    @cmd = ("bcftools view $outputDir/1000g/${id}/${id}_HGDP_subset_chr${chr}_filtered_phased.vcf.gz \\
      | bcftools query -f '[%GT]\\n' - | sed 's/|//g' > $outputDir/1000g/${id}/${id}_HGDP_subset_chr${chr}.alleles && \\
    bcftools view $outputDir/1000g/${id}/${id}_HGDP_subset_chr${chr}_filtered_phased.vcf.gz \\
        | bcftools query -f '%POS\\n' - | Rscript $outputDir/utilities/generate_markerLocations_file.R stdin $outputDir/genetic_map_GRCh37/genetic_map_chr${chr}_combined_b37.txt $outputDir/1000g/${id}/${id}_HGDP_subset_chr${chr}.locations && \\
    bcftools query -l $outputDir/1000g/${id}/${id}_HGDP_subset_chr${chr}_filtered_phased.vcf.gz | Rscript $outputDir/utilities/make_classes_file.R \\
        stdin $outputDir/1000g/${id}/${id}_HGDP_subset_chr${chr}.classes && \\
    sed 's/4/2/g' $outputDir/1000g/${id}/${id}_HGDP_subset_chr${chr}.classes | sed 's/5/3/g' > $outputDir/1000g/${id}/${id}_HGDP_subset_chr${chr}.classes.tmp && \\
        mv $outputDir/1000g/${id}/${id}_HGDP_subset_chr${chr}.classes.tmp $outputDir/1000g/${id}/${id}_HGDP_subset_chr${chr}.classes");
    makeJob($launchMethod, $tgt, $dep, @cmd);
}

######################
#2.4. run RFMix
######################
$inputFiles=""; #clean up
$inputFilesOK=""; #clean up
for my $chr (1..22)
{
    $inputFiles .= " ${outputDir}/1000g/${id}/${id}_HGDP_subset_chr${chr}.0.Viterbi.txt";
    $inputFiles .= " $outputDir/common_site/1000g_chr${chr}_HGDP_common.txt";
    $inputFilesOK .= " $outputDir/1000g/${id}/${id}_HGDP_subset_chr${chr}_RFMix_run.OK";
    $tgt = "$outputDir/1000g/${id}/${id}_HGDP_subset_chr${chr}_RFMix_run.OK";
    $dep = "$outputDir/1000g/${id}/${id}_HGDP_subset_chr${chr}_RFMix_files_prep.OK";
    @cmd = ("cd /net/snowwhite/home/khlin/tools/RFMix_v1.5.4 && \\
    python ./RunRFMix.py PopPhased $outputDir/1000g/${id}/${id}_HGDP_subset_chr${chr}.alleles \\
        ${outputDir}/1000g/${id}/${id}_HGDP_subset_chr${chr}.classes \\
        ${outputDir}/1000g/${id}/${id}_HGDP_subset_chr${chr}.locations \\
        -o ${outputDir}/1000g/${id}/${id}_HGDP_subset_chr${chr} --forward-backward --num-threads 1 && \\
    sed 's/2/4/g' ${outputDir}/1000g/${id}/${id}_HGDP_subset_chr${chr}.0.Viterbi.txt | sed 's/3/5/g' > ${outputDir}/1000g/${id}/${id}_HGDP_subset_chr${chr}.0.Viterbi.txt.tmp && \\
        mv ${outputDir}/1000g/${id}/${id}_HGDP_subset_chr${chr}.0.Viterbi.txt.tmp ${outputDir}/1000g/${id}/${id}_HGDP_subset_chr${chr}.0.Viterbi.txt");
    makeJob($launchMethod, $tgt, $dep, @cmd);
}

######################
#2.5. make local ancestry plot
######################
$outputFile = "${outputDir}/output/${id}_HGDP_subset_local_ancestry_RFMix.png";
$tgt = "${outputDir}/output/${id}_HGDP_subset_local_ancestry_RFMix.png.OK";
$dep = "$inputFilesOK";
@cmd = ("Rscript $outputDir/utilities/plot_local_ancestry_pipeline.R ${id}_HGDP_subset_RFMix $outputFile $inputFiles 40000");
makeJob($launchMethod, $tgt, $dep, @cmd);

######################
#2.6. log end time
######################
$tgt = "$outputDir/log/end.runmake.${id}_HGDP_subset_RFMix.OK";
$dep = "${outputDir}/output/${id}_HGDP_subset_local_ancestry_RFMix.png.OK";
@cmd = ("date | awk '{print \"\\nend: \"\$\$0}' >> $outputDir/log/runmake_${id}_HGDP_subset_RFMix_time.log");
makeJob("local", $tgt, $dep, @cmd);


#############################################
#local ancestry inference using LAMPLD
#############################################


######################
#3.0. create subset reference panel by ancestral populations
######################
$tgt = "$outputDir/HGDP_938/LAMPLD/create_LAMP_ref_hap.OK";
$dep = "";
@cmd = ("./HGDP_938/create_LAMP_ref_hap.sh");
makeJob("local", $tgt, $dep, @cmd);

######################
#3.1. log the start time
######################
$tgt = "$outputDir/log/start.runmake.${id}_HGDP_subset_LAMPLD.OK";
$dep = "";
@cmd = ("date | awk '{print \"Local ancestry pipeline using HGDP subset \\n\\nstart: \"\$\$0}' > $outputDir/log/runmake_${id}_HGDP_subset_LAMPLD_time.log");
makeJob("local", $tgt, $dep, @cmd);

######################
#3.2. convert sample vcf with common snp to HGDP to genotype dosage
######################
for my $chr (1..22)
{
    $tgt = "$outputDir/1000g/${id}/${id}_chr${chr}_lamp.012.OK";
    $dep = "";
    @cmd = ("vcftools --gzvcf $outputDir/1000g/${id}/${id}_filtered_phased_chr${chr}.vcf.gz --012 --out $outputDir/1000g/${id}/${id}_chr${chr}_lamp --positions $outputDir/common_site/1000g_chr${chr}_HGDP_common.txt && \\
        rm -f $outputDir/1000g/${id}/${id}_chr${chr}_lamp.012.indv $outputDir/1000g/${id}/${id}_chr${chr}_lamp.012.pos $outputDir/1000g/${id}/${id}_chr${chr}_lamp.log && \\
        cut -f 2- $outputDir/1000g/${id}/${id}_chr${chr}_lamp.012 | sed 's/\t//g' > $outputDir/1000g/${id}/${id}_chr${chr}_lamp.012.tmp && mv $outputDir/1000g/${id}/${id}_chr${chr}_lamp.012.tmp $outputDir/1000g/${id}/${id}_chr${chr}_lamp.012");
    makeJob($launchMethod, $tgt, $dep, @cmd);
}

######################
#3.3. run LAMPLD
######################
for my $chr (1..22)
{
    $tgt = "$outputDir/1000g/${id}/${id}_chr${chr}_lampped.out.OK";
    $dep = "$outputDir/HGDP_938/LAMPLD/create_LAMP_ref_hap.OK $outputDir/1000g/${id}/${id}_chr${chr}_lamp.012.OK";
    @cmd = ("perl /net/snowwhite/home/khlin/tools/LAMPLD-v1.1/run_LAMPLD.pl $outputDir/HGDP_938/LAMPLD/chr${chr}.pos \\
        $outputDir/HGDP_938/LAMPLD/HGDP_europe_chr${chr}.impute.hap \\
        $outputDir/HGDP_938/LAMPLD/HGDP_native_america_chr${chr}.impute.hap \\
        $outputDir/HGDP_938/LAMPLD/HGDP_africa_chr${chr}.impute.hap \\
        $outputDir/1000g/${id}/${id}_chr${chr}_lamp.012 \\
        $outputDir/1000g/${id}/${id}_chr${chr}_lampped.out");
    makeJob($launchMethod, $tgt, $dep, @cmd);
}

######################
#3.4. convert plot format 
######################
$inputFiles=""; #clean up
$inputFilesOK=""; #clean up
for my $chr (1..22)
{   $inputFiles .= " $outputDir/1000g/${id}/${id}_chr${chr}_lampped_plot.out";
    $inputFiles .= " $outputDir/common_site/1000g_chr${chr}_HGDP_common.txt";
    $inputFilesOK .= " $outputDir/1000g/${id}/${id}_chr${chr}_lampped_plot.out.OK";
    $tgt = "$outputDir/1000g/${id}/${id}_chr${chr}_lampped_plot.out.OK";
    $dep = "$outputDir/1000g/${id}/${id}_chr${chr}_lampped.out.OK";
    @cmd = ("python $outputDir/utilities/lamped_out_2_plot.py -i $outputDir/1000g/${id}/${id}_chr${chr}_lampped.out \\
        -p $outputDir/common_site/1000g_chr${chr}_HGDP_common.txt -o $outputDir/1000g/${id}/${id}_chr${chr}_lampped_plot.out");
    makeJob($launchMethod, $tgt, $dep, @cmd);
}

######################
#3.5. make local ancestry plot
######################
$outputFile = "${outputDir}/output/${id}_HGDP_subset_local_ancestry_LAMPLD.png";
$tgt = "${outputDir}/output/${id}_HGDP_subset_local_ancestry_LAMPLD.png.OK";
$dep = "$inputFilesOK";
@cmd = ("Rscript $outputDir/utilities/plot_local_ancestry_pipeline.R ${id}_HGDP_subset_LAMPLD $outputFile $inputFiles 40000");
makeJob($launchMethod, $tgt, $dep, @cmd);

######################
#3.6. log end time
######################
$tgt = "$outputDir/log/end.runmake.${id}_HGDP_subset_LAMPLD.OK";
$dep = "${outputDir}/output/${id}_HGDP_subset_local_ancestry_LAMPLD.png.OK";
@cmd = ("date | awk '{print \"\\nend: \"\$\$0}' >> $outputDir/log/runmake_${id}_HGDP_subset_LAMPLD_time.log");
makeJob("local", $tgt, $dep, @cmd);


#############################################
#local ancestry inference using Lanc-CSV
#############################################


######################
#4.0. log the start time
######################
$tgt = "$outputDir/log/start.runmake.${id}_1000g_subset_LancCSV.OK";
$dep = "";
@cmd = ("date | awk '{print \"Local ancestry pipeline using 1000g subset \\n\\nstart: \"\$\$0}' > $outputDir/log/runmake_${id}_1000g_subset_LancCSV_time.log");
makeJob("local", $tgt, $dep, @cmd);

######################
#4.1. convert sample vcf with all snps to to genotype dosage and create position file
######################
for my $chr (1..22)
{
    $tgt = "$outputDir/1000g/${id}/${id}_chr${chr}_LancCSV.012.OK";
    $dep = "$outputDir/1000g/${id}/${id}_filtered_phased_chr${chr}.vcf.gz.OK";
    @cmd = ("vcftools --gzvcf $outputDir/1000g/${id}/${id}_filtered_phased_chr${chr}.vcf.gz --012 --out $outputDir/1000g/${id}/${id}_chr${chr}_LancCSV && \\
        rm -f $outputDir/1000g/${id}/${id}_chr${chr}_LancCSV.012.indv $outputDir/1000g/${id}/${id}_chr${chr}_LancCSV.log && \\
        cut -f 2- $outputDir/1000g/${id}/${id}_chr${chr}_LancCSV.012 | sed 's/\t//g' > $outputDir/1000g/${id}/${id}_chr${chr}_LancCSV.012.tmp && mv $outputDir/1000g/${id}/${id}_chr${chr}_LancCSV.012.tmp $outputDir/1000g/${id}/${id}_chr${chr}_LancCSV.012 && \\
        cut -f 2 $outputDir/1000g/${id}/${id}_chr${chr}_LancCSV.012.pos > $outputDir/1000g/${id}/${id}_chr${chr}_LancCSV.012.pos.tmp && mv $outputDir/1000g/${id}/${id}_chr${chr}_LancCSV.012.pos.tmp $outputDir/1000g/${id}/${id}_chr${chr}_LancCSV.012.pos");
    makeJob($launchMethod, $tgt, $dep, @cmd);
}

######################
#4.2. run LancCSV
######################
for my $chr (1..22)
{
    $tgt = "$outputDir/1000g/${id}/${id}_chr${chr}_LancCSV.out.OK";
    $dep = "$outputDir/1000g/${id}/${id}_chr${chr}_LancCSV.012.OK";
    @cmd = ("$toolsDir/LancCSVRelease_beta_v0_1/LancCSV $toolsDir/LancCSVRelease_beta_v0_1/Database/DatabaseLancCSV_04302014_Chr${chr}.txt $outputDir/1000g/${id}/${id}_chr${chr}_LancCSV.012.pos $toolsDir/LancCSVRelease_beta_v0_1/AdmixtureFileASW.txt $outputDir/1000g/${id}/${id}_chr${chr}_LancCSV.012 $outputDir/1000g/${id}/${id}_chr${chr}_LancCSV.out");
    makeJob($launchMethod, $tgt, $dep, @cmd);
}

######################
#4.3. convert plot format and add chr column in position file
######################
$inputFiles=""; #clean up
$inputFilesOK=""; #clean up
for my $chr (1..22)
{   $inputFiles .= " $outputDir/1000g/${id}/${id}_chr${chr}_LancCSV_plot.out";
    $inputFiles .= " $outputDir/1000g/${id}/${id}_chr${chr}_LancCSV.012.pos";
    $inputFilesOK .= " $outputDir/1000g/${id}/${id}_chr${chr}_LancCSV_plot.out.OK";
    $tgt = "$outputDir/1000g/${id}/${id}_chr${chr}_LancCSV_plot.out.OK";
    $dep = "$outputDir/1000g/${id}/${id}_chr${chr}_LancCSV.out.OK";
    @cmd = ("python $outputDir/utilities/lancsv_out_2_plot.py -i $outputDir/1000g/${id}/${id}_chr${chr}_LancCSV.out -p $outputDir/1000g/${id}/${id}_chr${chr}_LancCSV.012.pos \\
        -o $outputDir/1000g/${id}/${id}_chr${chr}_LancCSV_plot.out && \\
        awk '{print \"${chr}\t\"\$\$1}' $outputDir/1000g/${id}/${id}_chr${chr}_LancCSV.012.pos > $outputDir/1000g/${id}/${id}_chr${chr}_LancCSV.012.pos.tmp && mv $outputDir/1000g/${id}/${id}_chr${chr}_LancCSV.012.pos.tmp $outputDir/1000g/${id}/${id}_chr${chr}_LancCSV.012.pos
");
    makeJob($launchMethod, $tgt, $dep, @cmd);
}

######################
#4.4. make local ancestry plot
######################
$outputFile = "${outputDir}/output/${id}_1000g_subset_local_ancestry_LancCSV.png";
$tgt = "${outputDir}/output/${id}_1000g_subset_local_ancestry_LancCSV.png.OK";
$dep = "$inputFilesOK";
@cmd = ("--mem-per-cpu=4096 Rscript $outputDir/utilities/plot_local_ancestry_pipeline.R ${id}_1000g_subset_LancCSV $outputFile $inputFiles 40000");
makeJob($launchMethod, $tgt, $dep, @cmd);

######################
#4.5. log end time
######################
$tgt = "$outputDir/log/end.runmake.${id}_1000g_subset_LancCSV.OK";
$dep = "${outputDir}/output/${id}_1000g_subset_local_ancestry_LancCSV.png.OK";
@cmd = ("date | awk '{print \"\\nend: \"\$\$0}' >> $outputDir/log/runmake_${id}_1000g_subset_LancCSV_time.log");
makeJob("local", $tgt, $dep, @cmd);


#############################################
#local ancestry inference using EILA
#############################################


######################
#5.0. log the start time
######################
$tgt = "$outputDir/log/start.runmake.${id}_HGDP_subset_EILA.OK";
$dep = "";
@cmd = ("date | awk '{print \"Local ancestry pipeline using HGDP subset \\n\\nstart: \"\$\$0}' > $outputDir/log/runmake_${id}_HGDP_subset_EILA_time.log");
makeJob("local", $tgt, $dep, @cmd);

######################
#5.1. convert sample vcf with all snps to to genotype dosage and create position file
######################
for my $chr (1..22)
{
    $tgt = "$outputDir/1000g/${id}/${id}_chr${chr}_EILA.012.OK";
    $dep = "$outputDir/1000g/${id}/${id}_filtered_phased_chr${chr}.vcf.gz.OK";
    @cmd = ("vcftools --gzvcf $outputDir/1000g/${id}/${id}_filtered_phased_chr${chr}.vcf.gz --012 --positions /net/snowwhite/home/khlin/topmed/common_site/1000g_chr${chr}_HGDP_common.txt --out $outputDir/1000g/${id}/${id}_chr${chr}_EILA && \\
            rm -f $outputDir/1000g/${id}/${id}_chr${chr}_EILA.012.indv $outputDir/1000g/${id}/${id}_chr${chr}_EILA.log && \\
            cut -f 2- $outputDir/1000g/${id}/${id}_chr${chr}_EILA.012 | awk -f /net/snowwhite/home/khlin/tools/transpose_withSpace.awk > $outputDir/1000g/${id}/${id}_chr${chr}_EILA.012.tmp && mv $outputDir/1000g/${id}/${id}_chr${chr}_EILA.012.tmp $outputDir/1000g/${id}/${id}_chr${chr}_EILA.012 && \\
            cut -f 2 $outputDir/1000g/${id}/${id}_chr${chr}_EILA.012.pos > $outputDir/1000g/${id}/${id}_chr${chr}_EILA.012.pos.tmp && mv $outputDir/1000g/${id}/${id}_chr${chr}_EILA.012.pos.tmp $outputDir/1000g/${id}/${id}_chr${chr}_EILA.012.pos");
    makeJob($launchMethod, $tgt, $dep, @cmd);
}

######################
#5.2. run EILA
######################
for my $chr (1..22)
{
    $tgt = "$outputDir/1000g/${id}/${id}_chr${chr}_EILA.out.OK";
    $dep = "$outputDir/1000g/${id}/${id}_chr${chr}_EILA.012.OK";
    @cmd = ("Rscript $outputDir/utilities/eila.R inp=\\\"$outputDir/1000g/${id}/${id}_chr${chr}_EILA.012\\\" pos=\\\"$outputDir/1000g/${id}/${id}_chr${chr}_EILA.012.pos\\\" eur=\\\"/net/snowwhite/home/khlin/topmed/HGDP_938/EILA/HGDP_europe_chr${chr}.012\\\" afr=\\\"/net/snowwhite/home/khlin/topmed/HGDP_938/EILA/HGDP_africa_chr${chr}.012\\\" nat=\\\"/net/snowwhite/home/khlin/topmed/HGDP_938/EILA/HGDP_native_america_chr${chr}.012\\\" out=\\\"/net/snowwhite/home/khlin/topmed/1000g/${id}/${id}_chr${chr}_EILA.out\\\"");
    makeJob($launchMethod, $tgt, $dep, @cmd);
}

######################
#5.3. convert plot format and add chr column in position file
######################
$inputFiles=""; #clean up
$inputFilesOK=""; #clean up
for my $chr (1..22)
{   $inputFiles .= " $outputDir/1000g/${id}/${id}_chr${chr}_EILA_plot.out";
    $inputFiles .= " $outputDir/1000g/${id}/${id}_chr${chr}_EILA.012.pos";
    $inputFilesOK .= " $outputDir/1000g/${id}/${id}_chr${chr}_EILA_plot.out.OK";
    $tgt = "$outputDir/1000g/${id}/${id}_chr${chr}_EILA_plot.out.OK";
    $dep = "$outputDir/1000g/${id}/${id}_chr${chr}_EILA.out.OK";
    @cmd = ("python $outputDir/utilities/eila_out_2_plot.py -i $outputDir/1000g/${id}/${id}_chr${chr}_EILA.out -p $outputDir/1000g/${id}/${id}_chr${chr}_EILA.012.pos \\
        -o $outputDir/1000g/${id}/${id}_chr${chr}_EILA_plot.out && \\
        awk '{print \"${chr}\t\"\$\$1}' $outputDir/1000g/${id}/${id}_chr${chr}_EILA.012.pos > $outputDir/1000g/${id}/${id}_chr${chr}_EILA.012.pos.tmp && mv $outputDir/1000g/${id}/${id}_chr${chr}_EILA.012.pos.tmp $outputDir/1000g/${id}/${id}_chr${chr}_EILA.012.pos
");
    makeJob($launchMethod, $tgt, $dep, @cmd);
}

######################
#5.4. make local ancestry plot
######################
$outputFile = "${outputDir}/output/${id}_HGDP_subset_local_ancestry_EILA.png";
$tgt = "${outputDir}/output/${id}_HGDP_subset_local_ancestry_EILA.png.OK";
$dep = "$inputFilesOK";
@cmd = ("Rscript $outputDir/utilities/plot_local_ancestry_pipeline.R ${id}_HGDP_subset_EILA $outputFile $inputFiles 40000");
makeJob($launchMethod, $tgt, $dep, @cmd);

######################
#5.5. log end time
######################
$tgt = "$outputDir/log/end.runmake.${id}_HGDP_subset_EILA.OK";
$dep = "${outputDir}/output/${id}_HGDP_subset_local_ancestry_EILA.png.OK";
@cmd = ("date | awk '{print \"\\nend: \"\$\$0}' >> $outputDir/log/runmake_${id}_HGDP_subset_EILA_time.log");
makeJob("local", $tgt, $dep, @cmd);

#*******************
#Write out make file
#*******************
open(MAK,">$makeFile") || die "Cannot open $makeFile\n";
print MAK ".DELETE_ON_ERROR:\n\n";
print MAK "all: @tgts\n\n";

#clean
push(@tgts, "clean");
push(@deps, "");
push(@cmds, "\t-rm -rf *NA* log/*.* output/*.OK $outputDir/1000g/NA*/*.* $slurmScriptsDir/*.*");

#cleam_ref
push(@tgts, "clean_ref");
push(@deps, "");
push(@cmds, "\t-rm -rf $outputDir/HGDP_938/*.OK $outputDir/common_site/*.*");

#clean_LAMPLD
push(@tgts, "clean_lampld");
push(@deps, "");
push(@cmds, "\t-rm -rf $outputDir/1000g/NA*/*LAMPLD* $outputDir/1000g/NA*/*lampped* $outputDir/1000g/NA*/*lamp*");

#clean_lanccsv
push(@tgts, "clean_lanccsv");
push(@deps, "");
push(@cmds, "\t-rm -rf $outputDir/1000g/NA*/*LancCSV* ");

#clean_jobs
push(@tgts, "clean_job");
push(@deps, "");
push(@cmds, "\tps xu | awk '{if(\$\$11==\"make\") print \$\$2}' | xargs kill; squeue -ukhlin | awk '{if(\$\$3 ~ /NA[0-9]/) print \$\$1}' | xargs -0 scancel");

 
for(my $i=0; $i < @tgts; ++$i)
{
    print MAK "$tgts[$i]: $deps[$i]\n";
    print MAK "$cmds[$i]\n";
}
close MAK;
 
##########
#functions
##########
 
#run a job either locally or by slurm
sub makeJob
{
    my ($method, $tgt, $dep, @cmd) = @_;
 
    if ($method eq "local")
    {
        makeLocalStep($tgt, $dep, @cmd);
    }
    elsif ($method eq "slurm")
    {
        makeSlurm($tgt, $dep, @cmd);
    }
}
 
#run slurm jobs
sub makeSlurm
{
    my ($tgt, $dep, @cmd) = @_;

    push(@tgts, $tgt);
    push(@deps, $dep);
    my $cmd = "";
    for my $c (@cmd)
    {
        #contains pipe or cd
        if ($c=~/\||^cd/)
        {
            ++$slurmScriptNo;
            my $slurmScriptFile = "$slurmScriptsDir/${slurmScriptNo}_${id}.sh";
            open(IN, ">$slurmScriptFile");
            print IN "#!/bin/bash\n"; 
            print IN "set pipefail; $c"; 
            close(IN);
            chmod(0755, $slurmScriptFile);
            
            # $cmd .= "echo '" . $c . "'\n";
            $cmd .= "\tsrun -p nomosix,main -J $id -D $outputDir $slurmScriptFile\n";
        }
        else
        {
            $cmd .= "\tsrun -p nomosix,main -J $id -D $outputDir " . $c . "\n";
        }
    }
    $cmd .= "\ttouch $tgt\n";
    push(@cmds, $cmd);
}
 
#run a local job
sub makeLocalStep
{
    my ($tgt, $dep, @cmd) = @_;
 
    push(@tgts, $tgt);
    push(@deps, $dep);
    my $cmd = "";
    for my $c (@cmd)
    {
        $cmd .= "\t" . $c . "\n";
    }
    $cmd .= "\ttouch $tgt\n";
    push(@cmds, $cmd);
}
