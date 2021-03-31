#!/usr/bin/env perl

use File::Basename;
use warnings;

#predict SS, SA for a single sequence, SS, SA is a combination of neural network and homology

#June 23, 2004, Jianlin Cheng

####################################
### MODIFIED to work along with QMEAN
### MODIFIED to eat our own MSA
####################################

#input 
	# 1: sequence file*
	# 2: output file
        # 3: msa file*

# * Important: filenames with a '.' in them make the whole thing fail!

#Output: pxml, casp, eva format contact map for both 12 and 8 A. 
#Notice: the white space or . in seq_file are replaced by _. 

#Oct 28, 2013, SB
# Modified to create query file for get[SS|ACC]Aligns.pl here

if (@ARGV != 3)
{
	die "Usage: seq_file(fasta) output_file alignment_file\n";
}

$seq_file = shift @ARGV;
$output_file = shift @ARGV;
$alignment_file = shift @ARGV;
$output_dir = "/tmp";

if (! -f $seq_file)
{
	die "can't find file: $seq_file.\n"; 
}

if (! -d $output_dir)
{
	die "the output directory doesn't exists.\n"; 
}

if ( substr($output_dir, length($output_dir) - 1, 1) ne "/" )
{
	$output_dir .= "/"; 
}

#extract sequence file name
$slash_pos = rindex($seq_file, "/");
if ($slash_pos != -1)
{
	$seq_filename = substr($seq_file, $slash_pos + 1, length($seq_file) - $slash_pos - 1); 
}
else
{
	$seq_filename = $seq_file; 
}
if (length($seq_filename) <= 0)
{
	die "sequence file name shouldn't be less or equal 0.\n"; 
}

#non-char and . is not allowed for ouput file name 
$seq_filename =~ s/\s/_/g; 
$seq_filename =~ s/\./_/g;  

#output prefix is used as the prefix name of output files
$output_prefix = $output_dir . $$ . $seq_filename; 
#set alignment file name
$seq_filename .= "alg";

open(SEQ_FILE, "$seq_file") || die "can't open sequence file.\n";
@content = <SEQ_FILE>;
close(SEQ_FILE);
$name = shift @content;
chomp $name; 
$sequence = shift @content; 

#remove unseen dos format (cause the program fail)
$name =~ s/\s//g;
$sequence =~ s/\s//g;

#check the sequence format
if (substr($name, 0, 1) ne ">") 
{
	die "sequence file: $seq_file is not in fasta format.\n"; 
}
$name = substr($name, 1, length($name) - 1); 
$target_name = $name; 
if (length($target_name) == 0)
{
	$target_name = "unknown"; 
}
if (length($sequence) < 1)
{
	die "sequence is empty. \n"; 
}
$install_dir = "/usr/local/sspro/";

# SB: create BLAST alignments here, hand down to follow up scripts
# this should prevent calling BLAST twice
$NCBI_PATH = "/usr/local/bin";
# SB: create query file
open(TEMP_QUERY,">$output_prefix.tmp.homo.$$")
                or die "Couldn't open $output_prefix.tmp.homo.$$ for writing\n";
#add > to convert to fasta format
print TEMP_QUERY ">$seq_filename\n$sequence\n";
close(TEMP_QUERY);
# run BLAST
my $blast_file = "${output_prefix}.tmp.homo.$$.foob.txt";
`$NCBI_PATH/blastall -i $output_prefix.tmp.homo.$$ -d ${install_dir}data/pdb_large/dataset -p blastp -F F -g F -o ${blast_file}`;

#Predict Solvent Accessibility
#create input file for accpro
open(TMP, ">$output_prefix.tmp") || die "can't create temporary file for accpro.\n";
#a little bit ugly here: for accpro: alignment_file = align_dir + seq_filename
print TMP "$seq_filename\n";
print TMP "$sequence\n";
close(TMP);

$pymod_dir = dirname(__FILE__);

print "Predict solvent accessibility...";
`${pymod_dir}/homology_sa.pl $install_dir $output_prefix.tmp $alignment_file  $output_prefix.accpro $blast_file`;
print "done\n"; 

`mv $output_prefix.accpro $output_file`;

#remove the intermediate files
`rm $output_prefix.tmp`; 
`rm $blast_file`;
`rm $output_prefix.tmp.homo.$$`;
