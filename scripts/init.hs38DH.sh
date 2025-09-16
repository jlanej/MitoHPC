#!/usr/bin/env bash 
#set -e

if [ -z $HP_SDIR ] ; then echo "Variable HP_SDIR not defined. Make sure you followed the SETUP ENVIRONMENT instructions" ;  fi

##############################################################################################################

# Program that setups the environmnet

# Variable HP_SDIR must be pre-set !!!

##############################################################################################################
#DIRECTORY PATHS

export HP_HDIR=`readlink -f $HP_SDIR/..`	#HP home directory
export HP_BDIR=$HP_HDIR/bin/			#bin directory
export HP_JDIR=$HP_HDIR/java/			#java directory

#Human
export HP_RDIR=$HP_HDIR/RefSeq/			#reference directory

#Mouse
#export HP_RDIR=$HP_HDIR/RefSeqMouse/           #Mouse reference directory

###############################################################
#SOFTWARE PATH

export PATH=$HP_SDIR:$HP_BDIR:$PATH

################################################################
#ALIGNMNET REFERENCE

#hs38DH(default)
export HP_RNAME=hs38DH
export HP_RMT=chrM
export HP_RNUMT="chr1:629084-634672 chr17:22521208-22521639"										                          # 150bp reads
#export HP_RNUMT="chr1:629080-634925 chr2:148881723-148881858 chr5:80651184-80651597 chr11:10508892-10509738 chr13:109424123-109424381 chr17:22521208-22521639"   # 100bp reads
export HP_RCOUNT=3366																		  # 195(hs38DH-no_alt); 194 (hs38DH-no_alt_EBV)
export HP_RURL=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa

###############################################################

export HP_E=300                  # extension(circularization)

################################################################
#GENOME REFERENCES

export HP_O=Human		 # organism: Human, Mouse...
export HP_MT=chrM                # chrM, rCRS or RSRS, FASTA file available under $HP_RDIR
export HP_MTC=chrMC
export HP_MTR=chrMR
export HP_MTLEN=16569
export HP_NUMT=NUMT              # NUMT FASTA file under $HP_RDIR

################################################################
#OTHER

export HP_CN=1			 # do compute mtDNA copy number
export HP_L=222000               # number of MT reads to subsample; empty: no subsampling; 222000 150bp reads => ~2000x MT coverage
export HP_FOPT="-q 15 -e 0"      # FASTP options: Ex: " -q 20 -e 30 "; -q: min base quality; -e: avg quality thold
export HP_DOPT="--removeDups"    # samblaster option; leave empty if no deduplication should be done
export HP_GOPT=                  # gatk mutect2 additional options : Ex "-max-reads-per-alignment-start 50" , "--mitochondria-mode"

export HP_M=mutect2 	         # SNV caller: mutect2,mutserve or freebayes
export HP_I=2		         # number of SNV iterations : 0,1,2
				 #  0: compute read counts,mtDNA-CN
                                 #  1:1 iteration (mutect2,mutserve)
                                 #  2:2 iterations (mutect2)


export HP_T1=03                  # heteroplasmy tholds
export HP_T2=05
export HP_T3=10

export HP_V=                     # SV caller: gridss
export HP_DP=100                 # minimum coverage: Ex 100
export HP_FRULE="perl -ane 'print unless(/strict_strand|strand_bias|base_qual|map_qual|weak_evidence|slippage|position|Homopolymer/ and /:0\.[01234]\d+$/);'"   # filter rule

# Automatically detect available CPU cores, with user override capability
if [ -z "$HP_P" ]; then
    # Default to number of available processors
    HP_P_AUTO=$(nproc 2>/dev/null || echo 1)
    export HP_P=$HP_P_AUTO
else
    export HP_P=$HP_P  # Use user-provided value
fi

# Set memory to 2G per core
HP_MM_CALC=$((HP_P * 2))
export HP_MM="${HP_MM_CALC}G"                                        # maximum memory (2G per core)
export HP_JOPT="-Xms$HP_MM -Xmx$HP_MM -XX:ParallelGCThreads=$HP_P"  # JAVA options
################################################################
#INPUT/OUTPUT

PWD=`pwd -P`
export HP_FDIR=$PWD/fastq/      # fastq input file directory ; .fq or .fq.gz file extension
export HP_ADIR=$PWD/bams/	# bams or crams input file directory
export HP_ODIR=$PWD/out/        # output dir
export HP_IN=$PWD/in.txt        # input file to be generated

if [ -d $HP_ADIR ] ; then
  if [ ! -s $HP_IN ] ; then
    find $HP_ADIR/  -name "*.bam" -o -name "*.cram" -readable | ls2in.pl -out $HP_ODIR | sort -V > $HP_IN
  fi
fi

###############################################################
#JOB SCHEDULING

export HP_SH="bash" ;                                                                        export HP_SHS="$HP_SH"                     # bash
#export HP_SH="sbatch -J HP_$$ --cpus-per-task=$HP_P --nodes=1 --mem=$HP_MM --time=20:00" ;  export HP_SHS="$HP_SH -d singleton"        # SLURM
#export HP_SH="qsub -V -N HP_$$ -l mem_free=$HP_MM,h_vmem=$HP_MM -pe local $HP_P -cwd" ;     export HP_SHS="$HP_SH -hold_jid HP_$$"     # SGE

