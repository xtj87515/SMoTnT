README file for Stochastic Model of Transcription and Translation (SMoTnT)

*************************************************************************************

VERSION: SMoTnT
DATE: Apr 28, 2021

gcc sourceCode/SMOTNT.c -g -lm -lgsl -lgslcblas -mtune=generic -O3 -o sourceCode/SMOTNT

./sourceCode/SMOTNT -Tt 1500 -Tb 1000 -R 200000 -t 3300000 -N 4839 -dc 0.2 -r 0.4 -F ../../publicInput/S.cer.genom -C ../../publicInput/S.cer.tRNA -D simInput/allGenesDecrEqualsSynrWithScaling_dc_0.2_r_0.4_S.cer.mRNA.ini.abndc.syn.dec.tsv -s 16629 -O simOutput/output.seed.16629_dc_0.2_r_0.4_allGenesDecrEqualsSynrScaling -p1 -p2 -p3 -p4 -p6 -p7 


Options:

	-V <value>    Volume of the cell in m^3/s. The minimum volume of the cell is set to
			contain at least 1000 ribosomes and tRNAs.
			[DEFAULT]  -V 4.2E-17 (volume of yeast cell)

	-T[CHAR]    Specify times for various events

			-Tt:    Total simulation time in seconds.
			    [DEFAULT]  -Tt 1500

			-Tb:    Burn-in/threshold time. Time spent by the cell to reach equilibrium.
			    Only calculations after this time will be included in the analyses.
			    [DEFAULT]  -Tb 1000

	-R <value>    Total number of ribosomes in the cell.
			[DEFAULT]  -R 200000

	-t <value>    Total number of tRNAs in the cell.
			[DEFAULT]  -t 3300000

	-N <value>    Total number of genes. This needs to be specified by the user.
			[DEFAULT]  -N 1

	-dc <value>    The ratio of 5'-3' co-translational mRNA decay for all mRNAs [0,1].
			[DEFAULT]  -dc 0

	-r <value>    The ribosome protection index [0,1]
			[DEFAULT]  -r 0

	-F <FILE>    File containing processed fasta file into a numeric sequence.
			This file is an output of the code utilities/convert.fasta.to.genom.pl
			It contains the information regarding initiation probability, mRNA
			abundance and codon sequence of each gene.
			[DEFAULT]  -F ../publicInput/truncated_S.cer.genom

	-D <FILE>    File containing the information about gene specific initiation
			probability, mRNA level, gene specific synthesis rate and decay rate.This
			file is an output of the code utilities/createinput.R
			[DEFAULT]  -D input/allGenesDecrEqualsSynrWithScaling_r_0.2_w_0.4_S.cer.mRNA.ini.abndc.syn.dec

	-C <FILE>    File containing the information about codon, tRNA, tRNA abundance and wobble.
			This file is an output of the code utilities/create.Scer.cod.anticod.numeric.pl
			[DEFAULT]  -C ../publicInput/S.cer.tRNA

	-J <FILE>    File containing the initial state of the system to begin simulations from.
			This file is an output of this simulation code containing '*_ribo_pos_*'
			[DEFAULT]  -C output/output_final_ribo_pos.out

	-s <value>    Random number seed. *MUST SETUP*
			[DEFAULT]  -s 1

	-O <prefix>    Specifies the prefix for the output files.
			[DEFAULT]  -O output


	-p[INTEGER]    Specify which output files to print

			-p1:    Generates a file of average elongation times
			    of all codons.

			-p2:    Generates a file of total average elongation
			    time of each gene.

			-p3:    Generates a file of average time between initiation
			    events on mRNAs of each gene.

			-p4:    Generates a file of average number of free ribosomes,
			    and free tRNAs of each type.

			-p5:    Generates a file of the final state of all mRNAs in a cell.
			    It contains the poistions of all bound ribosomes on mRNAs.

			-p6:    Generates 10 files:
			    Two files each containing the averaged number of mRNAs being synthsized/decayed
			    wihin a minute for each gene;
			    Two files each containing the mean/variance of total translating time among 
			    all the ribosomes that finsihed translating within a minute for each gene;
			    One file containing the averaged number of proteins being generated per minute for each gene;
			    One file containing the averaged number of mRNAs (being the result of synthesis and decay)
			    present in the cell during each minute for a particular gene.
			    One file containing the averaged number of mRNAs marked for decay present in the cell during
			    each minute for a particular gene.
			    One file containing the free ribosome number and free tRNA number per tRNA type at each minute mark
			    One file containing the number of elongatable ribosomes per codon type at each minute mark.
			    One file containing variance in the number of bound ribosomes among all mRNAs (whether marked 
			    or not for decay) for each gene at each minute mark

			-p7:    Generates one file containing the up to the first 10,000 mRNA life times for each gene
