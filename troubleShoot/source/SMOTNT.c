/// Code to simulate the transcription and translation process based on explicit diffusion and mass action properties of tRNAs, ribosome and mRNAs,
/// while the mRNA supply is in dynamic equilibrium

/*
To compile and run the code:

gcc source/SMOTNT.c -g -lm -lgsl -lgslcblas -mtune=generic -O3 -o source/SMOTNT

./source/SMOTNT -Tt 1500 -Tb 1000 -R 200000 -t 3300000 -N 4839 -dc 0.2 -r 0.4 -F ../publicInput/truncated_S.cer.genom -C ../publicInput/S.cer.tRNA -D input/allGenesDecrEqualsSynrWithScaling_dc_0.2_r_0.4_S.cer.mRNA.ini.abndc.syn.dec -s 16629 -O output/output.seed.16629_dc_0.2_r_0.4_allGenesDecrEqualsSynrScaling -p1 -p2 -p3 -p4 -p6 -p7 
 

 
*/
// Declaring Header Files
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <stdlib.h>

// Fixed parameters
#define MAX_GENES 5000                // Maximum number of genes supported
#define MAX_GENE_LEN_ALLW 5000        // Maximum allowed gene length in codons
#define MAX_TIME 2400000            // Maximum time of the simulation
#define char_len_tRNA 1.5e-8        // Characteristic length of tRNA
#define char_len_ribo 3e-8            // Characteristic length of ribosome
#define char_time_tRNA 5.719e-4        // Characteristic time of movement for tRNA (4.45e-7 * 1285.1)
#define char_time_ribo 5e-4            // Characteristic time of movement for ribosome


// Default global variables
int seed = 0;                        // Seed for RNG
int n_genes = 1;                    // Number of genes
int tot_ribo = 2e5;                    // Total ribosomes
int tot_tRNA = 3.3e6;                // Total tRNAs
int tot_mRNA = 0;                    // Total mRNA (at initialization)
int dbl_tot_mRNA=0;                 // Double totle mRNA
double tot_space = 4.2e-17;            // Total space within a cell
double avail_space_t = 1.24e7;        // Available space for tRNAs
double avail_space_r = 1.56e6;        // Available space for ribosomes
double tot_time = 1500;                // Total time for simulation
double thresh_time = 1000;            // Threshold time for analysis of e_times
double ratio_cotrans_decay = 0.0;           // Ratio of the 5'-3' co-translational mRNA decay pathway [0,1]
// 1 = exclusively co-translational/delayed decay pathway(5'-3'), 0= exclusively immidiate decay pathway(3'-5')
double ribo_prot_index = 0.0;           // The ribosome protection index [0,1]
// 0 = when the number of bound ribos doesn't at all affect the choice of an mRNA for decay (when the gene is already selected)

// Run options
int printOpt[11] = {0,0,0,0,0,0,0,0,0,0,0};
char out_prefix[150] = "output";    // Prefix for output file names
char out_file[150];
char code_file[150] = "../publicInput/S.cer.tRNA";
char fasta_file[150] = "../publicInput/truncated_S.cer.genom";// original S.cer.genom also contains initiation   
char m_decsyn_file[150] = "input/allGenesDecrEqualsSynrWithScaling_r_0.2_w_0.4_S.cer.mRNA.ini.abndc.syn.dec";
char state_file[150] = "";



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Create ribosome and mRNA based structures
typedef struct
{    int mRNA;                    // Bound to which mRNA
    int pos;                        // Position on mRNA
    double t_trans_ini;            // Time of translation initiation
    double t_elong_ini;            // Time of arrival at current codon
    int elng_cod_list;            // The id of the list of elongatable codons
    int elng_pos_list;            // Position in the list of elongatable codons
    int inhbtr_bound;            // Is the ribosome bound with CHX or Harr?
} ribosome;

typedef struct
{ int gene;                    // Gene id
    double last_ini;            // Time of last initiation event   
    int mdec_delayed;               /// whether the mRNA will be decayed immidiately or delayed in the elongation
    int num_boundr;                  /// number of bound ribosomes
    double msp_mdec_prob;              /// mRNA specific decay probability, depending on the # of bound ribos
    double t_mbirth;                  /// the time when an mRNA is newly synthesized (after burn-in)
} transcript;

typedef struct
{   int seq[MAX_GENE_LEN_ALLW];    // Codon sequence of the gene
    int len;                    // Length of the gene
    int exp;                    /// Gene expression level in real time (does not include mRNAs marked for decay)
    double ini_prob;            // Initiation probability of the mRNAT
    double mdec_r;       // mRNA decay rate
    double msyn_r;       // mRNA synthesis rate
    int ini_exp;          // original expression level
} gene;

typedef struct
{   char codon[4];                // Codon
    int tid;                    // tRNA id
    int gcn;                    // tRNA gene copy number
    double wobble;                // Wobble parmeter
} trna;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Reading the processed sequence file
int Read_FASTA_File(char *filename, gene *Gene)
{   FILE *fh;
    int c1=0, c2=1;
    char curr_char;
    
    fh=fopen(filename, "r");
    
    if(!fh)                    // Check if file exists
    {    printf("\nModified FASTA/Sequence File Doesn't Exist\n");
        fflush(stdout);
        Help_out();
        exit(1);
    }
    
    fscanf(fh,"%d",&Gene[c1].seq[0]);
    
    do
    {   c2=1;
        do
        {   fscanf(fh,"%d",&Gene[c1].seq[c2]);
            c2++;
            curr_char = fgetc(fh);
        }while(curr_char != '\n');
        Gene[c1].len = c2;
        
        c1++;
    }while(fscanf(fh,"%d",&Gene[c1].seq[0]) ==1);
    
}

// Read in tRNA file
int Read_tRNA_File(char *filename, trna *cTRNA)
{    FILE *fh;
    int c1=0,c2=0;
    char curr_char;
    
    fh=fopen(filename, "r");
    
    if(!fh)                    // Check if file exists
    {    printf("\ntRNA File Doesn't Exist\n");
        fflush(stdout);
        Help_out();
        exit(1);
    }
    c1 = 0;
    c2 = 0;
    while(fscanf(fh,"%s",&cTRNA[c1].codon)==1)
    {    fscanf(fh,"%d%d%lf",&cTRNA[c1].tid,&cTRNA[c1].gcn,&cTRNA[c1].wobble);
        c1++;
    }
}


// Note: Read in mRNA decay and synthesis rates file, sampled rates generated during simulation
int Read_DecSynRate_File(char *filename, gene *Gene)
{   FILE *fh;
    int c1=0;
    char curr_char;
    
    fh=fopen(filename, "r");
    
    if(!fh)                 // Check if file exists
    {   printf("\nmRNA Decay and Synthesis Rates File Doesn't Exist\n");
        fflush(stdout);
        Help_out();
        exit(1);
    }
    
    fscanf(fh,"%lf",&Gene[c1].ini_prob);
    
    do
    {   fscanf(fh,"%d %lf %lf",&Gene[c1].exp, &Gene[c1].msyn_r, &Gene[c1].mdec_r);  /// msyn_r here is the mean synthesis rates
        c1++;
    }while(fscanf(fh,"%lf",&Gene[c1].ini_prob)==1);
    
}


// Reading the state of the system
int Read_STATE_File(char *filename, int **R_grid)
{    FILE *fh;
    int c1=0,c2=0;
    char curr_char;
    
    fh=fopen(filename, "r");
    
    if(!fh)                    // Check if file exists
    {    printf("\nState File Doesn't Exist\n");
        fflush(stdout);
        Help_out();
        exit(1);
    }
    
    c1 = 0;
    c2 = -1;
    do
    {    c2++;
        do
        {    fscanf(fh,"%d",&R_grid[c1][c2]);
            c2++;
            curr_char = fgetc(fh);
        }while(curr_char != '\n');
        c1++;
        c2 = 0;
    }while(fscanf(fh,"%d",&R_grid[c1][c2]) ==1);
}




// Help output
int Help_out()
{    printf("\nUsage:\n");
    printf("\t./bin/SMoPT [options]\n\n");
    printf("Options:\n");
    printf("\n\t-V <value>    Volume of the cell in m^3/s. The minimum volume of the cell is set to\n");
    printf("\t\t\tcontain at least 1000 ribosomes and tRNAs.\n");
    printf("\t\t\t[DEFAULT]  -V 4.2E-17 (volume of yeast cell)\n");
    printf("\n");
    printf("\t-T[CHAR]    Specify times for various events\n");
    printf("\n");
    printf("\t\t\t-Tt:    Total simulation time in seconds.\n");
    printf("\t\t\t    [DEFAULT]  -Tt 1500\n");
    printf("\n");
    printf("\t\t\t-Tb:    Burn-in/threshold time. Time spent by the cell to reach equilibrium.\n");
    printf("\t\t\t    Only calculations after this time will be included in the analyses.\n");
    printf("\t\t\t    [DEFAULT]  -Tb 1000\n");
    printf("\n");
    printf("\t-R <value>    Total number of ribosomes in the cell.\n");
    printf("\t\t\t[DEFAULT]  -R 200000\n");
    printf("\n");
    printf("\t-t <value>    Total number of tRNAs in the cell.\n");
    printf("\t\t\t[DEFAULT]  -t 3300000\n");
    printf("\n");
    printf("\t-N <value>    Total number of genes. This needs to be specified by the user.\n");
    printf("\t\t\t[DEFAULT]  -N 1\n");
    printf("\n");
    printf("\t-dc <value>    The ratio of 5'-3' co-translational mRNA decay for all mRNAs [0,1].\n");
    printf("\t\t\t[DEFAULT]  -dc 0.0\n");
    printf("\n");
    printf("\t-r <value>    The ribosome protection index [0,1]\n");
    printf("\t\t\t[DEFAULT]  -r 0\n");
    printf("\n");
    printf("\t-F <FILE>    File containing processed fasta file into a numeric sequence.\n");
    printf("\t\t\tThis file is an output of the code utilities/convert.fasta.to.genom.pl\n");
    printf("\t\t\tIt contains the information regarding initiation probability, mRNA\n");
    printf("\t\t\tabundance and codon sequence of each gene.\n");
    printf("\t\t\t[DEFAULT]  -F ../publicInput/truncated_S.cer.genom\n");
    printf("\n");
    printf("\t-D <FILE>    File containing the information about gene specific initiation\n");
    printf("\t\t\tprobability, mRNA level, gene specific synthesis rate and decay rate.This\n");
    printf("\t\t\tfile is an output of the code utilities/createinput.R\n");
    printf("\t\t\t[DEFAULT]  -D input/allGenesDecrEqualsSynrWithScaling_r_0.2_w_0.4_S.cer.mRNA.ini.abndc.syn.dec\n");
    printf("\n");
    printf("\t-C <FILE>    File containing the information about codon, tRNA, tRNA abundance and wobble.\n");
    printf("\t\t\tThis file is an output of the code utilities/create.Scer.cod.anticod.numeric.pl\n");
    printf("\t\t\t[DEFAULT]  -C ../publicInput/S.cer.tRNA\n");
    printf("\n");
    printf("\t-J <FILE>    File containing the initial state of the system to begin simulations from.\n");
    printf("\t\t\tThis file is an output of this simulation code containing '*_ribo_pos_*'\n");
    printf("\t\t\t[DEFAULT]  -C output/output_final_ribo_pos.out\n");
    printf("\n");
    printf("\t-s <value>    Random number seed. *MUST SETUP*\n");
    printf("\t\t\t[DEFAULT]  -s 1\n");
    printf("\n");
    printf("\t-O <prefix>    Specifies the prefix for the output files.\n");
    printf("\t\t\t[DEFAULT]  -O output\n");
    printf("\n");
    printf("\n");
    printf("\t-p[INTEGER]    Specify which output files to print\n");
    printf("\n");
    printf("\t\t\t-p1:    Generates a file of average elongation times\n");
    printf("\t\t\t    of all codons.\n");
    printf("\n");
    printf("\t\t\t-p2:    Generates a file of total average elongation\n");
    printf("\t\t\t    time of each gene.\n");
    printf("\n");
    printf("\t\t\t-p3:    Generates a file of average time between initiation\n");
    printf("\t\t\t    events on mRNAs of each gene.\n");
    printf("\n");
    printf("\t\t\t-p4:    Generates a file of average number of free ribosomes,\n");
    printf("\t\t\t    and free tRNAs of each type.\n");
    printf("\n");
    printf("\t\t\t-p5:    Generates a file of the final state of all mRNAs in a cell.\n");
    printf("\t\t\t    It contains the poistions of all bound ribosomes on mRNAs.\n");
    printf("\n");
    printf("\t\t\t-p6:    Generates 10 files:\n");
    printf("\t\t\t    Two files each containing the averaged number of mRNAs being synthsized/decayed\n");
    printf("\t\t\t    wihin a minute for each gene;\n");
    printf("\t\t\t    Two files each containing the mean/variance of total translating time among \n");
    printf("\t\t\t    all the ribosomes that finsihed translating within a minute for each gene;\n");
    printf("\t\t\t    One file containing the averaged number of proteins being generated per minute for each gene;\n");
    printf("\t\t\t    One file containing the averaged number of mRNAs (being the result of synthesis and decay)\n");
    printf("\t\t\t    present in the cell during each minute for a particular gene.\n");
    printf("\t\t\t    One file containing the averaged number of mRNAs marked for decay present in the cell during\n");
    printf("\t\t\t    each minute for a particular gene.\n");
    printf("\t\t\t    One file containing the free ribosome number and free tRNA number per tRNA type at each minute mark\n");
    printf("\t\t\t    One file containing the number of elongatable ribosomes per codon type at each minute mark.\n");
    printf("\t\t\t    One file containing variance in the number of bound ribosomes among all mRNAs (whether marked \n");
    printf("\t\t\t    or not for decay) for each gene at each minute mark\n");
    printf("\n");
    printf("\t\t\t-p7:    Generates one file containing the up to the first 10,000 mRNA life times for each gene\n");
    printf("\n");
}

// Read in commandline arguments
void Read_Commandline_Args(int argc, char *argv[])
{  int i;
    
    for(i=1;i<argc;i++)
    {    if(argv[i][0] == '-')
    {    switch(argv[i][1])
    {    case 'h':
    case '-':
        fflush(stdout);
        Help_out();
        exit(1);
        break;
    case 'V':
        if((argv[i][2] != '\0') || (i==argc-1))
        {    printf("\nTotal space not specified or Incorrect usage\n");
            fflush(stdout);
            Help_out();
            exit(1);
        }
        else
        {    tot_space = atof(argv[++i]);
            
            avail_space_t = (double)floor(tot_space/pow(char_len_tRNA,3));
            avail_space_r = (double)floor(tot_space/pow(char_len_ribo,3));
            
            if(avail_space_r<1e3 || avail_space_t < 1e3)
            {    printf("\nAvailable cytoplasmic space is too small\n\n");
                fflush(stdout);
                Help_out();
                exit(1);
            }
            break;
        }
    case 's':
        if((argv[i][2] != '\0') || (i==argc-1))
        {    printf("\nSeed for RNG not specified or Incorrect usage\n");
            fflush(stdout);
            Help_out();
            exit(1);
        }
        else
        {    seed = atoi(argv[++i]);
            break;
        }
    case 't':
        if((argv[i][2] != '\0') || (i==argc-1))
        {   printf("\nTotal # of tRNAs not specified or Incorrect usage\n");
            fflush(stdout);
            Help_out();
            exit(1);
        }
        else
        {    tot_tRNA = atoi(argv[++i]);
            break;
        }
    case 'R':
        if((argv[i][2] != '\0') || (i==argc-1))
        {    printf("\nTotal # of ribosomes not specified or Incorrect usage\n");
            fflush(stdout);
            Help_out();
            exit(1);
        }
        else
        {    tot_ribo = atoi(argv[++i]);
            break;
        }
    case 'N':
        if((argv[i][2] != '\0') || (i==argc-1))
        {    printf("\nTotal # of genes not specified or Incorrect usage\n");
            fflush(stdout);
            Help_out();
            exit(1);
        }
        else
        {    n_genes = atoi(argv[++i]);
            if(n_genes>MAX_GENES)
            {    printf("\nTotal # of genes for simulation exceeds maximum genes = %d\n", MAX_GENES);
                fflush(stdout);
                Help_out();
                exit(1);
            }
            break;
        }
        
        
    case 'r':
        if((argv[i][2] != '\0') || (i==argc-1))
        {    printf("\nThe ribosome protection index not SPecified or incorrect usage\n");
            fflush(stdout);
            Help_out();
            exit(1);
        }
        else
        {    ribo_prot_index = atof(argv[++i]);
            if(ribo_prot_index >1 || ribo_prot_index <0)
            {    printf("\nThe ribosome protection index is outside of [0,1] range\n");
                fflush(stdout);
                Help_out();
                exit(1);
            }
            break;
        }
        
    case 'O':
        if((argv[i][2] != '\0') || (i==argc-1))
        {    printf("\nOutput prefix not specified or Incorrect usage\n");
            fflush(stdout);
            Help_out();
            exit(1);
        }
        else
        {    strcpy(out_prefix,argv[++i]);
            break;
        }
    case 'F':
        if((argv[i][2] != '\0') || (i==argc-1))
        {    printf("\nFasta file not specified or Incorrect usage\n");
            fflush(stdout);
            Help_out();
            exit(1);
        }
        else
        {    strcpy(fasta_file,argv[++i]);
            break;
        }
    case 'D':
        if((argv[i][2] != '\0') || (i==argc-1))
        {    printf("\nmRNA decay and synthesis rate file not specified or Incorrect usage\n");
            fflush(stdout);
            Help_out();
            exit(1);
        }
        else
        {    strcpy(m_decsyn_file,argv[++i]);
            break;
        }
    case 'J':
        if((argv[i][2] != '\0') || (i==argc-1))
        {    printf("\nState file not specified or Incorrect usage\n");
            fflush(stdout);
            Help_out();
            exit(1);
        }
        else
        {    strcpy(state_file,argv[++i]);
            break;
        }
    case 'C':
        if((argv[i][2] != '\0') || (i==argc-1))
        {    printf("\nS.cer Code file not specified or Incorrect usage\n");
            fflush(stdout);
            Help_out();
            exit(1);
        }
        else
        {    strcpy(code_file,argv[++i]);
            break;
        }
    //case 'P':
    case 'p':
        switch(argv[i][2])
        {   
        case '1':
            printOpt[0]=1;        // Elongation times of all codons
            break;
        case '2':
            printOpt[1]=1;        // Average total elongation time of all genes
            break;
        case '3':
            printOpt[2]=1;        // Average time between initiation events of all genes
            break;
        case '4':
            printOpt[3]=1;        // Average number of free ribosomes and tRNAs of each type
            break;
        case '5':
            printOpt[4]=1;        // Final state of the system - positions on mRNAs bound by ribosomes
            break;
        case '6':
            printOpt[5]=1;        // Generate ten files that are new to SMoTnT
            break;
        case '7':
            printOpt[6]=1;        // New to SMoTnT: a file containing up to 10,000 of mRNA lifetimes for all genes
            break;
        default:
            printf("\nInvalid print options\n");
        fflush(stdout);
        Help_out();
        exit(1);
        break;
        }
    case 'T':
        switch(argv[i][2])
        {    case 't':
            tot_time = atof(argv[++i]);
            if(tot_time>MAX_TIME)
            {    printf("\nTotal time for simulation exceeds maximum allowed time = %d\n",MAX_TIME);
                fflush(stdout);
                Help_out();
                exit(1);
            }
            break;
        case 'b':
            thresh_time = atof(argv[++i]);
            if(thresh_time<0)
            {    printf("\nTime threshold %g should be > 0 and < Total time %g\n",thresh_time,tot_time);
                fflush(stdout);
                Help_out();
                exit(1);
            }
            break;
        }
        
    case 'd':
        switch(argv[i][2])
        {   case 'c':
            if((argv[i][3] != '\0') || (i==argc-1))
            {    printf("\nThe ratio of 5'3' co-translational mRNA decay not specified or Incorrect usage\n");
                fflush(stdout);
                Help_out();
                exit(1);
            }
            else
            {    ratio_cotrans_decay = atof(argv[++i]);
                if(ratio_cotrans_decay>1 || ratio_cotrans_decay <0)
                {    printf("\nThe ratio of 5'-3' co-translational mRNA decay is outside of [0,1] range\n");
                    fflush(stdout);
                    Help_out();
                    exit(1);
                }
                break;
            }
        }
    }
    }
    }
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Main Function
int main(int argc, char *argv[])
{    int c1, c2, c3, c4;
    FILE *f0, *f1, *f2, *f3, *f4, *f5, *f6, *f7, *f8, *f9, *f10, *f11, *f12;
    
    // Read in arguments from the commandline
    Read_Commandline_Args(argc, argv);
    
    /// 11 files that are generated in real time during the simulation, that are new to SMoTnT
    // Note: in the comment below, "during each minute" or "in that minute" means all events happened within that minute was accounted; 
    // "at each minute mark" means only the snapshot of the last time step of that minute is accounted
    if(printOpt[5]==1)
    {   ///print out N lines, each line represents a minute(after burn-in) and has n elements(n=n_genes), each number means how many mRNAs has been decayed/synthesized in that minute
        strcpy(out_file,out_prefix);
        f0 = fopen(strcat(out_file,".gene_sp_decCount_perMin.out"),"w");
        
        strcpy(out_file,out_prefix);
        f1 = fopen(strcat(out_file,".gene_sp_synCount_perMin.out"),"w");
        
        /// print out N (=simuTime in min) lines, each line represents a minute(after burn-in) and has n elements(n=n_genes),each number means the average elongation time it takes to translate that gene in each minute
        /// (translation efficiency) how long does each ribo spend on translating an mRNA for a particular gene
        strcpy(out_file,out_prefix);
        f2 = fopen(strcat(out_file,".gene_sp_meanElngTime_perMin.out"),"w");
        
        strcpy(out_file,out_prefix);
        f3 = fopen(strcat(out_file,".gene_sp_varianceElngTime_perMin.out"),"w");
        
        /// （translation rate）how many ribosomes hop off (= how many proteins have been translated) per minute for a particular gene
        strcpy(out_file,out_prefix);
        f4 = fopen(strcat(out_file,".gene_sp_countElng_perMin.out"),"w");
        
        /// the average number of mRNAs (being the result of synthesis and decay) within each minute for a particular gene
        strcpy(out_file,out_prefix);
        f5 = fopen(strcat(out_file,".gene_sp_mRNAcount_perMin.out"),"w");
        
        /// the average number of mRNAs marked for decay within each minute for a particular gene
        strcpy(out_file,out_prefix);
        f6 = fopen(strcat(out_file,".gene_sp_mRNAmarkedDecay_perMin.out"),"w");
        
        /// the free ribosome number and free tRNA number per tRNA type at each minute mark 
        /// the output file contains (N=simuTime in min) lines, each line= time mark,free ribo, free tRNA1, free tRNA2...Free tRNA41
        strcpy(out_file,out_prefix);
        f7 = fopen(strcat(out_file,".free_ribo_tRNA_perMin.out"),"w");
        
        /// the number of elongatable ribosomes per codon type at each minute mark 
        /// the output file contains (N=simuTime in min) lines, each line= time mark, codon 1, codon2 ...codon 61
        strcpy(out_file,out_prefix);
        f8 = fopen(strcat(out_file,".codon_sp_elngRibo_perMin.out"),"w");
        
        /// the number of total bound ribosomes for each gene (on all the mRNAs of this gene, whether marked or not for decay) at each minute mark
        strcpy(out_file,out_prefix);
        f9 = fopen(strcat(out_file,".gene_sp_totBoundRibo_perMin.out"),"w");
        
        /// the variance in the number of bound ribosomes for all mRNAs (whether marked or not for decay) for each gene at each minute mark
        strcpy(out_file,out_prefix);
        f10 = fopen(strcat(out_file,".gene_sp_varianceBoundRibo_perMin.out"),"w");
        
    }
    
    
    // Random number generation setup
    gsl_rng * r;
    gsl_rng_env_setup();
    
    r = gsl_rng_alloc (gsl_rng_mt19937);
    gsl_rng_set (r, (unsigned long) seed);
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // User specified parameters for quick test
    gene *Gene = (gene *)malloc(sizeof(gene) * n_genes);
    if(Gene == NULL)
    {    printf("Too many genes\nOut of memory\n");fflush(stdout);
    exit(1);
    }
    
    // Initialize the various structures
    ribosome *Ribo = (ribosome *)malloc(sizeof(ribosome) * tot_ribo);
    if(Ribo == NULL)
    {    printf("Too many ribosomes\nOut of memory\n");fflush(stdout);
    exit(1);
    }
    
    trna *cTRNA = (trna *)malloc(sizeof(trna) * 61);    
    if(cTRNA == NULL)
    {    printf("Too many tRNAs\nOut of memory\n");fflush(stdout);
    exit(1);
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // Read in the numeric seq
    Read_FASTA_File(fasta_file, Gene);
    
    // Read in the tRNA code file
    Read_tRNA_File(code_file, cTRNA);
    
    // Read in the mRNA decay and synthesis rate file
    Read_DecSynRate_File(m_decsyn_file, Gene);
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    for(c1=0;c1<n_genes;c1++)
    {   Gene[c1].ini_exp = Gene[c1].exp;
        tot_mRNA += Gene[c1].exp;   // Gene[c1].exp does not include the mRNAmarkedForDecay
    }
    
    dbl_tot_mRNA = 2*tot_mRNA;
    
    transcript *mRNA = (transcript *)malloc(sizeof(transcript) * dbl_tot_mRNA);  // Doubled the memory for mRNAs
    if(mRNA == NULL)
    {    printf("Too many mRNAs\nOut of memory\n");fflush(stdout);
    exit(1);
    }
    
    c3=0;
    for(c1=0;c1<n_genes;c1++)                    // Fill the mRNA struct with
    {    for(c2=0;c2<Gene[c1].exp;c2++)            // length and ini_rate info.
    {   mRNA[c3].gene = c1;
        mRNA[c3].last_ini = 0.0;
        mRNA[c3].mdec_delayed = 0;
        mRNA[c3].num_boundr = 0;
        mRNA[c3].msp_mdec_prob = 0.0;
        mRNA[c3].t_mbirth = 0.0;
        c3++;
    }
    }
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    //     Variables to track the translation process
    int Tf[61];                                    // Number of free tRNAs of each type
    int Rf = tot_ribo;                            // Total number of free ribosomes (at start all ribosomes are free)
    int Rb[61];                                    // Number of ribosomes bound to each codon
    int Mf[n_genes];                            // Number of initiable mRNAs of each gene (first 10 codons unbound by ribosomes)
    int M_map[n_genes];                         // Number of all mRNAs of each gene in the gene_mRNA_map
    int M_dec = 0;                              // Number of decayed mRNAs in real time
    int M_decSlt[n_genes];                      // Number of mRNAs for each gene that can be selected for decay
    int count_mlifetime[n_genes];                      // Number of recorded life time for each gene
    int tot_Mf;                                    // Total number of initiable mRNAs
    int x;
    int n_Rb_e[61];                                // Number of elongable bound ribosomes to each codon
    int r_id;
    int m_id;
    int next_m_id;
    int c_id;
    int c2_id;
    int g_id;
    int tot_gcn=0;
    int obs_max_len = 0;                        // Observed max gene length among all genes
    int obs_max_exp = 0;                        // Observed max gene expression among all genes
    int triple_max_exp = 0;
    int next_avail_ribo = 0;
    int termtn_now=0;
    int n_e_times[61];                            // Number of times a codon type is elongated
    
    double e_times[61];                          /// Total time length of a codon type being elongated
    double current_etimes = 0.0;
    int del_n_etimes[61];                        /// Number of times a codon type is elongated for all deleted mRNAs
    double del_etimes[61];
    double scld_Mf[n_genes];
    double tot_scld_Mf = 0.0;
    double r_ini;
    double r_elng[61];
    double avg_tRNA_abndc[61];                    // Average number of free tRNAs of each type (averaged by time)
    double avg_Rf = {0.0};                        // Average number of free ribosomes (averaged by time)
    
    
    // Begin simulation of the translation process
    double t = 0.0;
    double prob_ini;
    double tmp_elng_prob;
    double prob_e[61];
    double prob_g;
    
    double tot_rate = 0.0;
    double inv_rate;
    double prob_mdec = 0.0;
    double prob_msyn = 0.0;
    double tot_mdec_r = 0.0;
    double tot_msyn_r = 0.0;
    int decCount[n_genes];
    int synCount[n_genes];
    double meanElngTime[n_genes];                  // Average elongation time among all the ribosomes that finished translating a gene within a minute
    double varElngTime[n_genes];                   // Corresponding var for above meanElngTime
    double sdElngTime[n_genes];  
    int ctElngTime[n_genes];                       // Count of how many times a gene has been translated during each time slot
    int mRNAcount[n_genes];                       // Total mRNA count for each gene, including the mRNAs that are marked for decay
    int M_decTrans[n_genes];                     // Number of mRNAs for each gene that have been selected for decay but are still translating (delayed decay) during each time slot
    int M_decTransPlus[n_genes];                     
    int M_decTransMinus[n_genes];  
    double elngTime = 0.0;                      // Used to update the running mean and variance for elongation time per ribosome
    double old_mean = 0.0;    
    double meanBoundRibo[n_genes];              // Average number of bound tibos among all the mRNAs of a particular gene at the last time step of each minute
    double intermedBoundRibo[n_genes];          // Intermediate for calculating running standard deviation
    double varBoundRibo[n_genes];               // Variance in the number of bound ribos among all the mRNAs of a particular gene at the last time step of each minute
    double m_lifeTime = 0.0;                    // The time between when a new mRNA is synthesized (after burn-in) and when it is completely decayed (note: not when it's marked for decay)
    int tot_decCount[n_genes];                  // The total number of decayed mRNAs for each gene
    int tot_bound_ribo[n_genes];                // The number of total bound ribosomes for all the mRNAs of a gene

    int n_trans[n_genes];                // Number of translation events per gene
    int n_ini[n_genes];                 // Number of initiation events per gene
    double avg_time_to_ini[n_genes];    // Gene-specific mean time between initn events
    double avg_time_to_trans[n_genes];  // Gene-specific mean time of total translation time (time taken from start to stop codon)
    double var_time_to_ini[n_genes];    // Corresponding variance for the time between initn events
    double var_time_to_trans[n_genes];  // Corresponding variance for the time of total translation time
    double sd_time_to_ini[n_genes];    // Corresponding standard deviation, need this to calculate the running variance
    double sd_time_to_trans[n_genes];  // Corresponding standard deviation
    double current_ini = 0.0;
    double current_trans = 0.0;
    
    double coin;
    int t_print = floor(thresh_time)+60;  // The first time the output is printed to file is already 1 min after the thresh_time
    
    // Note: when delete or synthesize new mRNA, Rb_e needs to be updated
    int **Rb_e;                                    // Bound ribosomes to each codon that can be elongated
    Rb_e = malloc(sizeof(int *) * 62);
    for(c1=0;c1<62;c1++)
    {    Rb_e[c1] = malloc(sizeof(int *) * (0.25*tot_ribo));
    }
    
    
    for(c1=0;c1<n_genes;c1++)
    { if(Gene[c1].exp > obs_max_exp)
        {    obs_max_exp = Gene[c1].exp;
        }
        
        if(Gene[c1].len > obs_max_len)
        {    obs_max_len = Gene[c1].len;
        }
        
        Mf[c1] = Gene[c1].exp;                              // Mf[n_genes] --Number of initiable mRNAs of each gene
        M_map[c1] = Gene[c1].exp;                           // Number of all present mRNAs of each gene
        M_decSlt[c1] = Gene[c1].exp;                        // Number of all mRNAs that can be selected for mRNA decay of each gene
        scld_Mf[c1] = Mf[c1]*Gene[c1].ini_prob;             // Indicates the total capability of a gene being initiated
        tot_scld_Mf += scld_Mf[c1];                         // Scale gene_exp with ini_prob for a particular gene
        count_mlifetime[c1] = 0;
    }
    
    triple_max_exp = 3*obs_max_exp;
    
    int **free_mRNA;                     // List to figure out which mRNAs can be initiated based on no bound ribosomes from pos=0->pos=10
    int **gene_mRNA_map;
    int **mdec_select_map;               // List to keep track of which mRNAs have never been (and therefore can be) selected for decay (whether it has bound ribos or not)
    double **gene_sp_mlifetime;          // List to keep track of the life time of first 10,000 decay events for each gene
    
    free_mRNA = malloc(sizeof(int *) * n_genes);
    gene_mRNA_map = malloc(sizeof(int *) * n_genes);
    mdec_select_map = malloc(sizeof(int *) * n_genes);
    gene_sp_mlifetime = malloc(sizeof(double *) *n_genes);
    
    int expmax;
    for(c1=0;c1<n_genes;c1++)
    {   
        
        if(Gene[c1].exp<100)
        {
            expmax = 1000;   // *3 this number here can be changed to larger when the model runs for environment stimulation experiments (e.g. transcriptional bursting)
        }
        else
        {
            expmax = Gene[c1].exp*5;   /// *3 this number here can be changed to larger when the model runs for environment stimulation experiments 
        }
        free_mRNA[c1] = malloc(sizeof(int *) * expmax);    
        gene_mRNA_map[c1] = malloc(sizeof(int *) * expmax); 
        mdec_select_map[c1] = malloc(sizeof(int *) * expmax);
        
        gene_sp_mlifetime[c1] = malloc(sizeof(double *) * 10000); 
    }
    
    // Initialize free mRNAs and gene_mRNA_map
    c3=0;
    for(c1=0;c1<n_genes;c1++)
    {   for(c2=0;c2<Gene[c1].exp;c2++)
    {   free_mRNA[c1][c2] = c3;                            // At the begininning all mRNAs can be initiated (store their ids here)
        gene_mRNA_map[c1][c2] = c3;                        // Updated once an mRNA has actually been decayed (whether the decay was delayed or not)
        mdec_select_map[c1][c2] = c3;                      // Updated once the decay of an mRNA has been delayed (whether the mRNA has already been 
        // Selected for decay or not, this structure stores the ones that can still be selected for decay)
        c3++;
    }
    }
    
    next_m_id = free_mRNA[n_genes-1][Gene[n_genes-1].exp-1] +1 ;  // A new m_id that can be assigned to newly synthesized mRNA
    // The first next_m_id=60000, b/c m_id starts from 0,so if the total initial mRNA is 60000, the last m_id=59999.
    
    
    // Initialize the output variables
    for(c1=0;c1<n_genes;c1++)
    {   decCount[c1] = 0;
        synCount[c1] = 0;
        meanElngTime[c1] = 0.0;
        varElngTime[c1] = 0.0;
        sdElngTime[c1] = 0.0;
        ctElngTime[c1] = 0;
        mRNAcount[c1] = Gene[c1].ini_exp;
        M_decTrans[c1] = 0;
        M_decTransPlus[c1] = 0;
        M_decTransMinus[c1] = 0;
        meanBoundRibo[c1] = 0.0;
        intermedBoundRibo[c1] = 0.0;
        varBoundRibo[c1] = 0.0;
        
        
        n_trans[c1] = 0;
        n_ini[c1] = 0;
        avg_time_to_ini[c1] = 0.0;
        avg_time_to_trans[c1] = 0.0;
        var_time_to_ini[c1] = 0.0;
        var_time_to_trans[c1] = 0.0;
        sd_time_to_ini[c1] = 0.0;
        sd_time_to_trans[c1] = 0.0;
        
    }
    
    // Initialize free tRNAs and elongation time vectors
    for(c1=0;c1<61;c1++)
    {   if(cTRNA[c1].wobble==1.0)
    {    tot_gcn += cTRNA[c1].gcn;
    }
    e_times[c1] = 0.0;
    n_e_times[c1] = 0;
    
    Tf[c1]=0;
    avg_tRNA_abndc[c1]=0.0;
    }
    c2 = 0;
    
    for(c1=0;c1<61;c1++)
    {    if(cTRNA[c1].wobble==1.0)
    {    Tf[cTRNA[c1].tid] = floor((double)cTRNA[c1].gcn*tot_tRNA/tot_gcn);    // Number of tRNAs for tRNA-type cTRNA[c1].tid
    }
    n_Rb_e[c1] = 0;                                                            // Initialize number of elongatable bound ribosomes to each codon
    cTRNA[c1].wobble = cTRNA[c1].wobble/(char_time_tRNA*avail_space_t);        // For faster computation
    }
    
    int **R_grid;                            // The state of the system with respect to mRNAs and bound ribosomes

    R_grid = malloc(sizeof(int *) * dbl_tot_mRNA);   // Initialize R_grid as a larger grid - 2X larger than initial tot_mRNA
    for(c1=0;c1<dbl_tot_mRNA;c1++)
    {   R_grid[c1] = malloc(sizeof(int *) * obs_max_len);
        for(c2=0;c2<Gene[mRNA[c1].gene].len;c2++)
        {    R_grid[c1][c2] = tot_ribo;
        }
    }
    
    for (c1=0;c1<n_genes;c1++)
    {   tot_mdec_r += Gene[c1].mdec_r*Gene[c1].exp;               // Assuming all mRNAs (free or bound) can be decayed
        tot_msyn_r += Gene[c1].msyn_r;
    }
    
    int ii=0;
    
    printf("\n\n-Tt %lf -Tb %lf -R %d -t %d -N %d -dc %lf -r %lf -F %s -C %s -D %s -s %d -O %s\n",tot_time,thresh_time,tot_ribo,tot_tRNA,n_genes,ratio_cotrans_decay,ribo_prot_index,fasta_file,code_file,m_decsyn_file,seed,out_prefix);
    
    
    /////////////////////////////////////////////////
    // Begin the actual continuous time simulation
    /////////////////////////////////////////////////
    while(t<tot_time)                                                                // Till current time is less than max time
    {
        r_ini = tot_scld_Mf*Rf/(char_time_ribo*avail_space_r);                        // Initiation rate        
        tot_rate = r_ini;
        
        if(t>thresh_time)
        {
            tot_rate += tot_mdec_r;
            tot_rate += tot_msyn_r;
        }
        
        for(c1=0;c1<61;c1++)
        {    r_elng[c1] = Tf[cTRNA[c1].tid]*cTRNA[c1].wobble;        // Elongation rate of codon c1; 
            r_elng[c1] *= (double)n_Rb_e[c1];                        // n_Rb_e = number of elongatable bound ribosomes to each codon
            // The more ribosome bound to a particular codon type, the higher chance this codon will gets elongated
            tot_rate += r_elng[c1];
        }
        
        if(tot_rate>0)
        {    inv_rate = 1/tot_rate;
        }
        else
        {    printf("\nNo further events to process.\nSimulation stopped at time %g\n\n",t);fflush(stdout);
        t=tot_time;
        continue;
        }
        
        // Increment time
        if(t>thresh_time)
        {   if(printOpt[3]==1)                                                        // Tracking free ribosomes and tRNAs
        {   for(c1=0;c1<61;c1++)
        {    avg_tRNA_abndc[c1] += (double)Tf[c1]*inv_rate;
        }
        avg_Rf += (double)Rf*inv_rate; // Total of (Rf that is normalized by time)
        }
        }
        
        t+=inv_rate;
        
        // Calculate prob of events
        prob_ini = r_ini*inv_rate;                                                    // Probability of initiation in current t
        
        if(t>thresh_time)
        {
            prob_mdec = prob_ini + tot_mdec_r*inv_rate;
            prob_msyn = prob_mdec + tot_msyn_r*inv_rate;
        }
        
        time_t rawtime;
        struct tm * timeinfo;
        
        if (ii%1000000000==0)
        {   time ( &rawtime );
            timeinfo = localtime ( &rawtime );
            printf("simulation time=%lf, system time=%s\n",t,asctime (timeinfo));fflush(stdout);

        };
        
        coin = gsl_rng_uniform(r);                                                    // Pick a random uniform to pick an event
        // Translation initiation
        if(coin<prob_ini)
        {
            // Ribosomes are picked sequentially
            r_id = next_avail_ribo;
            
            
            Ribo[r_id].pos = 0;
            Ribo[r_id].t_trans_ini = t;
            Ribo[r_id].t_elong_ini = t;
            next_avail_ribo++;
            
            // Pick a random mRNA for initiation
            coin = gsl_rng_uniform(r);
            prob_g = 0.0;
            
            for(c1=0;c1<n_genes;c1++)
            {   prob_g += scld_Mf[c1]/tot_scld_Mf;                                    // Pick a gene randomly first as they may differ in ini_prob
                if(coin<prob_g)
                {    c2 = gsl_rng_uniform_int(r, (unsigned long)Mf[c1]);                // Once a gene is selected pick a random mRNA
                    /// This function returns a random integer from 0 to n-1 --> [0,n-1] , n=Mf[c1]
                    m_id = free_mRNA[c1][c2];      
                    g_id = mRNA[m_id].gene;
                    
                    Mf[c1]--;
                    scld_Mf[c1]-=Gene[c1].ini_prob;
                    tot_scld_Mf-=Gene[c1].ini_prob;
                    if(c2!=Mf[c1])
                    {    free_mRNA[c1][c2] = free_mRNA[c1][Mf[c1]];        // Accounting of free mRNAs for next round, swap the info with last mRNA for this gene
                    }
                    break;
                }
            }
            
            
            if(t>thresh_time && printOpt[2]==1)
            {   
                current_ini = t-mRNA[m_id].last_ini;        // Time between initn
                old_mean = avg_time_to_ini[g_id];
                n_ini[g_id] ++;                         // Store the total # of initn events for this gene
                
                // Update the mean
                avg_time_to_ini[g_id] = avg_time_to_ini[g_id] + (current_ini - avg_time_to_ini[g_id])/n_ini[g_id];
                
                // Update the variance
                if(n_ini[g_id] >= 2)
                {
                    sd_time_to_ini[g_id] = sd_time_to_ini[g_id] + ((current_ini - old_mean)*(current_ini - avg_time_to_ini[g_id]));
                    var_time_to_ini[g_id] = sd_time_to_ini[g_id]/(n_ini[g_id]-1);
                }
                
            }
            
            mRNA[m_id].last_ini = t;
            
            c_id = Gene[mRNA[m_id].gene].seq[0];                                    // Codon at position 1
            
            Ribo[r_id].mRNA = m_id;
            
            if(R_grid[m_id][10]==tot_ribo)                                        // Check if the current ribosome can be elongated
            {    Rb_e[c_id][n_Rb_e[c_id]] = r_id;                                // If no ribosome at pos+10 then it can be elongated
                Ribo[r_id].elng_cod_list = c_id;
                Ribo[r_id].elng_pos_list = n_Rb_e[c_id];
                n_Rb_e[c_id]++;
            }
            Ribo[r_id].inhbtr_bound = 0;
            
            R_grid[m_id][0] = r_id;                                                   // Update the ribosome grid upon initiation
            
            Rf--;                                                                    // Update number of free ribosomes

            mRNA[m_id].num_boundr ++;
   
        }
        
        // mRNA decay
        else if(coin<prob_mdec && t>thresh_time)
        {   // Pick a random gene to process degradation
            
            coin = gsl_rng_uniform(r);
            double prob_dec_g = 0.0;
            
            for(c1=0;c1<n_genes;c1++)
            {   prob_dec_g += (Gene[c1].mdec_r*Gene[c1].exp)/tot_mdec_r;   // Gene[c1].exp does not include the mRNAmarkedForDecay
                
                if(coin<prob_dec_g)                                                  // First pick a gene
                {   
                    coin = gsl_rng_uniform(r);   // Next pick which mRNA of the gene will be decayed, depending on the number of bound ribosomes on it

                    double prob_dec_m = 0.0;    
                    double gtot_mdec_r = 0.0;   

                    // 1st: calculate the gtot_mdec_r (total decay rate for all the mRNAs belonging to this gene)
                    for(c2=0;c2<M_decSlt[c1];c2++)
                    {
                        mRNA[mdec_select_map[c1][c2]].msp_mdec_prob = exp(-mRNA[mdec_select_map[c1][c2]].num_boundr*ribo_prot_index);
                        gtot_mdec_r += mRNA[mdec_select_map[c1][c2]].msp_mdec_prob;

                    }
                    
                    
                    // 2nd: calculate the probability of each mRNA of the selected gene for decay, based on their number of bound ribosomes
                    for(c2=0;c2<M_decSlt[c1];c2++)
                    {
                        prob_dec_m += (mRNA[mdec_select_map[c1][c2]].msp_mdec_prob/gtot_mdec_r);
                        if(coin < prob_dec_m)
                        {
                            m_id = mdec_select_map[c1][c2];
                            break;
                        }
                    }
                    

                    g_id = c1;
                    
                    coin = gsl_rng_uniform(r);  // Pick which way of decay will proceed: either the 5'-3' Xrn1 or the 3'-5' exosome pathway
                    if(coin < (1-ratio_cotrans_decay))  // This is when the exosome 3'-5' immediate decay is selected to proceed
                    {
                        // Below is to update the moments of mRNA life times for each gene, life time = death time(decayed, not just marked for decay)- birth time for any newly synthesized mRNAs after burn-in
                        if(t>=thresh_time && printOpt[6]==1 && mRNA[m_id].t_mbirth >= thresh_time)    // If this selected mRNA is not one of the original 60,000 mRNAs
                        {
                            tot_decCount[g_id] ++;
                            
                            if(count_mlifetime[g_id] < 10000)
                            {
                                m_lifeTime = t - mRNA[m_id].t_mbirth;
                                
                                gene_sp_mlifetime[g_id][count_mlifetime[g_id]] = m_lifeTime;
                                count_mlifetime[g_id]++;
                            }
                        }
                        
                        int c3 = 0;          // Pos on Gene[c1]
                        int c4 = 0;          // To check if there's any bound ribo on pos0-pos10
                        int c5 = 0;          // Count of positions in R_grid
                        int c6 = 0;         
                        
                        // 1st for loop- go through R_grid to find any bound ribo
                        // 2nd for loop- go through Rb_e to find any elongatable ribo
                        for (c3=0;c3<Gene[c1].len;c3++)  // Go through R_grid to see if any ribo is bound， if there is
                        {   if(R_grid[m_id][c3]!= tot_ribo)
                        {
                            Rf++;           // Free a ribosome
                            next_avail_ribo--;      // Add a ribosome that can be used for initiation
                            
                            r_id = R_grid[m_id][c3];
                            c_id = Gene[mRNA[m_id].gene].seq[c3];
                            
                            
                            // Release the tRNA bound at the earlier position, if the c_id is not the first codon
                            if(Ribo[r_id].pos>0)
                            {    Tf[cTRNA[Gene[mRNA[m_id].gene].seq[Ribo[r_id].pos-1]].tid]++;
                            }
                            
                            // Check if next_avail_ribo is stalled
                            int c_id_nextribo = Gene[mRNA[Ribo[next_avail_ribo].mRNA].gene].seq[Ribo[next_avail_ribo].pos];
                            int judge_nextribo = 0;
                            for (int i=0;i<n_Rb_e[c_id_nextribo];i++)
                            {   if(Rb_e[c_id_nextribo][i]==next_avail_ribo)
                            {   judge_nextribo = 1;  // next_ribo is NOT stalled
                            }
                            }
                            
                            
                            // Check if selected ribo is stalled or not
                            int judge=0; // judge=0--->current r_id isn't in the Rb_e[c_id] array; judge=1--->current r_id is in the Rb_e[c_id] array=is elongatable
                            for(int i=0;i<n_Rb_e[c_id];i++)
                            {   if(Rb_e[c_id][i]==r_id)
                            {   judge=1;
                                // When there is 10 bound ribos for codon X, then n_Rb_e=10, but the last index in Rb_e should be 9 because it starts from 0.
                                n_Rb_e[c_id]--;         // Update number of elongatable ribosomes for current codon
                                
                                // below is to update Ribo struct and Rb_e array
                                x = Ribo[r_id].elng_pos_list;
                                
                                
                                // Needs two steps to update the ribo that's gone off the decayed mRNA
                                // step 1: put the last r_id* in the array of Rb_e[c_id] into current block, and change the Ribo[r_id*].elng_pos_list to current position
                                if(x!=n_Rb_e[c_id])
                                {   Ribo[Rb_e[c_id][n_Rb_e[c_id]]].elng_pos_list = x;// Update info for the last r_id in the Rb_e list for this codon
                                    Rb_e[c_id][x] = Rb_e[c_id][n_Rb_e[c_id]];  // Put the last r_id in the list for this codon into current position in Rb_e[c_id]
                                }
                                
                                
                                
                                // Update ribo struct, copy the info for ribo with the largest r_id to the r_id that's just been released
                                // But the info for current Ribo position is lost
                                // Step 2: copy all the contents/features of Ribo[next_avail_ribo] onto Ribo[current selected r_id]
                                // Four conditions: 1) the selected_rid stalled while the next_avail_ribo is unstalled; 
                                                 // 2) selected_rid is unstalled while the next_avail_ribo is stalled;
                                                 // 3) both selected_rid and next_avail_ribo are unstalled;
                                                 // 4) both selected_rid and next_avail_ribo are stalled;
                                // When current r_id is elongatable (next_r_id can be either elongatable OR NOT)
                                if(r_id != next_avail_ribo)
                                {   if(judge_nextribo==1)// When next_r_id is also elongatable, otherwise no need to do anything with Rb_e
                                {
                                    Rb_e[Ribo[next_avail_ribo].elng_cod_list][Ribo[next_avail_ribo].elng_pos_list]=r_id;
                                }
                                
                                // When either the next_avail_ribo is elongatable OR not
                                R_grid[Ribo[next_avail_ribo].mRNA][Ribo[next_avail_ribo].pos] = r_id;
                                
                                Ribo[r_id].mRNA = Ribo[next_avail_ribo].mRNA;
                                Ribo[r_id].pos = Ribo[next_avail_ribo].pos;
                                Ribo[r_id].t_trans_ini = Ribo[next_avail_ribo].t_trans_ini;
                                Ribo[r_id].t_elong_ini = Ribo[next_avail_ribo].t_elong_ini;
                                Ribo[r_id].elng_cod_list = Ribo[next_avail_ribo].elng_cod_list;
                                Ribo[r_id].elng_pos_list = Ribo[next_avail_ribo].elng_pos_list;
                                Ribo[r_id].inhbtr_bound = Ribo[next_avail_ribo].inhbtr_bound;
                                }
                                
                                
                                // Update the elng rate
                                r_elng[c_id] -= Tf[cTRNA[c_id].tid]*cTRNA[c_id].wobble;        // Elongation rate of codon c1; cTRNA[c1].tid is the tid of corresponding tRNA
                                
                                break;  // Once the ribo is recycled, no longer need to loop through the Rb_e[current c_id]
                            }
                            }
                            
                            // When current r_id is stalled
                            if(judge==0)
                            {   
                                // Copy the contents of Ribo[next_avail_ribo] onto Ribo[current_selected_rid]
                                if(r_id != next_avail_ribo)
                                {   if(judge_nextribo==1) // Current r_id is stalled, next_avail_ribo is NOT stalled
                                {
                                    Rb_e[Ribo[next_avail_ribo].elng_cod_list][Ribo[next_avail_ribo].elng_pos_list]=r_id;
                                }
                                
                                // Current r_id is stalled, next_avail_ribo is either stalled OR NOT
                                R_grid[Ribo[next_avail_ribo].mRNA][Ribo[next_avail_ribo].pos] = r_id;
                                
                                Ribo[r_id].mRNA = Ribo[next_avail_ribo].mRNA;
                                Ribo[r_id].pos = Ribo[next_avail_ribo].pos;
                                Ribo[r_id].t_trans_ini = Ribo[next_avail_ribo].t_trans_ini;
                                Ribo[r_id].t_elong_ini = Ribo[next_avail_ribo].t_elong_ini;
                                Ribo[r_id].elng_cod_list = Ribo[next_avail_ribo].elng_cod_list;
                                Ribo[r_id].elng_pos_list = Ribo[next_avail_ribo].elng_pos_list;
                                Ribo[r_id].inhbtr_bound = Ribo[next_avail_ribo].inhbtr_bound;
                                }
                            }
                            
                        }
                        }
                        
                        // If this mRNA is in the free_mRNA list
                        for(c3=0;c3<Mf[c1];c3++)
                        {   if(free_mRNA[c1][c3]==m_id)
                        {
                            Mf[c1]--;
                            if(c3!=Mf[c1])
                            {   
                                free_mRNA[c1][c3] = free_mRNA[c1][Mf[c1]];
                            }        // Accounting of free mRNAs for next round
                            scld_Mf[c1]-=Gene[c1].ini_prob;
                            tot_scld_Mf-=Gene[c1].ini_prob;
                            
                        }
                        }
                        
                        M_decSlt[g_id] --;
                        if(c2 != M_decSlt[g_id])  
                        {   mdec_select_map[g_id][c2] = mdec_select_map[g_id][M_decSlt[g_id]];
                        }
                        
                        // Once the mRNA is degraded, update mRNA struct, R_grid, gene_mRNA_map, free the dynamic memory,
                        // Decay step 1- put the last m_id in the gene_mRNA_map[current gene] into selected c2 position, no change of mRNA struct; 
                        // change the positions of mRNA[next-m_id].gene in gene_mRNA_map and free_mRNA
                        // Decay step 2- copy the contents in mRNA[next_m_id-1/current largest m_id] into the mRNA struct of current dedayable m_id, this equals to the largest m_id has been deleted.
                        
                        M_map[c1]--;   // Counts of total mRNA for each gene in the gene_mRNA_map
                        
                        int imap;
                        for (imap=0;imap<M_map[c1];imap++)  // imap is the index of selected m_id in gene_mRNA_map[c1]
                        {   if(gene_mRNA_map[c1][imap]==m_id)
                        {   break;
                        }
                        }
                        
                        // Step 1: Update gene_mRNA_map and free_mRNA
                        if(imap!=M_map[c1])
                        {   // Move the last m_id for gene c1 into current c2 pos
                            gene_mRNA_map[c1][imap] = gene_mRNA_map[c1][M_map[c1]];
                        }
                        
                        next_m_id--;
                        for (c3=0;c3<M_map[mRNA[next_m_id].gene];c3++)
                        {   if(gene_mRNA_map[mRNA[next_m_id].gene][c3]==next_m_id)
                        {   // Change the largest m_id to current selected m_id in gene_mRNA_map
                            gene_mRNA_map[mRNA[next_m_id].gene][c3]=m_id;
                        }
                        }
                        
                        for (c3=0;c3<Mf[mRNA[next_m_id].gene];c3++)
                        {   if(free_mRNA[mRNA[next_m_id].gene][c3]==next_m_id)
                        {   // Change the largest m_id to current selected m_id in free_mRNA, if the next_m_id does exist in free_mRNA
                            free_mRNA[mRNA[next_m_id].gene][c3]=m_id;
                        }
                        }
                        
                        for (c3=0;c3<M_decSlt[mRNA[next_m_id].gene];c3++)
                        { if(mdec_select_map[mRNA[next_m_id].gene][c3]==next_m_id)
                        {   // Change the largest m_id to current selected m_id in mdec_select_map, if the next_m_id does exist in free_mRNA
                            mdec_select_map[mRNA[next_m_id].gene][c3]=m_id;
                            break;
                        }
                        }
                        
                        // Step 2: Update mRNA struct and R_grid
                        if(m_id!=next_m_id)
                        {   // Copy the info from the mRNA of largest m_id to current m_id
                            mRNA[m_id].gene = mRNA[next_m_id].gene;
                            mRNA[m_id].last_ini = mRNA[next_m_id].last_ini;
                            mRNA[m_id].mdec_delayed = mRNA[next_m_id].mdec_delayed;
                            mRNA[m_id].num_boundr = mRNA[next_m_id].num_boundr;
                            mRNA[m_id].msp_mdec_prob = mRNA[next_m_id].msp_mdec_prob;
                            mRNA[m_id].t_mbirth = mRNA[next_m_id].t_mbirth;
                            
                            // Update R_grid
                            for(c5=0;c5 < obs_max_len;c5++)
                            {   if (c5 < Gene[mRNA[next_m_id].gene].len)
                            {   R_grid[m_id][c5] = R_grid[next_m_id][c5];
                                R_grid[next_m_id][c5]=tot_ribo;  // Clean the info on R_grid[next_m_id]
                            }
                            else if(c5 >=Gene[mRNA[next_m_id].gene].len)
                            {   R_grid[m_id][c5] = tot_ribo;
                            }
                            
                            if(R_grid[m_id][c5]!=tot_ribo)
                            {   // Update the m_id associted with any bound ribosomes
                                Ribo[R_grid[m_id][c5]].mRNA=m_id;
                            }
                            }
                        }
                        
                        tot_mRNA--;
                        Gene[c1].exp --;    // This does not include the mRNAmarkedForDecay
                        tot_mdec_r -= Gene[c1].mdec_r;
                        
                        decCount[c1]++;
                    }

                    else   // This is when the Xrn1 5'-3' co-translational/delayed decay is selected to proceed
                    {
                        // mRNA numbers and decay rate decreases immidiately, whether the selected mRNA is decayed immidiately or decayed later.
                        tot_mRNA--;
                        Gene[g_id].exp --;    // this does not include the mRNAmarkedForDecay
                        tot_mdec_r -= Gene[g_id].mdec_r;  

                        // If this mRNA is in the free_mRNA list, remove it so to make it non-initiable, whether it has bound ribos on or not
                        for(c3=0;c3<Mf[g_id];c3++)
                        { if(free_mRNA[g_id][c3]==m_id)
                        {
                            Mf[g_id]--;
                            if(c3!=Mf[g_id])
                            {   
                                free_mRNA[g_id][c3] = free_mRNA[g_id][Mf[g_id]];
                            }        // Accounting of free mRNAs for next round
                            scld_Mf[g_id]-=Gene[g_id].ini_prob;
                            tot_scld_Mf-=Gene[g_id].ini_prob;
                        }
                        }
                        
                        // Update the M_decSlt so that current selected mRNA won't be selected again for decay either
                        // When the decay is delayed; OR when current mRNA is empty and is deleted immediately
                        M_decSlt[g_id] --;
                        if(c2 != M_decSlt[g_id])  
                        {   mdec_select_map[g_id][c2] = mdec_select_map[g_id][M_decSlt[g_id]];
                        }
                        
                        int c4 = 0;          // To check if there's any bound ribo on pos0-pos10
                        int c5 = 0;          // Count of positions in R_grid
                        
                        // Go through R_grid to find any bound ribo
                        for (c3=0;c3<Gene[c1].len;c3++)  /// Go through R_grid to see if any ribo is bound， if there is then delay the degradation
                        { if(R_grid[m_id][c3]!= tot_ribo)   // if at least one ribo is bound                             
                        {
                            mRNA[m_id].mdec_delayed = 1;   // 1 means mRNA decay is delayed (1=True, 0=False)
                            M_decTransPlus[g_id] ++;           // Number of decay-delayed transcripts for a particular gene
                            break;
                        }
                        }
                        
                        
                        // If there is 0 bound ribosome on current selected mRNA=empty mRNA, then proceed to immediate decay
                        if(mRNA[m_id].mdec_delayed == 0)
                        { 
                            // For empty mRNA (0 bound ribos) to be degraded, update mRNA struct, R_grid, gene_mRNA_map, free the dynamic memory,
                            // Decay step 1- put the last m_id in the gene_mRNA_map into where the current selected c2 is; don't change the mRNA struct contents; change the position of mRNA[next_m_id].gene in gene_mRNA_map and free_mRNA
                            // Decay step 2- copy the mRNA struct of current next_m_id-1 (i.e. the mRNA with the largest m_id) into the mRNA struct of the currently selected mRNA (that's gonna be decayed). This is equivalent to deleting the mRNA with the largest m_id

                            // Below is to update the moments of mRNA life times for each gene, life time = death time(decayed, not just marked for decay)- birth time for any newly synthesized mRNAs after burn-in
                            if(t>=thresh_time && printOpt[6]==1 && mRNA[m_id].t_mbirth >= thresh_time)    // If this selected mRNA is not one of the original 60,000 mRNAs
                            {
                                tot_decCount[g_id] ++;
                                
                                if(count_mlifetime[g_id] < 10000)
                                {
                                    m_lifeTime = t - mRNA[m_id].t_mbirth;

                                    gene_sp_mlifetime[g_id][count_mlifetime[g_id]] = m_lifeTime;
                                    count_mlifetime[g_id]++;
                                }
                            }

                            // Step 1: Update gene_mRNA_map and free_mRNA and mdec_select_map
                            M_map[g_id]--;   // Counts of total mRNA for each gene in the gene_mRNA_map
                            int imap;
                            for(imap=0;imap<M_map[g_id];imap++)
                            {
                                if(gene_mRNA_map[g_id][imap]==m_id)
                                {
                                    if(imap!=M_map[g_id])
                                    {
                                        gene_mRNA_map[g_id][imap] = gene_mRNA_map[g_id][M_map[g_id]];
                                    }
                                    break;
                                }
                            }
                            
                            
                            next_m_id--;
                            for (imap=0;imap<M_map[mRNA[next_m_id].gene];imap++)
                            { if(gene_mRNA_map[mRNA[next_m_id].gene][imap]==next_m_id)
                            {   // Change the largest m_id to current selected m_id in gene_mRNA_map
                                gene_mRNA_map[mRNA[next_m_id].gene][imap]=m_id;
                                break;
                            }
                            }
                            
                            for (c3=0;c3<Mf[mRNA[next_m_id].gene];c3++)
                            { if(free_mRNA[mRNA[next_m_id].gene][c3]==next_m_id)
                            {   // Change the largest m_id to current selected m_id in free_mRNA, if the next_m_id does exist in free_mRNA
                                free_mRNA[mRNA[next_m_id].gene][c3]=m_id;
                                break;
                            }
                            }
                            
                            for (c3=0;c3<M_decSlt[mRNA[next_m_id].gene];c3++)
                            { if(mdec_select_map[mRNA[next_m_id].gene][c3]==next_m_id)
                            {   // Change the largest m_id to current selected m_id in mdec_select_map, if the next_m_id does exist in free_mRNA
                                mdec_select_map[mRNA[next_m_id].gene][c3]=m_id;
                                break;
                            }
                            }
                            
                            // Step 2: Update mRNA struct and R_grid
                            if(m_id!=next_m_id)
                            {   // Copy the info from the mRNA of largest m_id to current m_id
                                mRNA[m_id].gene = mRNA[next_m_id].gene;
                                mRNA[m_id].last_ini = mRNA[next_m_id].last_ini;
                                mRNA[m_id].mdec_delayed = mRNA[next_m_id].mdec_delayed;
                                mRNA[m_id].num_boundr = mRNA[next_m_id].num_boundr;
                                mRNA[m_id].msp_mdec_prob = mRNA[next_m_id].msp_mdec_prob;
                                mRNA[m_id].t_mbirth = mRNA[next_m_id].t_mbirth;
                                
                                // Update R_grid
                                for(c5=0;c5 < obs_max_len;c5++)
                                { 
                                    if (c5 < Gene[mRNA[next_m_id].gene].len)
                                    {   R_grid[m_id][c5] = R_grid[next_m_id][c5];
                                        R_grid[next_m_id][c5]=tot_ribo;  // Clean the info on R_grid[next_m_id]
                                    }
                                    else if(c5 >=Gene[mRNA[next_m_id].gene].len)
                                    {   R_grid[m_id][c5] = tot_ribo;
                                    }
                                    
                                    if(R_grid[m_id][c5]!=tot_ribo)
                                    {   // Update the m_id associted with any bound ribosomes
                                        Ribo[R_grid[m_id][c5]].mRNA=m_id;
                                    }
                                }
                            }
                            
                            decCount[c1]++;
                        }
                        
                    }
                    break;  // Once an mRNA is decayed or delayed for decay, no need to continue the loop
                }
            }
        }
        
        //       mRNA synthesis
        else if (coin < prob_msyn && t>thresh_time)
        {   // Pick a random gene to synthesize
            
            coin = gsl_rng_uniform(r);
            double prob_syn_g = 0.0;
            
            for(c1=0;c1<n_genes;c1++)
            {   prob_syn_g += Gene[c1].msyn_r/tot_msyn_r;
                if(coin < prob_syn_g)           // First pick a gene
                {  m_id = next_m_id;
                    
                    if (m_id>dbl_tot_mRNA)
                    {   printf("Total mRNA number exceeds two times of initial total mRNAs\nOut of memory\n");fflush(stdout);
                    exit(1);
                    }
                    
                    Gene[c1].exp++;       // This does not include the mRNAmarkedForDecay
                    if (Gene[c1].exp>triple_max_exp)
                    {   printf("Gene expression exceeds three times of max mRNAs expressions\nOut of memory\n");fflush(stdout);
                    exit(1);
                    }
                    
                    tot_mRNA++;
                    
                    gene_mRNA_map[c1][M_map[c1]] = m_id;
                    free_mRNA[c1][Mf[c1]] = m_id;
                    mdec_select_map[c1][M_decSlt[c1]] = m_id;
                    
                    M_map[c1]++;
                    Mf[c1]++;
                    M_decSlt[c1] ++;
                    next_m_id ++;
                    
                    // Accounting of initiation probability for next round
                    scld_Mf[c1]+=Gene[c1].ini_prob;
                    tot_scld_Mf+=Gene[c1].ini_prob;
                    
                    
                    // Update R_grid because the current m_id may be used before and may contain old info
                    for(int c5=0;c5<obs_max_len;c5++)
                    {   R_grid[m_id][c5] = tot_ribo;
                    }
                    
                    // Update mRNA struct
                    mRNA[m_id].gene = c1;
                    mRNA[m_id].last_ini = t; 
                    mRNA[m_id].mdec_delayed = 0;
                    mRNA[m_id].num_boundr = 0;
                    mRNA[m_id].msp_mdec_prob = 0.0;
                    mRNA[m_id].t_mbirth = t;                    // Updating the mRNA birth time for the newly synthsized mRNAs after burn-in.
                    
                    tot_mdec_r += Gene[c1].mdec_r;    
                    
                    synCount[c1]++;
                    
                    break;    // Once the selected gene is synthesized, no longer need to continue the loop
                }
            }
        }
        
        
        else    // Elongation cycle
        { 
            
            if(t>thresh_time)
            {
                tmp_elng_prob = prob_msyn;
            }
            
            else
            {
                tmp_elng_prob = prob_ini;
            }
            
            for(c1=0;c1<61;c1++)
            {    prob_e[c1] = tmp_elng_prob + r_elng[c1]*inv_rate;               // Prob of elongn of codon 1
                // r_elng depends on 1)Tf[c1]  2)wobble para   3)n_Rb_e[c1]
                if(coin<prob_e[c1])
                {   c_id = c1;                                      // Codon with number=c1 as its c_id selected
                    break;
                }
                tmp_elng_prob = prob_e[c1];
            }
            
            x = gsl_rng_uniform_int(r, (unsigned long)n_Rb_e[c_id]);                // Randomly pick an elongatable ribosome bound to codon c_id
            // c_id can be 0:61
            r_id = Rb_e[c_id][x];               // Selected ribosome that's bound to the selected codon
            m_id = Ribo[r_id].mRNA;             // mRNA that's bound to the selected ribosome with the r_id
            g_id = mRNA[m_id].gene;
            
            if(t>thresh_time && n_e_times[c_id] < 1e9)
            {   if(printOpt[0]==1)
            {   
                current_etimes = t-Ribo[r_id].t_elong_ini;                    // For estimation of avg elongation times of codons
                n_e_times[c_id]++;
                
                // Update the mean of e_times
                e_times[c_id] = e_times[c_id] + (current_etimes - e_times[c_id])/(double)n_e_times[c_id];
            }
            }

            // Normal elongation cycle
            if(Ribo[r_id].pos>0)                                                    // Release the tRNA bound at the earlier position
            {
                Tf[cTRNA[Gene[mRNA[m_id].gene].seq[Ribo[r_id].pos-1]].tid]++;
                
            }
            
            if(Ribo[r_id].pos==(Gene[mRNA[m_id].gene].len-1))     // Check if the current elongation has led to termination
                // If it has led to termination
            {

                mRNA[m_id].num_boundr --;                         // Update the number of bound ribosomes on the mRNA
                
                R_grid[m_id][Ribo[r_id].pos] = tot_ribo;                            // Update the position of ribosomes on the mRNA
                
                Ribo[r_id].pos++;
                Rf++;                                                                // Free a ribosome upon termination
                
                n_Rb_e[c_id]--;
                if(x!=n_Rb_e[c_id])
                {   Ribo[Rb_e[c_id][n_Rb_e[c_id]]].elng_pos_list = x;                // Update the ids and number of elongatable ribosomes
                    Rb_e[c_id][x] = Rb_e[c_id][n_Rb_e[c_id]];                       // This is doing a swap
                }

                if(t>thresh_time && printOpt[1]==1)
                {
                    current_trans = t - Ribo[r_id].t_trans_ini;
                    old_mean = avg_time_to_trans[g_id];
                    n_trans[g_id] ++;
                    
                    // Update the mean
                    avg_time_to_trans[g_id] = avg_time_to_trans[g_id] + (current_trans - avg_time_to_trans[g_id])/n_trans[g_id];
                    
                    // Update the variance
                    if(n_trans[g_id] >= 2)
                    {
                        sd_time_to_trans[g_id] = sd_time_to_trans[g_id] + ((current_trans - old_mean)*(current_trans - avg_time_to_trans[g_id]));
                        var_time_to_trans[g_id] = sd_time_to_trans[g_id]/(n_trans[g_id]-1);
                    }
                }

                // Update any previously unelongatable ribosomes
                if(R_grid[m_id][Ribo[r_id].pos-11]<tot_ribo)                        // When the ribosome moves, a previously unelongatable ribosome
                {    c2_id = Gene[mRNA[m_id].gene].seq[Ribo[r_id].pos-11];            // can now be elongatable on the same mRNA if its 11 codon behind
                    Rb_e[c2_id][n_Rb_e[c2_id]] = R_grid[m_id][Ribo[r_id].pos-11];
                    
                    Ribo[R_grid[m_id][Ribo[r_id].pos-11]].elng_cod_list = c2_id;
                    Ribo[R_grid[m_id][Ribo[r_id].pos-11]].elng_pos_list = n_Rb_e[c2_id];
                    n_Rb_e[c2_id]++;
                }
                
                // Update the running mean and var elongation time of a particular gene
                // Used formula by B.P. Welford www.johndcook.com/blog/standard_deviation/
                if(t>=thresh_time && printOpt[5]==1)
                {
                    elngTime = t - Ribo[r_id].t_trans_ini;
                    old_mean = meanElngTime[g_id];  // For the convience of calculating new variance
                    
                    ctElngTime[g_id]++;
                    
                    // Update the mean
                    meanElngTime[g_id] = meanElngTime[g_id] + (elngTime - meanElngTime[g_id])/ctElngTime[g_id];
                    
                    // Update the variance
                    if(ctElngTime[g_id] >= 2)
                    {
                        sdElngTime[g_id] = sdElngTime[g_id] + ((elngTime - old_mean)*(elngTime - meanElngTime[g_id]));
                        varElngTime[g_id] = sdElngTime[g_id]/(ctElngTime[g_id]-1);
                    }
                }
                
                
                // Update the ribosomes avail for initiation
                next_avail_ribo--;
                
                
                if(r_id!=next_avail_ribo)                                            // The last initiated ribosome's data is swapped with
                {    Ribo[r_id].mRNA = Ribo[next_avail_ribo].mRNA;                    // the ribosome that just finished translation (terminated)
                    Ribo[r_id].pos = Ribo[next_avail_ribo].pos;                        // This is done to primarily keep track of ONLY ribosomes
                    Ribo[r_id].t_trans_ini = Ribo[next_avail_ribo].t_trans_ini;        // that are currently bound for faster computation
                    Ribo[r_id].t_elong_ini = Ribo[next_avail_ribo].t_elong_ini;
                    Ribo[r_id].elng_cod_list = Gene[mRNA[Ribo[next_avail_ribo].mRNA].gene].seq[Ribo[next_avail_ribo].pos];
                    Ribo[r_id].elng_pos_list = Ribo[next_avail_ribo].elng_pos_list;
                    R_grid[Ribo[r_id].mRNA][Ribo[r_id].pos] = r_id;
                    
                    if(R_grid[Ribo[r_id].mRNA][Ribo[r_id].pos+10]==tot_ribo || (Ribo[r_id].pos+10)>=Gene[mRNA[Ribo[r_id].mRNA].gene].len)
                    {    Rb_e[Ribo[r_id].elng_cod_list][Ribo[r_id].elng_pos_list] = r_id;
                    }
                }
                termtn_now = 1;
                
                
                if(t>thresh_time)
                {
                    int check_last_ribo = 0;  // 0=when current ribo is the last ribo on the mRNA, 1=when there are other ribos on current mRNA
                    for (c3=0;c3<Gene[g_id].len;c3++)  // Go through R_grid to see if any ribo is bound,if there is
                    { if(R_grid[m_id][c3]!= tot_ribo)
                    {
                        check_last_ribo = 1;
                        break;
                    }
                    }
                    
                    
                    if(check_last_ribo==0 && mRNA[m_id].mdec_delayed==1) // If current ribo is the last one, and if current mRNA was marked degradable
                    {
                        // Below is to update the running moments of mRNA life times for each gene, life time = death time(decayed, not just marked for decay)- birth time for any newly synthesized mRNAs after burn-in
                        if(t>=thresh_time && printOpt[6]==1 && mRNA[m_id].t_mbirth >= thresh_time)    // if this selected mRNA is not one of the original 60,000 mRNAs
                        {
                            tot_decCount[g_id] ++;
                            
                            if(count_mlifetime[g_id] < 10000)
                            {
                                m_lifeTime = t - mRNA[m_id].t_mbirth;

                                gene_sp_mlifetime[g_id][count_mlifetime[g_id]] = m_lifeTime;
                                count_mlifetime[g_id]++;
                            }
                        }
      
                        M_decTransMinus[g_id] ++;   // Update the number of mRNAs that's marked for deletion
                        
                        // For empty mRNA (0 bound ribos) to be degraded, update mRNA struct, R_grid, gene_mRNA_map, free the dynamic memory, then similar to as above the mRNA decay
                        
                        M_map[g_id]--;   // Counts of total mRNA for each gene in the gene_mRNA_map
                        int imap;
                        for (imap=0;imap<M_map[g_id];imap++)  // imap is the index of selected m_id in gene_mRNA_map[c1]
                        { if(gene_mRNA_map[g_id][imap]==m_id)
                        {   
                            // Step 1: Update gene_mRNA_map and free_mRNA
                            if(imap!=M_map[g_id])
                            {   // Move the last m_id for gene c1 into current imap pos
                                gene_mRNA_map[g_id][imap] = gene_mRNA_map[g_id][M_map[g_id]];
                            }
                            
                            break;
                        }
                        }

                        next_m_id--;
                        for (c3=0;c3<M_map[mRNA[next_m_id].gene];c3++)
                        { if(gene_mRNA_map[mRNA[next_m_id].gene][c3]==next_m_id)
                        {   // Change the largest m_id to current selected m_id in gene_mRNA_map
                            gene_mRNA_map[mRNA[next_m_id].gene][c3]=m_id;
                            break;
                        }
                        }
                        
                        for (c3=0;c3<Mf[mRNA[next_m_id].gene];c3++)
                        { if(free_mRNA[mRNA[next_m_id].gene][c3]==next_m_id)
                        {   // Change the largest m_id to current selected m_id in free_mRNA, if the next_m_id does exist in free_mRNA
                            free_mRNA[mRNA[next_m_id].gene][c3]=m_id;
                            break;
                        }
                        }
                        
                        for (c3=0;c3<M_decSlt[mRNA[next_m_id].gene];c3++)
                        { if(mdec_select_map[mRNA[next_m_id].gene][c3]==next_m_id)
                        {   // Change the largest m_id to current selected m_id in mdec_select_map, if the next_m_id does exist in mdec_select_map
                            mdec_select_map[mRNA[next_m_id].gene][c3]=m_id;
                            break;
                        }
                        }
                        
                        // Step 2: Update mRNA struct and R_grid
                        if(m_id!=next_m_id)
                        {   // Copy the info from the mRNA of largest m_id to current m_id
                            mRNA[m_id].gene = mRNA[next_m_id].gene;
                            mRNA[m_id].last_ini = mRNA[next_m_id].last_ini;
                            mRNA[m_id].mdec_delayed = mRNA[next_m_id].mdec_delayed;
                            mRNA[m_id].num_boundr = mRNA[next_m_id].num_boundr;
                            mRNA[m_id].msp_mdec_prob = mRNA[next_m_id].msp_mdec_prob;
                            mRNA[m_id].t_mbirth = mRNA[next_m_id].t_mbirth;
                            
                            // Update R_grid
                            for(int c5=0;c5 < obs_max_len;c5++)
                            { if (c5 < Gene[mRNA[next_m_id].gene].len)
                            {     
                                R_grid[m_id][c5] = R_grid[next_m_id][c5];
                                R_grid[next_m_id][c5]=tot_ribo;  // Clean the info on R_grid[next_m_id]
                            }
                            else if(c5 >=Gene[mRNA[next_m_id].gene].len)
                            { R_grid[m_id][c5] = tot_ribo;
                            }
                            
                            if(R_grid[m_id][c5]!=tot_ribo)
                            { // Update the m_id associted with any bound ribosomes
                                Ribo[R_grid[m_id][c5]].mRNA=m_id;
                            }
                            }
                        }
                        
                        decCount[g_id]++;
                    }
                }
            }
            
            
            
            else if((Ribo[r_id].pos+11)>=Gene[mRNA[m_id].gene].len || R_grid[m_id][Ribo[r_id].pos+11]==tot_ribo)    // Check if the ribosome is still elongatable and, 
                // below is when if it's still elongatable (when current elongation hasn't led to termination or stalled ribo)
            {
                
                R_grid[m_id][Ribo[r_id].pos] = tot_ribo;                    // Update the position of ribosomes on the mRNA
                Ribo[r_id].pos++;
                R_grid[m_id][Ribo[r_id].pos] = r_id;
                c2_id = Gene[mRNA[m_id].gene].seq[Ribo[r_id].pos];
                
                if(c2_id!=c_id)                                            // If the codon has changed shift the elongatable ribosome
                {                                                          // to the other codon
                    Rb_e[c2_id][n_Rb_e[c2_id]] = r_id;  
                    Ribo[r_id].elng_cod_list = c2_id;
                    Ribo[r_id].elng_pos_list = n_Rb_e[c2_id];
                    n_Rb_e[c2_id]++;
                    
                    n_Rb_e[c_id]--;
                    if(x!=n_Rb_e[c_id])
                    {    Ribo[Rb_e[c_id][n_Rb_e[c_id]]].elng_pos_list = x;                // Update the ids and number of elongatable ribosomes
                        
                        Rb_e[c_id][x] = Rb_e[c_id][n_Rb_e[c_id]];
                    }
                }
                
                
                Tf[cTRNA[c_id].tid]--;
                
            }
            else           // When ribosome is not elongatable anymore because it'll lead to stalled ribo(has another ribo at +10 position)
            { 
                
                R_grid[m_id][Ribo[r_id].pos] = tot_ribo;                     // Update the position of ribosomes on the mRNA
                
                Ribo[r_id].pos++;
                R_grid[m_id][Ribo[r_id].pos] = r_id;
                Ribo[r_id].elng_cod_list = Gene[mRNA[m_id].gene].seq[Ribo[r_id].pos];
                
                n_Rb_e[c_id]--;
                if(x!=n_Rb_e[c_id])
                {    Ribo[Rb_e[c_id][n_Rb_e[c_id]]].elng_pos_list = x;                    // Update the ids and number of elongatable ribosomes
                    Rb_e[c_id][x] = Rb_e[c_id][n_Rb_e[c_id]];
                }
                
                Tf[cTRNA[c_id].tid]--;
                
            }
            // End of normal elongation cycle
            
            if(termtn_now==0)                                                            // When elongation does not lead to termination
            {
                if(Ribo[r_id].pos>10 && R_grid[m_id][Ribo[r_id].pos-11]<tot_ribo)  // When the ribosome moves, a previously unelongatable ribosome
                { if(Ribo[R_grid[m_id][Ribo[r_id].pos-11]].inhbtr_bound==0)    // can now be elongated on the same mRNA if its 11 codon behind
                { c2_id = Gene[mRNA[m_id].gene].seq[Ribo[r_id].pos-11];
                    Rb_e[c2_id][n_Rb_e[c2_id]] = R_grid[m_id][Ribo[r_id].pos-11];
                    
                    Ribo[R_grid[m_id][Ribo[r_id].pos-11]].elng_cod_list = c2_id;
                    Ribo[R_grid[m_id][Ribo[r_id].pos-11]].elng_pos_list = n_Rb_e[c2_id];
                    n_Rb_e[c2_id]++;
                }
                }
                if(Ribo[r_id].pos==10 && mRNA[m_id].mdec_delayed==0)                    // If the current position is at codon 11 and if the mRNA is not marked for degradation
                    // make the current mRNA initiable
                {   
                    g_id = mRNA[m_id].gene;
                    free_mRNA[g_id][Mf[g_id]] = m_id;
                    Mf[g_id]++;
                    scld_Mf[g_id] += Gene[g_id].ini_prob;
                    tot_scld_Mf += Gene[g_id].ini_prob;
                }
                
                Ribo[r_id].t_elong_ini = t;                        // Upon elongation, update the elong ini time for the next event
            }
            else
            { termtn_now = 0;
            }
        }
        
        
        if(t>((double)t_print) && printOpt[5]==1)    // Output the results every minute
        {
            fprintf(f0,"%d ",t_print);
            fprintf(f1,"%d ",t_print);
            fprintf(f2,"%d ",t_print);
            fprintf(f3,"%d ",t_print);
            fprintf(f4,"%d ",t_print);
            fprintf(f5,"%d ",t_print);
            fprintf(f6,"%d ",t_print);
            fprintf(f7,"%d %d ",t_print,Rf);
            fprintf(f8,"%d ",t_print);
            fprintf(f9,"%d ",t_print);
            fprintf(f10,"%d ",t_print);
            
            
            
            for(c1=0;c1<n_genes;c1++)
            {   fprintf(f0,"%d ",decCount[c1]);
                fprintf(f1,"%d ",synCount[c1]);
                fprintf(f2,"%g ",meanElngTime[c1]);
                fprintf(f3,"%g ",varElngTime[c1]);
                fprintf(f4,"%d ",ctElngTime[c1]);  // How many times a gene has been translated during each time slot
                
                mRNAcount[c1] = mRNAcount[c1] + synCount[c1] - decCount[c1];
                fprintf(f5,"%d ",mRNAcount[c1]);
                
                M_decTrans[c1] = M_decTrans[c1] + M_decTransPlus[c1] - M_decTransMinus[c1];
                fprintf(f6,"%d ",M_decTrans[c1]); // How many mRNAs have been marked for decay during each time slot
                
            }
            
            fflush(stdout);
            
            for(c3=0;c3<61;c3++)
            {
                if(Tf[c3]>0)
                {
                    fprintf(f7,"%d ",Tf[c3]);
                }
                
                fprintf(f8,"%d ",n_Rb_e[c3]);
            }
            
            for(c1=0;c1<n_genes;c1++)
            {
                for(c2=0;c2<M_map[c1];c2++)  // The running mean and variances are accounting for the snapshot of the end of each minute (one time step),
                                             // but has to iterate through all the mRNAs of gene c1
                {
                    tot_bound_ribo[c1] += mRNA[gene_mRNA_map[c1][c2]].num_boundr;
                    old_mean = meanBoundRibo[c1];  // old_mean is the mean before accounting the current mRNA
                    meanBoundRibo[c1] = old_mean + (mRNA[gene_mRNA_map[c1][c2]].num_boundr-old_mean)/(c2+1);
                    
                    if(c2>=1)  // When c2>=1, that means there are 2 mRNAs for this gene
                    {
                        intermedBoundRibo[c1] =  intermedBoundRibo[c1] + ((mRNA[gene_mRNA_map[c1][c2]].num_boundr - old_mean)*(mRNA[gene_mRNA_map[c1][c2]].num_boundr - meanBoundRibo[c1]));
                        varBoundRibo[c1] = intermedBoundRibo[c1]/c2;
                    }
                }

                fprintf(f9,"%d ",tot_bound_ribo[c1]);
                fprintf(f10,"%lf ",varBoundRibo[c1]);
            }
            fflush(stdout);
            
            fprintf(f0,"\n");
            fprintf(f1,"\n");
            fprintf(f2,"\n");
            fprintf(f3,"\n");
            fprintf(f4,"\n");
            fprintf(f5,"\n");
            fprintf(f6,"\n");
            fprintf(f7,"\n");
            fprintf(f8,"\n");
            fprintf(f9,"\n");
            fprintf(f10,"\n");
            
            
            for(c1=0;c1<n_genes;c1++)
            {   decCount[c1] = 0;
                synCount[c1] = 0;
                meanElngTime[c1] = 0.0;
                sdElngTime[c1] = 0.0;
                varElngTime[c1] = 0.0;
                ctElngTime[c1] = 0;
                M_decTransPlus[c1] = 0;
                M_decTransMinus[c1] = 0;
                tot_bound_ribo[c1] = 0;
                meanBoundRibo[c1] = 0.0;
                intermedBoundRibo[c1] = 0.0;
                varBoundRibo[c1] = 0.0;
            }
            fflush(stdout);
            t_print = t_print + 60;
        }
        ii ++;
        }
    
    fclose(f0);
    fclose(f1);
    fclose(f2);
    fclose(f3);
    fclose(f4);
    fclose(f5);
    fclose(f6);
    fclose(f7);
    fclose(f8);
    fclose(f9);
    fclose(f10);
    
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    // Process output for printing
    
    // Print the output
    // Elongation times of all codons
    if(printOpt[0]==1)
    {    strcpy(out_file,out_prefix);
        f2 = fopen(strcat(out_file,".etimes.out"),"w");
        
        fprintf(f2,"Codon\tNum_of_events\tAvg_elong_time(sec)\n");
        for(c1=0;c1<61;c1++)
        {    
            fprintf(f2,"%d\t%d\t%g\n",c1,n_e_times[c1],e_times[c1]);
            
        }
        fclose(f2);
    }
    
    // Average total elongation times of all genes
    if(printOpt[1]==1)
    {    strcpy(out_file,out_prefix);
        f3 = fopen(strcat(out_file,".gene_totetimes.out"),"w");

        fprintf(f3,"Gene\tNum_of_events\tAvg_trans_time(sec)\tVariance_trans_time\n");
        for(c1=0;c1<n_genes;c1++)
        {    
            fprintf(f3,"%d\t%d\t%g\t%g\n",c1,n_trans[c1],avg_time_to_trans[c1],var_time_to_trans[c1]);
        }
        fclose(f3);
    }
    
    // Average time between initiation of all genes
    if(printOpt[2]==1)
    {    strcpy(out_file,out_prefix);
        f4 = fopen(strcat(out_file,".gene_initimes.out"),"w");

        fprintf(f4,"Gene\tNum_of_events\tAvg_initiation_time(sec)\tVariance_initiation_time\n");
        for(c1=0;c1<n_genes;c1++)
        {   
            fprintf(f4,"%d\t%d\t%g\t%g\n",c1,n_ini[c1],avg_time_to_ini[c1],var_time_to_ini[c1]);
        }
        
        fclose(f4);
    }
    
    // Average number of free ribosomes and tRNAs at equilibrium
    if(printOpt[3]==1)
    {    strcpy(out_file,out_prefix);
        f5 = fopen(strcat(out_file,".avg_ribo_tRNA.out"),"w");
        
        avg_Rf = avg_Rf/(tot_time-thresh_time);
        fprintf(f5,"Free_ribo\t%g\n",avg_Rf);
        for(c1=0;c1<61;c1++)
        {    if(avg_tRNA_abndc[c1]>0)
        {    
            avg_tRNA_abndc[c1] = avg_tRNA_abndc[c1]/(tot_time-thresh_time);
            fprintf(f5,"Free_tRNA%d\t%g\n",c1,avg_tRNA_abndc[c1]);
        }
        }
        fclose(f5);
    }
    
    // The final state of the system - positions of bound ribosomes on mRNAs
    if(printOpt[4]==1)
    {   strcpy(out_file,out_prefix);
        f6 = fopen(strcat(out_file,".final_ribo_pos.out"),"w");
        
        // Print final state for individual mRNAs
        for(c1=0;c1<tot_mRNA;c1++)
        { if(R_grid[c1][0]==tot_ribo)
        {    fprintf(f6,"0");
        }
        else
        {    fprintf(f6,"1");
        }
        
        for(c2=1;c2<Gene[mRNA[c1].gene].len;c2++)
        { if(R_grid[c1][c2]==tot_ribo)
        {    fprintf(f6," 0");
        }
        else
        {    fprintf(f6," 1");
        }
        }
        fprintf(f6,"\n");
        }
        fclose(f6);
    }
    
    
    if(printOpt[6]==1)
    {   strcpy(out_file,out_prefix);
        f7 = fopen(strcat(out_file,".mRNA_lifeTime_10000.out"),"w");

        for(c1=0;c1<n_genes;c1++)
        {
            fprintf(f7, "%d %d ",c1+1,count_mlifetime[c1]);  // print out the g_id and the count recorded for this gene
            for(c2=0;c2<count_mlifetime[c1];c2++)
            {
                fprintf(f7, "%lf ",gene_sp_mlifetime[c1][c2]);  // print out the first 10,000 life times of this gene
            }
            fprintf(f7,"\n");
        }
        fclose(f7);

    }
    // Free the malloc structures and arrays
    free(Gene);
    free(Ribo);
    free(mRNA);
    free(cTRNA);
    free(Rb_e);
    free(free_mRNA);
    free(R_grid);
    free(gene_mRNA_map);
    free(mdec_select_map);
    free(gene_sp_mlifetime);
    gsl_rng_free(r);
    
    time_t rawtime;
    struct tm * timeinfo;
    
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    printf("\nEnd of program,system time=%s\n",asctime (timeinfo));fflush(stdout);
    
    }

