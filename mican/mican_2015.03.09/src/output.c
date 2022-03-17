#include "mican.h"

void output(int naa_q, int naa_t, int natom_q, int natom_t, RESDAT *resdat_q,
            RESDAT *resdat_t, ALLATM *allatm_q, ALLATM *allatm_t, ALIGN *align,
            ALIGN ***align_chain, PDBDAT pdbdat_q, PDBDAT pdbdat_t) {
  int i, j;
  FILE *fp;
  int ichain, jchain;

  char filename1[STRSIZE], filename2[STRSIZE];
  int naa1, nchain1, nchain2;
  char *chainID1, *chainID2;
  float TMscore1, TMscore2;
  float coverage1, coverage2;

  /*************************/
  /**    main output    ****/
  /*************************/
  /****  output protein data  ****/
  if (input.qtchange == 0) {
    strcpy(filename1, input.file_t);
    strcpy(filename2, input.file_q);
    naa1 = naa_t;
    nchain1 = pdbdat_t.nchain;
    nchain2 = pdbdat_q.nchain;
  } else {
    strcpy(filename1, input.file_q);
    strcpy(filename2, input.file_t);
    naa1 = naa_q;
    nchain1 = pdbdat_q.nchain;
    nchain2 = pdbdat_t.nchain;
  }

  /*****************************/
  /**    brief description    **/
  /*****************************/
  if(input.silent == OFF){
    printf(" ===== Brief description of top %d alignments =====\n", input.nsub);
    printf(" Rank   MICAN   TMscore  Dali.Z   SPscore  Length   RMSD   Seq.Id.\n");
    for (i = 0; i <= input.nsub - 1; i++) {
      if (align[i].check == OFF) {
	break;
      }
      if (input.qtchange == 0) {
	TMscore1 = align[i].TMscore_t;
      } else {
	TMscore1 = align[i].TMscore_q;
      }
      printf(" %4d  %6.3f   %6.3f   %6.2f   %6.3f    %4d  %6.3f   %5.1f\n",
	     i + 1, align[i].TMscore_mod, TMscore1, align[i].DaliZ,
	     align[i].SPscore, align[i].naa_align,
	     align[i].rmsd, align[i].seqID);
    }
    printf(" (TMscore was normalized by size of Protein1. 'Dali.Z'=Dali"
	   " Zscore.\n  'Length'=number of aligned residues. 'Seq.Id.'=Sequence Identity.)\n");
    if (input.nsub_org > input.nsub) {
      fprintf(stderr,
	      "!! Warning: Mican detected only %d alignments, though\n",
	      input.nsub);
      fprintf(stderr,
	      "!! maximum number of solutions were set to %d.\n",
	      input.nsub_org);
      fprintf(stderr,
	      "!! For more alignments, please set a number of GH candidates\n");
      fprintf(stderr, "!! with -g option (E.g. '-g 100', '-g 500', or more.)\n");
    }
    printf("\n");

  }

  if (align[input.isub_out - 1].check != TRUE) {
    printf(" Warning: %d-th solution can not be detected.\n\n",
	   input.isub_out);
    exit(0);
  }

  /********************************/
  /**    detailed description    **/
  /********************************/
  if(input.silent == OFF){
    printf(" ===== Detailed description of alignment %d =====\n", input.isub_out);
    
    if (input.qtchange == 0) {
      TMscore1 = align[input.isub_out - 1].TMscore_t;
      TMscore2 = align[input.isub_out - 1].TMscore_q;
      coverage1 = align[input.isub_out - 1].cover_t;
      coverage2 = align[input.isub_out - 1].cover_q;
    } else {
      TMscore1 = align[input.isub_out - 1].TMscore_q;
      TMscore2 = align[input.isub_out - 1].TMscore_t;
      coverage1 = align[input.isub_out - 1].cover_q;
      coverage2 = align[input.isub_out - 1].cover_t;
    }

    printf(
	   " TM-score=%5.3f, Coverage=%5.1f%% (if normalized by size of Protein1)\n",
	   TMscore1, coverage1);
    printf(
	   " TM-score=%5.3f, Coverage=%5.1f%% (if normalized by size of Protein2)\n",
	   TMscore2, coverage2);
    printf(" Dali-score=%7.2f, Dali Zscore=%5.2f\n",
	   align[input.isub_out-1].Daliscore, align[input.isub_out-1].DaliZ);
    
    /****  chain results  ****/
    printf(" ----- Results for each chain -----\n");
    for (jchain = 1; jchain <= nchain1; jchain++) {
      if (input.qtchange == OFF) {
	coverage1 = align_chain[input.isub_out - 1][0][jchain].cover_t;
      } else {
	coverage1 = align_chain[input.isub_out - 1][jchain][0].cover_q;
      }

      if (coverage1 > 0.0) {
	if (input.qtchange == OFF) {
	  chainID1 = pdbdat_t.chainID_org[jchain];
	  naa1 = pdbdat_t.naa[jchain];
	} else {
	  chainID1 = pdbdat_q.chainID_org[jchain];
	  naa1 = pdbdat_q.naa[jchain];
	}
	printf(" [P1:%1s (size %4d)] Coverage=%5.1f%% (",
	       chainID1, naa1, coverage1);
	
	for (ichain = 1; ichain <= nchain2; ichain++) {
	  if (input.qtchange == OFF) {
	    coverage2 = align_chain[input.isub_out - 1][ichain][jchain].cover_t;
	  } else {
	    coverage2 = align_chain[input.isub_out - 1][jchain][ichain].cover_q;
	  }

	  if (coverage2 > 0.0) {
	    if (input.qtchange == OFF) {
	      chainID2 = pdbdat_q.chainID_org[ichain];
	    } else {
	      chainID2 = pdbdat_t.chainID_org[ichain];
	    }
	    printf("P2:%1s %5.1f%% ", chainID2, coverage2);
	  }
	}
	printf(")\n");
      }
    }
  }

  /****  translation matrix  ****/
  if (input.qtchange == ON) {
    for (i = 0; i <= input.nsub - 1; i++) {
      inverse_mat(align[i].rot, align[i].vec);
    }
  }

  /****  small value  ****/
  for (i = 1; i <= 3; i++) {
    for (j = 1; j <= 3; j++) {
      if (fabs(align[input.isub_out - 1].rot[i][j]) < EPS) {
        align[input.isub_out - 1].rot[i][j] = 0.0;
      }
    }
    if (fabs(align[input.isub_out - 1].vec[i]) < 50 * EPS) {
      align[input.isub_out - 1].vec[i] = 0.0;
    }
  }

  if(input.silent == OFF){
    printf(" ----- Rotation matrix to rotate Protein1 to Protein2 -----\n");
    printf(" i       t(i)         u(i,1)       u(i,2)       u(i,3)\n");
    printf(" 1   %9.4f      %9.6f    %9.6f    %9.6f\n",
	   align[input.isub_out - 1].vec[1], align[input.isub_out - 1].rot[1][1],
	   align[input.isub_out - 1].rot[1][2],
	   align[input.isub_out - 1].rot[1][3]);
    printf(" 2   %9.4f      %9.6f    %9.6f    %9.6f\n",
	   align[input.isub_out - 1].vec[2], align[input.isub_out - 1].rot[2][1],
	   align[input.isub_out - 1].rot[2][2],
	   align[input.isub_out - 1].rot[2][3]);
    printf(" 3   %9.4f      %9.6f    %9.6f    %9.6f\n",
	   align[input.isub_out - 1].vec[3], align[input.isub_out - 1].rot[3][1],
	   align[input.isub_out - 1].rot[3][2],
	   align[input.isub_out - 1].rot[3][3]);
    printf(" ----------------------------------------------------------\n");
    printf(" t(i):Translation vector, ");
    printf(" u(i,j):Rotation matrix\n");
    printf("\n");
  }

  /****  filename  ****/
  if(input.silent == OFF){
    if (strcmp(input.aliout, "OFF") != 0) {
      printf(" alignment file          = %s\n", input.aliout);
    }
    if (strcmp(input.pdbout, "OFF") != 0) {
      printf(" superposition pdb file  = %s\n", input.pdbout);
    }
    if (strcmp(input.matout, "OFF") != 0) {
      printf(" translation matrix file = %s\n", input.matout);
    }
  }

  /***************************************/
  /**       output alignment file       **/
  /***************************************/
  if (strncmp(input.aliout, "OFF", 3) != 0) {
    fp = fopen(input.aliout, "w");
    printali(fp, naa_q, naa_t, resdat_q, resdat_t, align);
    fclose(fp);
  }

  /****************************************/
  /**    output superposition pdbfile    **/
  /****************************************/
  if (strncmp(input.pdbout, "OFF", 3) != 0) {
    fp = fopen(input.pdbout, "w");
    printsup(fp, naa_q, naa_t, natom_q, natom_t, resdat_q, resdat_t, allatm_q,
             allatm_t, align);
    fclose(fp);
  }

  /****************************************/
  /**   output translation matrix file   **/
  /****************************************/
  if (strncmp(input.matout, "OFF", 3) != 0) {
    fp = fopen(input.matout, "w");
    printmat(fp, align);
    fclose(fp);
  }

  return;
}
