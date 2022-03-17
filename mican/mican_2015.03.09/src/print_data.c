#include "mican.h"

/**  function print INPUT protein data  **/
void printinp(int naa_q, int naa_t, PDBDAT pdbdat_q, PDBDAT pdbdat_t) {
  char filename1[STRSIZE], filename2[STRSIZE];
  int naa1, naa2, nchain1, nchain2;
  char chainuse1[STRSIZE], chainuse2[STRSIZE];

  if (input.qtchange == 0) {
    strcpy(filename1, input.file_t);
    strcpy(filename2, input.file_q);
    naa1 = naa_t;
    naa2 = naa_q;
    nchain1 = pdbdat_t.nchain;
    nchain2 = pdbdat_q.nchain;
    strcpy(chainuse1, input.chain_t);
    strcpy(chainuse2, input.chain_q);
  } else {
    strcpy(filename1, input.file_q);
    strcpy(filename2, input.file_t);
    naa1 = naa_q;
    naa2 = naa_t;
    nchain1 = pdbdat_q.nchain;
    nchain2 = pdbdat_t.nchain;
    strcpy(chainuse1, input.chain_q);
    strcpy(chainuse2, input.chain_t);
  }

  /****  output input protein data  ****/
  if (strncmp(chainuse1, "0", 1) == 0) {
    printf(" Protein1 (%4d residues,%2d chain) = %s\n", naa1, nchain1,
           filename1);
  } else {
    printf(" Protein1 (%4d residues,%2d chain) = %s (chain %s)\n", naa1,
           nchain1, filename1, chainuse1);
  }
  if (strncmp(chainuse2, "0", 1) == 0) {
    printf(" Protein2 (%4d residues,%2d chain) = %s\n", naa2, nchain2,
           filename2);
  } else {
    printf(" Protein2 (%4d residues,%2d chain) = %s (chain %s)\n", naa2,
           nchain2, filename2, chainuse2);
  }
  printf("\n");

  return;
}

/**  function print superposition pdb  **/
void printsup(FILE *fp, int naa_q, int naa_t, int natom_q, int natom_t,
              RESDAT *resdat_q, RESDAT *resdat_t, ALLATM *allatm_q,
              ALLATM *allatm_t, ALIGN *align) {
  int i, iatom, isub_out;
  float co[4];

  isub_out = input.isub_out - 1;

  /**********************************/
  /**   output all atom  PDB file  **/
  /**********************************/
  /****  rotation  ****/
  if (input.qtchange == ON) {
    for (iatom = 1; iatom <= natom_q; iatom++) {
      for (i = 1; i <= 3; i++) {
        co[i] = align[isub_out].rot[i][1] * allatm_q[iatom].coord[1] +
                align[isub_out].rot[i][2] * allatm_q[iatom].coord[2] +
                align[isub_out].rot[i][3] * allatm_q[iatom].coord[3];
      }
      allatm_q[iatom].coord[1] = co[1] + align[isub_out].vec[1];
      allatm_q[iatom].coord[2] = co[2] + align[isub_out].vec[2];
      allatm_q[iatom].coord[3] = co[3] + align[isub_out].vec[3];
    }
  } else {
    for (iatom = 1; iatom <= natom_t; iatom++) {
      for (i = 1; i <= 3; i++) {
        co[i] = align[isub_out].rot[i][1] * allatm_t[iatom].coord[1] +
                align[isub_out].rot[i][2] * allatm_t[iatom].coord[2] +
                align[isub_out].rot[i][3] * allatm_t[iatom].coord[3];
      }
      allatm_t[iatom].coord[1] = co[1] + align[isub_out].vec[1];
      allatm_t[iatom].coord[2] = co[2] + align[isub_out].vec[2];
      allatm_t[iatom].coord[3] = co[3] + align[isub_out].vec[3];
    }
  }
  /****  output  ****/
  if (input.qtchange == OFF) {
    printrascript(fp, naa_t, naa_q, resdat_t, resdat_q, align[isub_out].align_t,
                  align[isub_out].align_q);
    printallatm(fp, natom_t, allatm_t, 1);
    printallatm(fp, natom_q, allatm_q, 2);
  } else {
    printrascript(fp, naa_q, naa_t, resdat_q, resdat_t, align[isub_out].align_q,
                  align[isub_out].align_t);
    printallatm(fp, natom_q, allatm_q, 1);
    printallatm(fp, natom_t, allatm_t, 2);
  }

  return;
}

/**  function print alignment  **/
void printali(FILE *fp, int naa_q, int naa_t, RESDAT *resdat_q,
                 RESDAT *resdat_t, ALIGN *align) {
  int iaa, jaa;
  int isub;

  if (input.qtchange == OFF) {
    fprintf(fp, "# PDB file1 = %s\n", input.file_t);
    fprintf(fp, "# PDB file2 = %s\n", input.file_q);
  } else {
    fprintf(fp, "# PDB file1 = %s\n", input.file_q);
    fprintf(fp, "# PDB file2 = %s\n", input.file_t);
  }

  fprintf(fp, "#\n");
  fprintf(fp, "# Alignment No.");
  for (isub = 1; isub <= input.nsub; isub++) {
    fprintf(fp, "       %02d  ", isub);
  }
  fprintf(fp, "\n");
  fprintf(fp, "#\n");
  /** RMSD **/
  fprintf(fp, "# RMSD        ");
  for (isub = 0; isub <= input.nsub - 1; isub++) {
    fprintf(fp, "      %5.3f", align[isub].rmsd);
  }
  fprintf(fp, "\n");
  /** Aligned length **/
  fprintf(fp, "# Length      ");
  for (isub = 0; isub <= input.nsub - 1; isub++) {
    fprintf(fp, "       %4d", align[isub].naa_align);
  }
  fprintf(fp, "\n");
  /** mTM-score **/
  fprintf(fp, "# mTM-score   ");
  for (isub = 0; isub <= input.nsub - 1; isub++) {
    fprintf(fp, "     %6.4f", align[isub].TMscore_mod);
  }
  fprintf(fp, "\n");
  /** SP-score **/
  fprintf(fp, "# SP-score    ");
  for (isub = 0; isub <= input.nsub - 1; isub++) {
    fprintf(fp, "     %6.4f", align[isub].SPscore);
  }
  fprintf(fp, "\n");
  /** Dali Zscore **/
  fprintf(fp, "# Dali Zscore ");
  for (isub = 0; isub <= input.nsub - 1; isub++) {
    fprintf(fp, "     %6.3f", align[isub].DaliZ);
  }
  fprintf(fp, "\n");
  /** TM-score **/
  if (input.qtchange == OFF) {
    fprintf(fp, "# TM-score: P1 ");
    for (isub = 0; isub <= input.nsub - 1; isub++) {
      fprintf(fp, "    %6.4f ", align[isub].TMscore_t);
    }
    fprintf(fp, "\n# TM-score: P2 ");
    for (isub = 0; isub <= input.nsub - 1; isub++) {
      fprintf(fp, "    %6.4f ", align[isub].TMscore_q);
    }
    fprintf(fp, "\n#\n");
    fprintf(fp, "# Protein1    ");
    for (isub = 0; isub <= input.nsub - 1; isub++) {
      fprintf(fp, "   Protein2");
    }
    fprintf(fp, "\n#\n");

    for (iaa = 1; iaa <= naa_t; iaa++) {
      fprintf(fp, "  %4d %1c %1s    ", resdat_t[iaa].iaa_org,
	      resdat_t[iaa].reschar, resdat_t[iaa].chainID_org);
      for (isub = 0; isub <= input.nsub - 1; isub++) {
        jaa = align[isub].align_t[iaa];
        if (jaa == 0) {
          fprintf(fp, "      . . .");
        } else {
	  fprintf(fp, "   %4d %1c %1s", resdat_q[jaa].iaa_org,
                  resdat_q[jaa].reschar, resdat_q[jaa].chainID_org);
	}
      }
      fprintf(fp, "\n");
    }
  } else {
    fprintf(fp, "# TM-score: P1 ");
    for (isub = 0; isub <= input.nsub - 1; isub++) {
      fprintf(fp, "     %5.3f ", align[isub].TMscore_q);
    }
    fprintf(fp, "\n# TM-score: P2 ");
    for (isub = 0; isub <= input.nsub - 1; isub++) {
      fprintf(fp, "     %5.3f ", align[isub].TMscore_t);
    }
    fprintf(fp, "\n#\n");
    fprintf(fp, "# Protein1    ");
    for (isub = 0; isub <= input.nsub - 1; isub++) {
      fprintf(fp, "   Protein2");
    }
    fprintf(fp, "\n#\n");

    for (iaa = 1; iaa <= naa_q; iaa++) {
      fprintf(fp, "  %4d %1c %1s    ", resdat_q[iaa].iaa_org,
              resdat_q[iaa].reschar, resdat_q[iaa].chainID_org);
      for (isub = 0; isub <= input.nsub - 1; isub++) {
        jaa = align[isub].align_q[iaa];
        if (jaa == 0) {
          fprintf(fp, "      . . .");
        } else {
	  fprintf(fp, "   %4d %1c %1s", resdat_t[jaa].iaa_org,
		  resdat_t[jaa].reschar, resdat_t[jaa].chainID_org);
        }
      }
      fprintf(fp, "\n");
    }
  }

  return;
}

/**  function print translation matrix  **/
void printmat(FILE *fp, ALIGN *align) {
  int i, j, isub;

  for (isub = 0; isub <= input.nsub - 1; isub++) {

    /**  small value  **/
    for (i = 1; i <= 3; i++) {
      for (j = 1; j <= 3; j++) {
        if (fabs(align[isub].rot[i][j]) < EPS) {
          align[isub].rot[i][j] = 0.0;
        }
      }
      if (fabs(align[isub].vec[i]) < EPS) {
        align[isub].vec[i] = 0.0;
      }
    }

    /**  output  **/
    if (input.qtchange == ON) {
      fprintf(fp, " Matrix%02d  ( TM-score=%5.3f, SP-score=%5.3f )\n", isub + 1,
              align[isub].TMscore_q, align[isub].SPscore);
    } else {
      fprintf(fp, " Matrix%02d  ( TM-score=%5.3f, SP-score=%5.3f )\n", isub + 1,
              align[isub].TMscore_t, align[isub].SPscore);
    }
    fprintf(fp, " ----- Rotation matrix to rotate Chain1 to Chain2 -----\n");
    fprintf(fp, " i       t(i)           u(i,1)      u(i,2)      u(i,3)\n");
    fprintf(fp, " 1   %11.6f      %9.6f   %9.6f   %9.6f\n", align[isub].vec[1],
            align[isub].rot[1][1], align[isub].rot[1][2],
            align[isub].rot[1][3]);
    fprintf(fp, " 2   %11.6f      %9.6f   %9.6f   %9.6f\n", align[isub].vec[2],
            align[isub].rot[2][1], align[isub].rot[2][2],
            align[isub].rot[2][3]);
    fprintf(fp, " 3   %11.6f      %9.6f   %9.6f   %9.6f\n", align[isub].vec[3],
            align[isub].rot[3][1], align[isub].rot[3][2],
            align[isub].rot[3][3]);
    fprintf(fp, " ------------------------------------------------------\n");
    fprintf(fp, "\n");
  }

  return;
}

void inverse_mat(float **rot, float *vec) {
  int i, j;
  float rot_tmp[4][4];
  float vec_tmp[4];

  for (i = 1; i <= 3; i++) {
    vec_tmp[i] = vec[i];
    for (j = 1; j <= 3; j++) {
      rot_tmp[i][j] = rot[j][i];
    }
  }

  for (i = 1; i <= 3; i++) {
    vec[i] = 0;
    for (j = 1; j <= 3; j++) {
      rot[i][j] = rot_tmp[i][j];
      vec[i] -= rot_tmp[i][j] * vec_tmp[j];
    }
  }

  return;
}
