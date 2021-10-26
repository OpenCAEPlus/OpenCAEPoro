#include "LinearSolver.hpp"


void LinearSolver::decoupling(dBSRmat* Absr, dvector* b, int scal_type, dBSRmat* Asc,
    dvector* fsc, ivector* order, double* Dmatvec, int decouple_type) {
    int nrow = Absr->ROW;
    int nb = Absr->nb;
    double* Dmat = Dmatvec;
    precond_diag_bsr diagA;

    switch (decouple_type) {
    case 2:
        decouple_anl(Absr, Dmat, Asc); break;
    case 3:
        decouple_quasi(Absr, Dmat, Asc); break;
    case 4:
        decouple_trueabf(Absr, Dmat, Asc); break;
    case 5:
        decouple_abftrue(Absr, Dmat, Asc); break;
    case 6:
        decouple_abftrue(Absr, Dmat, Asc); break;
    case 7:
        decouple_truetrans_alg(Absr, Dmat, Asc); break;
    case 8:
        decouple_truetrans(Absr, Dmat, Asc); break;
    case 9:
        decouple_true_scale(Absr, Dmat, Asc); break;
    case 10:
        decouple_rotate(Absr, Dmat, Asc); break;
    default: // case 1:
        decouple_abf(Absr, Dmat, Asc);
    }

    diagA.diag.row = nrow * nb * nb;
    diagA.diag.val = Dmat;
    diagA.nb = nb;
    fasp_precond_dbsr_diag(b->val, fsc->val, &diagA);

    // ordering
    for (int i = 0; i < nrow; ++i) order->val[i] = i;
}

// ABF decoupling method
void LinearSolver::decouple_abf(dBSRmat* A,
    REAL* diaginv,
    dBSRmat* B) {
    // members of A
    const INT ROW = A->ROW;
    const INT NNZ = A->NNZ;
    const INT nb = A->nb;
    const INT nb2 = nb * nb;
    const INT* IA = A->IA;
    const INT* JA = A->JA;
    REAL* val = A->val;

    INT* IAb = NULL;
    INT* JAb = NULL;
    REAL* valb = NULL;

    INT i, k, m, j;

    // Create a link to dBSRmat 'B'
    IAb = B->IA;
    JAb = B->JA;
    valb = B->val;
    memcpy(IAb, IA, (ROW + 1) * sizeof(INT));
    memcpy(JAb, JA, NNZ * sizeof(INT));

    for (i = 0; i < ROW; ++i) {
        // get the diagonal sub-blocks D and their inverse
        for (k = IA[i]; k < IA[i + 1]; ++k) {
            if (JA[k] == i) {
                memcpy(diaginv + i * nb2, val + k * nb2, nb2 * sizeof(REAL));
                fasp_smat_inv(diaginv + i * nb2, nb);
            }
        }
        // compute D^{-1}*A
        for (k = IA[i]; k < IA[i + 1]; ++k) {
            m = k * nb2;
            j = JA[k];
            fasp_blas_smat_mul(diaginv + i * nb2, val + m, valb + m, nb);
        }
    } // end of main loop
}

// Analytical decoupling method
void LinearSolver::decouple_anl(dBSRmat* A,
    REAL* diaginv,
    dBSRmat* B) {
    // members of A
    const INT ROW = A->ROW;
    const INT NNZ = A->NNZ;
    const INT nb = A->nb;
    const INT nb2 = nb * nb;
    const INT* IA = A->IA;
    const INT* JA = A->JA;
    REAL* val = A->val;

    INT* IAb = NULL;
    INT* JAb = NULL;
    REAL* valb = NULL;

    INT i, k, m, j;

    // Create a dBSRmat 'B'
    IAb = B->IA;
    JAb = B->JA;
    valb = B->val;
    memcpy(IAb, IA, (ROW + 1) * sizeof(INT));
    memcpy(JAb, JA, NNZ * sizeof(INT));

    for (i = 0; i < ROW; ++i) {
        // form the diagonal sub-blocks for analytical decoupling
        for (k = IA[i]; k < IA[i + 1]; ++k) {
            if (JA[k] == i) {
                m = k * nb2;
                fasp_smat_identity(diaginv + i * nb2, nb, nb2);
                for (int l = 0; l < nb - 1; l++)
                    diaginv[i * nb2 + 1 + l] = -val[m + 1 + l];
            }
        }

        // compute D^{-1}*A
        for (k = IA[i]; k < IA[i + 1]; ++k) {
            m = k * nb2;
            j = JA[k];
            fasp_blas_smat_mul(diaginv + i * nb2, val + m, valb + m, nb);
        }
    } // end of main loop
}

void LinearSolver::decouple_truetrans(dBSRmat* A,
    REAL* diaginv,
    dBSRmat* B) {
    // members of A
    const INT ROW = A->ROW;
    const INT NNZ = A->NNZ;
    const INT nb = A->nb;
    const INT nb2 = nb * nb;
    const INT* IA = A->IA;
    const INT* JA = A->JA;
    REAL* val = A->val;

    INT* IAb = NULL;
    INT* JAb = NULL;
    REAL* valb = NULL;

    INT i, k, m, j;

    // Create a dBSRmat 'B'
    IAb = B->IA;
    JAb = B->JA;
    valb = B->val;
    memcpy(IAb, IA, (ROW + 1) * sizeof(INT));
    memcpy(JAb, JA, NNZ * sizeof(INT));

    double* mat1, * mat2;
    mat1 = new double[nb2];
    mat2 = new double[nb2];

    //Dset0(nb2*ROW, diaginv);
    for (i = 0; i < ROW; ++i) {
        double Tt = 0;
        // get the diagonal sub-blocks
        //mat2 is true-IMPES matrix.
        fasp_smat_identity(mat2, nb, nb2);
        //mat1 is the mobility ratio matrix
        fasp_smat_identity(mat1, nb, nb2);

        for (int l = 0; l < nb - 1; l++) {
            mat2[1 + l] = diaginv[i * nb2 + 1 + l];
            Tt += diaginv[i * nb2 + 1 + l] * diaginv[i * nb2 + (l + 1) * nb];
        }
        for (k = IA[i]; k < IA[i + 1]; ++k) {
            if (JA[k] == i) {
                break;
            }
        }
        for (int l = 0; l < nb - 1; l++) {
            if (val[k * nb2 + (l + 1) * nb] > 0)
                mat1[(l + 1) * nb] = -diaginv[i * nb2 + (l + 1) * nb] / Tt;//diag(l)
        }
        DaABpbC(nb, nb, nb, 1, mat1, mat2, 0, diaginv + i * nb2);

        // compute D^{-1}*A
        for (k = IA[i]; k < IA[i + 1]; ++k) {
            m = k * nb2;
            j = JA[k];
            fasp_blas_smat_mul(diaginv + i * nb2, val + m, valb + m, nb);
        }
    } // end of main loop
    delete mat1;
    delete mat2;

}

void LinearSolver::decouple_truetrans_alg(dBSRmat* A,
    REAL* diaginv,
    dBSRmat* B) {
    // members of A
    const INT ROW = A->ROW;
    const INT NNZ = A->NNZ;
    const INT nb = A->nb;
    const INT nb2 = nb * nb;
    const INT* IA = A->IA;
    const INT* JA = A->JA;
    REAL* val = A->val;

    INT* IAb = NULL;
    INT* JAb = NULL;
    REAL* valb = NULL;

    INT i, k, m, j;

    // Create a link to dBSRmat 'B'
    IAb = B->IA;
    JAb = B->JA;
    valb = B->val;

    memcpy(IAb, IA, (ROW + 1) * sizeof(INT));
    memcpy(JAb, JA, NNZ * sizeof(INT));
    double* mat1, * mat2;
    mat1 = new double[nb2];
    mat2 = new double[nb2];

    //Dset0(nb2*ROW, diaginv);
    for (i = 0; i < ROW; ++i) {
        for (k = IA[i]; k < IA[i + 1]; ++k) {
            double Tt = 0;
            if (JA[k] == i) {
                m = k * nb2;
                fasp_smat_identity(mat2, nb, nb2);
                for (int l = 0; l < nb - 1; l++) {
                    mat2[1 + l] = -val[m + 1 + l];
                    if (val[m + (l + 1) * nb] > 0) {
                        Tt += -val[m + 1 + l] * val[m + (l + 1) * nb];
                    }
                }

                fasp_smat_identity(mat1, nb, nb2);
                for (int l = 0; l < nb - 1; l++) {
                    if (val[m + (l + 1) * nb] > 0) {
                        mat1[(l + 1) * nb] = -val[m + (l + 1) * nb] / Tt;//diag(l)
                    }
                }

                DaABpbC(nb, nb, nb, 1, mat1, mat2, 0, diaginv + i * nb2);
            }
        }

        // compute D^{-1}*A
        for (k = IA[i]; k < IA[i + 1]; ++k) {
            m = k * nb2;
            j = JA[k];
            fasp_blas_smat_mul(diaginv + i * nb2, val + m, valb + m, nb);
        }
    } // end of main loop
    delete mat1;
    delete mat2;

}

// Semi-analytical decoupling method
void LinearSolver::decouple_abftrue(dBSRmat* A,
    REAL* diaginv,
    dBSRmat* B) {
    // members of A
    const INT ROW = A->ROW;
    const INT NNZ = A->NNZ;
    const INT nb = A->nb;
    const INT nb2 = nb * nb;
    const INT* IA = A->IA;
    const INT* JA = A->JA;
    REAL* val = A->val;

    INT* IAb = NULL;
    INT* JAb = NULL;
    REAL* valb = NULL;

    INT i, k, m, j;

    // Create a dBSRmat 'B'
    IAb = B->IA;
    JAb = B->JA;
    valb = B->val;
    memcpy(IAb, IA, (ROW + 1) * sizeof(INT));
    memcpy(JAb, JA, NNZ * sizeof(INT));

    for (i = 0; i < ROW; ++i) {
        // get the diagonal sub-blocks
        for (k = IA[i]; k < IA[i + 1]; ++k) {
            if (JA[k] == i) {
                // Form the ABF part first
                memcpy(diaginv + i * nb2, val + k * nb2, nb2 * sizeof(REAL));
                fasp_smat_inv(diaginv + i * nb2, nb);
                // Replace the first line with analytical decoupling
                m = k * nb2;
                diaginv[i * nb2] = 1;
                for (int l = 0; l < nb - 1; l++)
                    diaginv[i * nb2 + 1 + l] = -val[m + 1 + l];
            }
        }

        // compute D^{-1}*A
        for (k = IA[i]; k < IA[i + 1]; ++k) {
            m = k * nb2;
            j = JA[k];
            fasp_blas_smat_mul(diaginv + i * nb2, val + m, valb + m, nb);
        }
    } // end of main loop
}

void LinearSolver::decouple_true_scale(dBSRmat* A,
    REAL* diaginv,
    dBSRmat* B) {
    // members of A
    const INT ROW = A->ROW;
    const INT NNZ = A->NNZ;
    const INT nb = A->nb;
    const INT nb2 = nb * nb;
    const INT* IA = A->IA;
    const INT* JA = A->JA;
    REAL* val = A->val;

    INT* IAb = NULL;
    INT* JAb = NULL;
    REAL* valb = NULL;

    INT i, k, m, j;

    // Variables for OpenMP

    // Create a dBSRmat 'B'

    IAb = B->IA;
    JAb = B->JA;
    valb = B->val;

    memcpy(IAb, IA, (ROW + 1) * sizeof(INT));
    memcpy(JAb, JA, NNZ * sizeof(INT));

    for (i = 0; i < ROW; ++i) {
        // get the diagonal sub-blocks
        for (k = IA[i]; k < IA[i + 1]; ++k) {
            if (JA[k] == i) {
                m = k * nb2;
                fasp_smat_identity(diaginv + i * nb2, nb, nb2);
                diaginv[i * nb2] = 1 / val[m];
                for (int l = 0; l < nb - 1; l++)
                    diaginv[i * nb2 + 1 + l] = -val[m + 1 + l] / val[m];
            }
        }

        // compute D^{-1}*A
        for (k = IA[i]; k < IA[i + 1]; ++k) {
            m = k * nb2;
            j = JA[k];
            fasp_blas_smat_mul(diaginv + i * nb2, val + m, valb + m, nb);
        }
    } // end of main loop
}

void LinearSolver::decouple_rotate(dBSRmat* A,
    REAL* diaginv,
    dBSRmat* B) {
    // members of A
    const INT ROW = A->ROW;
    const INT NNZ = A->NNZ;
    const INT nb = A->nb;
    const INT nb2 = nb * nb;
    const INT* IA = A->IA;
    const INT* JA = A->JA;
    REAL* val = A->val;

    INT* IAb = NULL;
    INT* JAb = NULL;
    REAL* valb = NULL;

    INT i, k, m, j;

    // Variables for OpenMP

    // Create a dBSRmat 'B'

    IAb = B->IA;
    JAb = B->JA;
    valb = B->val;

    memcpy(IAb, IA, (ROW + 1) * sizeof(INT));
    memcpy(JAb, JA, NNZ * sizeof(INT));

    for (i = 0; i < ROW; ++i) {
        // get the diagonal sub-blocks
        for (k = IA[i]; k < IA[i + 1]; ++k) {
            if (JA[k] == i) {
                m = k * nb2;
                for (int li = 0; li < nb; li++) {
                    for (int lj = 0; lj < nb; lj++) {
                        diaginv[i * nb2 + li * nb + lj] = 0;
                        if (lj - li == 1) {
                            diaginv[i * nb2 + li * nb + lj] = 1;
                        }
                        if (lj == 0 && li == nb - 1) {
                            diaginv[i * nb2 + li * nb + lj] = 1;
                        }
                    }
                }
            }
        }

        // compute D^{-1}*A
        for (k = IA[i]; k < IA[i + 1]; ++k) {
            m = k * nb2;
            j = JA[k];
            fasp_blas_smat_mul(diaginv + i * nb2, val + m, valb + m, nb);
        }
    } // end of main loop
}

// Quasi-IMPES decoupling method
void LinearSolver::decouple_quasi(dBSRmat* A,
    REAL* diaginv,
    dBSRmat* B) {
    // members of A
    const INT ROW = A->ROW;
    const INT NNZ = A->NNZ;
    const INT nb = A->nb;
    const INT nb2 = nb * nb;
    const INT* IA = A->IA;
    const INT* JA = A->JA;
    REAL* val = A->val;

    INT* IAb = NULL;
    INT* JAb = NULL;
    REAL* valb = NULL;

    INT i, k, m, j;

    IAb = B->IA;
    JAb = B->JA;
    valb = B->val;
    memcpy(IAb, IA, (ROW + 1) * sizeof(INT));
    memcpy(JAb, JA, NNZ * sizeof(INT));

    for (i = 0; i < ROW; ++i) {
        // get the diagonal sub-blocks
        for (k = IA[i]; k < IA[i + 1]; ++k) {
            if (JA[k] == i) {
                m = k * nb2;
                fasp_smat_identity(diaginv + i * nb2, nb, nb2);
                for (int l = 0; l < nb - 1; l++)
                    diaginv[i * nb2 + 1 + l] = -val[m + 1 + l] / val[m + (l + 1) * nb + l + 1];
            }
        }

        // compute D^{-1}*A
        for (k = IA[i]; k < IA[i + 1]; ++k) {
            m = k * nb2;
            j = JA[k];
            fasp_blas_smat_mul(diaginv + i * nb2, val + m, valb + m, nb);
        }
    } // end of main loop
}

// What's the difference with decouple_anl?
void LinearSolver::decouple_trueabf(dBSRmat* A,
    REAL* diaginv,
    dBSRmat* B) {
    // members of A
    const INT ROW = A->ROW;
    const INT NNZ = A->NNZ;
    const INT nb = A->nb;
    const INT nb2 = nb * nb;
    const INT* IA = A->IA;
    const INT* JA = A->JA;
    REAL* val = A->val;

    INT* IAb = NULL;
    INT* JAb = NULL;
    REAL* valb = NULL;

    INT i, k, m, j;

    IAb = B->IA;
    JAb = B->JA;
    valb = B->val;
    memcpy(IAb, IA, (ROW + 1) * sizeof(INT));
    memcpy(JAb, JA, NNZ * sizeof(INT));

    for (i = 0; i < ROW; ++i) {
        // get the diagonal sub-blocks
        for (k = IA[i]; k < IA[i + 1]; ++k) {
            if (JA[k] == i) {
                m = k * nb2;
                fasp_smat_identity(diaginv + i * nb2, nb, nb2);
                for (int l = 0; l < nb; l++) diaginv[i * nb2 + l] = val[m + l];
                fasp_smat_inv(diaginv + i * nb2, nb);
            }
        }

        // compute D^{-1}*A
        for (k = IA[i]; k < IA[i + 1]; ++k) {
            m = k * nb2;
            j = JA[k];
            fasp_blas_smat_mul(diaginv + i * nb2, val + m, valb + m, nb);
        }
    } // end of main loop
}

void LinearSolver::decouple_rowsum(dBSRmat* A,
    REAL* diaginv,
    dBSRmat* B) {
    // members of A
    const INT ROW = A->ROW;
    const INT NNZ = A->NNZ;
    const INT nb = A->nb;
    const INT nb2 = nb * nb;
    const INT* IA = A->IA;
    const INT* JA = A->JA;
    REAL* val = A->val;

    INT* IAb = NULL;
    INT* JAb = NULL;
    REAL* valb = NULL;

    INT i, k, m, j;

    // Variables for OpenMP

    IAb = B->IA;
    JAb = B->JA;
    valb = B->val;

    memcpy(IAb, IA, (ROW + 1) * sizeof(INT));
    memcpy(JAb, JA, NNZ * sizeof(INT));

    for (i = 0; i < ROW; ++i) {
        // get the diagonal sub-blocks
        for (k = IA[i]; k < IA[i + 1]; ++k) {
            if (JA[k] == i) {
                m = k * nb2;
                fasp_smat_identity(diaginv + i * nb2, nb, nb2);
                for (int l = 0; l < nb - 1; l++)
                    diaginv[i * nb2 + 1 + l] = 1;
            }
        }

        // compute D^{-1}*A
        for (k = IA[i]; k < IA[i + 1]; ++k) {
            m = k * nb2;
            j = JA[k];
            fasp_blas_smat_mul(diaginv + i * nb2, val + m, valb + m, nb);
        }
    } // end of main loop
}