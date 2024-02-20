#ifndef CM_PBC_HPP
#define CM_PBC_HPP

#include <iostream>

#define DB
#ifdef DB
#define OUT(x) std::cout << x << std::endl
#else
#define OUT(x)
#endif

static void pr_indent(FILE *fp, int indent)
{
    for (int i = 0; i < indent; i++)
    {
        fprintf(fp, " ");
    }
}
bool available(FILE* fp, const void* p, int indent, const char* title)
{
    if (!p)
    {
        if (indent > 0)
        {
            pr_indent(fp, indent);
        }
        fprintf(fp, "%s: not available\n", title);
    }
    return (p != nullptr);
}


int pr_title_nxn(FILE* fp, int indent, const char* title, int n1, int n2)
{
    pr_indent(fp, indent);
    fprintf(fp, "%s (%dx%d):\n", title, n1, n2);
    return (indent );
}
void pr_rvecs(FILE* fp, int indent, const char* title, const rvec vec[], int n)
{
    const char* fshort = "%12.5e";
    const char* flong  = "%15.8e";
    const char* format;
    int         i, j;

    if (getenv("GMX_PRINT_LONGFORMAT") != nullptr)
    {
        format = flong;
    }
    else
    {
        format = fshort;
    }

    if (available(fp, vec, indent, title))
    {
        indent = pr_title_nxn(fp, indent, title, n, DIM);
        for (i = 0; i < n; i++)
        {
            pr_indent(fp, indent);
            fprintf(fp, "%s[%5d]={", title, i);
            for (j = 0; j < DIM; j++)
            {
                if (j != 0)
                {
                    fprintf(fp, ", ");
                }
                fprintf(fp, format, vec[i][j]);
            }
            fprintf(fp, "}\n");
        }
    }
}
enum
{
    epbcdxRECTANGULAR = 1,
    epbcdxTRICLINIC,
    epbcdx2D_RECT,
    epbcdx2D_TRIC,
    epbcdx1D_RECT,
    epbcdx1D_TRIC,
    epbcdxSCREW_RECT,
    epbcdxSCREW_TRIC,
    epbcdxNOPBC,
    epbcdxUNSUPPORTED
};

FILE *debug = stderr;

static float sc_skewnessMargin = 1.1;

const char* check_box(PbcType pbcType, const matrix box)
{
    const char* ptr;
    const float sc_boxSkewnessMarginForWarning = 1.1;

    if (pbcType == PbcType::Unset)
    {
        pbcType = guessPbcType(box);
    }

    if (pbcType == PbcType::No)
    {
        return nullptr;
    }

    GMX_ASSERT(box != nullptr, "check_box requires a valid box unless pbcType is No");

    if (pbcType == PbcType::Xyz && box[XX][XX] == 0 && box[YY][YY] == 0 && box[ZZ][ZZ] == 0)
    {
        ptr = "Empty diagonal for a 3-dimensional periodic box";
    }
    else if (pbcType == PbcType::XY && box[XX][XX] == 0 && box[YY][YY] == 0)
    {
        ptr = "Empty diagonal for a 2-dimensional periodic box";
    }
    else if ((box[XX][YY] != 0) || (box[XX][ZZ] != 0) || (box[YY][ZZ] != 0))
    {
        ptr = "Only triclinic boxes with the first vector parallel to the x-axis and the second "
              "vector in the xy-plane are supported.";
    }
    else if (pbcType == PbcType::Screw && (box[YY][XX] != 0 || box[ZZ][XX] != 0))
    {
        ptr = "The unit cell can not have off-diagonal x-components with screw pbc";
    }
    else if (std::fabs(box[YY][XX]) > sc_boxSkewnessMarginForWarning * 0.5_real * box[XX][XX]
             || (pbcType != PbcType::XY
                 && (std::fabs(box[ZZ][XX]) > sc_boxSkewnessMarginForWarning * 0.5_real * box[XX][XX]
                     || std::fabs(box[ZZ][YY]) > sc_boxSkewnessMarginForWarning * 0.5_real * box[YY][YY])))
    {
        ptr = "Triclinic box is too skewed.";
    }
    else
    {
        ptr = nullptr;
    }

    return ptr;
}

template <typename tp, typename t, typename iv, typename mat>
void set_porco(tp* pbc, t pbcType, const iv dd_pbc, const mat box)
{
    std::cout << "set_porco" << std::endl;
    std::cout << "dd_pbc: " << dd_pbc[0] << std::endl;
    int         order[3] = { 0, -1, 1 };
    iv        bPBC;
    const char* ptr;

    pbc->pbcType   = pbcType;
    pbc->ndim_ePBC = numPbcDimensions(pbcType);

    if (pbc->pbcType == t::No)
    {
        pbc->pbcTypeDX = epbcdxNOPBC;

        return;
    }

    copy_mat(box, pbc->box);
    pbc->max_cutoff2 = 0;
    pbc->dim         = -1;
    pbc->ntric_vec   = 0;

    for (int i = 0; (i < DIM); i++)
    {
        pbc->fbox_diag[i]  = box[i][i];
        pbc->hbox_diag[i]  = pbc->fbox_diag[i] * 0.5_real;
        pbc->mhbox_diag[i] = -pbc->hbox_diag[i];
    }
    OUT("checkpoint 1");
    ptr = check_box(pbcType, box);
    if (!ptr) std::cout << "ptr is null" << std::endl;
    if (ptr)
    {
        OUT("checkbox ptr is set");
        fprintf(stderr, "Warning: %s\n", ptr);
        pr_rvecs(stderr, 0, "         Box", box, DIM);
        fprintf(stderr, "         Can not fix pbc.\n\n");
        pbc->pbcTypeDX = epbcdxUNSUPPORTED;
    }
    else
    {
        if (pbcType == t::Screw && nullptr != dd_pbc)
        {
            /* This combinated should never appear here */
            printf("low_set_pbc called with screw pbc and dd_nc != NULL");
        }

        int npbcdim = 0;
        if (bPBC == nullptr) std::cout << "bPBC is null" << std::endl;
        else std::cout << "bPBC is not null" << std::endl;
        
        for (int i = 0; i < DIM; i++)
        {
            if ((dd_pbc && dd_pbc[i] == 0) || (pbcType == t::XY && i == ZZ))
            {
                bPBC[i] = 0;
            }
            else
            {
                bPBC[i] = 1;
                npbcdim++;
            }
        }
        OUT("MORTACI");
        switch (npbcdim)
        {
            case 1:
                OUT("npbcdim is 1");
                /* 1D pbc is not an mdp option and it is therefore only used
                 * with single shifts.
                 */
                pbc->pbcTypeDX = epbcdx1D_RECT;
                for (int i = 0; i < DIM; i++)
                {
                    if (bPBC[i])
                    {
                        pbc->dim = i;
                    }
                }
                GMX_ASSERT(pbc->dim < DIM, "Dimension for PBC incorrect");
                for (int i = 0; i < pbc->dim; i++)
                {
                    if (pbc->box[pbc->dim][i] != 0)
                    {
                        pbc->pbcTypeDX = epbcdx1D_TRIC;
                    }
                }
                break;
            case 2:
                OUT("npbcdim is 2");
                pbc->pbcTypeDX = epbcdx2D_RECT;
                for (int i = 0; i < DIM; i++)
                {
                    if (!bPBC[i])
                    {
                        pbc->dim = i;
                    }
                }
                for (int i = 0; i < DIM; i++)
                {
                    if (bPBC[i])
                    {
                        for (int j = 0; j < i; j++)
                        {
                            if (pbc->box[i][j] != 0)
                            {
                                pbc->pbcTypeDX = epbcdx2D_TRIC;
                            }
                        }
                    }
                }
                break;
            case 3:
                OUT("npbcdim is 3");
                if (pbcType != t::Screw)
                {
                    if (TRICLINIC(box))
                    {
                        pbc->pbcTypeDX = epbcdxTRICLINIC;
                    }
                    else
                    {
                        pbc->pbcTypeDX = epbcdxRECTANGULAR;
                    }
                }
                else
                {
                    pbc->pbcTypeDX = (box[ZZ][YY] == 0 ? epbcdxSCREW_RECT : epbcdxSCREW_TRIC);
                    if (pbc->pbcTypeDX == epbcdxSCREW_TRIC)
                    {
                        fprintf(stderr,
                                "Screw pbc is not yet implemented for triclinic boxes.\n"
                                "Can not fix pbc.\n");
                        pbc->pbcTypeDX = epbcdxUNSUPPORTED;
                    }
                }
                break;
            default: printf("Incorrect number of pbc dimensions with DD: %d", npbcdim); exit(1);
        }
        pbc->max_cutoff2 = max_cutoff2(pbcType, box);
        OUT("max_cutoff2 is set");
        if (pbc->pbcTypeDX == epbcdxTRICLINIC || pbc->pbcTypeDX == epbcdx2D_TRIC
            || pbc->pbcTypeDX == epbcdxSCREW_TRIC)
        {
            if (debug)
            {
                pr_rvecs(debug, 0, "Box", box, DIM);
                fprintf(debug, "max cutoff %.3f\n", sqrt(pbc->max_cutoff2));
            }
            /* We will only need single shifts here */
            for (int kk = 0; kk < 3; kk++)
            {
                int k = order[kk];
                if (!bPBC[ZZ] && k != 0)
                {
                    continue;
                }
                for (int jj = 0; jj < 3; jj++)
                {
                    int j = order[jj];
                    if (!bPBC[YY] && j != 0)
                    {
                        continue;
                    }
                    for (int ii = 0; ii < 3; ii++)
                    {
                        int i = order[ii];
                        if (!bPBC[XX] && i != 0)
                        {
                            continue;
                        }
                        /* A shift is only useful when it is trilinic */
                        if (j != 0 || k != 0)
                        {
                            rvec trial;
                            rvec pos;
                            real d2old = 0;
                            real d2new = 0;

                            for (int d = 0; d < DIM; d++)
                            {
                                trial[d] = i * box[XX][d] + j * box[YY][d] + k * box[ZZ][d];
                                /* Choose the vector within the brick around 0,0,0 that
                                 * will become the shortest due to shift try.
                                 */
                                if (d == pbc->dim)
                                {
                                    trial[d] = 0;
                                    pos[d]   = 0;
                                }
                                else
                                {
                                    if (trial[d] < 0)
                                    {
                                        pos[d] = std::min(pbc->hbox_diag[d], -trial[d]);
                                    }
                                    else
                                    {
                                        pos[d] = std::max(-pbc->hbox_diag[d], -trial[d]);
                                    }
                                }
                                d2old += gmx::square(pos[d]);
                                d2new += gmx::square(pos[d] + trial[d]);
                            }
                            if (sc_skewnessMargin * d2new < d2old)
                            {
                                /* Check if shifts with one box vector less do better */
                                gmx_bool bUse = TRUE;
                                for (int dd = 0; dd < DIM; dd++)
                                {
                                    int shift = (dd == 0 ? i : (dd == 1 ? j : k));
                                    if (shift)
                                    {
                                        real d2new_c = 0;
                                        for (int d = 0; d < DIM; d++)
                                        {
                                            d2new_c += gmx::square(pos[d] + trial[d] - shift * box[dd][d]);
                                        }
                                        if (d2new_c <= sc_skewnessMargin * d2new)
                                        {
                                            bUse = FALSE;
                                        }
                                    }
                                }
                                if (bUse)
                                {
                                    /* Accept this shift vector. */
                                    if (pbc->ntric_vec >= MAX_NTRICVEC)
                                    {
                                        fprintf(stderr,
                                                "\nWARNING: Found more than %d triclinic "
                                                "correction vectors, ignoring some.\n"
                                                "  There is probably something wrong with your "
                                                "box.\n",
                                                MAX_NTRICVEC);
                                        pr_rvecs(stderr, 0, "         Box", box, DIM);
                                    }
                                    else
                                    {
                                        copy_rvec(trial, pbc->tric_vec[pbc->ntric_vec]);
                                        pbc->tric_shift[pbc->ntric_vec][XX] = i;
                                        pbc->tric_shift[pbc->ntric_vec][YY] = j;
                                        pbc->tric_shift[pbc->ntric_vec][ZZ] = k;
                                        pbc->ntric_vec++;

                                        if (debug)
                                        {
                                            fprintf(debug,
                                                    "  tricvec %2d = %2d %2d %2d  %5.2f %5.2f  "
                                                    "%5.2f %5.2f %5.2f  %5.2f %5.2f %5.2f\n",
                                                    pbc->ntric_vec,
                                                    i,
                                                    j,
                                                    k,
                                                    sqrt(d2old),
                                                    sqrt(d2new),
                                                    trial[XX],
                                                    trial[YY],
                                                    trial[ZZ],
                                                    pos[XX],
                                                    pos[YY],
                                                    pos[ZZ]);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    printf("done diocane\n");
}

#endif // CM_PBC_HPP