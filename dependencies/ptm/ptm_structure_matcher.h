#ifndef PTM_STRUCTURE_MATCHER_H
#define PTM_STRUCTURE_MATCHER_H

#include "ptm_initialize_data.h"
#include "ptm_constants.h"


namespace ptm {

typedef struct
{
    double rmsd;
    double scale;
    double q[4];                //rotation in quaternion form (rigid body transformation)
    int8_t mapping[PTM_MAX_POINTS];
    const refdata_t* ref_struct;
} result_t;

int match_general(const refdata_t* s, double (*ch_points)[3], double (*points)[3], convexhull_t* ch, result_t* res);
int match_fcc_hcp_ico(double (*ch_points)[3], double (*points)[3], int32_t flags, convexhull_t* ch, result_t* res);
int match_dcub_dhex(double (*ch_points)[3], double (*points)[3], int32_t flags, convexhull_t* ch, result_t* res);
int match_graphene(double (*points)[3], result_t* res);

// Scale-invariant QCP RMSD between an ideal template and a normalized environment.
// G1/G2/E0 are the precomputed Gram sums (G1=sum|ideal|^2, G2=sum|normalized|^2,
// E0=(G1+G2)/2). Writes the orientation quaternion q (w,x,y,z) and optimal scale.
double calc_rmsd(int num_points, const double (*ideal_points)[3], double (*normalized)[3], int8_t* mapping,
            double G1, double G2, double E0, double* q, double* p_scale);

}

#endif

