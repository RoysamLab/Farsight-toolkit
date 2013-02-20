#ifndef SPINECONSTS_H
#define SPINECONSTS_H
//SpineRingConsts.h
// Ring Consts 
#define spr_RADISCALE 1.1
#define spr_RADOSCALE 3.5
#define spr_RINGFITMIN 2
#define spr_SEQHITMIN  7
// Queue consts
#define spr_MAXQSIZE  7

// Spine Candidate Consts
#define spr_MINSPINESIZE 3
#define spr_LIKELIHOODTHRESH 0

//make max size 4k pixels
#define spr_MAXCANDSIZE  1<<12

// Spine Path Extraction Consts
#define spr_MINPATHSIZE  3
#define spr_PATHOVERLAPRATIO .4
#define spr_PATHVERBOSE 0
#define spr_MAXPATHITERINIT 50000
#define spr_PATHREGIONOFFSET 10
// MUOFFSET is the number of neighboring centerpoints
// (at each side) to be taken into consideration for best paths
#define MUOFFSET 10

#define spr_SPIMGDIMENSION 3
#define spr_EPSILON   1.0e-8
#define spr_CONSISTENCYTOL 1

#define spr_PXLMAXVAL  USHRT_MAX;

enum spr_SpineRingType {RingTypeCyl, RingTypeSE};
// cylindrical ring type not implemented yet
// When implementing, make USEDONUT an argument
#define spr_USEDONUT 0

// up_samplimng not implemented yet. Modify iteraions
// in RingSE ...  and make this an argument
const double spr_UP_SAMPLING[3] = {1.0,1.0,1.0};


enum spr_DebugType { no_DBG, feat_DBG, path_DBG, cand_DBG, cc_DBG, ring_DBG, im_DBG, xml_DBG} ;
#define spr_DEBUGGING no_DBG

#if spr_DEBUGGING==no_DBG
#define DEBUGSTMT(x) 
#define DEBUGSTMT2(x,y) 
#else
#define DEBUGSTMT(x) x
#define DEBUGSTMT2(x,y) x y
#endif

#endif