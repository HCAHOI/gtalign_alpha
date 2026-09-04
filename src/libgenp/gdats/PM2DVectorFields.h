/***************************************************************************
 *   Copyright (C) 2021-2025 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __PM2DVectorFields_h__
#define __PM2DVectorFields_h__

#include <stdlib.h>

// #include "libutil/mybase.h"

// for the case when FPTYPE==INTYPE
// #define FPTYPEisINTYPE
// #define LNTYPEisINTYPE
// #define CHTYPEisINTYPE

//NOTE: for alignment and proper data type conversion by CUDA in a device,
// make all types of the same size

#define FPTYPE      float
#define SZFPTYPE    (sizeof(FPTYPE))
#define FP0 0.0f

#define LNTYPE      unsigned int
#define SZLNTYPE    (sizeof(LNTYPE))

#define INTYPE      int //short
#define SZINTYPE    (sizeof(INTYPE))

#ifdef CHTYPEisINTYPE
#   define CHTYPE   int
#else
#   define CHTYPE   char
#endif
#define SZCHTYPE    (sizeof(CHTYPE))

//make integer variable by combining 1-byte residue insertion code, 
// chain id and serial number (order):
#define PM2D_MAKEINT_Ins_Ch_Ord(Ins, Ch, Ord) \
    ( (((Ins)&0xff)<<16) | (((Ch)&0xff)<<8) | ((Ord)&0xff) )
//decompose:
#define PM2D_GET_INSCODE(Ins_Ch_Ord)  ( ((Ins_Ch_Ord)>>16)&0xff )
#define PM2D_GET_CHID(Ins_Ch_Ord)  ( ((Ins_Ch_Ord)>>8)&0xff )
#define PM2D_GET_CHORD(Ins_Ch_Ord)  ( (Ins_Ch_Ord)&0xff )

//secondary structure assignment macros
#define pmvLOOP ' '
#define pmvSTRND 'e'
#define pmvHELIX 'h'
#define pmvTURN 't'

//nucleic acid secondary structure assignment macros
#define pmnasOPEN '('
#define pmnasCLOSE ')'
#define pmnasUNPAIRED '.'

#if defined(__CUDA_ARCH__)
__host__ __device__ __forceinline__
#else
/**
 * @brief 在结构数据读取与布局中处理 `isNASS` 对应的数据。
 * @param sss 供该函数读取或更新的 `sss` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline
#endif
bool isNASS(const char sss) {return (sss >= pmnasOPEN) && (sss <= pmnasUNPAIRED);}

#if defined(__CUDA_ARCH__)
__host__ __device__ __forceinline__
#else
/**
 * @brief 在结构数据读取与布局中处理 `isloop` 对应的数据。
 * @param sss 供该函数读取或更新的 `sss` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline
#endif
bool isloop(const char sss) {return (sss == pmvTURN) || (sss == pmvLOOP);}

#if defined(__CUDA_ARCH__)
__host__ __device__ __forceinline__
#else
/**
 * @brief 在结构数据读取与布局中处理 `helix_strnd` 对应的数据。
 * @param ss1 供该函数读取或更新的 `ss1` 参数。
 * @param ss2 供该函数读取或更新的 `ss2` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline
#endif
bool helix_strnd(const char ss1, const char ss2) {
    return (ss1 == pmvHELIX && ss2 == pmvSTRND) ||
        (ss2 == pmvHELIX && ss1 == pmvSTRND);
}


//molecule types
enum TGTMoleculeTypes {
        gtmtProtein,//protein
        gtmtNA,//nucleic acid
        gtmtNMoleculeTypes
};

// GetMoleculeType: get molecule type 
#if defined(__CUDA_ARCH__)
__host__ __device__ __forceinline__
#else
/**
 * @brief 在结构数据读取与布局中读取 `GetMoleculeType` 对应的数据。
 * @param typefieldvalue 供该函数读取或更新的 `typefieldvalue` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline
#endif
int GetMoleculeType(int typefieldvalue) {
    if(0 < typefieldvalue) return gtmtNA;
    return gtmtProtein;
}

//{{GetWatermarkedOnType/[Cmb]WatermarkOnType: write/return a 
// watermark based on molecule type
/**
 * @brief 在结构数据读取与布局中读取 `GetWatermarkedOnType` 对应的数据。
 * @param moltype 供该函数读取或更新的 `moltype` 参数。
 * @param value 需要读取、写入或转换的值。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline
int GetWatermarkedOnType(int moltype, int value) {
    if(moltype == gtmtProtein) return -value;
    // if(moltype == gtmtNA) return value;
    return value;
}

/**
 * @brief 在结构数据读取与布局中处理 `WatermarkOnType` 对应的数据。
 * @param moltype 供该函数读取或更新的 `moltype` 参数。
 * @param value 需要读取、写入或转换的值。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
inline
void WatermarkOnType(int moltype, int& value) {
    value = abs(value);
    if(moltype == gtmtProtein) value = -value;
    // if(moltype == gtmtNA);
}

/**
 * @brief 在结构数据读取与布局中处理 `CmbWatermarkOnType` 对应的数据。
 * @param moltype 供该函数读取或更新的 `moltype` 参数。
 * @param addedvalue 供该函数读取或更新的 `addedvalue` 参数。
 * @param value 需要读取、写入或转换的值。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
inline
void CmbWatermarkOnType(int moltype, int addedvalue, int* value) {
    if(moltype == gtmtProtein) *value += -addedvalue;
    // if(moltype == gtmtNA);
    else *value += addedvalue;
}
//}}

// GetMoleculeType: get molecule type 
#if defined(__CUDA_ARCH__)
__host__ __device__ __forceinline__
#else
/**
 * @brief 在结构数据读取与布局中处理 `MoleculeTypesCompatible` 对应的数据。
 * @param typefieldvalue1 供该函数读取或更新的 `typefieldvalue1` 参数。
 * @param typefieldvalue2 供该函数读取或更新的 `typefieldvalue2` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline
#endif
int MoleculeTypesCompatible(int typefieldvalue1, int typefieldvalue2) {
    return GetMoleculeType(typefieldvalue1) == GetMoleculeType(typefieldvalue2);
}

// GetComplexTypeLength: get the length of dominating residues in a complex
#if defined(__CUDA_ARCH__)
__host__ __device__ __forceinline__
#else
/**
 * @brief 在结构数据读取与布局中读取 `GetComplexTypeLength` 对应的数据。
 * @param typefieldvalue 供该函数读取或更新的 `typefieldvalue` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline
#endif
int GetComplexTypeLength(int typefieldvalue) {
    return abs(typefieldvalue);
}

// GetComplexResultantType: get complex resultant type and calculate
// resultant length based on current data and component molecule's
// type and #aligned residues
#if defined(__CUDA_ARCH__)
__host__ __device__ __forceinline__
#else
/**
 * @brief 在结构数据读取与布局中读取 `GetComplexResultantType` 对应的数据。
 * @param cpxtype 供该函数读取或更新的 `cpxtype` 参数。
 * @param cpxlen 控制当前步骤范围或规模的 `cpxlen`。
 * @param moltype 供该函数读取或更新的 `moltype` 参数。
 * @param mollen 控制当前步骤范围或规模的 `mollen`。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
inline
#endif
void GetComplexResultantType(int& cpxtype, int& cpxlen, int moltype, int mollen) {
    if(moltype == gtmtProtein) cpxlen -= mollen;
    // if(moltype == gtmtNA);
    else cpxlen += mollen;
    cpxtype = GetMoleculeType(cpxlen);
}

/**
 * @brief 在结构数据读取与布局中读取 `GetMoleculeTypeStr` 对应的数据。
 * @param moltype 供该函数读取或更新的 `moltype` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline
const char* GetMoleculeTypeStr(int moltype) {
    if(moltype == gtmtProtein) return "Protein";
    if(moltype == gtmtNA) return "NA";
    return "Other";
}

//nucleic acid molecule types
enum TGTNAAtomTypes {
        gtnaatC3p,
        gtnaatC4p,
        gtnaatC5p,
        gtnaatO3p,
        gtnaatO5p,
        gtnaatP,
        gtnaatOther,
        gtnaatNAtomTypes
};

//bounds for distances between paired nucleic acid atoms of given type
#define gtnaatC3p_LOWER (12.5f)
#define gtnaatC4p_LOWER (14.0f)
#define gtnaatC5p_LOWER (16.0f)
#define gtnaatO3p_LOWER (13.5f)
#define gtnaatO5p_LOWER (15.5f)
#define gtnaatP_LOWER (16.5f)
#define gtnaatC3p_UPPER (15.0f)
#define gtnaatC4p_UPPER (16.0f)
#define gtnaatC5p_UPPER (18.0f)
#define gtnaatO3p_UPPER (16.5f)
#define gtnaatO5p_UPPER (18.5f)
#define gtnaatP_UPPER (21.0f)
//statistical average of distances between paired nucleic acid atoms of given type
#define gtnaatC3p_LUB_AVG ((gtnaatC3p_UPPER + gtnaatC3p_LOWER) * 0.5f)
#define gtnaatC4p_LUB_AVG ((gtnaatC4p_UPPER + gtnaatC4p_LOWER) * 0.5f)
#define gtnaatC5p_LUB_AVG ((gtnaatC5p_UPPER + gtnaatC5p_LOWER) * 0.5f)
#define gtnaatO3p_LUB_AVG ((gtnaatO3p_UPPER + gtnaatO3p_LOWER) * 0.5f)
#define gtnaatO5p_LUB_AVG ((gtnaatO5p_UPPER + gtnaatO5p_LOWER) * 0.5f)
#define gtnaatP_LUB_AVG ((gtnaatP_UPPER + gtnaatP_LOWER) * 0.5f)
//max deviation from statistical average of distances
#define gtnaatC3p_LUB_DLT ((gtnaatC3p_UPPER - gtnaatC3p_LOWER) * 0.5f)
#define gtnaatC4p_LUB_DLT ((gtnaatC4p_UPPER - gtnaatC4p_LOWER) * 0.5f)
#define gtnaatC5p_LUB_DLT ((gtnaatC5p_UPPER - gtnaatC5p_LOWER) * 0.5f)
#define gtnaatO3p_LUB_DLT ((gtnaatO3p_UPPER - gtnaatO3p_LOWER) * 0.5f)
#define gtnaatO5p_LUB_DLT ((gtnaatO5p_UPPER - gtnaatO5p_LOWER) * 0.5f)
#define gtnaatP_LUB_DLT ((gtnaatP_UPPER - gtnaatP_LOWER) * 0.5f)
#define gtnaat_LUB_AVG_ERROR (1.5f)


//constants
enum TPSVectorConsts {
        pmv2DX,
        pmv2DY,
        pmv2DZ,
        pmv2DNoElems,//number of coordinates (dimensions)
};

//NOTE: data alignment is NOT required, as a separate buffer is 
// allocated for each field

enum TPM2DVectorFields {
//fields and their indices of the 2D vector, representing structure-specific data for a 
//given structure:
//NOTE: pps2DLen and pps2DDist have to be adjacent for efficient reading in the device!
        pps2DLen = 0,//structure length; [SHORT INT]
        pps2DDist = pps2DLen + 1,//total number of positions up to this structure; [UINT]
        pps2DType = pps2DDist + 1,//reserved for clustering; used to be structure type (protein, NA); [SHORT INT]
//{{Fields and their indices of 2D vector representing one position of 
//structure for scoring;
//NOTE: structure address of type unsigned short is reasonable as long as the 
// max number of structures is unlikely to be processed in parallel simultaneously (on both CPU and GPU)
        pmv2DCoords = pps2DType + 1,//coordinates; pmv2DNoElems fields in total
        //NOTE: composite field of TYPE (the sign bit), INSERTION CODE, CHAIN ID, and CHAIN ORDER, 
        //NOTE: all wrapped in 4 bytes of an int variable:
        pmv2D_Ins_Ch_Ord = pmv2DCoords + pmv2DNoElems,//reserved for potential structure output; [SIGNED INT]
        pmv2DResNumber = pmv2D_Ins_Ch_Ord + 1,//residue serial number [SHORT INT]
        pmv2Drsd = pmv2DResNumber + 1,//residue letter; [CHAR]
        pmv2Dss = pmv2Drsd + 1,//secondary structure (SS); [CHAR]; [*OPTIONAL*]
//}}
        pmv2DTotFlds//total number of fields
};
enum {
        pps2DStrFlds = pps2DType + 1,//total number of structure-specific fields
};

struct TPM2DVectorFieldSize {
    constexpr static size_t szvfs_[pmv2DTotFlds] = {
        SZINTYPE,//pps2DLen
        SZLNTYPE,//pps2DDist
        SZINTYPE,//pps2DType
        SZFPTYPE,//pmv2DCoords+0
        SZFPTYPE,//pmv2DCoords+1
        SZFPTYPE,//pmv2DCoords+2
        SZINTYPE,//pmv2D_Ins_Ch_Ord
        SZINTYPE,//pmv2DResNumber
        SZCHTYPE,//pmv2Drsd
        SZCHTYPE//pmv2Dss
    };
};

enum TPM2DIndexFields {
//{{Fields and their indices of 2D vector representing one position of indexed structure;
        pmv2DNdxCoords = 0,//coordinates; pmv2DNoElems fields in total
        pmv2DNdxLeft = pmv2DNdxCoords + pmv2DNoElems,//pointer (index) to left node in the tree [SHORT INT]
        pmv2DNdxRight = pmv2DNdxLeft + 1,//pointer (index) to right node in the tree [SHORT INT]
        pmv2DNdxOrgndx = pmv2DNdxRight + 1,//original position (node) index [SHORT INT]
//         pmv2DNdxCutDim = pmv2DNdxOrgndx + 1,//cutting dimension [CHAR]
//}}
        pmv2DTotIndexFlds//total number of fields
};

struct TPM2DIndexFieldSize {
    constexpr static size_t szvfs_[pmv2DTotIndexFlds] = {
        SZFPTYPE,//pmv2DNdxCoords+0
        SZFPTYPE,//pmv2DNdxCoords+1
        SZFPTYPE,//pmv2DNdxCoords+2
        SZINTYPE,//pmv2DNdxLeft
        SZINTYPE,//pmv2DNdxRight
        SZINTYPE//pmv2DNdxOrgndx
//         ,SZCHTYPE//pmv2DNdxCutDim
    };
};

#endif//__PM2DVectorFields_h__
