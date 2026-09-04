/***************************************************************************
 *   Copyright (C) 2021-2026 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __MpSecStr_h__
#define __MpSecStr_h__

#include "libutil/macros.h"
#include "libutil/CLOptions.h"
#include "libgenp/gdats/PM2DVectorFields.h"
#include "libmymp/mpstages/transformbase.h"
#include "MpSecStrBase.h"

// -------------------------------------------------------------------------
// class MpSecStr for calculating secondary structures
//
class MpSecStr
{
    // indices for distances along str. positions: two, three, four residues apart
    enum
    {
        css2RESdst,
        css3RESdst,
        css4RESdst,
        cssTotal
    };

public:
    /**
     * @brief 构造 `MpSecStr`，初始化其负责的CPU 二级结构计算状态。
     * @param tmpdpdiagbuffers 供当前步骤读取或更新的 `tmpdpdiagbuffers` 缓冲区。
     * @param querypmbeg 查询结构打包字段的起始指针数组。
     * @param querypmend 查询结构打包字段的结束指针数组。
     * @param bdbCpmbeg 参考结构打包字段的起始指针数组。
     * @param bdbCpmend 参考结构打包字段的结束指针数组。
     * @param nqystrs 当前批次中的查询结构数量。
     * @param ndbCstrs 当前批次中的参考结构数量。
     * @param nqyposs 当前批次中查询结构的总位置数。
     * @param ndbCposs 当前批次中参考结构的总位置数。
     * @param qystr1len 批次中最长查询结构的长度。
     * @param dbstr1len 批次中最长参考结构的长度。
     * @param qystrnlen 批次中最短查询结构的长度。
     * @param dbstrnlen 批次中最短参考结构的长度。
     * @param dbxpad 参考数据行末用于对齐访问的填充长度。
     * @return 无返回值；完成对象构造与初始状态设置。
     */
    MpSecStr(
        float *tmpdpdiagbuffers,
        char **querypmbeg, char **querypmend,
        char **bdbCpmbeg, char **bdbCpmend,
        uint nqystrs, uint ndbCstrs,
        uint nqyposs, uint ndbCposs,
        uint qystr1len, uint dbstr1len,
        uint qystrnlen, uint dbstrnlen,
        uint dbxpad)
        : tmpdpdiagbuffers_(tmpdpdiagbuffers),
          querypmbeg_(querypmbeg), querypmend_(querypmend),
          bdbCpmbeg_(bdbCpmbeg), bdbCpmend_(bdbCpmend),
          nqystrs_(nqystrs), ndbCstrs_(ndbCstrs),
          nqyposs_(nqyposs), ndbCposs_(ndbCposs),
          qystr1len_(qystr1len), dbstr1len_(dbstr1len),
          qystrnlen_(qystrnlen), dbstrnlen_(dbstrnlen),
          dbxpad_(dbxpad)
    {
    }

    /**
     * @brief 在CPU 二级结构计算中运行 `Run` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void Run()
    {
        ssk_kernel_helper(querypmbeg_, querypmend_, nqystrs_, qystr1len_);
        ssk_kernel_helper(bdbCpmbeg_, bdbCpmend_, ndbCstrs_, dbstr1len_);
        //
        nass_helper(tmpdpdiagbuffers_, querypmbeg_, querypmend_, nqystrs_, nqyposs_, qystr1len_);
        nass_helper(tmpdpdiagbuffers_, bdbCpmbeg_, bdbCpmend_, ndbCstrs_, ndbCposs_, dbstr1len_);
    }

private:
    /**
     * @brief 在CPU 二级结构计算中并行计算 `ssk_kernel_helper` 对应的数据。
     * @param pmbeg 供该函数读取或更新的 `pmbeg` 参数。
     * @param pmend 供该函数读取或更新的 `pmend` 参数。
     * @param nstrs 控制当前步骤范围或规模的 `nstrs`。
     * @param str1len 控制当前步骤范围或规模的 `str1len`。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void ssk_kernel_helper(
        char *const *const pmbeg, char *const *const pmend,
        int nstrs, int str1len);


    /**
     * @brief 在CPU 二级结构计算中并行计算 `nass_initialize_kernel_helper` 对应的数据。
     * @param tmpdpdiagbuffers 供当前步骤读取或更新的 `tmpdpdiagbuffers` 缓冲区。
     * @param nposs 控制当前步骤范围或规模的 `nposs`。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void nass_initialize_kernel_helper(
        float *const __RESTRICT__ tmpdpdiagbuffers,
        const int nposs);

    /**
     * @brief 在CPU 二级结构计算中并行计算 `nass_calcdistances_kernel_helper` 对应的数据。
     * @param tmpdpdiagbuffers 供当前步骤读取或更新的 `tmpdpdiagbuffers` 缓冲区。
     * @param pmbeg 供该函数读取或更新的 `pmbeg` 参数。
     * @param pmend 供该函数读取或更新的 `pmend` 参数。
     * @param nstrs 控制当前步骤范围或规模的 `nstrs`。
     * @param nposs 控制当前步骤范围或规模的 `nposs`。
     * @param str1len 控制当前步骤范围或规模的 `str1len`。
     * @param atomtype 供该函数读取或更新的 `atomtype` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void nass_calcdistances_kernel_helper(
        float *const __RESTRICT__ tmpdpdiagbuffers,
        const char *const *const pmbeg, const char *const *const pmend,
        const int nstrs, const int nposs, const int str1len, const int atomtype);

    /**
     * @brief 在CPU 二级结构计算中并行计算 `nass_calcsecstrs_kernel_helper` 对应的数据。
     * @param tmpdpdiagbuffers 供当前步骤读取或更新的 `tmpdpdiagbuffers` 缓冲区。
     * @param pmbeg 供该函数读取或更新的 `pmbeg` 参数。
     * @param pmend 供该函数读取或更新的 `pmend` 参数。
     * @param nstrs 控制当前步骤范围或规模的 `nstrs`。
     * @param nposs 控制当前步骤范围或规模的 `nposs`。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void nass_calcsecstrs_kernel_helper(
        const float *const __RESTRICT__ tmpdpdiagbuffers,
        char *const *const pmbeg, char *const *const pmend,
        const int nstrs, const int nposs);


    /**
     * @brief 在CPU 二级结构计算中处理 `nass_helper` 对应的数据。
     * @param tmpdpdiagbuffers 供当前步骤读取或更新的 `tmpdpdiagbuffers` 缓冲区。
     * @param pmbeg 供该函数读取或更新的 `pmbeg` 参数。
     * @param pmend 供该函数读取或更新的 `pmend` 参数。
     * @param nstrs 控制当前步骤范围或规模的 `nstrs`。
     * @param nposs 控制当前步骤范围或规模的 `nposs`。
     * @param str1len 控制当前步骤范围或规模的 `str1len`。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void nass_helper(
        float *const __RESTRICT__ tmpdpdiagbuffers,
        char *const *const pmbeg, char *const *const pmend,
        const int nstrs, const int nposs, const int str1len)
    {
        static const int atomtype = CLOptions::GetI_NATOM_type();
        nass_initialize_kernel_helper(tmpdpdiagbuffers, nposs);
        nass_calcdistances_kernel_helper(tmpdpdiagbuffers,  pmbeg, pmend,  nstrs, nposs, str1len, atomtype);
        nass_calcsecstrs_kernel_helper(tmpdpdiagbuffers,  pmbeg, pmend,  nstrs, nposs);
    }

private:
    template <int S>
    /**
     * @brief 在CPU 二级结构计算中读取 `GetDistance` 对应的数据。
     * @param crds 供该函数读取或更新的 `crds` 参数。
     * @param pi 供该函数读取或更新的 `pi` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    float GetDistance(
        const float (*crds)[pmv2DNoElems], int pi);

    /**
     * @brief 在CPU 二级结构计算中读取 `GetDistance_yx` 对应的数据。
     * @param crdsy 供该函数读取或更新的 `crdsy` 参数。
     * @param yi 供该函数读取或更新的 `yi` 参数。
     * @param crdsx 供该函数读取或更新的 `crdsx` 参数。
     * @param xi 供该函数读取或更新的 `xi` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    float GetDistance_yx(
        const float (*crdsy)[pmv2DNoElems], int yi,
        const float (*crdsx)[pmv2DNoElems], int xi);

    /**
     * @brief 在CPU 二级结构计算中读取 `GetNAPairingCondition` 对应的数据。
     * @param rsdy 供该函数读取或更新的 `rsdy` 参数。
     * @param rsdx 供该函数读取或更新的 `rsdx` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    bool GetNAPairingCondition(const char rsdy, const char* rsdx);

    /**
     * @brief 在CPU 二级结构计算中读取 `GetNADstDeviation` 对应的数据。
     * @param atype 供该函数读取或更新的 `atype` 参数。
     * @param dst 接收目标数据的 `dst`。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    float GetNADstDeviation(const int atype, const float* dst);

    /**
     * @brief 在CPU 二级结构计算中处理 `AassignSecStr` 对应的数据。
     * @param dsts 接收目标数据的 `dsts`。
     * @param pi 供该函数读取或更新的 `pi` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    char AassignSecStr(
        const float (*dsts)[cssTotal], int pi);

private:
    float *const tmpdpdiagbuffers_;
    char *const *const querypmbeg_, *const *const querypmend_;
    char *const *const bdbCpmbeg_, *const *const bdbCpmend_;
    const uint nqystrs_, ndbCstrs_;
    const uint nqyposs_, ndbCposs_;
    const uint qystr1len_, dbstr1len_;
    const uint qystrnlen_, dbstrnlen_;
    const uint dbxpad_;
};

// -------------------------------------------------------------------------
// INLINES ...
// -------------------------------------------------------------------------
// SSKGetDistance: calculate distances between two residues in the same
// sequence;
// crds, cached coordinates;
// pi, destination position;
// S, step of several residues apart;
//
#if defined(OS_MS_WINDOWS)
#define OMPDECLARE_MpSecStr_GetDistance
#else
#pragma omp declare simd linear(pi: 1) uniform(crds) notinbranch
#endif
template <int S>
/**
 * @brief 在CPU 二级结构计算中读取 `MpSecStr::GetDistance` 对应的数据。
 * @param crds 供该函数读取或更新的 `crds` 参数。
 * @param pi 供该函数读取或更新的 `pi` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline float MpSecStr::GetDistance(
    const float (*crds)[pmv2DNoElems], int pi)
{
    return sqrtf(distance2(
        crds[pi][pmv2DX], crds[pi][pmv2DY], crds[pi][pmv2DZ],
        crds[pi + S][pmv2DX], crds[pi + S][pmv2DY], crds[pi + S][pmv2DZ]));
}

// -------------------------------------------------------------------------
// GetDistance_yx: calculate distances between two residues in the same
// structure;
// crdsy, crdsx, cached coordinates;
// yi, position for crdsy;
// xi, position for crdsx;
//
#if defined(OS_MS_WINDOWS)
#define OMPDECLARE_MpSecStr_GetDistance_yx
#else
#pragma omp declare simd linear(xi: 1) uniform(crdsy,crdsx, yi) inbranch
#endif
/**
 * @brief 在CPU 二级结构计算中读取 `MpSecStr::GetDistance_yx` 对应的数据。
 * @param crdsy 供该函数读取或更新的 `crdsy` 参数。
 * @param yi 供该函数读取或更新的 `yi` 参数。
 * @param crdsx 供该函数读取或更新的 `crdsx` 参数。
 * @param xi 供该函数读取或更新的 `xi` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline float MpSecStr::GetDistance_yx(
    const float (*crdsy)[pmv2DNoElems], int yi,
    const float (*crdsx)[pmv2DNoElems], int xi)
{
    return sqrtf(distance2(
        crdsy[yi][pmv2DX], crdsy[yi][pmv2DY], crdsy[yi][pmv2DZ],
        crdsx[xi][pmv2DX], crdsx[xi][pmv2DY], crdsx[xi][pmv2DZ]));
}

// -------------------------------------------------------------------------
// GetNAPairingCondition: get the condition for pairing nucleotides;
// rsdy, rsdx, nucleotides;
//
#if defined(OS_MS_WINDOWS)
#define OMPDECLARE_MpSecStr_GetNAPairingCondition
#else
#pragma omp declare simd linear(rsdx: 1) uniform(rsdy) notinbranch
#endif
/**
 * @brief 在CPU 二级结构计算中读取 `MpSecStr::GetNAPairingCondition` 对应的数据。
 * @param rsdy 供该函数读取或更新的 `rsdy` 参数。
 * @param rsdx 供该函数读取或更新的 `rsdx` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline bool MpSecStr::GetNAPairingCondition(const char rsdy, const char* rsdx)
{
    return
        ((rsdy == 'T' || rsdy == 'U') && (*rsdx == 'A')) ||
        ((rsdy == 'A') && (*rsdx == 'T' || *rsdx == 'U')) ||
        ((rsdy == 'G') && (*rsdx == 'C' || *rsdx == 'U')) ||
        ((rsdy == 'C' || rsdy == 'U') && (*rsdx == 'G'));
}

// -------------------------------------------------------------------------
// GetNADstDeviation: get distance deviation from the statistical 
// average of distances between paired nucleic acid atoms of given type;
// atype, nucleic acid atom type;
// dst, observed distance between atoms of given type;
//
#if defined(OS_MS_WINDOWS)
#define OMPDECLARE_MpSecStr_GetNADstDeviation
#else
#pragma omp declare simd linear(dst: 1) uniform(atype) notinbranch
#endif
/**
 * @brief 在CPU 二级结构计算中读取 `MpSecStr::GetNADstDeviation` 对应的数据。
 * @param atype 供该函数读取或更新的 `atype` 参数。
 * @param dst 接收目标数据的 `dst`。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline float MpSecStr::GetNADstDeviation(const int atype, const float* dst)
{
    float dev;
    const float devl = 9999.9f;
    if(atype == gtnaatC3p) {dev = fabsf(*dst - gtnaatC3p_LUB_AVG); return (gtnaatC3p_LUB_DLT < dev)? devl: dev;}
    if(atype == gtnaatC4p) {dev = fabsf(*dst - gtnaatC4p_LUB_AVG); return (gtnaatC4p_LUB_DLT < dev)? devl: dev;}
    if(atype == gtnaatC5p) {dev = fabsf(*dst - gtnaatC5p_LUB_AVG); return (gtnaatC5p_LUB_DLT < dev)? devl: dev;}
    if(atype == gtnaatO3p) {dev = fabsf(*dst - gtnaatO3p_LUB_AVG); return (gtnaatO3p_LUB_DLT < dev)? devl: dev;}
    if(atype == gtnaatO5p) {dev = fabsf(*dst - gtnaatO5p_LUB_AVG); return (gtnaatO5p_LUB_DLT < dev)? devl: dev;}
    if(atype == gtnaatP) {dev = fabsf(*dst - gtnaatP_LUB_AVG); return (gtnaatP_LUB_DLT < dev)? devl: dev;}
    return devl;
}

// -------------------------------------------------------------------------
// AassignSecStr: assign secondary structure for given position once
// the required distances have been calculated;
// dsts, cache of calculated distances;
// pi, position to calculate secondary structure for;
// NOTE: pi corresponds to the structure position, hence
// NOTE: dsts position is nMRG0 greater;
// NOTE: pi assumed to be valid;
//
#if defined(OS_MS_WINDOWS)
#define OMPDECLARE_MpSecStr_AassignSecStr
#else
#pragma omp declare simd linear(pi : 1) uniform(dsts) notinbranch
#endif
/**
 * @brief 在CPU 二级结构计算中处理 `MpSecStr::AassignSecStr` 对应的数据。
 * @param dsts 接收目标数据的 `dsts`。
 * @param pi 供该函数读取或更新的 `pi` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline char MpSecStr::AassignSecStr(
    const float (*dsts)[cssTotal], int pi)
{
    // NOTE: calculating for position 3 of 5 consecutive position
    float d13 = dsts[pi][css2RESdst];
    float d14 = dsts[pi][css3RESdst];
    float d15 = dsts[pi][css4RESdst];
    float d24 = dsts[pi + 1][css2RESdst];
    float d25 = dsts[pi + 1][css3RESdst];
    float d35 = dsts[pi + 2][css2RESdst];

    float error = 2.1f;

    if (fabsf(d15 - 6.37f) < error && fabsf(d14 - 5.18f) < error &&
        fabsf(d25 - 5.18f) < error && fabsf(d13 - 5.45f) < error &&
        fabsf(d24 - 5.45f) < error && fabsf(d35 - 5.45f) < error)
        return pmvHELIX; // helix

    error = 1.42f;

    if (fabsf(d15 - 13.0f) < error && fabsf(d14 - 10.4f) < error &&
        fabsf(d25 - 10.4f) < error && fabsf(d13 - 6.1f) < error &&
        fabsf(d24 - 6.1f) < error && fabsf(d35 - 6.1f) < error)
        return pmvSTRND; // strand

    if (d15 < 8.0f)
        return pmvTURN; // turn

    return pmvLOOP; // loop/coil
}

#endif //__MpSecStr_h__
