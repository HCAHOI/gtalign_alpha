/***************************************************************************
 *   Copyright (C) 2021-2023 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __MpReform_h__
#define __MpReform_h__

#include <math.h>
#include "libutil/mybase.h"
#include "libgenp/gproc/gproc.h"
#include "libgenp/gdats/PM2DVectorFields.h"
#include "libgenp/gdats/PMBatchStrData.h"
#include "libgenp/gdats/PMBatchStrDataIndex.h"
#include "libmymp/mpproc/mpprocconf.h"
#include "libmymp/mputil/simdscan.h"
#include "libmycu/custages/stagecnsts.cuh"
#include "libmycu/custages/fragment.cuh"
#include "libmycu/cucom/cudef.h"

// -------------------------------------------------------------------------
// class MpReform for reformatting data of structure pairs that 
// proceed to the next stage
//
class MpReform {
public:
    /**
     * @brief 构造 `MpReform`，初始化其负责的CPU 对齐流水线状态。
     * @param maxnsteps 每对结构保留的候选搜索步数。
     * @param querypmbeg 查询结构打包字段的起始指针数组。
     * @param querypmend 查询结构打包字段的结束指针数组。
     * @param bdbCpmbeg 参考结构打包字段的起始指针数组。
     * @param bdbCpmend 参考结构打包字段的结束指针数组。
     * @param queryndxpmbeg 描述查询结构的 `queryndxpmbeg`。
     * @param queryndxpmend 描述查询结构的 `queryndxpmend`。
     * @param bdbCndxpmbeg 描述参考结构的 `bdbCndxpmbeg`。
     * @param bdbCndxpmend 描述参考结构的 `bdbCndxpmend`。
     * @param nqystrs 当前批次中的查询结构数量。
     * @param ndbCstrs 当前批次中的参考结构数量。
     * @param nqyposs 当前批次中查询结构的总位置数。
     * @param ndbCposs 当前批次中参考结构的总位置数。
     * @param tmpdpdiagbuffers 供当前步骤读取或更新的 `tmpdpdiagbuffers` 缓冲区。
     * @param wrkmemaux 保存分数、收敛标记等辅助状态的工作缓冲区。
     * @param tfmmem 保存最终刚体变换矩阵的缓冲区。
     * @param globvarsbuf 供当前步骤读取或更新的 `globvarsbuf` 缓冲区。
     * @param filterdata 供该函数读取或更新的 `filterdata` 参数。
     * @return 无返回值；完成对象构造与初始状态设置。
     */
    MpReform(
        const uint maxnsteps,
        char** querypmbeg, char** querypmend,
        char** bdbCpmbeg, char** bdbCpmend,
        char** queryndxpmbeg, char** queryndxpmend,
        char** bdbCndxpmbeg, char** bdbCndxpmend,
        uint nqystrs, uint ndbCstrs,
        uint nqyposs, uint ndbCposs,
        float* tmpdpdiagbuffers,
        float* wrkmemaux, float* tfmmem, uint* globvarsbuf, uint* filterdata)
    :
        maxnsteps_(maxnsteps),
        querypmbeg_(querypmbeg), querypmend_(querypmend),
        bdbCpmbeg_(bdbCpmbeg), bdbCpmend_(bdbCpmend),
        queryndxpmbeg_(queryndxpmbeg), queryndxpmend_(queryndxpmend),
        bdbCndxpmbeg_(bdbCndxpmbeg), bdbCndxpmend_(bdbCndxpmend),
        nqystrs_(nqystrs), ndbCstrs_(ndbCstrs),
        nqyposs_(nqyposs), ndbCposs_(ndbCposs),
        tmpdpdiagbuffers_(tmpdpdiagbuffers),
        wrkmemaux_(wrkmemaux), tfmmem_(tfmmem), globvarsbuf_(globvarsbuf),
        filterdata_(filterdata)
    {}

    /**
     * @brief 在CPU 对齐流水线中构造 `MakeDbCandidateList` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void MakeDbCandidateList() {
        constexpr int memalignment = 
            mycemin((size_t)PMBSdatalignment, CuMemoryBase::GetMinMemAlignment());
        MakeDbCandidateListHelper<memalignment>(
            nqystrs_, ndbCstrs_, maxnsteps_,
            querypmbeg_, bdbCpmbeg_, wrkmemaux_, filterdata_);
    }

    /**
     * @brief 在CPU 对齐流水线中处理 `SelectAndReformat` 对应的数据。
     * @param ndbCstrs2 控制当前步骤范围或规模的 `ndbCstrs2`。
     * @param maxndbCposs 描述参考结构的 `maxndbCposs`。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void SelectAndReformat(const int ndbCstrs2, const int maxndbCposs) {
        SelectAndReformatKernel(
            ndbCstrs2, maxndbCposs,  filterdata_,
            querypmbeg_, bdbCpmbeg_, tfmmem_, wrkmemaux_, tmpdpdiagbuffers_);
    }


protected:
    template<int DATALN>
    /**
     * @brief 在CPU 对齐流水线中构造 `MakeDbCandidateListHelper` 对应的数据。
     * @param nqystrs 当前批次中的查询结构数量。
     * @param ndbCstrs 当前批次中的参考结构数量。
     * @param maxnsteps 每对结构保留的候选搜索步数。
     * @param querypmbeg 查询结构打包字段的起始指针数组。
     * @param bdbCpmbeg 参考结构打包字段的起始指针数组。
     * @param wrkmemaux 保存分数、收敛标记等辅助状态的工作缓冲区。
     * @param filterdata 供该函数读取或更新的 `filterdata` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void MakeDbCandidateListHelper(
        const int nqystrs, const int ndbCstrs, const int maxnsteps,
        const char* const * const __RESTRICT__ /* querypmbeg */,
        const char* const * const __RESTRICT__ bdbCpmbeg,
        const float* const __RESTRICT__ wrkmemaux,
        uint* const __RESTRICT__ filterdata);

    /**
     * @brief 在CPU 对齐流水线中并行计算 `SelectAndReformatKernel` 对应的数据。
     * @param ndbCstrs2 控制当前步骤范围或规模的 `ndbCstrs2`。
     * @param maxndbCposs 描述参考结构的 `maxndbCposs`。
     * @param filterdata 供该函数读取或更新的 `filterdata` 参数。
     * @param querypmbeg 查询结构打包字段的起始指针数组。
     * @param bdbCpmbeg 参考结构打包字段的起始指针数组。
     * @param tfmmem 保存最终刚体变换矩阵的缓冲区。
     * @param wrkmemaux 保存分数、收敛标记等辅助状态的工作缓冲区。
     * @param tmpdpdiagbuffers 供当前步骤读取或更新的 `tmpdpdiagbuffers` 缓冲区。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void SelectAndReformatKernel(
        const int ndbCstrs2,
        const int maxndbCposs,
        const uint* const __RESTRICT__ filterdata,
        const char* const * const __RESTRICT__ querypmbeg,
        const char* const * const __RESTRICT__ bdbCpmbeg,
        float* const __RESTRICT__ tfmmem,
        float* const __RESTRICT__ wrkmemaux,
        float* const __RESTRICT__ tmpdpdiagbuffers);


protected:
    const uint maxnsteps_;
    char* const * const querypmbeg_, * const * const querypmend_;
    char* const * const bdbCpmbeg_, * const *const bdbCpmend_;
    char* const * const queryndxpmbeg_, * const * const queryndxpmend_;
    char* const * const bdbCndxpmbeg_, * const *const bdbCndxpmend_;
    const uint nqystrs_, ndbCstrs_;
    const uint nqyposs_, ndbCposs_;
    float* const tmpdpdiagbuffers_;
    float* const wrkmemaux_, *const tfmmem_;
    uint* const globvarsbuf_, *const filterdata_;
};



// -------------------------------------------------------------------------
// INLINES ...
// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
// MakeDbCandidateListHelper: make list of reference structure (database)
// candidates proceeding to stages of more detailed superposition search and
// refinement;
// nqystrs, total number of query structures in the chunk;
// ndbCstrs, total number of reference structures in the chunk;
// maxnsteps, max number of steps to perform for each reference structure;
// NOTE: memory pointers should be aligned!
// wrkmemaux, auxiliary working memory;
// filterdata, memory of new indices and addresses of passing references;
// NOTE: processes the reference structures over all queries for flags;
//
template<int DATALN>
/**
 * @brief 在CPU 对齐流水线中构造 `MpReform::MakeDbCandidateListHelper` 对应的数据。
 * @param nqystrs 当前批次中的查询结构数量。
 * @param ndbCstrs 当前批次中的参考结构数量。
 * @param maxnsteps 每对结构保留的候选搜索步数。
 * @param querypmbeg 查询结构打包字段的起始指针数组。
 * @param bdbCpmbeg 参考结构打包字段的起始指针数组。
 * @param wrkmemaux 保存分数、收敛标记等辅助状态的工作缓冲区。
 * @param filterdata 供该函数读取或更新的 `filterdata` 参数。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
inline
void MpReform::MakeDbCandidateListHelper(
    const int nqystrs, const int ndbCstrs, const int maxnsteps,
    const char* const * const __RESTRICT__ /* querypmbeg */,
    const char* const * const __RESTRICT__ bdbCpmbeg,
    const float* const __RESTRICT__ wrkmemaux,
    uint* const __RESTRICT__ filterdata)
{
    enum {
        PREPXS = 0,//index for the previous prefix sum and padding
        pad = 1,//padding
        lXFLG = fdNewReferenceIndex,//index for reference structure convergence flags/new indices
        lXLEN = fdNewReferenceAddress,//index for reference structure convergence lengths/new addresses
        lYDIM = nTFilterData,//MPFL_MAKECANDIDATELIST_YDIM,
        lXDIM = MPFL_MAKECANDIDATELIST_XDIM,
        lXDIM1 = MPFL_MAKECANDIDATELIST_XDIM + pad
    };

    //indices and lengths of selected structures
    int strd[lYDIM][lXDIM1];
    int tmp[lXDIM];

    strd[lXFLG][PREPXS] = 0;
    strd[lXLEN][PREPXS] = 0;

    for(int ri0 = 0; ri0 < ndbCstrs; ri0 += lXDIM)
    {
        //update the prefix sums originated from processing the last data block:
        if(ri0) {
            strd[lXFLG][PREPXS] = strd[lXFLG][lXDIM1 - 1];
            strd[lXLEN][PREPXS] = strd[lXLEN][lXDIM1 - 1];
        }

        // uint dbstrndx = dbstrndx0 + threadIdx.x;//reference index
        // int value = 0;//convflag for lXFLG and dbstrlen for lXLEN

        const int riend = mymin(ndbCstrs, ri0 + lXDIM);

        #pragma omp simd
        for(int ri = ri0; ri < riend; ri++) {
            int ii = ri - ri0 + pad;
            strd[lXFLG][ii] = 0;
        }

        //get convergence flags (over all queries)
        for(int qi/*qryndx*/ = 0; qi < nqystrs; qi++) {
            int mloc0 = ((qi * maxnsteps + 0) * nTAuxWorkingMemoryVars) * ndbCstrs;
            #pragma omp simd aligned(wrkmemaux:DATALN)
            for(int ri = ri0; ri < riend; ri++) {
                int ii = ri - ri0 + pad;
                int lconv = wrkmemaux[mloc0 + tawmvConverged * ndbCstrs + ri];//float->int
                strd[lXFLG][ii] += ((lconv & CONVERGED_LOWTMSC_bitval) != 0);
            }
        }

        //selected structures have no convergence flags set for all queries;
        //get reference lengths too;
        #pragma omp simd aligned(bdbCpmbeg:DATALN)
        for(int ri = ri0; ri < riend; ri++) {
            int ii = ri - ri0 + pad;
            int value = strd[lXFLG][ii];
            strd[lXFLG][ii] = (value < nqystrs);
            strd[lXLEN][ii] = PMBatchStrData::GetLengthAt(bdbCpmbeg, ri);
        }

        //set reference lengths to 0 where convflag is set
        #pragma omp simd
        for(int ri = ri0; ri < riend; ri++) {
            int ii = ri - ri0 + pad;
            if(strd[lXFLG][ii] == 0)
                strd[lXLEN][ii] = 0;
        }

        //calculate inclusive (!) prefix sums for both flags, which then give indices,
        //and lengths for addresses:
        //TODO: change 1st agrument to the real size;
        mysimdincprefixsum<lXDIM>(lXDIM, &strd[lXFLG][1], tmp);
        mysimdincprefixsum<lXDIM>(lXDIM, &strd[lXLEN][1], tmp);

        //correct the prefix sums by adding the previously obtained values:
        #pragma omp simd
        for(int ri = ri0; ri < riend; ri++) {
            int ii = ri - ri0 + pad;
            strd[lXFLG][ii] += strd[lXFLG][PREPXS];
            strd[lXLEN][ii] += strd[lXLEN][PREPXS];
        }

        //write to output:
        #pragma omp simd
        for(int ri = ri0; ri < riend; ri++) {
            int ii = ri - ri0 + pad;
            int mloc = lXFLG * ndbCstrs + ri;
            int valueprev = strd[lXFLG][ii - 1];
            int value = strd[lXFLG][ii];
            //set to 0 for filtered-out structures:
            if(value == valueprev) value = 0;
            filterdata[mloc] = value;//WRITE indices
            //same for addresses:
            mloc = lXLEN * ndbCstrs + ri;
            valueprev = strd[lXLEN][ii - 1];
            value = strd[lXLEN][ii];
            if(value == valueprev) valueprev = 0;
            filterdata[mloc] = valueprev;//WRITE adjusted addresses
        }
    }
}

// -------------------------------------------------------------------------

#endif//__MpReform_h__
