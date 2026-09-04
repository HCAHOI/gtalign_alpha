/***************************************************************************
 *   Copyright (C) 2021-2023 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include <string>
#include <vector>
#include <map>

#include "libutil/cnsts.h"
#include "libutil/macros.h"
#include "libutil/CLOptions.h"
#include "libgenp/gproc/gproc.h"
#include "libgenp/gproc/dputils.h"
#include "libgenp/gdats/PM2DVectorFields.h"

#include "libmycu/cucom/cucommon.h"
#include "libmycu/cucom/warpscan.cuh"
#include "libmycu/cucom/cugraphs.cuh"
#include "libmycu/cuproc/cuprocconf.h"
#include "libmycu/culayout/cuconstant.cuh"

#include "libmycu/custages/stagecnsts.cuh"
#include "libmycu/custages/scoring.cuh"
#include "libmycu/cudp/dpw_btck.cuh"
#include "libmycu/cudp/dpssw_btck.cuh"
#include "libmycu/cudp/btck2match.cuh"
#include "libmycu/custage1/custage1.cuh"
#include "custage_ssrr.cuh"

// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
// run_stage_ssrr: search for superposition between multiple molecules
// simultaneoulsy and identify alignments by DP using secondary 
// structure information and sequence similarity criteria;
// USESEQSCORING, template parameter, flag of using sequence similarity scoring;
// qystr1len, length of the largest query;
// dbstr1len, length of the largest reference;
// qystrnlen, length of the smallest query;
// dbstrnlen, length of the smallest reference;
//
template<bool USESEQSCORING>
/**
 * @brief 在 GPU 上运行二级结构及可选序列相似性引导的候选搜索与精修阶段。
 * @param stgraphs 供该函数读取或更新的 `stgraphs` 参数。
 * @param streamproc 执行当前计算阶段的 CUDA 流。
 * @param scorethld 当前步骤使用或写回的 `scorethld` 分数。
 * @param prescore 进入后续精修前要求的预筛选分数阈值。
 * @param maxndpiters 动态规划精修允许的最大迭代次数。
 * @param maxnsteps 每对结构保留的候选搜索步数。
 * @param minfraglen 参与初始叠合的最短片段长度。
 * @param nqystrs 当前批次中的查询结构数量。
 * @param ndbCstrs 当前批次中的参考结构数量。
 * @param nqyposs 当前批次中查询结构的总位置数。
 * @param ndbCposs 当前批次中参考结构的总位置数。
 * @param qystr1len 批次中最长查询结构的长度。
 * @param dbstr1len 批次中最长参考结构的长度。
 * @param qystrnlen 批次中最短查询结构的长度。
 * @param dbstrnlen 批次中最短参考结构的长度。
 * @param dbxpad 参考数据行末用于对齐访问的填充长度。
 * @param scores 保存或读取对齐分数的缓冲区。
 * @param tmpdpdiagbuffers 供当前步骤读取或更新的 `tmpdpdiagbuffers` 缓冲区。
 * @param tmpdpbotbuffer 供当前步骤读取或更新的 `tmpdpbotbuffer` 缓冲区。
 * @param tmpdpalnpossbuffer 供当前步骤读取或更新的 `tmpdpalnpossbuffer` 缓冲区。
 * @param maxscoordsbuf 保存动态规划最大分数坐标的缓冲区。
 * @param btckdata 保存动态规划回溯方向的缓冲区。
 * @param wrkmem 当前计算阶段的主工作缓冲区。
 * @param wrkmemccd 供当前步骤读取或更新的 `wrkmemccd` 缓冲区。
 * @param wrkmemtm 保存候选刚体变换的工作缓冲区。
 * @param wrkmemtmibest 保存各候选当前最佳刚体变换的缓冲区。
 * @param wrkmemaux 保存分数、收敛标记等辅助状态的工作缓冲区。
 * @param wrkmem2 供当前步骤读取或更新的 `wrkmem2` 缓冲区。
 * @param tfmmem 保存最终刚体变换矩阵的缓冲区。
 * @param globvarsbuf 供当前步骤读取或更新的 `globvarsbuf` 缓冲区。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
void stage_ssrr::run_stage_ssrr(
    std::map<CGKey,MyCuGraph>& stgraphs,
    cudaStream_t streamproc,
    const float /*scorethld*/,
    const float prescore,
    const int maxndpiters,
    const uint maxnsteps,
    const uint minfraglen,
    uint nqystrs, uint ndbCstrs,
    uint nqyposs, uint ndbCposs,
    uint qystr1len, uint dbstr1len,
    uint qystrnlen, uint dbstrnlen,
    uint dbxpad,
    float* __restrict__ /*scores*/, 
    float* __restrict__ tmpdpdiagbuffers,
    float* __restrict__ tmpdpbotbuffer,
    float* __restrict__ tmpdpalnpossbuffer,
    uint* __restrict__ maxscoordsbuf,
    char* __restrict__ btckdata,
    float* __restrict__ wrkmem,
    float* __restrict__ wrkmemccd,
    float* __restrict__ wrkmemtm,
    float* __restrict__ wrkmemtmibest,
    float* __restrict__ wrkmemaux,
    float* __restrict__ wrkmem2,
    float* __restrict__ tfmmem,
    uint* __restrict__ /*globvarsbuf*/)
{
    //alignment based on ss information and sequence similarity:
    stage_ssrr_align<USESEQSCORING>(
        streamproc,
        maxnsteps,
        nqystrs, ndbCstrs,
        nqyposs, ndbCposs,
        qystr1len, dbstr1len,
        qystrnlen, dbstrnlen,
        dbxpad,
        tmpdpdiagbuffers,
        tmpdpbotbuffer,
        tmpdpalnpossbuffer,
        maxscoordsbuf,
        btckdata,
        wrkmemaux);

    //refine alignment boundaries to improve scores
    stage1::stage1_refinefrag<false/* CONDITIONAL */>(
        stgraphs,
        stg1REFINE_INITIAL_DP/*fragments identified by DP*/,
        FRAGREF_NMAXCONVIT/*#iterations until convergence*/,
        streamproc,
        maxnsteps, minfraglen,
        nqystrs, ndbCstrs,
        nqyposs, ndbCposs,
        qystr1len, dbstr1len,
        qystrnlen, dbstrnlen,
        dbxpad,
        tmpdpdiagbuffers,
        tmpdpalnpossbuffer,
        wrkmem, wrkmemccd, wrkmemtm, wrkmemtmibest,
        wrkmemaux, wrkmem2, tfmmem);

    //refine alignment boundaries identified in the previous
    //substage by applying DP;
    //1. With a gap cost:
    stage1_dprefine<false/*GAP0*/,false/*PRESCREEN*/>(
        stgraphs,
        streamproc,
        2/* maxndpiters */,
        prescore,
        maxnsteps, minfraglen,
        nqystrs, ndbCstrs,
        nqyposs, ndbCposs,
        qystr1len, dbstr1len,
        qystrnlen, dbstrnlen,
        dbxpad,
        tmpdpdiagbuffers,
        tmpdpbotbuffer,
        tmpdpalnpossbuffer,
        maxscoordsbuf,
        btckdata,
        wrkmem, wrkmemccd, wrkmemtm, wrkmemtmibest,
        wrkmemaux, wrkmem2, tfmmem/*out*/);

    //2. No gap cost:
    stage1_dprefine<true/*GAP0*/,false/*PRESCREEN*/>(
        stgraphs,
        streamproc,
        maxndpiters,
        prescore,
        maxnsteps, minfraglen,
        nqystrs, ndbCstrs,
        nqyposs, ndbCposs,
        qystr1len, dbstr1len,
        qystrnlen, dbstrnlen,
        dbxpad,
        tmpdpdiagbuffers,
        tmpdpbotbuffer,
        tmpdpalnpossbuffer,
        maxscoordsbuf,
        btckdata,
        wrkmem, wrkmemccd, wrkmemtm, wrkmemtmibest,
        wrkmemaux, wrkmem2, tfmmem/*out*/);
}

// Instantiations
// 
#define INSTANTIATE_stage_ssrr__run_stage_ssrr(USESEQSCORING) \
    template void stage_ssrr::run_stage_ssrr<USESEQSCORING>( \
    std::map<CGKey,MyCuGraph>& stgraphs, \
    cudaStream_t streamproc, \
    const float scorethld, \
    const float prescore, \
    const int maxndpiters, \
    const uint maxnsteps, \
    const uint minfraglen, \
    uint nqystrs, uint ndbCstrs, \
    uint nqyposs, uint ndbCposs, \
    uint qystr1len, uint dbstr1len, \
    uint qystrnlen, uint dbstrnlen, \
    uint dbxpad, \
    float* __restrict__ /*scores*/,  \
    float* __restrict__ tmpdpdiagbuffers, \
    float* __restrict__ tmpdpbotbuffer, \
    float* __restrict__ tmpdpalnpossbuffer, \
    uint* __restrict__ maxscoordsbuf, \
    char* __restrict__ btckdata, \
    float* __restrict__ wrkmem, \
    float* __restrict__ wrkmemccd, \
    float* __restrict__ wrkmemtm, \
    float* __restrict__ wrkmemtmibest, \
    float* __restrict__ wrkmemaux, \
    float* __restrict__ wrkmem2, \
    float* __restrict__ tfmmem, \
    uint* __restrict__ /*globvarsbuf*/);

INSTANTIATE_stage_ssrr__run_stage_ssrr(false);
INSTANTIATE_stage_ssrr__run_stage_ssrr(true);



// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
// stage_ssrr_align: get alignment based on secondary structure information
// and sequence similarity;
// qystr1len, length of the largest query;
// dbstr1len, length of the largest reference;
// qystrnlen, length of the smallest query;
// dbstrnlen, length of the smallest reference;
//
template<bool USESEQSCORING>
/**
 * @brief 在CUDA 候选搜索与精修中处理 `stage_ssrr::stage_ssrr_align` 对应的数据。
 * @param streamproc 执行当前计算阶段的 CUDA 流。
 * @param maxnsteps 每对结构保留的候选搜索步数。
 * @param nqystrs 当前批次中的查询结构数量。
 * @param ndbCstrs 当前批次中的参考结构数量。
 * @param nqyposs 当前批次中查询结构的总位置数。
 * @param ndbCposs 当前批次中参考结构的总位置数。
 * @param qystr1len 批次中最长查询结构的长度。
 * @param dbstr1len 批次中最长参考结构的长度。
 * @param qystrnlen 批次中最短查询结构的长度。
 * @param dbstrnlen 批次中最短参考结构的长度。
 * @param dbxpad 参考数据行末用于对齐访问的填充长度。
 * @param tmpdpdiagbuffers 供当前步骤读取或更新的 `tmpdpdiagbuffers` 缓冲区。
 * @param tmpdpbotbuffer 供当前步骤读取或更新的 `tmpdpbotbuffer` 缓冲区。
 * @param tmpdpalnpossbuffer 供当前步骤读取或更新的 `tmpdpalnpossbuffer` 缓冲区。
 * @param maxscoordsbuf 保存动态规划最大分数坐标的缓冲区。
 * @param btckdata 保存动态规划回溯方向的缓冲区。
 * @param wrkmemaux 保存分数、收敛标记等辅助状态的工作缓冲区。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
void stage_ssrr::stage_ssrr_align(
    cudaStream_t streamproc,
    const uint maxnsteps,
    uint nqystrs, uint ndbCstrs,
    uint /*nqyposs*/, uint ndbCposs,
    uint qystr1len, uint dbstr1len,
    uint /*qystrnlen*/, uint /*dbstrnlen*/,
    uint dbxpad,
    float* __restrict__ tmpdpdiagbuffers,
    float* __restrict__ tmpdpbotbuffer,
    float* __restrict__ tmpdpalnpossbuffer,
    uint* __restrict__ /*maxscoordsbuf*/,
    char* __restrict__ btckdata,
    float* __restrict__ wrkmemaux)
{
    constexpr float gcost = -1.0f;
    static const float wgt4ss = 1.0f;//weight for scoring ss
    static const float wgt4rr = 0.2f;//weight for pairwise residue scoring

    //execution configuration for DP:
    //1D thread block processes 2D DP matrix oblique block of dimension 
    //CUDP_2DCACHE_DIM_D x CUDP_2DCACHE_DIM_X;
    //NOTE: using block diagonals, where blocks share a common point 
    //NOTE: (corner) with a neighbour in a diagonal;
    const uint maxblkdiagelems = GetMaxBlockDiagonalElems(
            dbstr1len, qystr1len, CUDP_2DCACHE_DIM_D, CUDP_2DCACHE_DIM_X);
    dim3 nthrds_dp(CUDP_2DCACHE_DIM_D,1,1);
    dim3 nblcks_dp(maxblkdiagelems,ndbCstrs,nqystrs);
    uint nblkdiags = (uint)
        (((dbstr1len + qystr1len) + CUDP_2DCACHE_DIM_X-1) / CUDP_2DCACHE_DIM_X);
    nblkdiags += (uint)(qystr1len - 1) / CUDP_2DCACHE_DIM_D;

    //execution configuration for extracting matched positions
    //identified during DP:
    dim3 nthrds_mtch(CUDP_MATCHED_DIM_X,CUDP_MATCHED_DIM_Y,1);
    dim3 nblcks_mtch(ndbCstrs,nqystrs,1);

    //launch blocks along block diagonals to perform DP;
    //nblkdiags, total number of diagonals:
    for(uint d = 0; d < nblkdiags; d++)
    {
        ExecDPSSwBtck3264x<USESEQSCORING>
            <<<nblcks_dp,nthrds_dp,0,streamproc>>>(
                d, ndbCstrs, ndbCposs, dbxpad, maxnsteps,
                wgt4ss, wgt4rr, gcost,
                wrkmemaux, tmpdpdiagbuffers, tmpdpbotbuffer, btckdata);
        MYCUDACHECKLAST;
    }

    //process the result of DP
    BtckToMatched32x<false/*ANCHORRGN*/,false/*BANDED*/>
        <<<nblcks_mtch,nthrds_mtch,0,streamproc>>>(
            ndbCstrs, ndbCposs, dbxpad, maxnsteps, 0/*stepnumber*/,
            btckdata, wrkmemaux, tmpdpalnpossbuffer);
    MYCUDACHECKLAST;
}
