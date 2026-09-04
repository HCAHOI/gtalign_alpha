/***************************************************************************
 *   Copyright (C) 2021-2023 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __MpStageSSRR_h__
#define __MpStageSSRR_h__

#include "libutil/mybase.h"
#include "libgenp/gproc/gproc.h"
#include "libgenp/gdats/PM2DVectorFields.h"
#include "libmymp/mpproc/mpprocconfbase.h"
#include "libmymp/mpstages/transformbase.h"
#include "libmymp/mpstages/scoringbase.h"
#include "libmycu/custages/stagecnsts.cuh"
#include "libmycu/custages/fragment.cuh"
#include "libmymp/mpstage1/MpStageBase.h"
#include "libmymp/mpstage1/MpStage1.h"
#include "libmymp/mpdp/MpDPHub.h"
#include "libmycu/cucom/cudef.h"

// -------------------------------------------------------------------------
// class MpStageSSRR for implementing structure comparison by matching 
// secondary structure and sequence similarity
//
class MpStageSSRR: public MpStage1 {
public:
    /**
     * @brief 构造 `MpStageSSRR`，初始化其负责的CPU 候选搜索与精修状态。
     * @param maxndpiters 动态规划精修允许的最大迭代次数。
     * @param maxnsteps 每对结构保留的候选搜索步数。
     * @param minfraglen 参与初始叠合的最短片段长度。
     * @param prescore 进入后续精修前要求的预筛选分数阈值。
     * @param stepinit 供该函数读取或更新的 `stepinit` 参数。
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
     * @param alndatamem 保存最终对齐统计量的缓冲区。
     * @param tfmmem 保存最终刚体变换矩阵的缓冲区。
     * @param globvarsbuf 供当前步骤读取或更新的 `globvarsbuf` 缓冲区。
     * @return 无返回值；完成对象构造与初始状态设置。
     */
    MpStageSSRR(
        const int maxndpiters,
        const uint maxnsteps,
        const uint minfraglen,
        const float prescore,
        const int stepinit,
        char** querypmbeg, char** querypmend,
        char** bdbCpmbeg, char** bdbCpmend,
        uint nqystrs, uint ndbCstrs,
        uint nqyposs, uint ndbCposs,
        uint qystr1len, uint dbstr1len,
        uint qystrnlen, uint dbstrnlen,
        uint dbxpad,
        float* scores, 
        float* tmpdpdiagbuffers, float* tmpdpbotbuffer,
        float* tmpdpalnpossbuffer, uint* maxscoordsbuf, char* btckdata,
        float* wrkmem, float* wrkmemccd, float* wrkmemtm, float* wrkmemtmibest,
        float* wrkmemaux, float* wrkmem2, float* alndatamem, float* tfmmem,
        uint* globvarsbuf)
    :
        MpStage1(
            maxndpiters, maxnsteps, minfraglen, prescore, stepinit,
            querypmbeg, querypmend, bdbCpmbeg, bdbCpmend,
            nqystrs, ndbCstrs, nqyposs, ndbCposs,
            qystr1len, dbstr1len, qystrnlen, dbstrnlen, dbxpad,
            scores,
            tmpdpdiagbuffers, tmpdpbotbuffer, tmpdpalnpossbuffer,
            maxscoordsbuf, btckdata,
            wrkmem, wrkmemccd, wrkmemtm, wrkmemtmibest,
            wrkmemaux, wrkmem2, alndatamem, tfmmem,
            globvarsbuf
        )
    {}

    /**
     * @brief 在CPU 候选搜索与精修中运行 `Run` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    virtual void Run() {}

    template<bool USESEQSCORING>
    /**
     * @brief 按模板开关选择是否加入序列分数，运行二级结构引导的候选搜索与 DP 精修。
     * @param maxndpiters 动态规划精修允许的最大迭代次数。
     * @return 无返回值；二级结构/序列候选及其精修变换写入工作缓冲区。
     */
    void RunSpecialized(const int maxndpiters)
    {
        //draw alignment based on ss and sequence similarity:
        Align<USESEQSCORING>();
        //refine alignment boundaries to improve scores
        RefineFragDPKernelCaller(false/*readlocalconv*/, FRAGREF_NMAXCONVIT);
        //refine alignment boundaries identified by applying DP;
        //1. With a gap cost:
        DPRefine<false/*GAP0*/,false/*PRESCREEN*/,false/*WRKMEMTM1*/>(2/*maxndpiters_*/, prescore_);
        //2. No gap cost:
        DPRefine<true/*GAP0*/,false/*PRESCREEN*/,false/*WRKMEMTM1*/>(maxndpiters, prescore_);
    }

protected:
    template<bool USESEQSCORING>
    /**
     * @brief 用二级结构相容性以及可选的序列相似性建立一条 DP 对齐路径。
     * @par 参数
     * 无。
     * @return 无返回值；回溯方向与匹配位置写入 DP 缓冲区。
     */
    void Align();
};



// -------------------------------------------------------------------------
// INLINES ...
// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
// Align: get alignment based on secondary structure information and
// sequence similarity;
// USESEQSCORING, flag of using sequence similarity scoring;
//
template<bool USESEQSCORING>
/**
 * @brief 在CPU 候选搜索与精修中对齐 `MpStageSSRR::Align` 对应的数据。
 * @par 参数
 * 无。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
inline
void MpStageSSRR::Align()
{
    MYMSG("MpStageSSRR::Align", 5);
    // static std::string preamb = "MpStageSSRR::Align: ";

    constexpr float gcost = -1.0f;
    static const float wgt4ss = 1.0f;//weight for scoring ss
    static const float wgt4rr = 0.2f;//weight for pairwise residue scoring

    MpDPHub dphub(
        maxnsteps_,
        querypmbeg_, querypmend_, bdbCpmbeg_, bdbCpmend_,
        nqystrs_, ndbCstrs_, nqyposs_, ndbCposs_,  qystr1len_, dbstr1len_, dbxpad_,
        tmpdpdiagbuffers_, tmpdpbotbuffer_, tmpdpalnpossbuffer_, maxscoordsbuf_, btckdata_,
        wrkmem_, wrkmemccd_,  wrkmemtm_,  wrkmemtmibest_,
        wrkmemaux_, wrkmem2_, alndatamem_, tfmmem_, globvarsbuf_
    );

    dphub.ExecDPSSwBtck128xKernel<USESEQSCORING>(
        gcost, wgt4ss, wgt4rr,
        querypmbeg_, bdbCpmbeg_,
        wrkmemaux_, tmpdpdiagbuffers_, tmpdpbotbuffer_, btckdata_);

    dphub.BtckToMatched128xKernel<false/*ANCHORRGN*/,false/*BANDED*/>(
        0/*stepnumber*/,//slot 0
        querypmbeg_, bdbCpmbeg_, btckdata_, wrkmemaux_, tmpdpalnpossbuffer_);
}


// -------------------------------------------------------------------------

#endif//__MpStageSSRR_h__
