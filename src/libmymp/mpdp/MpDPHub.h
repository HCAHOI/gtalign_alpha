/***************************************************************************
 *   Copyright (C) 2021-2023 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __MpDPHub_h__
#define __MpDPHub_h__

#include <omp.h>

#include "libutil/mybase.h"
#include "libgenp/gproc/gproc.h"
#include "libgenp/gdats/PM2DVectorFields.h"
#include "libgenp/gdats/PMBatchStrData.h"
#include "libmymp/mpproc/mpprocconfbase.h"
#include "libmymp/mpstages/transformbase.h"
#include "libmymp/mpstages/scoringbase.h"
#include "libmycu/custages/stagecnsts.cuh"
#include "libmycu/custages/fragment.cuh"
#include "libmymp/mpdp/mpdpbase.h"
#include "libmycu/cucom/cudef.h"

// -------------------------------------------------------------------------
// class MpDPHub implementing DP variants for structure alignment
//
class MpDPHub {
public:
    /**
     * @brief 构造 `MpDPHub`，初始化其负责的CPU 动态规划状态。
     * @param maxnsteps 每对结构保留的候选搜索步数。
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
     * @param dbxpad 参考数据行末用于对齐访问的填充长度。
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
    MpDPHub(
        const uint maxnsteps,
        const char* const * const querypmbeg,
        const char* const * const querypmend,
        const char* const * const bdbCpmbeg,
        const char* const * const bdbCpmend,
        uint nqystrs, uint ndbCstrs,
        uint nqyposs, uint ndbCposs,
        uint qystr1len, uint dbstr1len,
        uint dbxpad,
        float* const tmpdpdiagbuffers, float* const tmpdpbotbuffer,
        float* const tmpdpalnpossbuffer, uint* const maxscoordsbuf, char* const btckdata,
        float* const wrkmem, float* const wrkmemccd, float* const wrkmemtm, float* const wrkmemtmibest,
        float* const wrkmemaux, float* const wrkmem2, float* const alndatamem, float* const tfmmem,
        uint* const globvarsbuf)
    :
        maxnsteps_(maxnsteps),
        querypmbeg_(querypmbeg), querypmend_(querypmend),
        bdbCpmbeg_(bdbCpmbeg), bdbCpmend_(bdbCpmend),
        nqystrs_(nqystrs), ndbCstrs_(ndbCstrs),
        nqyposs_(nqyposs), ndbCposs_(ndbCposs),
        qystr1len_(qystr1len), dbstr1len_(dbstr1len),
        dbxpad_(dbxpad),
        tmpdpdiagbuffers_(tmpdpdiagbuffers), tmpdpbotbuffer_(tmpdpbotbuffer),
        tmpdpalnpossbuffer_(tmpdpalnpossbuffer),
        maxscoordsbuf_(maxscoordsbuf), btckdata_(btckdata),
        wrkmem_(wrkmem), wrkmemccd_(wrkmemccd), wrkmemtm_(wrkmemtm), wrkmemtmibest_(wrkmemtmibest),
        wrkmemaux_(wrkmemaux), wrkmem2_(wrkmem2), alndatamem_(alndatamem), tfmmem_(tfmmem),
        globvarsbuf_(globvarsbuf)
    {}

public:
    template<bool ANCHOR, bool BANDED, bool GAP0, int D02IND, bool ALTSCTMS = false>
    /**
     * @brief 在CPU 动态规划中并行计算 `ExecDPwBtck128xKernel` 对应的数据。
     * @param gapopencost 供该函数读取或更新的 `gapopencost` 参数。
     * @param stepnumber 供该函数读取或更新的 `stepnumber` 参数。
     * @param querypmbeg 查询结构打包字段的起始指针数组。
     * @param bdbCpmbeg 参考结构打包字段的起始指针数组。
     * @param wrkmemtmibest 保存各候选当前最佳刚体变换的缓冲区。
     * @param wrkmemaux 保存分数、收敛标记等辅助状态的工作缓冲区。
     * @param tmpdpdiagbuffers 供当前步骤读取或更新的 `tmpdpdiagbuffers` 缓冲区。
     * @param tmpdpbotbuffer 供当前步骤读取或更新的 `tmpdpbotbuffer` 缓冲区。
     * @param btckdata 保存动态规划回溯方向的缓冲区。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void ExecDPwBtck128xKernel(
        const float gapopencost,
        const int stepnumber,
        const char* const * const __RESTRICT__ querypmbeg,
        const char* const * const __RESTRICT__ bdbCpmbeg,
        const float* const __RESTRICT__ wrkmemtmibest,
        const float* const __RESTRICT__ wrkmemaux,
        float* const __RESTRICT__ tmpdpdiagbuffers,
        float* const __RESTRICT__ tmpdpbotbuffer,
    //     uint* const __RESTRICT__ maxscoordsbuf,
        char* const __RESTRICT__ btckdata
    );

    template<bool GAP0>
    /**
     * @brief 在CPU 动态规划中并行计算 `ExecDPScore128xKernel` 对应的数据。
     * @param gapopencost 供该函数读取或更新的 `gapopencost` 参数。
     * @param querypmbeg 查询结构打包字段的起始指针数组。
     * @param bdbCpmbeg 参考结构打包字段的起始指针数组。
     * @param wrkmemtm 保存候选刚体变换的工作缓冲区。
     * @param wrkmemaux 保存分数、收敛标记等辅助状态的工作缓冲区。
     * @param tmpdpdiagbuffers 供当前步骤读取或更新的 `tmpdpdiagbuffers` 缓冲区。
     * @param tmpdpbotbuffer 供当前步骤读取或更新的 `tmpdpbotbuffer` 缓冲区。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void ExecDPScore128xKernel(
        const float gapopencost,
        const char* const * const __RESTRICT__ querypmbeg,
        const char* const * const __RESTRICT__ bdbCpmbeg,
        const float* const __RESTRICT__ wrkmemtm,
        const float* const __RESTRICT__ wrkmemaux,
        float* const __RESTRICT__ tmpdpdiagbuffers,
        float* const __RESTRICT__ tmpdpbotbuffer);

    template<bool GLOBTFM, bool GAP0, bool USESS, int D02IND>
    /**
     * @brief 在CPU 动态规划中并行计算 `ExecDPTFMSSwBtck128xKernel` 对应的数据。
     * @param gapopencost 供该函数读取或更新的 `gapopencost` 参数。
     * @param ssweight 供该函数读取或更新的 `ssweight` 参数。
     * @param stepnumber 供该函数读取或更新的 `stepnumber` 参数。
     * @param querypmbeg 查询结构打包字段的起始指针数组。
     * @param bdbCpmbeg 参考结构打包字段的起始指针数组。
     * @param tfmmem 保存最终刚体变换矩阵的缓冲区。
     * @param wrkmemaux 保存分数、收敛标记等辅助状态的工作缓冲区。
     * @param tmpdpdiagbuffers 供当前步骤读取或更新的 `tmpdpdiagbuffers` 缓冲区。
     * @param tmpdpbotbuffer 供当前步骤读取或更新的 `tmpdpbotbuffer` 缓冲区。
     * @param btckdata 保存动态规划回溯方向的缓冲区。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void ExecDPTFMSSwBtck128xKernel(
        const float gapopencost,
        const float ssweight,
        const int stepnumber,
        const char* const * const __RESTRICT__ querypmbeg,
        const char* const * const __RESTRICT__ bdbCpmbeg,
        const float* const __RESTRICT__ tfmmem,
        const float* const __RESTRICT__ wrkmemaux,
        float* const __RESTRICT__ tmpdpdiagbuffers,
        float* const __RESTRICT__ tmpdpbotbuffer,
    //     uint* const __RESTRICT__ maxscoordsbuf,
        char* const __RESTRICT__ btckdata
    );

    template<bool USESEQSCORING>
    /**
     * @brief 在CPU 动态规划中并行计算 `ExecDPSSwBtck128xKernel` 对应的数据。
     * @param gapopencost 供该函数读取或更新的 `gapopencost` 参数。
     * @param weight4ss 供该函数读取或更新的 `weight4ss` 参数。
     * @param weight4rr 供该函数读取或更新的 `weight4rr` 参数。
     * @param querypmbeg 查询结构打包字段的起始指针数组。
     * @param bdbCpmbeg 参考结构打包字段的起始指针数组。
     * @param wrkmemaux 保存分数、收敛标记等辅助状态的工作缓冲区。
     * @param tmpdpdiagbuffers 供当前步骤读取或更新的 `tmpdpdiagbuffers` 缓冲区。
     * @param tmpdpbotbuffer 供当前步骤读取或更新的 `tmpdpbotbuffer` 缓冲区。
     * @param btckdata 保存动态规划回溯方向的缓冲区。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void ExecDPSSwBtck128xKernel(
        const float gapopencost,
        const float weight4ss,
        const float weight4rr,
        const char* const * const __RESTRICT__ querypmbeg,
        const char* const * const __RESTRICT__ bdbCpmbeg,
        const float* const __RESTRICT__ wrkmemaux,
        float* const __RESTRICT__ tmpdpdiagbuffers,
        float* const __RESTRICT__ tmpdpbotbuffer,
        char* const __RESTRICT__ btckdata);

    template<int NASEQ = 0>
    /**
     * @brief 在CPU 动态规划中并行计算 `ExecDPSSLocal128xKernel` 对应的数据。
     * @param gapcost 供该函数读取或更新的 `gapcost` 参数。
     * @param querypmbeg 查询结构打包字段的起始指针数组。
     * @param bdbCpmbeg 参考结构打包字段的起始指针数组。
     * @param wrkmemaux 保存分数、收敛标记等辅助状态的工作缓冲区。
     * @param tmpdpdiagbuffers 供当前步骤读取或更新的 `tmpdpdiagbuffers` 缓冲区。
     * @param tmpdpbotbuffer 供当前步骤读取或更新的 `tmpdpbotbuffer` 缓冲区。
     * @param dpscoremtx 当前步骤使用或写回的 `dpscoremtx` 分数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void ExecDPSSLocal128xKernel(
        const float gapcost,
        const char* const * const __RESTRICT__ querypmbeg,
        const char* const * const __RESTRICT__ bdbCpmbeg,
        const float* const __RESTRICT__ wrkmemaux,
        float* const __RESTRICT__ tmpdpdiagbuffers,
        float* const __RESTRICT__ tmpdpbotbuffer,
        char* const __RESTRICT__ dpscoremtx);


    template<bool ANCHORRGN, bool BANDED>
    /**
     * @brief 在CPU 动态规划中并行计算 `BtckToMatched128xKernel` 对应的数据。
     * @param stepnumber 供该函数读取或更新的 `stepnumber` 参数。
     * @param querypmbeg 查询结构打包字段的起始指针数组。
     * @param bdbCpmbeg 参考结构打包字段的起始指针数组。
     * @param btckdata 保存动态规划回溯方向的缓冲区。
     * @param wrkmemaux 保存分数、收敛标记等辅助状态的工作缓冲区。
     * @param tmpdpalnpossbuffer 供当前步骤读取或更新的 `tmpdpalnpossbuffer` 缓冲区。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void BtckToMatched128xKernel(
        const uint stepnumber,
        const char* const * const __RESTRICT__ querypmbeg,
        const char* const * const __RESTRICT__ bdbCpmbeg,
        const char* const __RESTRICT__ btckdata,
        float* const __RESTRICT__ wrkmemaux,
        float* const __RESTRICT__ tmpdpalnpossbuffer);

    /**
     * @brief 在CPU 动态规划中并行计算 `ConstrainedBtckToMatched128xKernel` 对应的数据。
     * @param querypmbeg 查询结构打包字段的起始指针数组。
     * @param bdbCpmbeg 参考结构打包字段的起始指针数组。
     * @param btckdata 保存动态规划回溯方向的缓冲区。
     * @param tfmmem 保存最终刚体变换矩阵的缓冲区。
     * @param wrkmemaux 保存分数、收敛标记等辅助状态的工作缓冲区。
     * @param tmpdpalnpossbuffer 供当前步骤读取或更新的 `tmpdpalnpossbuffer` 缓冲区。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void ConstrainedBtckToMatched128xKernel(
        const char* const * const __RESTRICT__ querypmbeg,
        const char* const * const __RESTRICT__ bdbCpmbeg,
        const char* const __RESTRICT__ btckdata,
        const float* const __RESTRICT__ tfmmem,
        float* const __RESTRICT__ wrkmemaux,
        float* const __RESTRICT__ tmpdpalnpossbuffer);


    /**
     * @brief 在CPU 动态规划中并行计算 `ProductionMatchToAlignment128xKernel` 对应的数据。
     * @param nodeletions 控制当前步骤范围或规模的 `nodeletions`。
     * @param d2equiv 供该函数读取或更新的 `d2equiv` 参数。
     * @param querypmbeg 查询结构打包字段的起始指针数组。
     * @param bdbCpmbeg 参考结构打包字段的起始指针数组。
     * @param tmpdpalnpossbuffer 供当前步骤读取或更新的 `tmpdpalnpossbuffer` 缓冲区。
     * @param wrkmemaux 保存分数、收敛标记等辅助状态的工作缓冲区。
     * @param alndatamem 保存最终对齐统计量的缓冲区。
     * @param alnsmem 供当前步骤读取或更新的 `alnsmem` 缓冲区。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void ProductionMatchToAlignment128xKernel(
        const bool nodeletions, const float d2equiv,
        const char* const * const __RESTRICT__ querypmbeg,
        const char* const * const __RESTRICT__ bdbCpmbeg,
        const float* const __RESTRICT__ tmpdpalnpossbuffer,
        const float* const __RESTRICT__ wrkmemaux,
        float* const __RESTRICT__ alndatamem,
        char* const __RESTRICT__ alnsmem);

protected:
    // {{---------------------------------------------
    template<int DIMD, int DATALN, int DEFV = 0>
    /**
     * @brief 在CPU 动态规划中读取 `ReadQryRE` 对应的数据。
     * @param x 供该函数读取或更新的 `x` 参数。
     * @param y 供该函数读取或更新的 `y` 参数。
     * @param qrydst 描述查询结构的 `qrydst`。
     * @param qrylen 控制当前步骤范围或规模的 `qrylen`。
     * @param querypmbeg 查询结构打包字段的起始指针数组。
     * @param qryRE 描述查询结构的 `qryRE`。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void ReadQryRE(
        const int /*x*/, const int y,
        const int qrydst, const int qrylen,
        const char* const * const __RESTRICT__ querypmbeg,
        char* const __RESTRICT__ qryRE);

    template<int DIMD, int DATALN>
    /**
     * @brief 在CPU 动态规划中读取 `ReadQrySS` 对应的数据。
     * @param x 供该函数读取或更新的 `x` 参数。
     * @param y 供该函数读取或更新的 `y` 参数。
     * @param qrydst 描述查询结构的 `qrydst`。
     * @param qrylen 控制当前步骤范围或规模的 `qrylen`。
     * @param querypmbeg 查询结构打包字段的起始指针数组。
     * @param qrySS 描述查询结构的 `qrySS`。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void ReadQrySS(
        const int /*x*/, const int y,
        const int qrydst, const int qrylen,
        const char* const * const __RESTRICT__ querypmbeg,
        char* const __RESTRICT__ qrySS);

    template<int DIMDpX, int DATALN>
    /**
     * @brief 在CPU 动态规划中读取 `ReadRfnSS` 对应的数据。
     * @param x 供该函数读取或更新的 `x` 参数。
     * @param y 供该函数读取或更新的 `y` 参数。
     * @param dbstrdst 描述参考结构的 `dbstrdst`。
     * @param dbstrlen 控制当前步骤范围或规模的 `dbstrlen`。
     * @param bdbCpmbeg 参考结构打包字段的起始指针数组。
     * @param rfnSS 描述参考结构的 `rfnSS`。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void ReadRfnSS(
        const int x, const int /*y*/,
        const int dbstrdst, const int dbstrlen,
        const char* const * const __RESTRICT__ bdbCpmbeg,
        char* const __RESTRICT__ rfnSS);

    template<int DIMDpX, int DATALN>
    /**
     * @brief 在CPU 动态规划中读取 `ReadRfnRE` 对应的数据。
     * @param x 供该函数读取或更新的 `x` 参数。
     * @param y 供该函数读取或更新的 `y` 参数。
     * @param dbstrdst 描述参考结构的 `dbstrdst`。
     * @param dbstrlen 控制当前步骤范围或规模的 `dbstrlen`。
     * @param bdbCpmbeg 参考结构打包字段的起始指针数组。
     * @param rfnRE 描述参考结构的 `rfnRE`。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void ReadRfnRE(
        const int x, const int /*y*/,
        const int dbstrdst, const int dbstrlen,
        const char* const * const __RESTRICT__ bdbCpmbeg,
        char* const __RESTRICT__ rfnRE);
    // }}---------------------------------------------


    // {{---------------------------------------------
    template<int DIMD, int DATALN>
    /**
     * @brief 在CPU 动态规划中读取 `ReadAndTransformQryCoords` 对应的数据。
     * @param x 供该函数读取或更新的 `x` 参数。
     * @param y 供该函数读取或更新的 `y` 参数。
     * @param qrydst 描述查询结构的 `qrydst`。
     * @param qrylen 控制当前步骤范围或规模的 `qrylen`。
     * @param querypmbeg 查询结构打包字段的起始指针数组。
     * @param tfm 表示或保存刚体变换的 `tfm`。
     * @param qryCoords 描述查询结构的 `qryCoords`。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void ReadAndTransformQryCoords(
        const int x, const int y,
        const int qrydst, const int qrylen,
        const char* const * const __RESTRICT__ querypmbeg,
        const float* const __RESTRICT__ tfm,
        float* const __RESTRICT__ qryCoords);

    template<int DIMDpX, int DATALN>
    /**
     * @brief 在CPU 动态规划中读取 `ReadRfnCoords` 对应的数据。
     * @param x 供该函数读取或更新的 `x` 参数。
     * @param y 供该函数读取或更新的 `y` 参数。
     * @param dbstrdst 描述参考结构的 `dbstrdst`。
     * @param dbstrlen 控制当前步骤范围或规模的 `dbstrlen`。
     * @param bdbCpmbeg 参考结构打包字段的起始指针数组。
     * @param rfnCoords 描述参考结构的 `rfnCoords`。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void ReadRfnCoords(
        const int x, const int /*y*/,
        const int dbstrdst, const int dbstrlen,
        const char* const * const __RESTRICT__ bdbCpmbeg,
        float* const __RESTRICT__ rfnCoords);
    // }}---------------------------------------------


    // {{---------------------------------------------
    template<int DIMD, int DIMD1, int DATALN>
    /**
     * @brief 在CPU 动态规划中读取 `ReadTwoDiagonals` 对应的数据。
     * @param x 供该函数读取或更新的 `x` 参数。
     * @param y 供该函数读取或更新的 `y` 参数。
     * @param dbstrdst 描述参考结构的 `dbstrdst`。
     * @param dbstrlen 控制当前步骤范围或规模的 `dbstrlen`。
     * @param dblen 控制当前步骤范围或规模的 `dblen`。
     * @param yofff 供该函数读取或更新的 `yofff` 参数。
     * @param tmpdpdiagbuffers 供当前步骤读取或更新的 `tmpdpdiagbuffers` 缓冲区。
     * @param diag1 供该函数读取或更新的 `diag1` 参数。
     * @param diag2 供该函数读取或更新的 `diag2` 参数。
     * @param value 需要读取、写入或转换的值。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void ReadTwoDiagonals(
        const int x, const int /*y*/,
        const int dbstrdst, const int dbstrlen,
        const int dblen, const int yofff,
        const float* const __RESTRICT__ tmpdpdiagbuffers,
        float* const __RESTRICT__ diag1,
        float* const __RESTRICT__ diag2,
        const float value = 0.0f);

    template<int DIMD, int DIMD1, int DATALN>
    /**
     * @brief 在CPU 动态规划中写出 `WriteTwoDiagonals` 对应的数据。
     * @param DIMX 供该函数读取或更新的 `DIMX` 参数。
     * @param x 供该函数读取或更新的 `x` 参数。
     * @param y 供该函数读取或更新的 `y` 参数。
     * @param dbstrdst 描述参考结构的 `dbstrdst`。
     * @param dbstrlen 控制当前步骤范围或规模的 `dbstrlen`。
     * @param dblen 控制当前步骤范围或规模的 `dblen`。
     * @param yofff 供该函数读取或更新的 `yofff` 参数。
     * @param tmpdpdiagbuffers 供当前步骤读取或更新的 `tmpdpdiagbuffers` 缓冲区。
     * @param diag1 供该函数读取或更新的 `diag1` 参数。
     * @param diag2 供该函数读取或更新的 `diag2` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void WriteTwoDiagonals(
        const int DIMX,
        const int x, const int /*y*/,
        const int dbstrdst, const int dbstrlen,
        const int dblen, const int yofff,
        float* const __RESTRICT__ tmpdpdiagbuffers,
        const float* const __RESTRICT__ diag1,
        const float* const __RESTRICT__ diag2);

    template<int DIMD, int DIMX, int DATALN>
    /**
     * @brief 在CPU 动态规划中读取 `ReadBottomEdge` 对应的数据。
     * @param x 供该函数读取或更新的 `x` 参数。
     * @param y 供该函数读取或更新的 `y` 参数。
     * @param dbstrdst 描述参考结构的 `dbstrdst`。
     * @param dbstrlen 控制当前步骤范围或规模的 `dbstrlen`。
     * @param dblen 控制当前步骤范围或规模的 `dblen`。
     * @param yofff 供该函数读取或更新的 `yofff` 参数。
     * @param tmpdpbotbuffer 供当前步骤读取或更新的 `tmpdpbotbuffer` 缓冲区。
     * @param bottm 供该函数读取或更新的 `bottm` 参数。
     * @param value 需要读取、写入或转换的值。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void ReadBottomEdge(
        const int x, const int y,
        const int dbstrdst, const int dbstrlen,
        const int dblen, const int yofff,
        const float* const __RESTRICT__ tmpdpbotbuffer,
        float* const __RESTRICT__ bottm,
        const float value = 0.0f);

    template<int DIMX, int DATALN>
    /**
     * @brief 在CPU 动态规划中写出 `WriteBottomEdge` 对应的数据。
     * @param x 供该函数读取或更新的 `x` 参数。
     * @param y 供该函数读取或更新的 `y` 参数。
     * @param dbstrdst 描述参考结构的 `dbstrdst`。
     * @param dbstrlen 控制当前步骤范围或规模的 `dbstrlen`。
     * @param dblen 控制当前步骤范围或规模的 `dblen`。
     * @param yofff 供该函数读取或更新的 `yofff` 参数。
     * @param tmpdpbotbuffer 供当前步骤读取或更新的 `tmpdpbotbuffer` 缓冲区。
     * @param bottm 供该函数读取或更新的 `bottm` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void WriteBottomEdge(
        const int x, const int /*y*/,
        const int dbstrdst, const int dbstrlen,
        const int dblen, const int yofff,
        float* const __RESTRICT__ tmpdpbotbuffer,
        const float* const __RESTRICT__ bottm);

    template<int DIMD, int DIMX, int DATALN>
    /**
     * @brief 在CPU 动态规划中写出 `WriteBtckInfo` 对应的数据。
     * @param x 供该函数读取或更新的 `x` 参数。
     * @param y 供该函数读取或更新的 `y` 参数。
     * @param qrydst 描述查询结构的 `qrydst`。
     * @param qrylen 控制当前步骤范围或规模的 `qrylen`。
     * @param dbstrdst 描述参考结构的 `dbstrdst`。
     * @param dbstrlen 控制当前步骤范围或规模的 `dbstrlen`。
     * @param dblen 控制当前步骤范围或规模的 `dblen`。
     * @param btckdata 保存动态规划回溯方向的缓冲区。
     * @param DIMX 供该函数读取或更新的 `DIMX` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void WriteBtckInfo(
        const int x, const int y,
        const int qrydst, const int qrylen,
        const int dbstrdst, const int dbstrlen,
        const int dblen,
        char* const __RESTRICT__ btckdata,
        const char (* const __RESTRICT__ btck)[DIMX]);
    // }}---------------------------------------------

protected:
    /**
     * @brief 在CPU 动态规划中处理 `lfNdx` 对应的数据。
     * @param i 供该函数读取或更新的 `i` 参数。
     * @param j 供该函数读取或更新的 `j` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    int lfNdx(int i,int j) {return pmv2DNoElems * j + i;}
    /**
     * @brief 在CPU 动态规划中处理 `lfDgNdx` 对应的数据。
     * @param i 供该函数读取或更新的 `i` 参数。
     * @param j 供该函数读取或更新的 `j` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    template<int DIM> int lfDgNdx(int i,int j) {return DIM * i + j;}

protected:
    // {{---------------------------------------------
    template<int DIMX, int DATALN>
    /**
     * @brief 在CPU 动态规划中写出 `WriteAlignmentFragment` 对应的数据。
     * @param qrydst 描述查询结构的 `qrydst`。
     * @param dbstrdst 描述参考结构的 `dbstrdst`。
     * @param alnofff 供该函数读取或更新的 `alnofff` 参数。
     * @param dbalnlen 控制当前步骤范围或规模的 `dbalnlen`。
     * @param dbalnbeg 描述参考结构的 `dbalnbeg`。
     * @param written 供该函数读取或更新的 `written` 参数。
     * @param lentowrite 控制当前步骤范围或规模的 `lentowrite`。
     * @param lentocheck 控制当前步骤范围或规模的 `lentocheck`。
     * @param querypmbeg 查询结构打包字段的起始指针数组。
     * @param bdbCpmbeg 参考结构打包字段的起始指针数组。
     * @param tmp 供该函数读取或更新的 `tmp` 参数。
     * @param DIMX 供该函数读取或更新的 `DIMX` 参数。
     * @param alnsmem 供当前步骤读取或更新的 `alnsmem` 缓冲区。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    int WriteAlignmentFragment(
        const int qrydst, const int dbstrdst,
        const int alnofff, const int dbalnlen, const int dbalnbeg,
        const int written, const int lentowrite, const int lentocheck,
        const char* const * const __RESTRICT__ querypmbeg,
        const char* const * const __RESTRICT__ bdbCpmbeg,
        int* __RESTRICT__ tmp,
        int (* __RESTRICT__ outAln)[DIMX],
        char* const __RESTRICT__ alnsmem);
    // }}---------------------------------------------

protected:
    const uint maxnsteps_;
    const char* const * const querypmbeg_, * const * const querypmend_;
    const char* const * const bdbCpmbeg_, * const *const bdbCpmend_;
    const uint nqystrs_, ndbCstrs_;
    const uint nqyposs_, ndbCposs_;
    const uint qystr1len_, dbstr1len_;
    const uint dbxpad_;
    float* const tmpdpdiagbuffers_, *const tmpdpbotbuffer_, *const tmpdpalnpossbuffer_;
    uint* const maxscoordsbuf_;
    char* const btckdata_;
    float* const wrkmem_, *const wrkmemccd_, *const wrkmemtm_, *const wrkmemtmibest_;
    float* const wrkmemaux_, *const wrkmem2_, *const alndatamem_, *const tfmmem_;
    uint* const globvarsbuf_;
};



// -------------------------------------------------------------------------
// INLINES ...
// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
// ReadQryRE: read query residues;
// x, y, bottom-left corner of an oblique block under consideration;
// qrydst, distance to query data;
// qrylen, query length;
// querypmbeg, query chunk data;
// qryRE, query residues read;
//
template<int DIMD, int DATALN, int DEFV>
/**
 * @brief 在CPU 动态规划中读取 `MpDPHub::ReadQryRE` 对应的数据。
 * @param x 供该函数读取或更新的 `x` 参数。
 * @param y 供该函数读取或更新的 `y` 参数。
 * @param qrydst 描述查询结构的 `qrydst`。
 * @param qrylen 控制当前步骤范围或规模的 `qrylen`。
 * @param querypmbeg 查询结构打包字段的起始指针数组。
 * @param qryRE 描述查询结构的 `qryRE`。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
inline
void MpDPHub::ReadQryRE(
    const int /*x*/, const int y,
    const int qrydst, const int qrylen,
    const char* const * const __RESTRICT__ querypmbeg,
    char* const __RESTRICT__ qryRE)
{
    #pragma omp simd
    for(int pi = 0; pi < DIMD; pi++)
        qryRE[pi] = DEFV;

    int qposbeg = mymax(y - qrylen + 1, 0);
    int qposend = mymin(y, DIMD - 1);
    #pragma omp simd aligned(querypmbeg:DATALN)
    for(int pi = qposbeg; pi <= qposend; pi++)
    {   //(0 <= y - pi && y - pi < qrylen)
        qryRE[pi] = PMBatchStrData::GetFieldAt<char,pmv2Drsd>(querypmbeg, qrydst + y - pi);
    }
}

// -------------------------------------------------------------------------
// ReadQrySS: read query secondaryy structure;
// x, y, bottom-left corner of an oblique block under consideration;
// qrydst, distance to query data;
// qrylen, query length;
// querypmbeg, query chunk data;
// qrySS, read query secondary structure;
//
template<int DIMD, int DATALN>
/**
 * @brief 在CPU 动态规划中读取 `MpDPHub::ReadQrySS` 对应的数据。
 * @param x 供该函数读取或更新的 `x` 参数。
 * @param y 供该函数读取或更新的 `y` 参数。
 * @param qrydst 描述查询结构的 `qrydst`。
 * @param qrylen 控制当前步骤范围或规模的 `qrylen`。
 * @param querypmbeg 查询结构打包字段的起始指针数组。
 * @param qrySS 描述查询结构的 `qrySS`。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
inline
void MpDPHub::ReadQrySS(
    const int /*x*/, const int y,
    const int qrydst, const int qrylen,
    const char* const * const __RESTRICT__ querypmbeg,
    char* const __RESTRICT__ qrySS)
{
    #pragma omp simd
    for(int pi = 0; pi < DIMD; pi++)
        qrySS[pi] = pmvLOOP;

    int qposbeg = mymax(y - qrylen + 1, 0);
    int qposend = mymin(y, DIMD - 1);
    #pragma omp simd aligned(querypmbeg:DATALN)
    for(int pi = qposbeg; pi <= qposend; pi++)
    {   //(0 <= y - pi && y - pi < qrylen)
        qrySS[pi] = PMBatchStrData::GetFieldAt<char,pmv2Dss>(querypmbeg, qrydst + y - pi);
    }
}

// -------------------------------------------------------------------------
// ReadRfnSS: read reference secondary structure;
// x, y, bottom-left corner of an oblique block under consideration;
// dbstrdst, distance to reference data;
// dbstrlen, reference length;
// bdbCpmbeg, reference chunk data;
// rfnSS, reference secondary structure read;
//
template<int DIMDpX, int DATALN>
/**
 * @brief 在CPU 动态规划中读取 `MpDPHub::ReadRfnSS` 对应的数据。
 * @param x 供该函数读取或更新的 `x` 参数。
 * @param y 供该函数读取或更新的 `y` 参数。
 * @param dbstrdst 描述参考结构的 `dbstrdst`。
 * @param dbstrlen 控制当前步骤范围或规模的 `dbstrlen`。
 * @param bdbCpmbeg 参考结构打包字段的起始指针数组。
 * @param rfnSS 描述参考结构的 `rfnSS`。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
inline
void MpDPHub::ReadRfnSS(
    const int x, const int /*y*/,
    const int dbstrdst, const int dbstrlen,
    const char* const * const __RESTRICT__ bdbCpmbeg,
    char* const __RESTRICT__ rfnSS)
{
    //db reference structure position corresponding to the oblique block's
    //bottom-left corner in the buffers dc_pm2dvfields_ (score matrix) 
    //plus the offset determined by thread id:
    int dbpos = x + dbstrdst;//going right

    #pragma omp simd
    for(int pi = 0; pi < DIMDpX; pi++)
        rfnSS[pi] = 0;//value

    int dbposbeg = mymax(-x, 0);
    int dbposend = mymin(dbstrlen - 1 - x, DIMDpX - 1);
    #pragma omp simd aligned(bdbCpmbeg:DATALN)
    for(int pi = dbposbeg; pi <= dbposend; pi++)
    {   //(0 <= x + pi && x + pi < dbstrlen)
        rfnSS[pi] = PMBatchStrData::GetFieldAt<char,pmv2Dss>(bdbCpmbeg, dbpos + pi);
    }
}

// -------------------------------------------------------------------------
// ReadRfnRE: read reference residues;
// x, y, bottom-left corner of an oblique block under consideration;
// dbstrdst, distance to reference data;
// dbstrlen, reference length;
// bdbCpmbeg, reference chunk data;
// rfnRE, reference residues read;
//
template<int DIMDpX, int DATALN>
/**
 * @brief 在CPU 动态规划中读取 `MpDPHub::ReadRfnRE` 对应的数据。
 * @param x 供该函数读取或更新的 `x` 参数。
 * @param y 供该函数读取或更新的 `y` 参数。
 * @param dbstrdst 描述参考结构的 `dbstrdst`。
 * @param dbstrlen 控制当前步骤范围或规模的 `dbstrlen`。
 * @param bdbCpmbeg 参考结构打包字段的起始指针数组。
 * @param rfnRE 描述参考结构的 `rfnRE`。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
inline
void MpDPHub::ReadRfnRE(
    const int x, const int /*y*/,
    const int dbstrdst, const int dbstrlen,
    const char* const * const __RESTRICT__ bdbCpmbeg,
    char* const __RESTRICT__ rfnRE)
{
    //db reference structure position corresponding to the oblique block's
    //bottom-left corner in the buffers dc_pm2dvfields_ (score matrix) 
    //plus the offset determined by thread id:
    int dbpos = x + dbstrdst;//going right

    #pragma omp simd
    for(int pi = 0; pi < DIMDpX; pi++)
        rfnRE[pi] = 0;//value

    int dbposbeg = mymax(-x, 0);
    int dbposend = mymin(dbstrlen - 1 - x, DIMDpX - 1);
    #pragma omp simd aligned(bdbCpmbeg:DATALN)
    for(int pi = dbposbeg; pi <= dbposend; pi++)
    {   //(0 <= x + pi && x + pi < dbstrlen)
        rfnRE[pi] = PMBatchStrData::GetFieldAt<char,pmv2Drsd>(bdbCpmbeg, dbpos + pi);
    }
}



// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
// ReadAndTransformQryCoords: read and optionally transform query coordinates;
// x, y, bottom-left corner of an oblique block under consideration;
// qrydst, distance to query data;
// qrylen, query length;
// querypmbeg, query chunk data;
// tfm, transformation matrix;
// qryCoords, (transformed and) read query coordinates;
//
template<int DIMD, int DATALN>
/**
 * @brief 在CPU 动态规划中读取 `MpDPHub::ReadAndTransformQryCoords` 对应的数据。
 * @param x 供该函数读取或更新的 `x` 参数。
 * @param y 供该函数读取或更新的 `y` 参数。
 * @param qrydst 描述查询结构的 `qrydst`。
 * @param qrylen 控制当前步骤范围或规模的 `qrylen`。
 * @param querypmbeg 查询结构打包字段的起始指针数组。
 * @param tfm 表示或保存刚体变换的 `tfm`。
 * @param qryCoords 描述查询结构的 `qryCoords`。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
inline
void MpDPHub::ReadAndTransformQryCoords(
    const int /*x*/, const int y,
    const int qrydst, const int qrylen,
    const char* const * const __RESTRICT__ querypmbeg,
    const float* const __RESTRICT__ tfm,
    float* const __RESTRICT__ qryCoords)
{
    #pragma omp simd
    for(int pi = 0; pi < DIMD; pi++) {
        qryCoords[lfNdx(pmv2DX, pi)] = (float)(CUDP_DEFCOORD_QRY);
        qryCoords[lfNdx(pmv2DY, pi)] = (float)(CUDP_DEFCOORD_QRY);
        qryCoords[lfNdx(pmv2DZ, pi)] = (float)(CUDP_DEFCOORD_QRY);
    }
    int qposbeg = mymax(y - qrylen + 1, 0);
    int qposend = mymin(y, DIMD - 1);
    #pragma omp simd aligned(querypmbeg:DATALN)
    for(int pi = qposbeg; pi <= qposend; pi++)
    {   //(0 <= y - pi && y - pi < qrylen)
        qryCoords[lfNdx(pmv2DX, pi)] =
            PMBatchStrData::GetFieldAt<float,pmv2DCoords+pmv2DX>(querypmbeg, qrydst + y - pi);
        qryCoords[lfNdx(pmv2DY, pi)] =
            PMBatchStrData::GetFieldAt<float,pmv2DCoords+pmv2DY>(querypmbeg, qrydst + y - pi);
        qryCoords[lfNdx(pmv2DZ, pi)] =
            PMBatchStrData::GetFieldAt<float,pmv2DCoords+pmv2DZ>(querypmbeg, qrydst + y - pi);
    }
    if(tfm) {
        #pragma omp simd
        for(int pi = qposbeg; pi <= qposend; pi++)
        {   //(0 <= y - pi && y - pi < qrylen); transform the query fragment read
            transform_point(tfm,
                qryCoords[lfNdx(pmv2DX, pi)],
                qryCoords[lfNdx(pmv2DY, pi)],
                qryCoords[lfNdx(pmv2DZ, pi)]);
        }
    }
}

// -------------------------------------------------------------------------
// ReadRfnCoords: read reference coordinates;
// x, y, bottom-left corner of an oblique block under consideration;
// dbstrdst, distance to reference data;
// dbstrlen, reference length;
// bdbCpmbeg, reference chunk data;
// rfnCoords, read reference coordinates;
//
template<int DIMDpX, int DATALN>
/**
 * @brief 在CPU 动态规划中读取 `MpDPHub::ReadRfnCoords` 对应的数据。
 * @param x 供该函数读取或更新的 `x` 参数。
 * @param y 供该函数读取或更新的 `y` 参数。
 * @param dbstrdst 描述参考结构的 `dbstrdst`。
 * @param dbstrlen 控制当前步骤范围或规模的 `dbstrlen`。
 * @param bdbCpmbeg 参考结构打包字段的起始指针数组。
 * @param rfnCoords 描述参考结构的 `rfnCoords`。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
inline
void MpDPHub::ReadRfnCoords(
    const int x, const int /*y*/,
    const int dbstrdst, const int dbstrlen,
    const char* const * const __RESTRICT__ bdbCpmbeg,
    float* const __RESTRICT__ rfnCoords)
{
    //db reference structure position corresponding to the oblique block's
    //bottom-left corner in the buffers dc_pm2dvfields_ (score matrix) 
    //plus the offset determined by thread id:
    int dbpos = x + dbstrdst;//going right

    #pragma omp simd
    for(int pi = 0; pi < DIMDpX; pi++) {
        rfnCoords[lfNdx(pmv2DX, pi)] = (float)(CUDP_DEFCOORD_RFN);
        rfnCoords[lfNdx(pmv2DY, pi)] = (float)(CUDP_DEFCOORD_RFN);
        rfnCoords[lfNdx(pmv2DZ, pi)] = (float)(CUDP_DEFCOORD_RFN);
    }
    int dbposbeg = mymax(-x, 0);
    int dbposend = mymin(dbstrlen - 1 - x, DIMDpX - 1);
    #pragma omp simd aligned(bdbCpmbeg:DATALN)
    for(int pi = dbposbeg; pi <= dbposend; pi++)
    {   //(0 <= x + pi && x + pi < dbstrlen)
        rfnCoords[lfNdx(pmv2DX, pi)] =
            PMBatchStrData::GetFieldAt<float,pmv2DCoords+pmv2DX>(bdbCpmbeg, dbpos + pi);
        rfnCoords[lfNdx(pmv2DY, pi)] =
            PMBatchStrData::GetFieldAt<float,pmv2DCoords+pmv2DY>(bdbCpmbeg, dbpos + pi);
        rfnCoords[lfNdx(pmv2DZ, pi)] =
            PMBatchStrData::GetFieldAt<float,pmv2DCoords+pmv2DZ>(bdbCpmbeg, dbpos + pi);
    }
}



// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
// ReadTwoDiagonals: read two diagonal edges of oblique block under process;
// x, y, bottom-left corner of an oblique block under consideration;
// dbstrdst, distance to reference data;
// dbstrlen, reference length;
// dblen, total length in positions along the x-axis; 
// yofff, offset (w/o a factor) to the beginning of the data along the y axis;
// tmpdpdiagbuffers, memory region where diagonals are stored (1D, along the x axis);
// diag1, diag2, two diagonals of scores;
//
template<int DIMD, int DIMD1, int DATALN>
/**
 * @brief 在CPU 动态规划中读取 `MpDPHub::ReadTwoDiagonals` 对应的数据。
 * @param x 供该函数读取或更新的 `x` 参数。
 * @param y 供该函数读取或更新的 `y` 参数。
 * @param dbstrdst 描述参考结构的 `dbstrdst`。
 * @param dbstrlen 控制当前步骤范围或规模的 `dbstrlen`。
 * @param dblen 控制当前步骤范围或规模的 `dblen`。
 * @param yofff 供该函数读取或更新的 `yofff` 参数。
 * @param tmpdpdiagbuffers 供当前步骤读取或更新的 `tmpdpdiagbuffers` 缓冲区。
 * @param diag1 供该函数读取或更新的 `diag1` 参数。
 * @param diag2 供该函数读取或更新的 `diag2` 参数。
 * @param value 需要读取、写入或转换的值。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
inline
void MpDPHub::ReadTwoDiagonals(
    const int x, const int /*y*/,
    const int dbstrdst, const int dbstrlen,
    const int dblen, const int yofff,
    const float* const __RESTRICT__ tmpdpdiagbuffers,
    float* const __RESTRICT__ diag1,
    float* const __RESTRICT__ diag2,
    const float value)
{
    int dbpos = x + dbstrdst;
    int doffs = nTDPDiagScoreSections * nTDPDiagScoreSubsections * yofff;
    int doffsDelta = dpdssDiag2 * nTDPDiagScoreSubsections * dblen;

    #pragma omp simd collapse(2)
    for(int i = 0; i < nTDPDiagScoreSubsections; i++)
        for(int pi = 0; pi < DIMD; pi++) {
            diag1[lfDgNdx<DIMD1>(i, pi)] = value;
            diag2[lfDgNdx<DIMD1>(i, pi+1)] = value;
        }
    int dbposbeg = mymax(1 - x, 0);
    int dbposend = mymin(dbstrlen - x, DIMD - 1);
    #pragma omp simd collapse(2) aligned(tmpdpdiagbuffers:DATALN)
    for(int i = 0; i < nTDPDiagScoreSubsections; i++)
        for(int pi = dbposbeg; pi <= dbposend; pi++)
        {   //(0 <= x-1 + pi && x-1 + pi < dbstrlen)
            diag1[lfDgNdx<DIMD1>(i, pi)] = tmpdpdiagbuffers[doffs + i * dblen + dbpos + pi - 1];
            diag2[lfDgNdx<DIMD1>(i, pi+1)] =
                tmpdpdiagbuffers[(doffs + doffsDelta) + i * dblen + dbpos + pi - 1];
        }
}

// -------------------------------------------------------------------------
// WriteTwoDiagonals: write two diagonal edges of processed oblique block;
// x, y, bottom-left corner of an oblique block under consideration;
// dbstrdst, distance to reference data;
// dbstrlen, reference length;
// dblen, total length in positions along the x-axis; 
// yofff, offset (w/o a factor) to the beginning of the data along the y axis;
// tmpdpdiagbuffers, memory region where diagonals are stored (1D, along the x axis);
// diag1, diag2, two diagonals of scores;
//
template<int DIMD, int DIMD1, int DATALN>
/**
 * @brief 在CPU 动态规划中写出 `MpDPHub::WriteTwoDiagonals` 对应的数据。
 * @param DIMX 供该函数读取或更新的 `DIMX` 参数。
 * @param x 供该函数读取或更新的 `x` 参数。
 * @param y 供该函数读取或更新的 `y` 参数。
 * @param dbstrdst 描述参考结构的 `dbstrdst`。
 * @param dbstrlen 控制当前步骤范围或规模的 `dbstrlen`。
 * @param dblen 控制当前步骤范围或规模的 `dblen`。
 * @param yofff 供该函数读取或更新的 `yofff` 参数。
 * @param tmpdpdiagbuffers 供当前步骤读取或更新的 `tmpdpdiagbuffers` 缓冲区。
 * @param diag1 供该函数读取或更新的 `diag1` 参数。
 * @param diag2 供该函数读取或更新的 `diag2` 参数。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
inline
void MpDPHub::WriteTwoDiagonals(
    const int DIMX,
    const int x, const int /*y*/,
    const int dbstrdst, const int dbstrlen,
    const int dblen, const int yofff,
    float* const __RESTRICT__ tmpdpdiagbuffers,
    const float* const __RESTRICT__ diag1,
    const float* const __RESTRICT__ diag2)
{
    const int dbpos = x + dbstrdst;
    const int doffs = nTDPDiagScoreSections * nTDPDiagScoreSubsections * yofff;
    const int doffsDelta = dpdssDiag2 * nTDPDiagScoreSubsections * dblen;

    int dbposbeg = mymax(1 - x - DIMX, 0);
    int dbposend = mymin(dbstrlen - x - DIMX, DIMD - 1);
    #pragma omp simd collapse(2) aligned(tmpdpdiagbuffers:DATALN)
    for(int i = 0; i < nTDPDiagScoreSubsections; i++)
        for(int pi = dbposbeg; pi <= dbposend; pi++)
        {   //(0 <= x+DIMX-1 + pi && x+DIMX-1 + pi < dbstrlen)
            tmpdpdiagbuffers[doffs + i * dblen + dbpos + DIMX + pi - 1] = diag1[lfDgNdx<DIMD1>(i, pi)];
            tmpdpdiagbuffers[(doffs + doffsDelta) + i * dblen + dbpos + DIMX + pi - 1] =
                diag2[lfDgNdx<DIMD1>(i, pi+1)];
        }
}

// -------------------------------------------------------------------------
// ReadBottomEdge: read the bottom edge of the upper oblique blocks;
// x, y, bottom-left corner of an oblique block under consideration;
// dbstrdst, distance to reference data;
// dbstrlen, reference length;
// dblen, total length in positions along the x-axis; 
// yofff, offset (w/o a factor) to the beginning of the data along the y axis;
// tmpdpbotbuffer, memory region of stored oblique blocks' bottom lines (1D, along the x axis);
// bottm, bottom edge of scores;
//
template<int DIMD, int DIMX, int DATALN>
/**
 * @brief 在CPU 动态规划中读取 `MpDPHub::ReadBottomEdge` 对应的数据。
 * @param x 供该函数读取或更新的 `x` 参数。
 * @param y 供该函数读取或更新的 `y` 参数。
 * @param dbstrdst 描述参考结构的 `dbstrdst`。
 * @param dbstrlen 控制当前步骤范围或规模的 `dbstrlen`。
 * @param dblen 控制当前步骤范围或规模的 `dblen`。
 * @param yofff 供该函数读取或更新的 `yofff` 参数。
 * @param tmpdpbotbuffer 供当前步骤读取或更新的 `tmpdpbotbuffer` 缓冲区。
 * @param bottm 供该函数读取或更新的 `bottm` 参数。
 * @param value 需要读取、写入或转换的值。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
inline
    void MpDPHub::ReadBottomEdge(
    const int x, const int y,
    const int dbstrdst, const int dbstrlen,
    const int dblen, const int yofff,
    const float* const __RESTRICT__ tmpdpbotbuffer,
    float* const __RESTRICT__ bottm,
    const float value)
{
    //cache the bottom of the upper oblique blocks;
    int dbpos = x + dbstrdst;
    int doffs = nTDPBottomScoreSections * nTDPDiagScoreSubsections * yofff;
    #pragma omp simd collapse(2)
    for(int i = 0; i < nTDPDiagScoreSubsections; i++)
        for(int pi = 0; pi < DIMX; pi++) {
            bottm[lfDgNdx<DIMX>(i, pi)] = value;
        }
    if(DIMD <= y) {
        int dbposbeg = mymax(1 - DIMD - x, 0);
        int dbposend = mymin(dbstrlen - x - DIMD, DIMX - 1);
        #pragma omp simd collapse(2) aligned(tmpdpbotbuffer:DATALN)
        for(int i = 0; i < nTDPDiagScoreSubsections; i++)
            for(int pi = dbposbeg; pi <= dbposend; pi++)
            {   //(0 <= x+DIMD-1 + pi && x+DIMD-1 + pi < dbstrlen)
                bottm[lfDgNdx<DIMX>(i, pi)] = tmpdpbotbuffer[doffs + i * dblen + dbpos + DIMD - 1 + pi];
            }
    }
}

// -------------------------------------------------------------------------
// WriteBottomEdge: write the bottom edge of the oblique blocks;
// x, y, bottom-left corner of an oblique block under consideration;
// dbstrdst, distance to reference data;
// dbstrlen, reference length;
// dblen, total length in positions along the x-axis; 
// yofff, offset (w/o a factor) to the beginning of the data along the y axis;
// tmpdpbotbuffer, memory region of stored oblique blocks' bottom lines (1D, along the x axis);
// bottm, bottom edge of scores;
//
template<int DIMX, int DATALN>
/**
 * @brief 在CPU 动态规划中写出 `MpDPHub::WriteBottomEdge` 对应的数据。
 * @param x 供该函数读取或更新的 `x` 参数。
 * @param y 供该函数读取或更新的 `y` 参数。
 * @param dbstrdst 描述参考结构的 `dbstrdst`。
 * @param dbstrlen 控制当前步骤范围或规模的 `dbstrlen`。
 * @param dblen 控制当前步骤范围或规模的 `dblen`。
 * @param yofff 供该函数读取或更新的 `yofff` 参数。
 * @param tmpdpbotbuffer 供当前步骤读取或更新的 `tmpdpbotbuffer` 缓冲区。
 * @param bottm 供该函数读取或更新的 `bottm` 参数。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
inline
void MpDPHub::WriteBottomEdge(
    const int x, const int /*y*/,
    const int dbstrdst, const int dbstrlen,
    const int dblen, const int yofff,
    float* const __RESTRICT__ tmpdpbotbuffer,
    const float* const __RESTRICT__ bottm)
{
    const int dbpos = x + dbstrdst;
    const int doffs = nTDPBottomScoreSections * nTDPDiagScoreSubsections * yofff;

    int dbposbeg = mymax(-x, 0);
    int dbposend = mymin(dbstrlen - x - 1, DIMX - 1);
    #pragma omp simd collapse(2) aligned(tmpdpbotbuffer:DATALN)
    for(int i = 0; i < nTDPDiagScoreSubsections; i++)
        for(int pi = dbposbeg; pi <= dbposend; pi++)
        {   //(0 <= x + pi && x + pi < dbstrlen)
            tmpdpbotbuffer[doffs + i * dblen + dbpos + pi] = bottm[lfDgNdx<DIMX>(i, pi)];
        }
}

// -------------------------------------------------------------------------
// WriteBtckInfo: write backtracking information;
// x, y, bottom-left corner of an oblique block under consideration;
// qrydst, distance to query data;
// qrylen, query length;
// dbstrdst, distance to reference data;
// dbstrlen, reference length;
// dblen, total length in positions along the x-axis; 
// btckdata, memory region dedicated to backtracking;
// btck, backtracking information stored in cache;
//
template<int DIMD, int DIMX, int DATALN>
/**
 * @brief 在CPU 动态规划中写出 `MpDPHub::WriteBtckInfo` 对应的数据。
 * @param x 供该函数读取或更新的 `x` 参数。
 * @param y 供该函数读取或更新的 `y` 参数。
 * @param qrydst 描述查询结构的 `qrydst`。
 * @param qrylen 控制当前步骤范围或规模的 `qrylen`。
 * @param dbstrdst 描述参考结构的 `dbstrdst`。
 * @param dbstrlen 控制当前步骤范围或规模的 `dbstrlen`。
 * @param dblen 控制当前步骤范围或规模的 `dblen`。
 * @param btckdata 保存动态规划回溯方向的缓冲区。
 * @param DIMX 供该函数读取或更新的 `DIMX` 参数。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
inline
void MpDPHub::WriteBtckInfo(
    const int x, const int y,
    const int qrydst, const int qrylen,
    const int dbstrdst, const int dbstrlen,
    const int dblen,
    char* const __RESTRICT__ btckdata,
    const char (* const __RESTRICT__ btck)[DIMX])
{
    const int dbpos = x + dbstrdst;

    for(int i = 0; i < DIMD; i++) {
        if(0 <= y-i && y-i < qrylen)
        {   //going upwards
            //starting position of line i of the oblq. block in the matrix:
            int qpos = (qrydst + (y-i)) * dblen + i;
            int dbposbeg = mymax(-x - i, 0);
            int dbposend = mymin(dbstrlen - x - i - 1, DIMX - 1);
            #pragma omp simd aligned(btckdata:DATALN)
            for(int pi = dbposbeg; pi <= dbposend; pi++)
            {   //(0 <= x+i + pi && x+i + pi < dbstrlen)
                btckdata[qpos + dbpos + pi] = btck[i][pi];
            }
        }
    }
}



// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
// WriteAlignmentFragment: write alignment fragment to memory;
// written, alignment length written already;
// lentowrite, length of alignment fragment to write;
// lentocheck, alignment fragment length to check for modification;
// tmp, temporary array;
// outAln, alignment fragment cache;
// alnsmem, global memory of alignments;
// 
template<int DIMX, int DATALN>
/**
 * @brief 在CPU 动态规划中写出 `MpDPHub::WriteAlignmentFragment` 对应的数据。
 * @param qrydst 描述查询结构的 `qrydst`。
 * @param dbstrdst 描述参考结构的 `dbstrdst`。
 * @param alnofff 供该函数读取或更新的 `alnofff` 参数。
 * @param dbalnlen 控制当前步骤范围或规模的 `dbalnlen`。
 * @param dbalnbeg 描述参考结构的 `dbalnbeg`。
 * @param written 供该函数读取或更新的 `written` 参数。
 * @param lentowrite 控制当前步骤范围或规模的 `lentowrite`。
 * @param lentocheck 控制当前步骤范围或规模的 `lentocheck`。
 * @param querypmbeg 查询结构打包字段的起始指针数组。
 * @param bdbCpmbeg 参考结构打包字段的起始指针数组。
 * @param tmp 供该函数读取或更新的 `tmp` 参数。
 * @param DIMX 供该函数读取或更新的 `DIMX` 参数。
 * @param alnsmem 供当前步骤读取或更新的 `alnsmem` 缓冲区。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline
int MpDPHub::WriteAlignmentFragment(
    const int qrydst, const int dbstrdst,
    const int alnofff, const int dbalnlen, const int dbalnbeg,
    const int written, const int lentowrite, const int lentocheck,
    const char* const * const __RESTRICT__ querypmbeg,
    const char* const * const __RESTRICT__ bdbCpmbeg,
    int* __RESTRICT__ tmp,
    int (* __RESTRICT__ outAln)[DIMX],
    char* const __RESTRICT__ alnsmem)
{
    int wpos, idnts;

    #pragma omp simd aligned(querypmbeg:DATALN)
    for(int pi = 0; pi < lentocheck; pi++) {
        if((wpos = outAln[dp2oaQuery][pi]) <= 0) {
            outAln[dp2oaQuery][pi] = PMBatchStrData::GetFieldAt<char,pmv2Drsd>(querypmbeg, qrydst - wpos);
            outAln[dp2oaQuerySSS][pi] = PMBatchStrData::GetFieldAt<char,pmv2Dss>(querypmbeg, qrydst - wpos);
        }
    }
    #pragma omp simd aligned(bdbCpmbeg:DATALN)
    for(int pi = 0; pi < lentocheck; pi++) {
        if((wpos = outAln[dp2oaTarget][pi]) <= 0) {
            outAln[dp2oaTarget][pi] = PMBatchStrData::GetFieldAt<char,pmv2Drsd>(bdbCpmbeg, dbstrdst - wpos);
            outAln[dp2oaTargetSSS][pi] = PMBatchStrData::GetFieldAt<char,pmv2Dss>(bdbCpmbeg, dbstrdst - wpos);
        }
    }

    #pragma omp simd
    for(int pi = 0; pi < lentocheck; pi++) {
        int mc = outAln[dp2oaMiddle][pi];
        int r1 = outAln[dp2oaQuery][pi];
        int r2 = outAln[dp2oaTarget][pi];
        tmp[pi] = idnts = (/* r1 &&  */r1 != '-' && r1 == r2);
        if(idnts) outAln[dp2oaMiddle][pi] = (mc == '+')? r1: (r1|32)/*lower*/;
    }

    //write alignment fragment
    #pragma omp simd collapse(2) aligned(alnsmem:DATALN)
    for(int f = 0; f < nTDP2OutputAlignmentSSS; f++)
        for(int pi = 0; pi < lentowrite; pi++) {
            alnsmem[alnofff + dbalnbeg + pi + written + dbalnlen * f] =
                (char)(outAln[f][pi]);
        }

    //reduction for identities
    idnts = 0;
    #pragma omp simd reduction(+:idnts)
    for(int pi = 0; pi < lentocheck; pi++)
        idnts += tmp[pi];

    return idnts;
}

// -------------------------------------------------------------------------

#endif//__MpDPHub_h__
