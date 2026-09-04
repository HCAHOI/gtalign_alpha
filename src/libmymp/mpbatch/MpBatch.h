/***************************************************************************
 *   Copyright (C) 2021-2023 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __MpBatch_h__
#define __MpBatch_h__

#include "libutil/mybase.h"

#include <stdio.h>
#include <stdlib.h>

#include <memory>
#include <mutex>

#include "tsafety/TSCounterVar.h"
#include "libgenp/goutp/TdAlnWriter.h"
#include "libgenp/goutp/TdFinalizer.h"
#include "libmymp/mplayout/MpGlobalMemory.h"

////////////////////////////////////////////////////////////////////////////
// CLASS MpBatch
// multi-processing Batch computation of structure alignment
//
class MpBatch
{
public:
    /**
     * @brief 构造 `MpBatch`，初始化其负责的CPU 对齐流水线状态。
     * @param MpGlobalMemory 供当前步骤读取或更新的 `MpGlobalMemory` 缓冲区。
     * @param dareano 供该函数读取或更新的 `dareano` 参数。
     * @param TdAlnWriter 供该函数读取或更新的 `TdAlnWriter` 参数。
     * @return 无返回值；完成对象构造与初始状态设置。
     */
    MpBatch(MpGlobalMemory*, int dareano, TdAlnWriter*);

    /**
     * @brief 销毁 `MpBatch`，释放其管理的CPU 对齐流水线资源。
     * @par 参数
     * 无。
     * @return 无返回值；对象持有的资源在返回前完成释放。
     */
    ~MpBatch();

    /**
     * @brief 在CPU 对齐流水线中处理 `ProcessBlock` 对应的数据。
     * @param qrysernrbeg 描述查询结构的 `qrysernrbeg`。
     * @param queryndxpmbeg 描述查询结构的 `queryndxpmbeg`。
     * @param queryndxpmend 描述查询结构的 `queryndxpmend`。
     * @param querydesc 描述查询结构的 `querydesc`。
     * @param querypmbeg 查询结构打包字段的起始指针数组。
     * @param querypmend 查询结构打包字段的结束指针数组。
     * @param bdbCdesc 描述参考结构的 `bdbCdesc`。
     * @param bdbCpmbeg 参考结构打包字段的起始指针数组。
     * @param bdbCpmend 参考结构打包字段的结束指针数组。
     * @param bdbCndxpmbeg 描述参考结构的 `bdbCndxpmbeg`。
     * @param bdbCndxpmend 描述参考结构的 `bdbCndxpmend`。
     * @param qrscnt 供该函数读取或更新的 `qrscnt` 参数。
     * @param cnt 供该函数读取或更新的 `cnt` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void ProcessBlock(
        int qrysernrbeg,
        char** queryndxpmbeg,
        char** queryndxpmend,
        const char** querydesc,
        char** querypmbeg,
        char** querypmend,
        const char** bdbCdesc,
        char** bdbCpmbeg,
        char** bdbCpmend,
        char** bdbCndxpmbeg,
        char** bdbCndxpmend,
        TSCounterVar* qrscnt,
        TSCounterVar* cnt
    );

    // void WaitForIdleChilds() {
    //     if(cbpfin_)
    //         std::lock_guard<std::mutex> lck(cbpfin_->GetPrivateMutex());
    // }

    /**
     * @brief 在CPU 对齐流水线中读取 `GetCurrentMaxDbPos` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    size_t GetCurrentMaxDbPos() const { return gmem_->GetCurrentMaxDbPos(); }
    /**
     * @brief 在CPU 对齐流水线中读取 `GetCurrentMaxNDbStrs` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    size_t GetCurrentMaxNDbStrs() const { return gmem_->GetCurrentMaxNDbStrs(); }

    /**
     * @brief 在CPU 对齐流水线中读取 `GetCurrentMaxDbPosPass2` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    size_t GetCurrentMaxDbPosPass2() const { return gmem_->GetCurrentMaxDbPosPass2(); }
    /**
     * @brief 在CPU 对齐流水线中读取 `GetCurrentMaxNDbStrsPass2` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    size_t GetCurrentMaxNDbStrsPass2() const { return gmem_->GetCurrentMaxNDbStrsPass2(); }

protected:

    /**
     * @brief 在CPU 对齐流水线中读取 `GetHeapSectionOffset` 对应的数据。
     * @param sec 供该函数读取或更新的 `sec` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    size_t GetHeapSectionOffset(int sec) const {return gmem_->GetHeapSectionOffset(devareano_,sec);}
    /**
     * @brief 在CPU 对齐流水线中读取 `GetHeapSectionAddress` 对应的数据。
     * @param sec 供该函数读取或更新的 `sec` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    char* GetHeapSectionAddress(int sec) const {
        return gmem_->GetHeap() + gmem_->GetHeapSectionOffset(devareano_,sec);
    }

    /**
     * @brief 在CPU 对齐流水线中读取 `GetDeviceName` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    const std::string& GetDeviceName() const {return gmem_->GetDeviceName();}

    //Filter out flagged references
    /**
     * @brief 在CPU 对齐流水线中筛选 `FilteroutReferences` 对应的数据。
     * @param queryndxpmbeg 描述查询结构的 `queryndxpmbeg`。
     * @param queryndxpmend 描述查询结构的 `queryndxpmend`。
     * @param querypmbeg 查询结构打包字段的起始指针数组。
     * @param querypmend 查询结构打包字段的结束指针数组。
     * @param bdbCdesc 描述参考结构的 `bdbCdesc`。
     * @param bdbCpmbeg 参考结构打包字段的起始指针数组。
     * @param bdbCpmend 参考结构打包字段的结束指针数组。
     * @param bdbCndxpmbeg 描述参考结构的 `bdbCndxpmbeg`。
     * @param bdbCndxpmend 描述参考结构的 `bdbCndxpmend`。
     * @param nqyposs 当前批次中查询结构的总位置数。
     * @param nqystrs 当前批次中的查询结构数量。
     * @param ndbCposs 当前批次中参考结构的总位置数。
     * @param ndbCstrs 当前批次中的参考结构数量。
     * @param dbstr1len 批次中最长参考结构的长度。
     * @param dbstrnlen 批次中最短参考结构的长度。
     * @param dbxpad 参考数据行末用于对齐访问的填充长度。
     * @param maxnsteps 每对结构保留的候选搜索步数。
     * @param ndbCposs2 控制当前步骤范围或规模的 `ndbCposs2`。
     * @param ndbCstrs2 控制当前步骤范围或规模的 `ndbCstrs2`。
     * @param dbstr1len2 控制当前步骤范围或规模的 `dbstr1len2`。
     * @param dbstrnlen2 控制当前步骤范围或规模的 `dbstrnlen2`。
     * @param dbxpad2 描述参考结构的 `dbxpad2`。
     * @param tmpdpdiagbuffers 供当前步骤读取或更新的 `tmpdpdiagbuffers` 缓冲区。
     * @param tfmmemory 供当前步骤读取或更新的 `tfmmemory` 缓冲区。
     * @param auxwrkmemory 供当前步骤读取或更新的 `auxwrkmemory` 缓冲区。
     * @param globvarsbuf 供当前步骤读取或更新的 `globvarsbuf` 缓冲区。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void FilteroutReferences(
        char** queryndxpmbeg, char** queryndxpmend,
        char** querypmbeg, char** querypmend,
        const char** bdbCdesc, char** bdbCpmbeg, char** bdbCpmend,
        char** bdbCndxpmbeg, char** bdbCndxpmend,
        const size_t nqyposs, const size_t nqystrs,
        const size_t ndbCposs, const size_t ndbCstrs,
        const size_t dbstr1len, const size_t dbstrnlen,
        const size_t dbxpad, const size_t maxnsteps,
        size_t& ndbCposs2, size_t& ndbCstrs2,
        size_t& dbstr1len2, size_t& dbstrnlen2,
        size_t& dbxpad2,
        float* tmpdpdiagbuffers,
        float* tfmmemory,
        float* auxwrkmemory,
        unsigned int* globvarsbuf);

    //{{ results processing
    /**
     * @brief 在CPU 对齐流水线中处理 `TriggerFinalization` 对应的数据。
     * @param tdrtn 供该函数读取或更新的 `tdrtn` 参数。
     * @param qrysernrbeg 描述查询结构的 `qrysernrbeg`。
     * @param nqyposs 当前批次中查询结构的总位置数。
     * @param nqystrs 当前批次中的查询结构数量。
     * @param ndbCposs 当前批次中参考结构的总位置数。
     * @param ndbCstrs 当前批次中的参考结构数量。
     * @param querydesc 描述查询结构的 `querydesc`。
     * @param querypmbeg 查询结构打包字段的起始指针数组。
     * @param querypmend 查询结构打包字段的结束指针数组。
     * @param bdbCdesc 描述参考结构的 `bdbCdesc`。
     * @param bdbCpmbeg 参考结构打包字段的起始指针数组。
     * @param bdbCpmend 参考结构打包字段的结束指针数组。
     * @param qrscnt 供该函数读取或更新的 `qrscnt` 参数。
     * @param cnt 供该函数读取或更新的 `cnt` 参数。
     * @param szaligns2 供该函数读取或更新的 `szaligns2` 参数。
     * @param passedstats 供该函数读取或更新的 `passedstats` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void TriggerFinalization(
        double tdrtn,
        int qrysernrbeg,
        size_t nqyposs, size_t nqystrs,
        size_t ndbCposs, size_t ndbCstrs,
        const char** querydesc,
        char** querypmbeg,
        char** querypmend,
        const char** bdbCdesc,
        char** bdbCpmbeg,
        char** bdbCpmend,
        TSCounterVar* qrscnt,
        TSCounterVar* cnt,
        size_t szaligns2,
        unsigned int* passedstats
    );
    //}}

private:
    MpGlobalMemory* gmem_;//global memory configuration
    const int devareano_;//memory area number
    std::unique_ptr<unsigned int[]> filterdata_;//passed structures (globals)
    //statistics of structures passed to the next processing stages (globals):
    std::unique_ptr<unsigned int[]> passedstatscntrd_;
    std::unique_ptr<TdFinalizer> cbpfin_;//results finalizer
};

// -------------------------------------------------------------------------
// INLINES ...
//

#endif//__MpBatch_h__
