/***************************************************************************
 *   Copyright (C) 2021-2023 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __TdFinalizer_h__
#define __TdFinalizer_h__

#include "libutil/mybase.h"

#include <stdio.h>

#include <cmath>
#include <string>
#include <vector>
#include <memory>
#include <functional>
#include <mutex>
#include <condition_variable>
#include <thread>

#ifdef GPUINUSE
#   include <cuda.h>
#   include <cuda_runtime_api.h>
#   include "libmycu/cuproc/Devices.h"
#endif

#include "tsafety/TSCounterVar.h"
#include "libgenp/gproc/gproc.h"

#include "libgenp/gdats/PM2DVectorFields.h"
// #include "libgenp/gdats/PMBatchStrData.h"
#include "libgenp/goutp/TdAlnWriter.h"
#include "libmycu/cuproc/cuprocconf.h"

#define CUBPTHREAD_MSG_UNSET -1
#define CUBPTHREAD_MSG_ERROR -2

// _________________________________________________________________________
// Class TdFinalizer
//
// thread class for finalizing results calculated on a device
//
class TdFinalizer
{
public:
    enum TCUBPThreadMsg {
        cubpthreadmsgFinalize,
        cubpthreadmsgTerminate
    };
    enum TCUBPThreadResponse {
        cubptrespmsgFinalizing,
        cubptrespmsgTerminating
    };

public:
#ifdef GPUINUSE
    TdFinalizer(
        cudaStream_t& strcopyres,
        DeviceProperties dprop, 
        TdAlnWriter*
    );
#endif

    /**
     * @brief 构造 `TdFinalizer`，初始化其负责的对齐结果输出状态。
     * @param TdAlnWriter 供该函数读取或更新的 `TdAlnWriter` 参数。
     * @return 无返回值；完成对象构造与初始状态设置。
     */
    TdFinalizer(
        TdAlnWriter*
    );

    /**
     * @brief 销毁 `TdFinalizer`，释放其管理的对齐结果输出资源。
     * @par 参数
     * 无。
     * @return 无返回值；对象持有的资源在返回前完成释放。
     */
    ~TdFinalizer();

    /**
     * @brief 在对齐结果输出中读取 `GetPrivateMutex` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    std::mutex& GetPrivateMutex() {return mx_dataccess_;}

    //{{NOTE: messaging functions accessed from outside!
    /**
     * @brief 在对齐结果输出中通知 `Notify` 对应的数据。
     * @param msg 供该函数读取或更新的 `msg` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void Notify(int msg) {
        {//mutex must be unlocked before notifying
            std::unique_lock<std::mutex> lck(mx_dataccess_);
            if( req_msg_ != CUBPTHREAD_MSG_ERROR)
                req_msg_ = msg;
        }
        cv_msg_.notify_one();
    }
    /**
     * @brief 在对齐结果输出中等待 `waitForDataAccess` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void waitForDataAccess() {
        std::unique_lock<std::mutex> lck(mx_dataccess_);
        cv_msg_.wait(lck,
            [this] {return
                rsp_msg_ == cubptrespmsgTerminating || 
                rsp_msg_ == CUBPTHREAD_MSG_ERROR ||
               (rsp_msg_ == CUBPTHREAD_MSG_UNSET && req_msg_ == CUBPTHREAD_MSG_UNSET);
            }
        );
    }
    /**
     * @brief 在对齐结果输出中读取 `GetResponse` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    int GetResponse() const {
        std::lock_guard<std::mutex> lck(mx_dataccess_);
        return rsp_msg_;
    }
    /**
     * @brief 在对齐结果输出中重置 `ResetResponse` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void ResetResponse() {
        std::lock_guard<std::mutex> lck(mx_dataccess_);
        if( rsp_msg_!= CUBPTHREAD_MSG_ERROR )
            rsp_msg_ = CUBPTHREAD_MSG_UNSET;
    }
    //}}


    /**
     * @brief 在对齐结果输出中设置 `SetCuBPBDbdata` 对应的数据。
     * @param tdrtn 供该函数读取或更新的 `tdrtn` 参数。
     * @param scorethld 当前步骤使用或写回的 `scorethld` 分数。
     * @param qrysernrbeg 描述查询结构的 `qrysernrbeg`。
     * @param devname 供该函数读取或更新的 `devname` 参数。
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
     * @param passedstats 供该函数读取或更新的 `passedstats` 参数。
     * @param h_results 接收当前步骤输出的 `h_results`。
     * @param sz_alndata 供该函数读取或更新的 `sz_alndata` 参数。
     * @param sz_tfmmatrices 表示或保存刚体变换的 `sz_tfmmatrices`。
     * @param sz_alns 供该函数读取或更新的 `sz_alns` 参数。
     * @param dbalnlen2 控制当前步骤范围或规模的 `dbalnlen2`。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void SetCuBPBDbdata(
        double tdrtn,
        float scorethld,
        int qrysernrbeg,
        const std::string& devname,
        size_t nqyposs, size_t nqystrs,
        size_t ndbCposs, size_t ndbCstrs,
        const char** querydesc, char** querypmbeg, char** querypmend, 
        const char** bdbCdesc, char** bdbCpmbeg, char** bdbCpmend,
        TSCounterVar* qrscnt,
        TSCounterVar* cnt,
        unsigned int* passedstats,
        const char* h_results,
        size_t sz_alndata,
        size_t sz_tfmmatrices,
        size_t sz_alns,
        unsigned int dbalnlen2)
    {
        std::unique_lock<std::mutex> lck(mx_dataccess_);
        cv_msg_.wait(lck,
            [this]{
                return (req_msg_ == CUBPTHREAD_MSG_UNSET ||
                        req_msg_ == CUBPTHREAD_MSG_ERROR ||
                        rsp_msg_ == CUBPTHREAD_MSG_ERROR);
            }
        );
        if(req_msg_ == CUBPTHREAD_MSG_ERROR || rsp_msg_ == CUBPTHREAD_MSG_ERROR)
            return;
        cubp_set_duration_ = tdrtn;
        cubp_set_scorethld_ = scorethld;
        cubp_set_qrysernrbeg_ = qrysernrbeg;
        cubp_set_devname_ = devname;
        cubp_set_nqyposs_ = (int)nqyposs;
        cubp_set_nqystrs_ = nqystrs;
        cubp_set_ndbCposs_ = ndbCposs;
        cubp_set_ndbCstrs_ = ndbCstrs;
        cubp_set_querydesc_ = querydesc;
        if(querypmbeg && querypmend) {
            memcpy(cubp_set_querypmbeg_, querypmbeg, pmv2DTotFlds * sizeof(void*));
            memcpy(cubp_set_querypmend_, querypmend, pmv2DTotFlds * sizeof(void*));
        }
        else {
            memset(cubp_set_querypmbeg_, 0, pmv2DTotFlds * sizeof(void*));
            memset(cubp_set_querypmend_, 0, pmv2DTotFlds * sizeof(void*));
        }
        cubp_set_bdbCdesc_ = bdbCdesc;
        if(bdbCpmbeg && bdbCpmend) {
            memcpy(cubp_set_bdbCpmbeg_, bdbCpmbeg, pmv2DTotFlds * sizeof(void*));
            memcpy(cubp_set_bdbCpmend_, bdbCpmend, pmv2DTotFlds * sizeof(void*));
        }
        else {
            memset(cubp_set_bdbCpmbeg_, 0, pmv2DTotFlds * sizeof(void*));
            memset(cubp_set_bdbCpmend_, 0, pmv2DTotFlds * sizeof(void*));
        }
        cubp_set_nposits_.clear();
        cubp_set_nstrs_.clear();
        for(size_t i = 0; i < nqystrs; i++) {
            size_t loff = nDevGlobVariables * i;
            cubp_set_nposits_.push_back(passedstats[loff + dgvNPosits]);
            cubp_set_nstrs_.push_back((int)passedstats[loff + dgvNPassedStrs]);
        }
        cubp_set_qrscnt_ = qrscnt;
        cubp_set_cnt_ = cnt;
        cubp_set_h_results_ = h_results;
        cubp_set_sz_alndata_ = sz_alndata;
        cubp_set_sz_tfmmatrices_ = sz_tfmmatrices;
        cubp_set_sz_alns_ = sz_alns;
        cubp_set_sz_dbalnlen2_ = dbalnlen2;
    }


protected:
    /**
     * @brief 在对齐结果输出中处理 `Execute` 对应的数据。
     * @param args 供该函数读取或更新的 `args` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void Execute(void* args);

    /**
     * @brief 在对齐结果输出中设置 `SetResponseError` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void SetResponseError() {
        std::lock_guard<std::mutex> lck(mx_dataccess_);
        rsp_msg_ = CUBPTHREAD_MSG_ERROR;
    }

    /**
     * @brief 在对齐结果输出中处理 `FinalizeQueries` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void FinalizeQueries();
    /**
     * @brief 在对齐结果输出中排序 `SortCompressedResults` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void SortCompressedResults();
    /**
     * @brief 在对齐结果输出中处理 `PassResultsToWriter` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void PassResultsToWriter();
    /**
     * @brief 在对齐结果输出中格式化输出 `PrintCompressedResults` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void PrintCompressedResults() const;


    //{{formating methods for PLAIN format
    /**
     * @brief 在对齐结果输出中处理 `CompressResultsPlain` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void CompressResultsPlain();
    /**
     * @brief 在对齐结果输出中读取 `GetSizeOfCompressedResultsPlain` 对应的数据。
     * @param szannot 供该函数读取或更新的 `szannot` 参数。
     * @param szalns 供该函数读取或更新的 `szalns` 参数。
     * @param szalnswodesc 供该函数读取或更新的 `szalnswodesc` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void GetSizeOfCompressedResultsPlain(
        size_t* szannot, size_t* szalns, size_t* szalnswodesc) const;

    /**
     * @brief 在对齐结果输出中构造 `MakeAnnotationPlain` 对应的数据。
     * @param outptr 接收当前步骤输出的 `outptr`。
     * @param strndx 供该函数读取或更新的 `strndx` 参数。
     * @param orgstrndx 供该函数读取或更新的 `orgstrndx` 参数。
     * @param desc 供该函数读取或更新的 `desc` 参数。
     * @param alnlen 控制当前步骤范围或规模的 `alnlen`。
     * @param dbstrlen 控制当前步骤范围或规模的 `dbstrlen`。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void MakeAnnotationPlain(char*& outptr,
        const int strndx, const unsigned int orgstrndx, const char* desc,
        const unsigned int alnlen, const int dbstrlen) const;

    /**
     * @brief 在对齐结果输出中处理 `FormatScoresPlain` 对应的数据。
     * @param outptr 接收当前步骤输出的 `outptr`。
     * @param strndx 供该函数读取或更新的 `strndx` 参数。
     * @param orgstrndx 供该函数读取或更新的 `orgstrndx` 参数。
     * @param alnlen 控制当前步骤范围或规模的 `alnlen`。
     * @param score 当前步骤使用或写回的 `score` 分数。
     * @param qrystrlen 控制当前步骤范围或规模的 `qrystrlen`。
     * @param dbstrlen 控制当前步骤范围或规模的 `dbstrlen`。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void FormatScoresPlain(char*& outptr,
        int strndx, unsigned int orgstrndx, unsigned int alnlen,
        float score, int qrystrlen, int dbstrlen);

    /**
     * @brief 在对齐结果输出中处理 `FormatAlignmentPlain` 对应的数据。
     * @param outptr 接收当前步骤输出的 `outptr`。
     * @param strndx 供该函数读取或更新的 `strndx` 参数。
     * @param orgstrndx 供该函数读取或更新的 `orgstrndx` 参数。
     * @param dbstr2dst 描述参考结构的 `dbstr2dst`。
     * @param alnlen 控制当前步骤范围或规模的 `alnlen`。
     * @param qrystrlen 控制当前步骤范围或规模的 `qrystrlen`。
     * @param dbstrlen 控制当前步骤范围或规模的 `dbstrlen`。
     * @param width 供该函数读取或更新的 `width` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void FormatAlignmentPlain(char*& outptr,
        int strndx, unsigned int orgstrndx, unsigned int dbstr2dst,
        int alnlen, int qrystrlen, int dbstrlen, const int width);

    /**
     * @brief 在对齐结果输出中处理 `FormatFooterPlain` 对应的数据。
     * @param outptr 接收当前步骤输出的 `outptr`。
     * @param strndx 供该函数读取或更新的 `strndx` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void FormatFooterPlain(char*& outptr, int strndx);
    //}}
    //{{formating methods for JSON format
    /**
     * @brief 在对齐结果输出中处理 `CompressResultsJSON` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void CompressResultsJSON();
    /**
     * @brief 在对齐结果输出中读取 `GetSizeOfCompressedResultsJSON` 对应的数据。
     * @param szannot 供该函数读取或更新的 `szannot` 参数。
     * @param szalns 供该函数读取或更新的 `szalns` 参数。
     * @param szalnswodesc 供该函数读取或更新的 `szalnswodesc` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void GetSizeOfCompressedResultsJSON(
        size_t* szannot, size_t* szalns, size_t* szalnswodesc) const;

    /**
     * @brief 在对齐结果输出中构造 `MakeAnnotationJSON` 对应的数据。
     * @param outptr 接收当前步骤输出的 `outptr`。
     * @param strndx 供该函数读取或更新的 `strndx` 参数。
     * @param desc 供该函数读取或更新的 `desc` 参数。
     * @param alnlen 控制当前步骤范围或规模的 `alnlen`。
     * @param dbstrlen 控制当前步骤范围或规模的 `dbstrlen`。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void MakeAnnotationJSON(char*& outptr,
        const int strndx, const char* desc,
        const unsigned int alnlen, const int dbstrlen) const;

    /**
     * @brief 在对齐结果输出中处理 `FormatScoresJSON` 对应的数据。
     * @param outptr 接收当前步骤输出的 `outptr`。
     * @param strndx 供该函数读取或更新的 `strndx` 参数。
     * @param alnlen 控制当前步骤范围或规模的 `alnlen`。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void FormatScoresJSON(char*& outptr, int strndx, unsigned int alnlen);

    /**
     * @brief 在对齐结果输出中处理 `FormatAlignmentJSON` 对应的数据。
     * @param outptr 接收当前步骤输出的 `outptr`。
     * @param strndx 供该函数读取或更新的 `strndx` 参数。
     * @param orgstrndx 供该函数读取或更新的 `orgstrndx` 参数。
     * @param dbstr2dst 描述参考结构的 `dbstr2dst`。
     * @param alnlen 控制当前步骤范围或规模的 `alnlen`。
     * @param qrystrlen 控制当前步骤范围或规模的 `qrystrlen`。
     * @param dbstrlen 控制当前步骤范围或规模的 `dbstrlen`。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void FormatAlignmentJSON(char*& outptr,
        int strndx, unsigned int orgstrndx, unsigned int dbstr2dst,
        int alnlen, int qrystrlen, int dbstrlen);

    /**
     * @brief 在对齐结果输出中处理 `FormatFooterJSON` 对应的数据。
     * @param outptr 接收当前步骤输出的 `outptr`。
     * @param strndx 供该函数读取或更新的 `strndx` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void FormatFooterJSON(char*& outptr, int strndx);
    //}}

    
    /**
     * @brief 在对齐结果输出中读取 `GetBegOfAlns` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    const char* GetBegOfAlns() {
        return cubp_set_h_results_ + cubp_set_sz_alndata_ + cubp_set_sz_tfmmatrices_;
    }

    /**
     * @brief 在对齐结果输出中读取 `GetAlnSectionAt` 对应的数据。
     * @param ptr 供该函数读取或更新的 `ptr` 参数。
     * @param sctndx 供该函数读取或更新的 `sctndx` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    const char* GetAlnSectionAt(const char* ptr, const int sctndx) const {
        return ptr + sctndx * dbalnlen_;
    }

    template<typename T>
    /**
     * @brief 在对齐结果输出中读取 `GetOutputAlnDataField` 对应的数据。
     * @param strndx 供该函数读取或更新的 `strndx` 参数。
     * @param field 供该函数读取或更新的 `field` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    T GetOutputAlnDataField(int strndx, int field) const
        {
            return 
                *(T*)((float*)cubp_set_h_results_ + 
                    nTDP2OutputAlnData * (cumnstrs_ + strndx) + field);
        }

    /**
     * @brief 在对齐结果输出中读取 `GetOutputTfmMtxField` 对应的数据。
     * @param strndx 供该函数读取或更新的 `strndx` 参数。
     * @param field 供该函数读取或更新的 `field` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    float GetOutputTfmMtxField(int strndx, int field) const
        {
            return 
                *((float*)(cubp_set_h_results_ + cubp_set_sz_alndata_) + 
                    nTTranformMatrix * (cumnstrs_ + strndx) + field);
        }

    /**
     * @brief 在对齐结果输出中读取 `GetOutputTfmMtxAddress` 对应的数据。
     * @param strndx 供该函数读取或更新的 `strndx` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    float* GetOutputTfmMtxAddress(int strndx) const
        {
            return 
                (float*)(cubp_set_h_results_ + cubp_set_sz_alndata_) + 
                    nTTranformMatrix * (cumnstrs_ + strndx);
        }

    template<typename T>
    /**
     * @brief 在对齐结果输出中读取 `GetStructureField` 对应的数据。
     * @param pmbeg 供该函数读取或更新的 `pmbeg` 参数。
     * @param ndx 控制当前步骤范围或规模的 `ndx`。
     * @param field 供该函数读取或更新的 `field` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    T GetStructureField(
        char* const pmbeg[pmv2DTotFlds],
        unsigned int ndx, unsigned int field) const
        {
            return ((T*)(pmbeg[field]))[ndx];
        }

    template<typename T>
    /**
     * @brief 在对齐结果输出中读取 `GetQueryField` 对应的数据。
     * @param orgqryndx 描述查询结构的 `orgqryndx`。
     * @param field 供该函数读取或更新的 `field` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    T GetQueryField(unsigned int orgqryndx, unsigned int field) const;
    template<typename T>
    /**
     * @brief 在对齐结果输出中读取 `GetQueryFieldPos` 对应的数据。
     * @param field 供该函数读取或更新的 `field` 参数。
     * @param pos 供该函数读取或更新的 `pos` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    T GetQueryFieldPos(unsigned int field, unsigned int pos) const;
    template<typename T>
    /**
     * @brief 在对齐结果输出中读取 `GetDbStructureField` 对应的数据。
     * @param orgstrndx 供该函数读取或更新的 `orgstrndx` 参数。
     * @param field 供该函数读取或更新的 `field` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    T GetDbStructureField(unsigned int orgstrndx, unsigned int field) const;
    template<typename T>
    /**
     * @brief 在对齐结果输出中读取 `GetDbStructureFieldPos` 对应的数据。
     * @param orgstrndx 供该函数读取或更新的 `orgstrndx` 参数。
     * @param field 供该函数读取或更新的 `field` 参数。
     * @param pos 供该函数读取或更新的 `pos` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    T GetDbStructureFieldPos(unsigned int orgstrndx, unsigned int field, unsigned int pos) const;

    /**
     * @brief 在对齐结果输出中读取 `GetDbStructureDesc` 对应的数据。
     * @param desc 供该函数读取或更新的 `desc` 参数。
     * @param orgstrndx 供该函数读取或更新的 `orgstrndx` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void GetDbStructureDesc(const char*& desc, unsigned int orgstrndx) const;


    /**
     * @brief 在对齐结果输出中处理 `ReserveVectors` 对应的数据。
     * @param capacity 供该函数读取或更新的 `capacity` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void ReserveVectors(int capacity) {
        srtindxs_.reset(new std::vector<int>);
        scores_.reset(new std::vector<float>);
        alnptrs_.reset(new std::vector<char*>);
        annotptrs_.reset(new std::vector<char*>);
        if(capacity > 0) {
            if(srtindxs_) srtindxs_->reserve(capacity);
            if(scores_) scores_->reserve(capacity);
            if(alnptrs_) alnptrs_->reserve(capacity);
            if(annotptrs_) annotptrs_->reserve(capacity);
        }
    }

private:
    //thread section
    std::thread* tobj_;//thread object
private:
    //{{messaging
    std::condition_variable cv_msg_;//condition variable for messaging
    mutable std::mutex mx_dataccess_;//mutex for accessing class data
    int req_msg_;//request message issued for thread
    int rsp_msg_;//private response message
    //}}
    //
#ifdef GPUINUSE
    //properties of device the thread is associated with:
    cudaStream_t& strcopyres_;
    DeviceProperties dprop_;
#endif
    //results writer:
    TdAlnWriter* alnwriter_;
    //{{data arguments: 
    // cubp-set data/addresses:
    int qrysernr_;//serial number of the query under processing
    std::string qrydesc_;//description of query qrysernr_
    int qrystrlen_;//length of query qrysernr_
    int qrynstrs_;//number of target db structures for query qrysernr_
    unsigned int qrynposits_;//total number of target db structure positions for query qrysernr_
    int cumnstrs_;//cumulative number of structures over all queries up to qrysernr_
    int offsetalns_;//offset to the start of the alignments for query qrysernr_
    int dbalnlen_;//alignment length for query qrysernr_ across all references
    double cubp_set_duration_;//time duration
    float cubp_set_scorethld_;//tm-score threshold
    int cubp_set_qrysernrbeg_;//serial number of the first query in the chunk
    int cubp_set_nqyposs_;//total length of queries in the chunk
    size_t cubp_set_nqystrs_;
    size_t cubp_set_ndbCposs_;
    size_t cubp_set_ndbCstrs_;
    std::string cubp_set_devname_;//device name
    const char** cubp_set_querydesc_;
    char* cubp_set_querypmbeg_[pmv2DTotFlds];
    char* cubp_set_querypmend_[pmv2DTotFlds];
    const char** cubp_set_bdbCdesc_;
    char* cubp_set_bdbCpmbeg_[pmv2DTotFlds];
    char* cubp_set_bdbCpmend_[pmv2DTotFlds];
    TSCounterVar* cubp_set_qrscnt_;//counter for queries
    TSCounterVar* cubp_set_cnt_;//counter for references
    std::vector<unsigned int> cubp_set_nposits_;//total number of positions in the transfered results for each query
    std::vector<int> cubp_set_nstrs_;//number of structures in the transfered results for each query
    const char* cubp_set_h_results_;//cubp-set host-side results
    size_t cubp_set_sz_alndata_;//size of alignment data of results
    size_t cubp_set_sz_tfmmatrices_;//size of transformation matrices of results
    size_t cubp_set_sz_alns_;//size of alignments of results
    int cubp_set_sz_dbalnlen2_;
    //}}
    //{{formatted results for one particular query:
    std::unique_ptr<char,WritersDataDestroyer> annotations_;
    std::unique_ptr<char,WritersDataDestroyer> alignments_;
    std::unique_ptr<std::vector<int>> srtindxs_;//index vector of sorted scores
    std::unique_ptr<std::vector<float>> scores_;//vector of scores
    std::unique_ptr<std::vector<char*>> alnptrs_;//vector of alignments
    std::unique_ptr<std::vector<char*>> annotptrs_;//vector of annotations
    //}}
};

// -------------------------------------------------------------------------
// GetQueryField: get a field of the given query structure;
// NOTE: the pointer of the query structure under process is accessed!
template<typename T>
/**
 * @brief 在对齐结果输出中读取 `TdFinalizer::GetQueryField` 对应的数据。
 * @param orgqryndx 描述查询结构的 `orgqryndx`。
 * @param field 供该函数读取或更新的 `field` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline
T TdFinalizer::GetQueryField(
    unsigned int orgqryndx, unsigned int field) const
{
    return GetStructureField<T>(cubp_set_querypmbeg_, orgqryndx, field);
}
// GetQueryFieldPos: get a field of the query structure at the given 
// position;
// NOTE: the pointer of the query structure under process is accessed!
template<typename T>
/**
 * @brief 在对齐结果输出中读取 `TdFinalizer::GetQueryFieldPos` 对应的数据。
 * @param field 供该函数读取或更新的 `field` 参数。
 * @param pos 供该函数读取或更新的 `pos` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline
T TdFinalizer::GetQueryFieldPos(
    unsigned int field, unsigned int pos) const
{
    return GetStructureField<T>(cubp_set_querypmbeg_, pos, field);
}
// GetDbStructureField: get a field of the given Db structure;
// orgstrndx, index of the structure over all pm data structures;
template<typename T>
/**
 * @brief 在对齐结果输出中读取 `TdFinalizer::GetDbStructureField` 对应的数据。
 * @param orgstrndx 供该函数读取或更新的 `orgstrndx` 参数。
 * @param field 供该函数读取或更新的 `field` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline
T TdFinalizer::GetDbStructureField(
    unsigned int orgstrndx, unsigned int field) const
{
    return GetStructureField<T>(cubp_set_bdbCpmbeg_, orgstrndx, field);
}
// GetDbStructureFieldPos: get a field of the given Db structure at the 
// given position;
// orgstrndx, index of the structure over all pm data structures;
template<typename T>
/**
 * @brief 在对齐结果输出中读取 `TdFinalizer::GetDbStructureFieldPos` 对应的数据。
 * @param orgstrndx 供该函数读取或更新的 `orgstrndx` 参数。
 * @param field 供该函数读取或更新的 `field` 参数。
 * @param pos 供该函数读取或更新的 `pos` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline
T TdFinalizer::GetDbStructureFieldPos(
    unsigned int /*orgstrndx*/, unsigned int field, unsigned int pos) const
{
    return GetStructureField<T>(cubp_set_bdbCpmbeg_, pos, field);
}

// GetDbStructureDesc: get the db structure description;
// orgstrndx, index of the structure over all pm data structures;
/**
 * @brief 在对齐结果输出中读取 `TdFinalizer::GetDbStructureDesc` 对应的数据。
 * @param desc 供该函数读取或更新的 `desc` 参数。
 * @param orgstrndx 供该函数读取或更新的 `orgstrndx` 参数。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
inline
void TdFinalizer::GetDbStructureDesc(
    const char*& desc, unsigned int orgstrndx) const
{
    desc = cubp_set_bdbCdesc_[orgstrndx];
}

#endif//__TdFinalizer_h__
