/***************************************************************************
 *   Copyright (C) 2021-2023 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __TdCommutator_h__
#define __TdCommutator_h__

#include "libutil/mybase.h"

#include <stdio.h>

#include <memory>
#include <mutex>
#include <condition_variable>
#include <thread>

#include "tsafety/TSCounterVar.h"
#include "libgenp/gdats/PM2DVectorFields.h"
// #include "libgenp/gdats/PMBatchStrData.h"
#include "libgenp/goutp/TdClustWriter.h"
#include "libgenp/goutp/TdAlnWriter.h"
#include "libmycu/cuproc/Devices.h"
#include "libmycu/culayout/CuDeviceMemory.cuh"
#include "libmycu/cubatch/CuBatch.cuh"

#define THREAD_MSG_UNSET -1
#define THREAD_MSG_ERROR -2
#define THREAD_MSG_ADDRESSEE_NONE -1

// _________________________________________________________________________
// Class TdCommutator
//
// thread class for ensuring data flow to and from a device (GPU)
//
class TdCommutator
{
public:
    enum TThreadMsg {
        tthreadmsgGetDataChunkSize,
        tthreadmsgProcessNewData,
        tthreadmsgProbe,
        tthreadmsgTerminate
    };
    enum TThreadResponseMsg {
        ttrespmsgChunkSizeReady,
        ttrespmsgInProgress,
        ttrespmsgProbed,
        ttrespmsgTerminating
    };

public:
    /**
     * @brief 构造 `TdCommutator`，初始化其负责的CUDA 对齐流水线状态。
     * @param ringsize 控制当前步骤范围或规模的 `ringsize`。
     * @param tid 供该函数读取或更新的 `tid` 参数。
     * @param dmem 供当前步骤读取或更新的 `dmem` 缓冲区。
     * @param areano 供该函数读取或更新的 `areano` 参数。
     * @param TdAlnWriter 供该函数读取或更新的 `TdAlnWriter` 参数。
     * @param TdClustWriter 供该函数读取或更新的 `TdClustWriter` 参数。
     * @param chunkdatasize 控制当前步骤范围或规模的 `chunkdatasize`。
     * @param chunkdatalen 控制当前步骤范围或规模的 `chunkdatalen`。
     * @param chunknstrs 供该函数读取或更新的 `chunknstrs` 参数。
     * @return 无返回值；完成对象构造与初始状态设置。
     */
    TdCommutator(
        size_t ringsize,
        int tid, CuDeviceMemory* dmem, int areano,
        TdAlnWriter*, TdClustWriter*,
        size_t chunkdatasize, size_t chunkdatalen, size_t chunknstrs
    );

    /**
     * @brief 销毁 `TdCommutator`，释放其管理的CUDA 对齐流水线资源。
     * @par 参数
     * 无。
     * @return 无返回值；对象持有的资源在返回前完成释放。
     */
    ~TdCommutator();

//     std::mutex& GetPrivateMutex() {return mx_rsp_msg_;}
//     std::condition_variable& GetPrivateCV() {return cv_rsp_msg_;}

    //{{NOTE: messaging functions accessed from outside!
    /**
     * @brief 在CUDA 对齐流水线中通知 `Notify` 对应的数据。
     * @param msg 供该函数读取或更新的 `msg` 参数。
     * @param adr 供该函数读取或更新的 `adr` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void Notify(int msg, int adr) {
        {//mutex must be unlocked before notifying
            std::lock_guard<std::mutex> lck(mx_rsp_msg_);
            req_msg_ = msg;
            msg_addressee_ = adr;
        }
        cv_rsp_msg_.notify_one();
    }
    /**
     * @brief 在CUDA 对齐流水线中等待 `waitForDataAccess` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    int waitForDataAccess() {
        int rsp = GetResponseAsync();
        if(rsp == THREAD_MSG_ERROR)
            return rsp;
        {   std::unique_lock<std::mutex> lck(mx_dataccess_);
            cv_dataccess_.wait(lck, [this]{return GetMstrDataEmpty();});
            //mutex release
        }
        return GetResponseAsync();
    }
    /**
     * @brief 在CUDA 对齐流水线中等待 `Wait` 对应的数据。
     * @param rsp 供该函数读取或更新的 `rsp` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    int Wait(int rsp) {
        //wait for the response
        std::unique_lock<std::mutex> lck_msg(mx_rsp_msg_);
        cv_rsp_msg_.wait(lck_msg,
            [this,rsp]{
                return (rsp_msg_ == rsp ||
                        rsp_msg_ == THREAD_MSG_ERROR);
            }
        );
        //lock is back; unset the response
        int rspmsg = rsp_msg_;
        if(rsp_msg_ != THREAD_MSG_ERROR)
            rsp_msg_ = THREAD_MSG_UNSET;
        return rspmsg;
    }
    /**
     * @brief 在CUDA 对齐流水线中处理 `IsIdle` 对应的数据。
     * @param mstr_set_data_empty 供该函数读取或更新的 `mstr_set_data_empty` 参数。
     * @param wait 供该函数读取或更新的 `wait` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    bool IsIdle(bool* mstr_set_data_empty, bool wait = false) {
        std::unique_lock<std::mutex> lck_busy(mx_rsp_msg_, std::defer_lock);
        if(mstr_set_data_empty)
        {   std::lock_guard<std::mutex> lck(mx_dataccess_);
            *mstr_set_data_empty = GetMstrDataEmpty();
        }
        //NOTE: when recycling the loop, lock may be acquired, 
        // although the thread is not idle;
        // data emptiness should be checked jointly!
        bool lck_acquired = true;
        if(wait)
            lck_busy.lock();
        else
            lck_acquired = lck_busy.try_lock();
        return lck_acquired;
        //release of the mutex if locked
    }
    /**
     * @brief 在CUDA 对齐流水线中读取 `GetResponseAsync` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    int GetResponseAsync() const {
        //get a response if available
        std::unique_lock<std::mutex> lck_busy(mx_rsp_msg_, std::defer_lock);
        int rsp = THREAD_MSG_UNSET;
        if(lck_busy.try_lock())
            rsp = rsp_msg_;
        return rsp;
    }
    /**
     * @brief 在CUDA 对齐流水线中读取 `GetResponse` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    int GetResponse() const {
        std::lock_guard<std::mutex> lck(mx_rsp_msg_);
        return rsp_msg_;
    }
//     void ResetResponse() {
//         std::lock_guard<std::mutex> lck(mx_rsp_msg_);
//         if( rsp_msg_!= THREAD_MSG_ERROR )
//             rsp_msg_ = THREAD_MSG_UNSET;
//     }
    //}}


//     //{{NOTE: these functions should be called, and member variables accessed, 
//     // only under locked mutex mx_rsp_msg_!
//     void SetBcastMessageAndAddressee( int msg, int adr ) {
//         req_msg_ = msg; msg_addressee_ = adr;
//     }
//     int GetMessage() const {return req_msg_;}
//     int GetResponseMsg() const {return rsp_msg_;}
//     void ResetResponseMsg() {
//         if( rsp_msg_!= THREAD_MSG_ERROR ) 
//             rsp_msg_ = THREAD_MSG_UNSET;
//     }
//     //}}


    /**
     * @brief 在CUDA 对齐流水线中设置 `SetQueryLen` 对应的数据。
     * @param nqrsposs 控制当前步骤范围或规模的 `nqrsposs`。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void SetQueryLen(size_t nqrsposs)
    {
        //safely read locked data
        std::lock_guard<std::mutex> lck(mx_dataccess_);
        mstr_set_nqrsposs_ = nqrsposs;
    }


    /**
     * @brief 在CUDA 对齐流水线中设置 `SetMstrQueryBDbdata` 对应的数据。
     * @param chunkno 供该函数读取或更新的 `chunkno` 参数。
     * @param lastchunk 供该函数读取或更新的 `lastchunk` 参数。
     * @param newsetqrs 控制当前步骤范围或规模的 `newsetqrs`。
     * @param qrysernrbeg 描述查询结构的 `qrysernrbeg`。
     * @param scorethld 当前步骤使用或写回的 `scorethld` 分数。
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
     * @param qrstscnt 供该函数读取或更新的 `qrstscnt` 参数。
     * @param tscnt 供该函数读取或更新的 `tscnt` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void SetMstrQueryBDbdata(
        int chunkno,
        bool lastchunk,
        bool newsetqrs,
        int qrysernrbeg,
        float scorethld,
        char** queryndxpmbeg, char** queryndxpmend,
        const char** querydesc, char** querypmbeg, char** querypmend,
        const char** bdbCdesc, char** bdbCpmbeg, char** bdbCpmend,
        char** bdbCndxpmbeg, char** bdbCndxpmend,
        TSCounterVar* qrstscnt,
        TSCounterVar* tscnt)
    {
        std::lock_guard<std::mutex> lck(mx_dataccess_);
        mstr_set_chunkno_ = chunkno;
        mstr_set_lastchunk_ = lastchunk;
        mstr_set_newsetqrs_ = newsetqrs;
        mstr_set_qrysernrbeg_ = qrysernrbeg;
        mstr_set_scorethld_ = scorethld;
        //
        mstr_set_querydesc_ = querydesc;
        mstr_set_bdbCdesc_ = bdbCdesc;
        mstr_set_qrscnt_ = qrstscnt;
        mstr_set_cnt_ = tscnt;
        if(newsetqrs) {
            if(queryndxpmbeg && queryndxpmend) {
                memcpy(mstr_set_queryndxpmbeg_, queryndxpmbeg, pmv2DTotIndexFlds * sizeof(void*));
                memcpy(mstr_set_queryndxpmend_, queryndxpmend, pmv2DTotIndexFlds * sizeof(void*));
            } else {
                memset(mstr_set_queryndxpmbeg_, 0, pmv2DTotIndexFlds * sizeof(void*));
                memset(mstr_set_queryndxpmend_, 0, pmv2DTotIndexFlds * sizeof(void*));
            }
            if(querypmbeg && querypmend) {
                memcpy(mstr_set_querypmbeg_, querypmbeg, pmv2DTotFlds * sizeof(void*));
                memcpy(mstr_set_querypmend_, querypmend, pmv2DTotFlds * sizeof(void*));
            } else {
                memset(mstr_set_querypmbeg_, 0, pmv2DTotFlds * sizeof(void*));
                memset(mstr_set_querypmend_, 0, pmv2DTotFlds * sizeof(void*));
            }
        }
        if(bdbCndxpmbeg && bdbCndxpmend) {
            memcpy(mstr_set_bdbCndxpmbeg_, bdbCndxpmbeg, pmv2DTotIndexFlds * sizeof(void*));
            memcpy(mstr_set_bdbCndxpmend_, bdbCndxpmend, pmv2DTotIndexFlds * sizeof(void*));
        }
        if(bdbCpmbeg && bdbCpmend) {
            memcpy(mstr_set_bdbCpmbeg_, bdbCpmbeg, pmv2DTotFlds * sizeof(void*));
            memcpy(mstr_set_bdbCpmend_, bdbCpmend, pmv2DTotFlds * sizeof(void*));
        }
    }


    /**
     * @brief 在CUDA 对齐流水线中读取 `GetChunkDataAttributes` 对应的数据。
     * @param chunkdatasize 控制当前步骤范围或规模的 `chunkdatasize`。
     * @param chunkdatalen 控制当前步骤范围或规模的 `chunkdatalen`。
     * @param chunknstrs 供该函数读取或更新的 `chunknstrs` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void GetChunkDataAttributes(
        size_t* chunkdatasize, size_t* chunkdatalen, size_t* chunknstrs)
    {
        std::lock_guard<std::mutex> lck(mx_dataccess_);
        *chunkdatasize = chunkdatasize_;
        *chunkdatalen = chunkdatalen_;
        *chunknstrs = chunknstrs_;
    }


protected:
    /**
     * @brief 在CUDA 对齐流水线中处理 `Execute` 对应的数据。
     * @param args 供该函数读取或更新的 `args` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void Execute(void* args);


    // GetArgsOnMsgGetDataChunkSize: copy data set by the master on
    // acceptance of the message tthreadmsgGetDataChunkSize
    /**
     * @brief 在CUDA 对齐流水线中读取 `GetArgsOnMsgGetDataChunkSize` 对应的数据。
     * @param cbpc 供该函数读取或更新的 `cbpc` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void GetArgsOnMsgGetDataChunkSize(CuBatch& /*cbpc*/)
    {
        //safely read locked data
        std::lock_guard<std::mutex> lck(mx_dataccess_);
        nqrsposs_ = mstr_set_nqrsposs_;
//         cbpc.SetQueryLen(nqrsposs_);
        mstr_set_nqrsposs_ = 0;
    }

    /**
     * @brief 在CUDA 对齐流水线中计算 `CalculateMaxDbDataChunkSize` 对应的数据。
     * @param CuBatch 供该函数读取或更新的 `CuBatch` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void CalculateMaxDbDataChunkSize(CuBatch&);

    /**
     * @brief 在CUDA 对齐流水线中设置 `SetChunkDataAttributes` 对应的数据。
     * @param chunkdatasize 控制当前步骤范围或规模的 `chunkdatasize`。
     * @param chunkdatalen 控制当前步骤范围或规模的 `chunkdatalen`。
     * @param chunknstrs 供该函数读取或更新的 `chunknstrs` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void SetChunkDataAttributes( 
        size_t chunkdatasize, size_t chunkdatalen, size_t chunknstrs)
    {
        std::lock_guard<std::mutex> lck(mx_dataccess_);
        chunkdatasize_ = chunkdatasize;
        chunkdatalen_ = chunkdatalen;
        chunknstrs_ = chunknstrs;
    }


    // CopyDataOnMsgProcessNewData: copy data set by the master on
    // acceptance of the message tthreadmsgProcessNewData
    //
    /**
     * @brief 在CUDA 对齐流水线中复制 `CopyDataOnMsgProcessNewData` 对应的数据。
     * @param CuBatch 供该函数读取或更新的 `CuBatch` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void CopyDataOnMsgProcessNewData(CuBatch&)
    {
        //safely read and write addresses written by the master
        {   std::lock_guard<std::mutex> lck(mx_dataccess_);
            chunkno_ = mstr_set_chunkno_;
            lastchunk_ = mstr_set_lastchunk_;
            newsetqrs_ = mstr_set_newsetqrs_;
            qrysernrbeg_ = mstr_set_qrysernrbeg_;
            nqrsposs_ = mstr_set_nqrsposs_;/**/
            scorethld_ = mstr_set_scorethld_;/**/
            //
            querydesc_ = mstr_set_querydesc_;
            bdbCdesc_ = mstr_set_bdbCdesc_;
            qrscnt_ = mstr_set_qrscnt_;
            cnt_ = mstr_set_cnt_;
            if(newsetqrs_) {
                memcpy(queryndxpmbeg_, mstr_set_queryndxpmbeg_, pmv2DTotIndexFlds * sizeof(void*));
                memcpy(queryndxpmend_, mstr_set_queryndxpmend_, pmv2DTotIndexFlds * sizeof(void*));
                memcpy(querypmbeg_, mstr_set_querypmbeg_, pmv2DTotFlds * sizeof(void*));
                memcpy(querypmend_, mstr_set_querypmend_, pmv2DTotFlds * sizeof(void*));/**/
            }
            memcpy(bdbCndxpmbeg_, mstr_set_bdbCndxpmbeg_, pmv2DTotIndexFlds * sizeof(void*));
            memcpy(bdbCndxpmend_, mstr_set_bdbCndxpmend_, pmv2DTotIndexFlds * sizeof(void*));
            memcpy(bdbCpmbeg_, mstr_set_bdbCpmbeg_, pmv2DTotFlds * sizeof(void*));
            memcpy(bdbCpmend_, mstr_set_bdbCpmend_, pmv2DTotFlds * sizeof(void*));
        }

        ResetMasterData();
    }

    /**
     * @brief 在CUDA 对齐流水线中重置 `ResetMasterData` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void ResetMasterData()
    {
        {   std::lock_guard<std::mutex> lck(mx_dataccess_);
            //reset addresses so that master can write the addresses of next data
            mstr_set_chunkno_ = -1;
            mstr_set_lastchunk_ = false;
            mstr_set_newsetqrs_ = false;
            mstr_set_qrysernrbeg_ = -1;
            mstr_set_nqrsposs_ = 0;
            mstr_set_scorethld_ = 0.0f;
            mstr_set_querydesc_ = NULL;
            mstr_set_bdbCdesc_ = NULL;
            mstr_set_qrscnt_ = NULL;
            mstr_set_cnt_ = NULL;
            memset(mstr_set_queryndxpmbeg_, 0, pmv2DTotIndexFlds * sizeof(void*));
            memset(mstr_set_queryndxpmend_, 0, pmv2DTotIndexFlds * sizeof(void*));
            memset(mstr_set_querypmbeg_, 0, pmv2DTotFlds * sizeof(void*));
            memset(mstr_set_querypmend_, 0, pmv2DTotFlds * sizeof(void*));
            memset(mstr_set_bdbCndxpmbeg_, 0, pmv2DTotIndexFlds * sizeof(void*));
            memset(mstr_set_bdbCndxpmend_, 0, pmv2DTotIndexFlds * sizeof(void*));
            memset(mstr_set_bdbCpmbeg_, 0, pmv2DTotFlds * sizeof(void*));
            memset(mstr_set_bdbCpmend_, 0, pmv2DTotFlds * sizeof(void*));
        }
        cv_dataccess_.notify_one();
    }

    /**
     * @brief 在CUDA 对齐流水线中读取 `GetMstrDataEmpty` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    bool GetMstrDataEmpty() {
        return 
            mstr_set_querydesc_ == NULL && 
                mstr_set_queryndxpmbeg_[0] == NULL && mstr_set_queryndxpmend_[0] == NULL &&
                mstr_set_querypmbeg_[0] == NULL && mstr_set_querypmend_[0] == NULL &&
            mstr_set_bdbCdesc_ == NULL && 
                mstr_set_bdbCndxpmbeg_[0] == NULL && mstr_set_bdbCndxpmend_[0] == NULL &&
                mstr_set_bdbCpmbeg_[0] == NULL && mstr_set_bdbCpmend_[0] == NULL;
    }

    /**
     * @brief 在CUDA 对齐流水线中处理 `ProcessBlock` 对应的数据。
     * @param CuBatch 供该函数读取或更新的 `CuBatch` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void ProcessBlock(CuBatch&);

private:
    //thread section
    size_t ringsize_;//total number of workers assigned to GPUs
    int mytid_;//private thread id
    int myareano_;//device memory area number
    std::thread* tobj_;//thread object
private:
    //{{messaging
    // broadcasting attributes:
    int req_msg_;//request message issued for thread
    int msg_addressee_;//message addressee seen by thread
    // response attributes:
    mutable std::mutex mx_rsp_msg_;//mutex for messaging between the thread and the master
    std::condition_variable cv_rsp_msg_;//condition variable for messaging
    mutable std::mutex mx_dataccess_;//mutex for accessing class data
    std::condition_variable cv_dataccess_;//condition variable for checking for null data (ready)
    int rsp_msg_;//private response message to master
    //}}
    //
    //properties of the device the thread's communicating with:
    CuDeviceMemory* dmem_;
    //results writer:
    TdAlnWriter* alnwriter_;
    TdClustWriter* clustwriter_;
    //can be reinitialized on msg tthreadmsgGetDataChunkSize
    size_t chunkdatasize_;
    size_t chunkdatalen_;
    size_t chunknstrs_;
    //query-specific arguments: 
    size_t nqrsposs_;//total length of queries accessed under locked mutex mx_msg_!
    size_t mstr_set_nqrsposs_;//master-set length guarded by mutex mx_msg_!
    float mstr_set_scorethld_;
    //{{data arguments: 
    // master-set data/addresses:
    bool mstr_set_lastchunk_;
    int mstr_set_chunkno_;
    int mstr_set_qrysernrbeg_;
    bool mstr_set_newsetqrs_;//new set of queries
    char* mstr_set_queryndxpmbeg_[pmv2DTotIndexFlds];
    char* mstr_set_queryndxpmend_[pmv2DTotIndexFlds];
    const char** mstr_set_querydesc_;
    char* mstr_set_querypmbeg_[pmv2DTotFlds];
    char* mstr_set_querypmend_[pmv2DTotFlds];
    const char** mstr_set_bdbCdesc_;
    char* mstr_set_bdbCndxpmbeg_[pmv2DTotIndexFlds];
    char* mstr_set_bdbCndxpmend_[pmv2DTotIndexFlds];
    char* mstr_set_bdbCpmbeg_[pmv2DTotFlds];
    char* mstr_set_bdbCpmend_[pmv2DTotFlds];
    TSCounterVar* mstr_set_qrscnt_;//counter for queries
    TSCounterVar* mstr_set_cnt_;//counter for references
    // adresses saved by the thread
    float scorethld_;//score threshold
    bool lastchunk_;
    int chunkno_;//data chunk serial number
    int qrysernrbeg_;//serial number of the first query in the chunk
    bool newsetqrs_;//new set of queries
    char* queryndxpmbeg_[pmv2DTotIndexFlds];
    char* queryndxpmend_[pmv2DTotIndexFlds];
    const char** querydesc_;
    char* querypmbeg_[pmv2DTotFlds];
    char* querypmend_[pmv2DTotFlds];
    const char** bdbCdesc_;
    char* bdbCndxpmbeg_[pmv2DTotIndexFlds];
    char* bdbCndxpmend_[pmv2DTotIndexFlds];
    char* bdbCpmbeg_[pmv2DTotFlds];
    char* bdbCpmend_[pmv2DTotFlds];
    TSCounterVar* qrscnt_;
    TSCounterVar* cnt_;
    //}}
};

#endif//__TdCommutator_h__
