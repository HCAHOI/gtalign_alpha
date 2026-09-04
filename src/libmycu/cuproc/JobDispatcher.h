/***************************************************************************
 *   Copyright (C) 2021-2023 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __JobDispatcher_h__
#define __JobDispatcher_h__

#include "libutil/mybase.h"

#include <stdio.h>

#include <memory>
#include <vector>

#include "libutil/CLOptions.h"
#include "tsafety/TSCounterVar.h"

#include "libgenp/gdats/InputFilelist.h"
#include "libgenp/gdats/PMBatchStrData.h"
#include "libgenp/gdats/TdDataReader.h"
// #include "libmycu/culayout/CuDeviceMemory.cuh"
#include "libmycu/cubatch/TdCommutator.h"
#include "libgenp/goutp/TdClustWriter.h"
#include "libgenp/goutp/TdAlnWriter.h"
// #include "libmycu/cuproc/Devices.h"

#define JDSP_WORKER_BUSY -1
#define JDSP_WORKER_NONE -2

class JobDispatcher;

template<typename T>
using TWriteDataForWorker1 = 
        void (JobDispatcher::*)(
            int req, int addr,
            TdCommutator*, 
            const std::vector<T>&,
            std::unique_ptr<PMBatchStrData>);
template<typename T>
using TGetDataFromWorker1 = 
        void (JobDispatcher::*)(
            int req, int rsp, int addr,
            TdCommutator*,
            std::vector<T>&);

// _________________________________________________________________________
// Class JobDispatcher
//
// Implementation of distributing jobs over CPU and GPU threads
//
class JobDispatcher
{
public:
    /**
     * @brief 构造 `JobDispatcher`，初始化其负责的CUDA 设备管理状态。
     * @param inputlist 供该函数读取或更新的 `inputlist` 参数。
     * @param dnamelist 供该函数读取或更新的 `dnamelist` 参数。
     * @param sfxlst 供该函数读取或更新的 `sfxlst` 参数。
     * @param output 接收当前步骤输出的 `output`。
     * @param cachedir 供该函数读取或更新的 `cachedir` 参数。
     * @return 无返回值；完成对象构造与初始状态设置。
     */
    JobDispatcher(
            const std::vector<std::string>& inputlist,
            const std::vector<std::string>& dnamelist,
            const std::vector<std::string>& sfxlst, 
            const char* output,
            const char* cachedir
    );

    /**
     * @brief 构造 `JobDispatcher`，初始化其负责的CUDA 设备管理状态。
     * @param clustlist 供该函数读取或更新的 `clustlist` 参数。
     * @param sfxlst 供该函数读取或更新的 `sfxlst` 参数。
     * @param output 接收当前步骤输出的 `output`。
     * @param cachedir 供该函数读取或更新的 `cachedir` 参数。
     * @return 无返回值；完成对象构造与初始状态设置。
     */
    JobDispatcher(
            const std::vector<std::string>& clustlist,
            const std::vector<std::string>& sfxlst, 
            const char* output,
            const char* cachedir
    );

    /**
     * @brief 销毁 `JobDispatcher`，释放其管理的CUDA 设备管理资源。
     * @par 参数
     * 无。
     * @return 无返回值；对象持有的资源在返回前完成释放。
     */
    ~JobDispatcher();

    /**
     * @brief 在CUDA 设备管理中运行 `Run` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void Run();
    /**
     * @brief 在CUDA 设备管理中运行 `RunClust` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void RunClust();

    /**
     * @brief 在CUDA 设备管理中读取 `GetInputList` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    const std::vector<std::string>& GetInputList() const { return inputlist_; }
    /**
     * @brief 在CUDA 设备管理中读取 `GetDNamelist` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    const std::vector<std::string>& GetDNamelist() const { return dnamelist_; }
    /**
     * @brief 在CUDA 设备管理中读取 `GetClustList` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    const std::vector<std::string>& GetClustList() const { return clustlist_; }
    /**
     * @brief 在CUDA 设备管理中读取 `GetOutput` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    const char* GetOutput() const { return output_; }

protected:
    /**
     * @brief 在CUDA 设备管理中处理 `CreateReader` 对应的数据。
     * @param maxstrlen 控制当前步骤范围或规模的 `maxstrlen`。
     * @param mapped 供该函数读取或更新的 `mapped` 参数。
     * @param ndatbufs 控制当前步骤范围或规模的 `ndatbufs`。
     * @param nagents 控制当前步骤范围或规模的 `nagents`。
     * @param chunkdatasize 控制当前步骤范围或规模的 `chunkdatasize`。
     * @param chunkdatalen 控制当前步骤范围或规模的 `chunkdatalen`。
     * @param chunknstrs 供该函数读取或更新的 `chunknstrs` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void CreateReader( 
        int maxstrlen, bool mapped, int ndatbufs, int nagents,
        size_t chunkdatasize, size_t chunkdatalen, size_t chunknstrs
    );

    /**
     * @brief 在CUDA 设备管理中处理 `CreateQrsReader` 对应的数据。
     * @param maxstrlen 控制当前步骤范围或规模的 `maxstrlen`。
     * @param mapped 供该函数读取或更新的 `mapped` 参数。
     * @param ndatbufs 控制当前步骤范围或规模的 `ndatbufs`。
     * @param nagents 控制当前步骤范围或规模的 `nagents`。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void CreateQrsReader(
        int maxstrlen,
        bool mapped,
        int ndatbufs,
        int nagents
    );

    /**
     * @brief 在CUDA 设备管理中读取 `GetDataFromReader` 对应的数据。
     * @param reader 供该函数读取或更新的 `reader` 参数。
     * @param bdbCdesc 描述参考结构的 `bdbCdesc`。
     * @param bdbCpmbeg 参考结构打包字段的起始指针数组。
     * @param bdbCpmend 参考结构打包字段的结束指针数组。
     * @param bdbCNdxpmbeg 描述参考结构的 `bdbCNdxpmbeg`。
     * @param bdbCNdxpmend 描述参考结构的 `bdbCNdxpmend`。
     * @param tscnt 供该函数读取或更新的 `tscnt` 参数。
     * @param lastchunk 供该函数读取或更新的 `lastchunk` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    bool GetDataFromReader(
        TdDataReader* reader,
        char**& bdbCdesc, char**& bdbCpmbeg, char**& bdbCpmend,
        char**& bdbCNdxpmbeg, char**& bdbCNdxpmend,
        TSCounterVar*& tscnt,
        bool* lastchunk);

    /**
     * @brief 在CUDA 设备管理中读取 `GetDataFromReader` 对应的数据。
     * @param reader 供该函数读取或更新的 `reader` 参数。
     * @param bdbCdesc 描述参考结构的 `bdbCdesc`。
     * @param bdbCpmbeg 参考结构打包字段的起始指针数组。
     * @param bdbCpmend 参考结构打包字段的结束指针数组。
     * @param tscnt 供该函数读取或更新的 `tscnt` 参数。
     * @param lastchunk 供该函数读取或更新的 `lastchunk` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    bool GetDataFromReader(
        TdDataReader* reader,
        char**& bdbCdesc, char**& bdbCpmbeg, char**& bdbCpmend,
        TSCounterVar*& tscnt,
        bool* lastchunk);

    /**
     * @brief 在CUDA 设备管理中读取 `GetReferenceData` 对应的数据。
     * @param bdbCdescs 描述参考结构的 `bdbCdescs`。
     * @param bdbCpmbeg 参考结构打包字段的起始指针数组。
     * @param bdbCpmend 参考结构打包字段的结束指针数组。
     * @param bdbCNdxpmbeg 描述参考结构的 `bdbCNdxpmbeg`。
     * @param bdbCNdxpmend 描述参考结构的 `bdbCNdxpmend`。
     * @param tscnt 供该函数读取或更新的 `tscnt` 参数。
     * @param lastchunk 供该函数读取或更新的 `lastchunk` 参数。
     * @param rewind 供该函数读取或更新的 `rewind` 参数。
     * @param ntotqstrs 控制当前步骤范围或规模的 `ntotqstrs`。
     * @param queryblocks 描述查询结构的 `queryblocks`。
     * @param querypmbegs 描述查询结构的 `querypmbegs`。
     * @param querypmends 描述查询结构的 `querypmends`。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    bool GetReferenceData(
        char**& bdbCdescs, char**& bdbCpmbeg, char**& bdbCpmend,
        char**& bdbCNdxpmbeg, char**& bdbCNdxpmend,
        TSCounterVar*& tscnt, bool* lastchunk, bool rewind,
        size_t ntotqstrs, int queryblocks,
        char* const * const * const querypmbegs,
        char* const * const * const querypmends);

    /**
     * @brief 在CUDA 设备管理中处理 `CreateAlnWriter` 对应的数据。
     * @param outdirname 接收当前步骤输出的 `outdirname`。
     * @param dnamelist 供该函数读取或更新的 `dnamelist` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void CreateAlnWriter( 
        const char* outdirname,
        const std::vector<std::string>& dnamelist
    );

    /**
     * @brief 在CUDA 设备管理中处理 `CreateClustWriter` 对应的数据。
     * @param outdirname 接收当前步骤输出的 `outdirname`。
     * @param clustlist 供该函数读取或更新的 `clustlist` 参数。
     * @param devnames 供该函数读取或更新的 `devnames` 参数。
     * @param nmaxchunkqueries 控制当前步骤范围或规模的 `nmaxchunkqueries`。
     * @param nagents 控制当前步骤范围或规模的 `nagents`。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void CreateClustWriter( 
        const char* outdirname,
        const std::vector<std::string>& clustlist,
        const std::vector<std::string>& devnames,
        const int nmaxchunkqueries,
        const int nagents
    );

    /**
     * @brief 在CUDA 设备管理中通知 `NotifyAlnWriter` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void NotifyAlnWriter();
    /**
     * @brief 在CUDA 设备管理中等待 `WaitForAlnWriterToFinish` 对应的数据。
     * @param error 供该函数读取或更新的 `error` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void WaitForAlnWriterToFinish(bool error);

    /**
     * @brief 在CUDA 设备管理中通知 `NotifyClustWriter` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void NotifyClustWriter();
    /**
     * @brief 在CUDA 设备管理中等待 `WaitForClustWriterToFinish` 对应的数据。
     * @param error 供该函数读取或更新的 `error` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void WaitForClustWriterToFinish(bool error);


    /**
     * @brief 在CUDA 设备管理中处理 `CreateDevMemoryConfigs` 对应的数据。
     * @param nareasperdevice 控制当前步骤范围或规模的 `nareasperdevice`。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void CreateDevMemoryConfigs(size_t nareasperdevice);


    /**
     * @brief 在CUDA 设备管理中处理 `CreateWorkerThreads` 对应的数据。
     * @param chunkdatasize 控制当前步骤范围或规模的 `chunkdatasize`。
     * @param chunkdatalen 控制当前步骤范围或规模的 `chunkdatalen`。
     * @param chunknstrs 供该函数读取或更新的 `chunknstrs` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void CreateWorkerThreads(
        size_t chunkdatasize, size_t chunkdatalen, size_t chunknstrs);

    /**
     * @brief 在CUDA 设备管理中处理 `SubmitWorkerJob` 对应的数据。
     * @param tid 供该函数读取或更新的 `tid` 参数。
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
    void SubmitWorkerJob(
        int tid, int chunkno, bool lastchunk,
        bool newsetqrs,
        int qrysernrbeg,
        float scorethld,
        char** queryndxpmbeg, char** queryndxpmend,
        const char** querydesc, char** querypmbeg, char** querypmend,
        const char** bdbCdesc, char** bdbCpmbeg, char** bdbCpmend,
        char** bdbCndxpmbeg, char** bdbCndxpmend,
        TSCounterVar* qrstscnt, TSCounterVar* tscnt);

    /**
     * @brief 在CUDA 设备管理中等待 `WaitForWorker` 对应的数据。
     * @param tid 供该函数读取或更新的 `tid` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void WaitForWorker(int tid);

    /**
     * @brief 在CUDA 设备管理中处理 `ProbeWorker` 对应的数据。
     * @param tid 供该函数读取或更新的 `tid` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void ProbeWorker(int tid);

    /**
     * @brief 在CUDA 设备管理中处理 `TerminateWorker` 对应的数据。
     * @param tid 供该函数读取或更新的 `tid` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void TerminateWorker(int tid);

    /**
     * @brief 在CUDA 设备管理中处理 `TerminateAllWorkers` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void TerminateAllWorkers();

    /**
     * @brief 在CUDA 设备管理中等待 `WaitForAllWorkersToFinish` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void WaitForAllWorkersToFinish();

    /**
     * @brief 在CUDA 设备管理中读取 `GetAvailableWorker` 对应的数据。
     * @param tid 供该函数读取或更新的 `tid` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void GetAvailableWorker(int* tid);
    /**
     * @brief 在CUDA 设备管理中等待 `WaitForAvailableWorker` 对应的数据。
     * @param tid 供该函数读取或更新的 `tid` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void WaitForAvailableWorker(int* tid);
    /**
     * @brief 在CUDA 设备管理中读取 `GetNextWorker` 对应的数据。
     * @param tid 供该函数读取或更新的 `tid` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void GetNextWorker(int* tid);



    /**
     * @brief 在CUDA 设备管理中处理 `ProcessPart` 对应的数据。
     * @param tid 供该函数读取或更新的 `tid` 参数。
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
    void ProcessPart(
        int tid, int chunkno, bool lastchunk,
        bool newsetqrs,
        int qrysernrbeg,
        float scorethld,
        char** queryndxpmbeg, char** queryndxpmend,
        const char** querydesc, char** querypmbeg, char** querypmend,
        const char** bdbCdesc, char** bdbCpmbeg, char** bdbCpmend,
        char** bdbCndxpmbeg, char** bdbCndxpmend,
        TSCounterVar* qrstscnt, TSCounterVar* tscnt);



    /**
     * @brief 在CUDA 设备管理中读取 `GetAlnWriter` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    const TdAlnWriter* GetAlnWriter() const {return writer_;}
    /**
     * @brief 在CUDA 设备管理中读取 `GetAlnWriter` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    TdAlnWriter* GetAlnWriter() {return writer_;}

    /**
     * @brief 在CUDA 设备管理中读取 `GetClustWriter` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    const TdClustWriter* GetClustWriter() const {return clustwriter_;}
    /**
     * @brief 在CUDA 设备管理中读取 `GetClustWriter` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    TdClustWriter* GetClustWriter() {return clustwriter_;}

    /**
     * @brief 在CUDA 设备管理中读取 `GetQrsReader` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    const TdDataReader* GetQrsReader() const {return qrsreader_;}
    /**
     * @brief 在CUDA 设备管理中读取 `GetQrsReader` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    TdDataReader* GetQrsReader() {return qrsreader_;}

    //set chunk attributes for each reference data reader
    /**
     * @brief 在CUDA 设备管理中设置 `SetChunkDataAttributesReaders` 对应的数据。
     * @param chunkdatasize 控制当前步骤范围或规模的 `chunkdatasize`。
     * @param chunkdatalen 控制当前步骤范围或规模的 `chunkdatalen`。
     * @param chunknstrs 供该函数读取或更新的 `chunknstrs` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void SetChunkDataAttributesReaders(
        size_t chunkdatasize, size_t chunkdatalen, size_t chunknstrs)
    {
        std::for_each(readers_.begin(), readers_.end(),
            [chunkdatasize, chunkdatalen, chunknstrs](std::unique_ptr<TdDataReader>& p) {
                if(p) p->SetChunkDataAttributes(
                        chunkdatasize, chunkdatalen, chunknstrs);
            }
        );
    }

private:
    const char* output_;//pattern for output file (null=standard output)
    const char* cachedir_;//directory for cached data
    std::vector<std::string> inputlist_;//input files/databases
    std::vector<std::string> dnamelist_;//reference structure files/databases
    std::vector<std::string> clustlist_;//input files/databases for clustering
    std::vector<std::string> sfxlst_;//suffix list 
    std::vector<CuDeviceMemory*> memdevs_;//memory configurations for devices
    std::vector<TdCommutator*> hostworkers_;//host worker threads

    std::unique_ptr<InputFilelist> queries_;
    std::unique_ptr<InputFilelist> references_;

    TdAlnWriter* writer_;//alignment results writer
    TdClustWriter* clustwriter_;//clustering results writer
    TdDataReader* qrsreader_;//reader of query structure files and databases
    std::vector<std::unique_ptr<TdDataReader>> readers_;//structure files and database readers
};


////////////////////////////////////////////////////////////////////////////
// INLINES
//

#endif//__JobDispatcher_h__
