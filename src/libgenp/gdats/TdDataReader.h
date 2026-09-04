/***************************************************************************
 *   Copyright (C) 2021-2023 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __TdDataReader_h__
#define __TdDataReader_h__

#include "libutil/mybase.h"

#include <stdio.h>
#include <string.h>

#include <memory>
#include <mutex>
#include <condition_variable>
#include <thread>
#include <vector>

#include "tsafety/TSCounterVar.h"
#include "PM2DVectorFields.h"
#include "PMBatchStrData.h"
#include "PMBatchStrDataIndex.h"
#include "FlexDataRead.h"

#define TREADER_MSG_UNSET -1
#define TREADER_MSG_ERROR -2
#define TREADER_MSG_EXIT -3

class TdDataReader;

using TGetNextData = void (TdDataReader::*)(
        size_t chunkdatasize, size_t chunkdatalen, size_t chunknstrs,
        const size_t ntotqstrs,
        const int queryblocks,
        char* const * const * const querypmbegs,
        char* const * const * const querypmends);

using TReadDataChunk = bool (TdDataReader::*)(
        PMBatchStrDataIndex*,
        PMBatchStrData*, const PMBatchStrData*, 
        size_t chunkdatasize, size_t chunkdatalen, size_t chunknstrs,
        const size_t ntotqstrs,
        const int queryblocks,
        char* const * const * const querypmbegs,
        char* const * const * const querypmends);

// _________________________________________________________________________
// Class TdDataReader
//
// thread class responsible for reading profile data
//
class TdDataReader
{
public:
    enum {
        tdrcheNotSet,
        tdrcheNoCaching,
        tdrcheCacheAndRead,
        tdrcheReadCached
    };
    enum TDRMsg {
        tdrmsgGetSize,
        tdrmsgGetData,
        tdrmsgTerminate
    };
    enum TDRResponseMsg {
        tdrrespmsgSize,
        tdrrespmsgDataReady,
        tdrrespmsgNoData,
        tdrrespmsgTerminating
    };

public:
    /**
     * @brief 构造 `TdDataReader`，初始化其负责的结构数据读取与布局状态。
     * @param cachedir 供该函数读取或更新的 `cachedir` 参数。
     * @param inputlist 供该函数读取或更新的 `inputlist` 参数。
     * @param strfilelist 供该函数读取或更新的 `strfilelist` 参数。
     * @param pntfilelist 供该函数读取或更新的 `pntfilelist` 参数。
     * @param strfilepositionlist 供该函数读取或更新的 `strfilepositionlist` 参数。
     * @param strfilesizelist 控制当前步骤范围或规模的 `strfilesizelist`。
     * @param strparenttypelist 供该函数读取或更新的 `strparenttypelist` 参数。
     * @param strfiletypelist 供该函数读取或更新的 `strfiletypelist` 参数。
     * @param filendxlist 控制当前步骤范围或规模的 `filendxlist`。
     * @param globalids 供该函数读取或更新的 `globalids` 参数。
     * @param ndxstartwith 控制当前步骤范围或规模的 `ndxstartwith`。
     * @param ndxstep 控制当前步骤范围或规模的 `ndxstep`。
     * @param maxstrlen 控制当前步骤范围或规模的 `maxstrlen`。
     * @param mapped 供该函数读取或更新的 `mapped` 参数。
     * @param indexed 供该函数读取或更新的 `indexed` 参数。
     * @param ndatbufs 控制当前步骤范围或规模的 `ndatbufs`。
     * @param nagents 控制当前步骤范围或规模的 `nagents`。
     * @param clustering 供该函数读取或更新的 `clustering` 参数。
     * @param clustmaster 供该函数读取或更新的 `clustmaster` 参数。
     * @return 无返回值；完成对象构造与初始状态设置。
     */
    TdDataReader(
        const char* cachedir,
        const std::vector<std::string>& inputlist,
        const std::vector<std::string>& strfilelist,
        const std::vector<std::string>& pntfilelist,
        const std::vector<size_t>& strfilepositionlist,
        const std::vector<size_t>& strfilesizelist,
        const std::vector<int>& strparenttypelist,
        const std::vector<int>& strfiletypelist,
        const std::vector<size_t>& filendxlist,
        std::vector<std::vector<int>>& globalids,
        const size_t ndxstartwith,
        const size_t ndxstep,
        int maxstrlen,
        bool mapped,
        bool indexed,
        int ndatbufs,
        int nagents,
        const bool clustering = false,
        const bool clustmaster = false
    );
    /**
     * @brief 构造 `TdDataReader`，初始化其负责的结构数据读取与布局状态。
     * @par 参数
     * 无。
     * @return 无返回值；完成对象构造与初始状态设置。
     */
    TdDataReader();
    /**
     * @brief 销毁 `TdDataReader`，释放其管理的结构数据读取与布局资源。
     * @par 参数
     * 无。
     * @return 无返回值；对象持有的资源在返回前完成释放。
     */
    ~TdDataReader();

    /**
     * @brief 在结构数据读取与布局中读取 `GetDbsCacheFlag` 对应的数据。
     * @param nagents 控制当前步骤范围或规模的 `nagents`。
     * @param cachedir 供该函数读取或更新的 `cachedir` 参数。
     * @param inputlist 供该函数读取或更新的 `inputlist` 参数。
     * @param chunkdatasize 控制当前步骤范围或规模的 `chunkdatasize`。
     * @param chunkdatalen 控制当前步骤范围或规模的 `chunkdatalen`。
     * @param chunknstrs 供该函数读取或更新的 `chunknstrs` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    static int GetDbsCacheFlag(
        int nagents, 
        const char* cachedir,
        const std::vector<std::string>& inputlist,
        size_t chunkdatasize, size_t chunkdatalen, size_t chunknstrs);

    // std::mutex& GetPrivateMutex() {return mx_dataccess_;}

    //{{NOTE: messaging functions accessed from outside!
    /**
     * @brief 在结构数据读取与布局中通知 `Notify` 对应的数据。
     * @param msg 供该函数读取或更新的 `msg` 参数。
     * @param ntotqstrs 控制当前步骤范围或规模的 `ntotqstrs`。
     * @param queryblocks 描述查询结构的 `queryblocks`。
     * @param querypmbegs 描述查询结构的 `querypmbegs`。
     * @param querypmends 描述查询结构的 `querypmends`。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void Notify(
        int msg,
        size_t ntotqstrs = 0,
        int queryblocks = 0,
        char* const * const * const querypmbegs = NULL,
        char* const * const * const querypmends = NULL)
    {
        {//mutex must be unlocked before notifying
            std::lock_guard<std::mutex> lck(mx_dataccess_);
            req_msg_ = msg;
            ntotqstrs_ = ntotqstrs;
            queryblocks_ = queryblocks;
            querypmbegs_ = querypmbegs;
            querypmends_ = querypmends;
        }
        cv_msg_.notify_one();
    }
    /**
     * @brief 在结构数据读取与布局中等待 `Wait` 对应的数据。
     * @param rsp1 供该函数读取或更新的 `rsp1` 参数。
     * @param rsp2 供该函数读取或更新的 `rsp2` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    int Wait(int rsp1, int rsp2 = TREADER_MSG_ERROR) {
        //wait until a response arrives
        std::unique_lock<std::mutex> lck_msg(mx_dataccess_);
        cv_msg_.wait(lck_msg,
            [this,rsp1,rsp2]{
                return (rsp_msg_ == rsp1 || rsp_msg_ == rsp2 || 
                        rsp_msg_ == TREADER_MSG_ERROR ||
                        rsp_msg_ == TREADER_MSG_EXIT);
            }
        );
        //lock is back; unset the response
        int rspmsg = rsp_msg_;
        //NOTE: master may change rsp after the reader set it repeatedly!
        //NOTE: may lead to dead lock when chunks delivered non-continuously!
        if( rsp_msg_!= TREADER_MSG_ERROR )
            rsp_msg_ = TREADER_MSG_UNSET;
        return rspmsg;
    }
    /**
     * @brief 在结构数据读取与布局中处理 `Rewind` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void Rewind() {
        std::lock_guard<std::mutex> lck(mx_dataccess_);
        rewind_ = true;
    }
    /**
     * @brief 在结构数据读取与布局中读取 `GetResponseAsync` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    int GetResponseAsync() const {
        //get a response if available
        std::unique_lock<std::mutex> lck(mx_dataccess_, std::defer_lock);
        int rsp = TREADER_MSG_UNSET;
        if(lck.try_lock()) rsp = rsp_msg_;
        return rsp;
    }
    /**
     * @brief 在结构数据读取与布局中读取 `GetResponse` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    int GetResponse() const {
        std::lock_guard<std::mutex> lck(mx_dataccess_);
        return rsp_msg_;
    }
    /**
     * @brief 在结构数据读取与布局中重置 `ResetResponse` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void ResetResponse() {
        std::lock_guard<std::mutex> lck(mx_dataccess_);
        if( rsp_msg_!= TREADER_MSG_ERROR )
            rsp_msg_ = TREADER_MSG_UNSET;
    }
    //}}

    /**
     * @brief 在结构数据读取与布局中读取 `GetbdbCdata` 对应的数据。
     * @param bdbCdescs 描述参考结构的 `bdbCdescs`。
     * @param bdbCpmbeg 参考结构打包字段的起始指针数组。
     * @param bdbCpmend 参考结构打包字段的结束指针数组。
     * @param bdbCNdxpmbeg 描述参考结构的 `bdbCNdxpmbeg`。
     * @param bdbCNdxpmend 描述参考结构的 `bdbCNdxpmend`。
     * @param tscnt 供该函数读取或更新的 `tscnt` 参数。
     * @param lastchunk 供该函数读取或更新的 `lastchunk` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void GetbdbCdata(
        char**& bdbCdescs,
        char**& bdbCpmbeg, char**& bdbCpmend,
        char**& bdbCNdxpmbeg, char**& bdbCNdxpmend,
        TSCounterVar*& tscnt,
        bool* lastchunk )
    {
        std::lock_guard<std::mutex> lck(mx_dataccess_);
        bdbCNdxpmbeg = bdbCNdxpmbeg_;
        bdbCNdxpmend = bdbCNdxpmend_;
        bdbCdescs = bdbCptrdescs_;
        bdbCpmbeg = bdbCpmbeg_;
        bdbCpmend = bdbCpmend_;
        tscnt = ccnt_;
        *lastchunk = lastchunk_;
    }

    /**
     * @brief 在结构数据读取与布局中读取 `GetbdbCdata` 对应的数据。
     * @param bdbCdescs 描述参考结构的 `bdbCdescs`。
     * @param bdbCpmbeg 参考结构打包字段的起始指针数组。
     * @param bdbCpmend 参考结构打包字段的结束指针数组。
     * @param tscnt 供该函数读取或更新的 `tscnt` 参数。
     * @param lastchunk 供该函数读取或更新的 `lastchunk` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void GetbdbCdata(
        char**& bdbCdescs,
        char**& bdbCpmbeg, char**& bdbCpmend,
//         size_t*& szpm2dvfields,
        TSCounterVar*& tscnt,
        bool* lastchunk )
    {
        std::lock_guard<std::mutex> lck(mx_dataccess_);
        bdbCdescs = bdbCptrdescs_;
        bdbCpmbeg = bdbCpmbeg_;
        bdbCpmend = bdbCpmend_;
//         szpm2dvfields = szpm2dvfields_;
        tscnt = ccnt_;
        *lastchunk = lastchunk_;

//         if( bdbCpmbeg && bdbCpmend ) {
//             memcpy( bdbCpmbeg, bdbCpmbeg_, pmv2DTotFlds * sizeof(void*));
//             memcpy( bdbCpmend, bdbCpmend_, pmv2DTotFlds * sizeof(void*));
//         }
//         if( szpm2dvfields )
//             memcpy( szpm2dvfields, szpm2dvfields_, pmv2DTotFlds * sizeof(size_t));
    }

    /**
     * @brief 在结构数据读取与布局中设置 `SetChunkDataAttributes` 对应的数据。
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

protected:
    /**
     * @brief 在结构数据读取与布局中处理 `Execute` 对应的数据。
     * @param args 供该函数读取或更新的 `args` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void Execute(void* args);

    /**
     * @brief 在结构数据读取与布局中更新 `UpdateCacheFlag` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void UpdateCacheFlag();

    /**
     * @brief 在结构数据读取与布局中读取 `GetData` 对应的数据。
     * @param chunkdatasize 控制当前步骤范围或规模的 `chunkdatasize`。
     * @param chunkdatalen 控制当前步骤范围或规模的 `chunkdatalen`。
     * @param chunknstrs 供该函数读取或更新的 `chunknstrs` 参数。
     * @param ntotqstrs 控制当前步骤范围或规模的 `ntotqstrs`。
     * @param queryblocks 描述查询结构的 `queryblocks`。
     * @param querypmbegs 描述查询结构的 `querypmbegs`。
     * @param querypmends 描述查询结构的 `querypmends`。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    bool GetData(
        size_t chunkdatasize, size_t chunkdatalen, size_t chunknstrs,
        const size_t ntotqstrs,
        const int queryblocks,
        char* const * const * const querypmbegs,
        char* const * const * const querypmends);
    /**
     * @brief 在结构数据读取与布局中读取 `GetNextData` 对应的数据。
     * @param chunkdatasize 控制当前步骤范围或规模的 `chunkdatasize`。
     * @param chunkdatalen 控制当前步骤范围或规模的 `chunkdatalen`。
     * @param chunknstrs 供该函数读取或更新的 `chunknstrs` 参数。
     * @param ntotqstrs 控制当前步骤范围或规模的 `ntotqstrs`。
     * @param queryblocks 描述查询结构的 `queryblocks`。
     * @param querypmbegs 描述查询结构的 `querypmbegs`。
     * @param querypmends 描述查询结构的 `querypmends`。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void GetNextData(
        size_t chunkdatasize, size_t chunkdatalen, size_t chunknstrs,
        const size_t ntotqstrs,
        const int queryblocks,
        char* const * const * const querypmbegs,
        char* const * const * const querypmends);
    /**
     * @brief 在结构数据读取与布局中读取 `GetNextDataClustCache` 对应的数据。
     * @param chunkdatasize 控制当前步骤范围或规模的 `chunkdatasize`。
     * @param chunkdatalen 控制当前步骤范围或规模的 `chunkdatalen`。
     * @param chunknstrs 供该函数读取或更新的 `chunknstrs` 参数。
     * @param ntotqstrs 控制当前步骤范围或规模的 `ntotqstrs`。
     * @param queryblocks 描述查询结构的 `queryblocks`。
     * @param querypmbegs 描述查询结构的 `querypmbegs`。
     * @param querypmends 描述查询结构的 `querypmends`。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void GetNextDataClustCache(
        size_t chunkdatasize, size_t chunkdatalen, size_t chunknstrs,
        const size_t ntotqstrs,
        const int queryblocks,
        char* const * const * const querypmbegs,
        char* const * const * const querypmends);
    /**
     * @brief 在结构数据读取与布局中读取 `ReadDataChunk` 对应的数据。
     * @param PMBatchStrDataIndex 供该函数读取或更新的 `PMBatchStrDataIndex` 参数。
     * @param PMBatchStrData 供该函数读取或更新的 `PMBatchStrData` 参数。
     * @param PMBatchStrData 供该函数读取或更新的 `PMBatchStrData` 参数。
     * @param chunkdatasize 控制当前步骤范围或规模的 `chunkdatasize`。
     * @param chunkdatalen 控制当前步骤范围或规模的 `chunkdatalen`。
     * @param chunknstrs 供该函数读取或更新的 `chunknstrs` 参数。
     * @param ntotqstrs 控制当前步骤范围或规模的 `ntotqstrs`。
     * @param queryblocks 描述查询结构的 `queryblocks`。
     * @param querypmbegs 描述查询结构的 `querypmbegs`。
     * @param querypmends 描述查询结构的 `querypmends`。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    bool ReadDataChunk(
        PMBatchStrDataIndex*,
        PMBatchStrData*, const PMBatchStrData*, 
        size_t chunkdatasize, size_t chunkdatalen, size_t chunknstrs,
        const size_t ntotqstrs,
        const int queryblocks,
        char* const * const * const querypmbegs,
        char* const * const * const querypmends);
    /**
     * @brief 在结构数据读取与布局中读取 `ReadDataChunkClustCache` 对应的数据。
     * @param PMBatchStrDataIndex 供该函数读取或更新的 `PMBatchStrDataIndex` 参数。
     * @param PMBatchStrData 供该函数读取或更新的 `PMBatchStrData` 参数。
     * @param PMBatchStrData 供该函数读取或更新的 `PMBatchStrData` 参数。
     * @param chunkdatasize 控制当前步骤范围或规模的 `chunkdatasize`。
     * @param chunkdatalen 控制当前步骤范围或规模的 `chunkdatalen`。
     * @param chunknstrs 供该函数读取或更新的 `chunknstrs` 参数。
     * @param ntotqstrs 控制当前步骤范围或规模的 `ntotqstrs`。
     * @param queryblocks 描述查询结构的 `queryblocks`。
     * @param querypmbegs 描述查询结构的 `querypmbegs`。
     * @param querypmends 描述查询结构的 `querypmends`。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    bool ReadDataChunkClustCache(
        PMBatchStrDataIndex*,
        PMBatchStrData*, const PMBatchStrData*, 
        size_t chunkdatasize, size_t chunkdatalen, size_t chunknstrs,
        const size_t ntotqstrs,
        const int queryblocks,
        char* const * const * const querypmbegs,
        char* const * const * const querypmends);

private:
    //thread section
    std::thread* tobj_;//thread object
private:
    TGetNextData fgetnextdata;
    TReadDataChunk freaddatachunk;
    //{{caching
    int cacheflag_;
    //}}
    //{{messaging
    std::condition_variable cv_msg_;//condition variable for messaging
    mutable std::mutex mx_dataccess_;//mutex for accessing data
    int req_msg_;//request message issued for thread
    int rsp_msg_;//private response message
    //}}
    //{{variables that determine the size of data chunks being read
    size_t chunkdatasize_;
    size_t chunkdatalen_;
    size_t chunknstrs_;
    //}}
    //{{db name and configuration:
    std::unique_ptr<FlexDataRead> dbobj_;//ptr to database(s)
    const char* cachedir_;//cache directory
    const std::vector<std::string>& inputlist_;//input list
    const std::vector<std::string>& strfilelist_;//structure file list
    const std::vector<std::string>& pntfilelist_;//parent file list
    const std::vector<size_t>& strfilepositionlist_;//list of structure file positions within an archive
    const std::vector<size_t>& strfilesizelist_;//list of structure file sizes
    const std::vector<int>& strparenttypelist_;//parent file type list of structure files (e.g, tar)
    const std::vector<int>& strfiletypelist_;//list of structure file types
    const std::vector<size_t>& filendxlist_;//list of the indices of files sorted by filesize
    std::vector<std::vector<int>>& globalids_;//global ids for structures across all files
    const size_t ndxstartwith_;//file index to start with
    const size_t ndxstep_;//file index step
    size_t ndxpacketfile_;//index of a packet file for caching
    size_t ndxlastpacketfile_;//index of the last packet file analyzed for clustering with caching
    const int maxstrlen_;//max structure length
    bool mapped_;//database mapped
    bool indexed_;//database to be indexed
    //}}
    //{{data:
//     std::unique_ptr<int,DRDataDeleter> bdbClengths_;//profile lengths read
//     size_t bdbClengths_from_, bdbClengths_to_;//number of lengths read (in the number of profiles)
//     std::unique_ptr<size_t,DRDataDeleter> bdbCdesc_end_addrs_;//description end addresses read (numbers as above)
    //
    PMBatchStrData bdbC4CC;//farthest data read for clustering with caching
    std::vector<PMBatchStrData> bdbCstruct_;//read data
    std::vector<PMBatchStrDataIndex> bdbCindex_;//indexed data
    int pbdbCstruct_ndx_;//index of the data block to be read 
    //
//     size_t bdbCdata_from_, bdbCdata_to_;//profile data read in the number of profiles
//     size_t bdbCdata_poss_from_, bdbCdata_poss_to_;//profile data read in the number of positions
//     size_t addrdescproced_;//address of the last description processed
    bool eodclustcached_;//flag of final EOD for clustering with caching
    bool endofdata_;//flag to indicate that no data read on the last call
    bool recycled_;//starting over the full cycle of data and the first data chunk has been read again
    bool rewind_;//flag set by the master and this thread for rewinding database
    //}}
    //{{addresses of indexed data to return
    char** bdbCNdxpmbeg_;//addresses of the field beginnings
    char** bdbCNdxpmend_;//addresses of the field endings
    //}}
    //{{addresses of read data to return
    char** bdbCptrdescs_;//structure descriptions
    char** bdbCpmbeg_;//addresses of the beginnings of the fields
    char** bdbCpmend_;//addresses of the endings of the fields
//     size_t* szpm2dvfields_;//beginnings in bytes (sizes) of the fields
    TSCounterVar* ccnt_;//counter associated with data to be return
    bool lastchunk_;//flag of whether the last chunk has been read
    int nagents_;//number of agents accessing the data
    //
    const bool clustering_;//flag indicating clustering
    const bool clustmaster_;//when clustering is on, the master will create global ids
    //}}
    size_t ntotqstrs_;
    //{{pairwise sequence similarity verification:
    int queryblocks_;
    char* const * const * querypmbegs_;
    char* const * const * querypmends_;
    //}}
    static const char* signaturefile_;
    static const char* packetfilebasename_;
};

#endif//__TdDataReader_h__
