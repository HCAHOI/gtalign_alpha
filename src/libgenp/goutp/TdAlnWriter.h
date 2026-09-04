/***************************************************************************
 *   Copyright (C) 2021-2023 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __TdAlnWriter_h__
#define __TdAlnWriter_h__

#include "libutil/mybase.h"

#include <stdio.h>
#include <cmath>

#include <string>
#include <memory>
#include <utility>
#include <functional>
#include <vector>
#include <mutex>
#include <condition_variable>
#include <thread>

#include "libutil/CLOptions.h"

#define WRITERTHREAD_MSG_UNSET -1
#define WRITERTHREAD_MSG_ERROR -2

// -------------------------------------------------------------------------

struct WritersDataDestroyer {
    /**
     * @brief 在对齐结果输出中处理 `operator()` 对应的数据。
     * @param p 供该函数读取或更新的 `p` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void operator()(char* p) const {
        std::free(p);
    };
};

// _________________________________________________________________________
// Class TdAlnWriter
//
// alignment writer thread
//
class TdAlnWriter
{
    enum {
        szTmpBuffer = KBYTE,
        szWriterBuffer = TIMES4(KBYTE)
    };

public:
    enum TWriterThreadMsg {
        wrtthreadmsgWrite,
        wrtthreadmsgTerminate
    };
    enum TWriterThreadResponse {
        wrttrespmsgWriting,
        wrttrespmsgTerminating
    };

public:
    /**
     * @brief 构造 `TdAlnWriter`，初始化其负责的对齐结果输出状态。
     * @param outdirname 接收当前步骤输出的 `outdirname`。
     * @param rfilelist 供该函数读取或更新的 `rfilelist` 参数。
     * @return 无返回值；完成对象构造与初始状态设置。
     */
    TdAlnWriter( 
        const char* outdirname,
        const std::vector<std::string>& rfilelist
    );

    /**
     * @brief 销毁 `TdAlnWriter`，释放其管理的对齐结果输出资源。
     * @par 参数
     * 无。
     * @return 无返回值；对象持有的资源在返回前完成释放。
     */
    ~TdAlnWriter();

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
            std::lock_guard<std::mutex> lck(mx_dataccess_);
            req_msg_ = msg;
        }
        cv_msg_.notify_all();
    }
    /**
     * @brief 在对齐结果输出中等待 `Wait` 对应的数据。
     * @param rsp 供该函数读取或更新的 `rsp` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    int Wait(int rsp) {
        //wait until the response is received
        std::unique_lock<std::mutex> lck_msg(mx_dataccess_);
        cv_msg_.wait(lck_msg,
            [this,rsp]{
                return (rsp_msg_ == rsp ||
                        rsp_msg_ == WRITERTHREAD_MSG_ERROR);
            }
        );
        //lock is back; unset the response
        int rspmsg = rsp_msg_;
        if(rsp_msg_ != WRITERTHREAD_MSG_ERROR)
            rsp_msg_ = WRITERTHREAD_MSG_UNSET;
        return rspmsg;
    }
    /**
     * @brief 在对齐结果输出中等待 `WaitDone` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    int WaitDone() {
        for(size_t i = 0; i < parts_qrs_.size();)
        {
            //wait until all queries have been processed
            std::unique_lock<std::mutex> lck_msg(mx_dataccess_);
            cv_msg_.wait(lck_msg,
                [this,i]{
                    return (parts_qrs_[i] <= 0 ||
                        rsp_msg_ == WRITERTHREAD_MSG_ERROR);
                }
            );
            if(rsp_msg_ == WRITERTHREAD_MSG_ERROR)
                return rsp_msg_;
            if(req_msg_ != WRITERTHREAD_MSG_UNSET)
                continue;
            i++;
        }
        return rsp_msg_;
    }
    /**
     * @brief 在对齐结果输出中处理 `IncreaseQueryNParts` 对应的数据。
     * @param qrysernrfrom 描述查询结构的 `qrysernrfrom`。
     * @param qrysernrto 描述查询结构的 `qrysernrto`。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void IncreaseQueryNParts(int qrysernrfrom, int qrysernrto) {
        std::lock_guard<std::mutex> lck(mx_dataccess_);
        if(qrysernrto < qrysernrfrom)
            return;
        if((int)parts_qrs_.size() <= qrysernrto) {
            //let the grand master decide on triggering the write by
            //initializing #parts to 1
            parts_qrs_.resize(qrysernrto+1, 0);
            ResizeVectors(qrysernrto+1);
        }
        for(int qsn = qrysernrfrom; qsn <= qrysernrto; qsn++)
            parts_qrs_[qsn]++;
    }
//     void DereaseNPartsAndTrigger(int qrysernr) {
//         std::unique_lock<std::mutex> lck(mx_dataccess_);
//         cv_msg_.wait(lck,
//             [this]{return req_msg_ == WRITERTHREAD_MSG_UNSET;}
//         );
//         if( --parts_qrs_[qrysernr] <= 0 ) {
//             //this is the last part for the given query:
//             //trigger write to a file
//             qrysernr_ = qrysernr;
//             req_msg_ = wrtthreadmsgWrite;
//             lck.unlock();
//             cv_msg_.notify_all();
//         }
//     }
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
//     void ResetResponse() {
//         std::lock_guard<std::mutex> lck(mx_dataccess_);
//         if( rsp_msg_!= WRITERTHREAD_MSG_ERROR )
//             rsp_msg_ = WRITERTHREAD_MSG_UNSET;
//     }
    //}}

    /**
     * @brief 在对齐结果输出中处理 `PushPartOfResults` 对应的数据。
     * @param qrysernr 描述查询结构的 `qrysernr`。
     * @param nqyposs 当前批次中查询结构的总位置数。
     * @param qrydesc 描述查询结构的 `qrydesc`。
     * @param devanme 供该函数读取或更新的 `devanme` 参数。
     * @param nqystrs 当前批次中的查询结构数量。
     * @param duration 供该函数读取或更新的 `duration` 参数。
     * @param tmsthld 供该函数读取或更新的 `tmsthld` 参数。
     * @param ndbCposs 当前批次中参考结构的总位置数。
     * @param ndbCstrs 当前批次中的参考结构数量。
     * @param annotations 供该函数读取或更新的 `annotations` 参数。
     * @param alignments 供该函数读取或更新的 `alignments` 参数。
     * @param srtindxs 供该函数读取或更新的 `srtindxs` 参数。
     * @param tmscores 当前步骤使用或写回的 `tmscores` 分数。
     * @param alnptrs 供该函数读取或更新的 `alnptrs` 参数。
     * @param annotptrs 供该函数读取或更新的 `annotptrs` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void PushPartOfResults( 
        int qrysernr,
        int nqyposs,
        const std::string& qrydesc,
        const std::string& devanme,
        const size_t nqystrs,
        const double duration,
        const float tmsthld,
        const size_t ndbCposs,
        const size_t ndbCstrs,
        std::unique_ptr<char,WritersDataDestroyer> annotations,
        std::unique_ptr<char,WritersDataDestroyer> alignments,
        std::unique_ptr<std::vector<int>> srtindxs,
        std::unique_ptr<std::vector<float>> tmscores,
        std::unique_ptr<std::vector<char*>> alnptrs,
        std::unique_ptr<std::vector<char*>> annotptrs)
    {
        std::unique_lock<std::mutex> lck(mx_dataccess_);
        //{{NOTE: [inserted]
        cv_msg_.wait(lck,
            [this]{
                return req_msg_ == WRITERTHREAD_MSG_UNSET ||
                    rsp_msg_ == WRITERTHREAD_MSG_ERROR ||
                    rsp_msg_ == wrttrespmsgTerminating;
            }
        );
        //}}
        if(rsp_msg_ == WRITERTHREAD_MSG_ERROR ||
           rsp_msg_ == wrttrespmsgTerminating)
            return;
        //
        if((int)parts_qrs_.size() <= qrysernr || qrysernr < 0)
            throw MYRUNTIME_ERROR(
            "TdAlnWriter::PushPartOfResults: Invalid query serial number.");
        vec_duration_[qrysernr] += duration;
        vec_nposschd_[qrysernr] += ndbCposs;
        vec_nentries_[qrysernr] += ndbCstrs;
        vec_nqystrs_[qrysernr] = nqystrs;
        vec_nqyposs_[qrysernr] = nqyposs;
        vec_qrydesc_[qrysernr] = qrydesc;
        vec_devname_[qrysernr] = devanme;
        vec_tmsthld_[qrysernr] = tmsthld;
        //
        vec_annotations_[qrysernr].push_back(std::move(annotations));
        vec_alignments_[qrysernr].push_back(std::move(alignments));
        vec_srtindxs_[qrysernr].push_back(std::move(srtindxs));
        vec_tmscores_[qrysernr].push_back(std::move(tmscores));
        vec_alnptrs_[qrysernr].push_back(std::move(alnptrs));
        vec_annotptrs_[qrysernr].push_back(std::move(annotptrs));
        //{{NOTE: [commented out] the following statement must be the last
        //lck.unlock();
        //DereaseNPartsAndTrigger(qrysernr);
        //}}
        if( --parts_qrs_[qrysernr] <= 0 ) {
            //this is the last part for the given query:
            //trigger write to a file
            qrysernr_ = qrysernr;
            req_msg_ = wrtthreadmsgWrite;
            lck.unlock();
            cv_msg_.notify_all();
        }
    }


public:
    /**
     * @brief 在对齐结果输出中写出 `WritePrognamePlain` 对应的数据。
     * @param outptr 接收当前步骤输出的 `outptr`。
     * @param maxsize 控制当前步骤范围或规模的 `maxsize`。
     * @param width 供该函数读取或更新的 `width` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    static int WritePrognamePlain( char*& outptr, int maxsize, const int width );
    /**
     * @brief 在对齐结果输出中写出 `WriteCommandLinePlain` 对应的数据。
     * @param fp 供该函数读取或更新的 `fp` 参数。
     * @param buffer 供当前步骤读取或更新的 `buffer` 缓冲区。
     * @param szbuffer 供当前步骤读取或更新的 `szbuffer` 缓冲区。
     * @param outptr 接收当前步骤输出的 `outptr`。
     * @param offset 供该函数读取或更新的 `offset` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    static void WriteCommandLinePlain(FILE* fp,
        char* const buffer, const int szbuffer, char*& outptr, int& offset);
    /**
     * @brief 在对齐结果输出中写出 `WriteCommandLineJSON` 对应的数据。
     * @param fp 供该函数读取或更新的 `fp` 参数。
     * @param buffer 供当前步骤读取或更新的 `buffer` 缓冲区。
     * @param szbuffer 供当前步骤读取或更新的 `szbuffer` 缓冲区。
     * @param outptr 接收当前步骤输出的 `outptr`。
     * @param offset 供该函数读取或更新的 `offset` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    static void WriteCommandLineJSON(FILE* fp,
        char* const buffer, const int szbuffer, char*& outptr, int& offset);
    /**
     * @brief 在对齐结果输出中写出 `WriteSearchInformationPlain` 对应的数据。
     * @param fp 供该函数读取或更新的 `fp` 参数。
     * @param buffer 供当前步骤读取或更新的 `buffer` 缓冲区。
     * @param szbuffer 供当前步骤读取或更新的 `szbuffer` 缓冲区。
     * @param outptr 接收当前步骤输出的 `outptr`。
     * @param offset 供该函数读取或更新的 `offset` 参数。
     * @param tmpbuf 供当前步骤读取或更新的 `tmpbuf` 缓冲区。
     * @param sztmpbuf 供当前步骤读取或更新的 `sztmpbuf` 缓冲区。
     * @param rfilelist 供该函数读取或更新的 `rfilelist` 参数。
     * @param npossearched 控制当前步骤范围或规模的 `npossearched`。
     * @param nentries 控制当前步骤范围或规模的 `nentries`。
     * @param tmsthrld 供该函数读取或更新的 `tmsthrld` 参数。
     * @param indent 供该函数读取或更新的 `indent` 参数。
     * @param found 供该函数读取或更新的 `found` 参数。
     * @param clustering 供该函数读取或更新的 `clustering` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    static void WriteSearchInformationPlain(FILE* fp,
        char* const buffer, const int szbuffer, char*& outptr, int& offset,
        char* tmpbuf, int sztmpbuf, 
        const std::vector<std::string>& rfilelist,
        const size_t npossearched, const size_t nentries,
        const float tmsthrld, const int indent, const bool found,
        const bool clustering = false);

    /**
     * @brief 在对齐结果输出中处理 `BufferData` 对应的数据。
     * @param fp 供该函数读取或更新的 `fp` 参数。
     * @param buffer 供当前步骤读取或更新的 `buffer` 缓冲区。
     * @param szbuffer 供当前步骤读取或更新的 `szbuffer` 缓冲区。
     * @param outptr 接收当前步骤输出的 `outptr`。
     * @param offset 供该函数读取或更新的 `offset` 参数。
     * @param data 供该函数读取或更新的 `data` 参数。
     * @param szdata 供该函数读取或更新的 `szdata` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    static void BufferData( 
        FILE* fp, 
        char* const buffer, const int szbuffer, char*& outptr, int& offset, 
        const char* data, int szdata );
    /**
     * @brief 在对齐结果输出中写出 `WriteToFile` 对应的数据。
     * @param fp 供该函数读取或更新的 `fp` 参数。
     * @param data 供该函数读取或更新的 `data` 参数。
     * @param szdata 供该函数读取或更新的 `szdata` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    static void WriteToFile( FILE* fp, char* data, int szdata );

protected:
    /**
     * @brief 在对齐结果输出中处理 `Execute` 对应的数据。
     * @param args 供该函数读取或更新的 `args` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void Execute( void* args );

    /**
     * @brief 在对齐结果输出中设置 `SetResponseError` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void SetResponseError() {
        std::lock_guard<std::mutex> lck(mx_dataccess_);
        rsp_msg_ = WRITERTHREAD_MSG_ERROR;
    }

    /**
     * @brief 在对齐结果输出中处理 `MergeResults` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void MergeResults();
    /**
     * @brief 在对齐结果输出中写出 `WriteResults` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void WriteResults();
    /**
     * @brief 在对齐结果输出中读取 `GetOutputFilename` 对应的数据。
     * @param outfilename 控制当前步骤范围或规模的 `outfilename`。
     * @param outdirname 接收当前步骤输出的 `outdirname`。
     * @param qrydesc 描述查询结构的 `qrydesc`。
     * @param qrynr 描述查询结构的 `qrynr`。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void GetOutputFilename( 
        std::string& outfilename,
        const char* outdirname,
        const std::string& qrydesc,
        const int qrynr);
    /**
     * @brief 在对齐结果输出中写出 `WriteProgname` 对应的数据。
     * @param outptr 接收当前步骤输出的 `outptr`。
     * @param maxsize 控制当前步骤范围或规模的 `maxsize`。
     * @param width 供该函数读取或更新的 `width` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    int WriteProgname( char*& outptr, int maxsize, const int width );
    /**
     * @brief 在对齐结果输出中写出 `WriteQueryDescription` 对应的数据。
     * @param outptr 接收当前步骤输出的 `outptr`。
     * @param maxsize 控制当前步骤范围或规模的 `maxsize`。
     * @param qrylen 控制当前步骤范围或规模的 `qrylen`。
     * @param desc 供该函数读取或更新的 `desc` 参数。
     * @param width 供该函数读取或更新的 `width` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    int WriteQueryDescription(char*& outptr, int maxsize,
        const int qrylen, const char* desc, const int width );

    /**
     * @brief 在对齐结果输出中写出 `WriteResultsPlain` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void WriteResultsPlain();
    /**
     * @brief 在对齐结果输出中写出 `WriteQueryDescriptionPlain` 对应的数据。
     * @param outptr 接收当前步骤输出的 `outptr`。
     * @param maxsize 控制当前步骤范围或规模的 `maxsize`。
     * @param qrylen 控制当前步骤范围或规模的 `qrylen`。
     * @param desc 供该函数读取或更新的 `desc` 参数。
     * @param width 供该函数读取或更新的 `width` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    int WriteQueryDescriptionPlain(char*& outptr, int maxsize,
        const int qrylen, const std::string& desc, const int width );
    /**
     * @brief 在对齐结果输出中写出 `WriteSummaryPlain` 对应的数据。
     * @param outptr 接收当前步骤输出的 `outptr`。
     * @param qrylen 控制当前步骤范围或规模的 `qrylen`。
     * @param npossearched 控制当前步骤范围或规模的 `npossearched`。
     * @param nentries 控制当前步骤范围或规模的 `nentries`。
     * @param nqystrs 当前批次中的查询结构数量。
     * @param duration 供该函数读取或更新的 `duration` 参数。
     * @param string 供该函数读取或更新的 `string` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    int WriteSummaryPlain(char*& outptr,
        const int qrylen, const size_t npossearched, const size_t nentries,
        const int nqystrs, const double duration, const std::string&);

    /**
     * @brief 在对齐结果输出中写出 `WriteResultsJSON` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void WriteResultsJSON();
    /**
     * @brief 在对齐结果输出中写出 `WritePrognameJSON` 对应的数据。
     * @param outptr 接收当前步骤输出的 `outptr`。
     * @param maxsize 控制当前步骤范围或规模的 `maxsize`。
     * @param width 供该函数读取或更新的 `width` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    int WritePrognameJSON( char*& outptr, int maxsize, const int width );
    /**
     * @brief 在对齐结果输出中写出 `WriteQueryDescriptionJSON` 对应的数据。
     * @param outptr 接收当前步骤输出的 `outptr`。
     * @param maxsize 控制当前步骤范围或规模的 `maxsize`。
     * @param qrylen 控制当前步骤范围或规模的 `qrylen`。
     * @param desc 供该函数读取或更新的 `desc` 参数。
     * @param width 供该函数读取或更新的 `width` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    int WriteQueryDescriptionJSON(char*& outptr, int maxsize,
        const int qrylen, const std::string& desc, const int width );
    /**
     * @brief 在对齐结果输出中写出 `WriteSearchInformationJSON` 对应的数据。
     * @param fp 供该函数读取或更新的 `fp` 参数。
     * @param buffer 供当前步骤读取或更新的 `buffer` 缓冲区。
     * @param szbuffer 供当前步骤读取或更新的 `szbuffer` 缓冲区。
     * @param outptr 接收当前步骤输出的 `outptr`。
     * @param offset 供该函数读取或更新的 `offset` 参数。
     * @param tmpbuf 供当前步骤读取或更新的 `tmpbuf` 缓冲区。
     * @param sztmpbuf 供当前步骤读取或更新的 `sztmpbuf` 缓冲区。
     * @param rfilelist 供该函数读取或更新的 `rfilelist` 参数。
     * @param npossearched 控制当前步骤范围或规模的 `npossearched`。
     * @param nentries 控制当前步骤范围或规模的 `nentries`。
     * @param tmsthld 供该函数读取或更新的 `tmsthld` 参数。
     * @param found 供该函数读取或更新的 `found` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void WriteSearchInformationJSON(FILE* fp,
        char* const buffer, const int szbuffer, char*& outptr, int& offset,
        char* tmpbuf, int /* sztmpbuf */,
        const std::vector<std::string>& rfilelist,
        const size_t npossearched, const size_t nentries,
        const float tmsthld, const bool found);
    /**
     * @brief 在对齐结果输出中写出 `WriteSummaryJSON` 对应的数据。
     * @param outptr 接收当前步骤输出的 `outptr`。
     * @param qrylen 控制当前步骤范围或规模的 `qrylen`。
     * @param npossearched 控制当前步骤范围或规模的 `npossearched`。
     * @param nentries 控制当前步骤范围或规模的 `nentries`。
     * @param nqystrs 当前批次中的查询结构数量。
     * @param duration 供该函数读取或更新的 `duration` 参数。
     * @param devname 供该函数读取或更新的 `devname` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    int WriteSummaryJSON(char*& outptr,
        const int qrylen, const size_t npossearched, const size_t /* nentries */,
        const int nqystrs, const double duration, std::string devname);

    /**
     * @brief 在对齐结果输出中读取 `GetTotalNumberOfRecords` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    int GetTotalNumberOfRecords() const
    {
        if((int)parts_qrs_.size() <= qrysernr_ || qrysernr_ < 0 )
            throw MYRUNTIME_ERROR(
            "TdAlnWriter::GetTotalNumberOfRecords: Invalid query serial number.");
        int ntot = 0;
        for(size_t i = 0; i < vec_srtindxs_[qrysernr_].size(); i++) {
            if(vec_srtindxs_[qrysernr_][i])
                ntot += (int)vec_srtindxs_[qrysernr_][i]->size();
        }
        return ntot;
    }

    /**
     * @brief 在对齐结果输出中初始化 `InitializeVectors` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void InitializeVectors();

    /**
     * @brief 在对齐结果输出中处理 `ResizeVectors` 对应的数据。
     * @param newsize 控制当前步骤范围或规模的 `newsize`。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void ResizeVectors(int newsize) {
        vec_duration_.resize(newsize, 0.0);
        vec_nposschd_.resize(newsize, 0);
        vec_nentries_.resize(newsize, 0);
        vec_nqystrs_.resize(newsize);
        vec_nqyposs_.resize(newsize);
        vec_qrydesc_.resize(newsize);
        vec_devname_.resize(newsize);
        vec_tmsthld_.resize(newsize);
        vec_annotations_.resize(newsize);
        vec_alignments_.resize(newsize);
        vec_srtindxs_.resize(newsize);
        vec_tmscores_.resize(newsize);
        vec_alnptrs_.resize(newsize);
        vec_annotptrs_.resize(newsize);
    }

    /**
     * @brief 在对齐结果输出中处理 `ReleaseAllocations` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void ReleaseAllocations() {
        if((int)parts_qrs_.size() <= qrysernr_ || qrysernr_ < 0 )
            throw MYRUNTIME_ERROR(
            "TdAlnWriter::ReleaseAllocations: Invalid query serial number.");
        vec_annotations_[qrysernr_].clear();
        vec_alignments_[qrysernr_].clear();
        vec_srtindxs_[qrysernr_].clear();
        vec_tmscores_[qrysernr_].clear();
        vec_alnptrs_[qrysernr_].clear();
        vec_annotptrs_[qrysernr_].clear();
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
    //{{variables global to all queries:
    const char* mstr_set_outdirname_;//output directory name
    std::vector<std::string> mstr_set_rfilelist_;//filelist of reference (target) structures
    std::vector<double> vec_duration_;//accumulated duration for a query batch
    std::vector<size_t> vec_nposschd_;//total target db size in positions
    std::vector<size_t> vec_nentries_;//number of database entries
    //}}
    //{{query and summary data:
    std::vector<int> vec_nqystrs_;//batch size: #queries
    std::vector<int> vec_nqyposs_;//query lengths
    std::vector<std::string> vec_qrydesc_;//query descriptions
    std::vector<std::string> vec_devname_;//device names
    std::vector<float> vec_tmsthld_;//TM-score thresholds
    //}}
    //{{vectors of formatted results for QUERIES:
    std::vector<std::vector< std::unique_ptr<char,WritersDataDestroyer> >> vec_annotations_;
    std::vector<std::vector< std::unique_ptr<char,WritersDataDestroyer> >> vec_alignments_;
    std::vector<std::vector< std::unique_ptr<std::vector<int>> >> vec_srtindxs_;//index vectors of sorted TM-scores
    std::vector<std::vector< std::unique_ptr<std::vector<float>> >> vec_tmscores_;//2D vector of TM-scores for queries
    std::vector<std::vector< std::unique_ptr<std::vector<char*>> >> vec_alnptrs_;//2D vector of alignments for queries
    std::vector<std::vector< std::unique_ptr<std::vector<char*>> >> vec_annotptrs_;//2D vector of annotations for queries
    //}}
    //{{sorted indices over all parts (vectors) of results for a query:
    std::vector<int> allsrtindxs_;//indices along all vectors
    std::vector<int> allsrtvecs_;//corresponding vector indices (part numbers)
    std::vector<int> finalsrtindxs_;//globally (over all parts) sorted indices 
    std::vector<int> finalsrtindxs_dup_;//duplicate of globally sorted indices (for efficient memory management)
    std::vector<int>* p_finalindxs_;//pointer to the final vector of sorted indices
    //}}
    //buffer for writing to file:
    char buffer_[szWriterBuffer];
    //vector of the number of parts to be processed for each query serial number:
    std::vector<int> parts_qrs_;
    int qrysernr_;//query serial number
};

////////////////////////////////////////////////////////////////////////////
// TdAlnWriter INLINES
//
/**
 * @brief 在对齐结果输出中写出 `TdAlnWriter::WriteResults` 对应的数据。
 * @par 参数
 * 无。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
inline
void TdAlnWriter::WriteResults()
{
    static const int outfmt = CLOptions::GetO_OUTFMT();

    if(outfmt == CLOptions::oofJSON) {
        WriteResultsJSON();
        return;
    }

    WriteResultsPlain();
    return;
}

// -------------------------------------------------------------------------
//
/**
 * @brief 在对齐结果输出中写出 `TdAlnWriter::WriteProgname` 对应的数据。
 * @param outptr 接收当前步骤输出的 `outptr`。
 * @param maxsize 控制当前步骤范围或规模的 `maxsize`。
 * @param width 供该函数读取或更新的 `width` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline
int TdAlnWriter::WriteProgname( char*& outptr, int maxsize, const int width )
{
    static const int outfmt = CLOptions::GetO_OUTFMT();

    if(outfmt == CLOptions::oofJSON)
        return WritePrognameJSON(outptr, maxsize, width);

    return WritePrognamePlain(outptr, maxsize, width);
}

// -------------------------------------------------------------------------
//
/**
 * @brief 在对齐结果输出中写出 `TdAlnWriter::WriteQueryDescription` 对应的数据。
 * @param outptr 接收当前步骤输出的 `outptr`。
 * @param maxsize 控制当前步骤范围或规模的 `maxsize`。
 * @param qrylen 控制当前步骤范围或规模的 `qrylen`。
 * @param desc 供该函数读取或更新的 `desc` 参数。
 * @param width 供该函数读取或更新的 `width` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline
int TdAlnWriter::WriteQueryDescription(
    char*& outptr, int maxsize,
    const int qrylen, const char* desc, const int width)
{
    static const int outfmt = CLOptions::GetO_OUTFMT();

    if(outfmt == CLOptions::oofJSON)
        return WriteQueryDescriptionJSON(outptr, maxsize, qrylen, desc, width);

    return 
        WriteQueryDescriptionPlain(
            outptr, maxsize, qrylen, desc, width);
}

// -------------------------------------------------------------------------
#endif//__TdAlnWriter_h__
