/***************************************************************************
 *   Copyright (C) 2021-2026 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __PMBatchStrData_h__
#define __PMBatchStrData_h__

#include "libutil/mybase.h"

#include <string.h>

#include <vector>
#include <memory>
#include <algorithm>

#include "tsafety/TSCounterVar.h"
#include "libutil/CLOptions.h"
#include "PM2DVectorFields.h"
#include "DRDataDeleter.h"
#include "gdconst.h"

class PMBatchStrDataIndex;

// -------------------------------------------------------------------------
//
enum PMBatchStrDataCnsts {
    PMBSdatalignment = CACHECLINESIZE,
    PMBSdatDEFDESCLEN = DEFAULT_DESCRIPTION_LENGTH,//default description length (>4)
    //max structure length, i.e. allowed number of support atoms:
    //NOTE: encoding of structure positions is currently limited to one word (2 bytes):
    //NOTE: assign max length to 2^16-2 when recording DP cell coordinates: value of 
    //NOTE: 2^16-1 is left for the stop marker;
    PMBSmaxonestructurelength = 65535 //100 * ONEK
};

enum PMBatchStrDataMolType {
    PMBSMTProtein = -1,//Protein
    PMBSMTNA = 0//NA
};

// -------------------------------------------------------------------------
// PMBatchStrData: Complete batch structure data for parallel processing
//
class PMBatchStrData {
public:
    enum TPMBSDFinRetCode {
        pmbsdfEmpty=1,//no structure data (0 for compatibility)
        pmbsdfShort,//structure too short
        pmbsdfAbandoned,//structure is to large to be kept
        pmbsdfLowSimilarity,//low sequence similarity of structures
        pmbsdfLimits,//structure in the temporary buffers due to space limits
        pmbsdfWritten//structure has been written in the buffers
    };
    enum {
        PMBatchStrData_TMPBUFFSIZE = 4096
    };
public:
    /**
     * @brief 构造 `PMBatchStrData`，初始化其负责的结构数据读取与布局状态。
     * @par 参数
     * 无。
     * @return 无返回值；完成对象构造与初始状态设置。
     */
    PMBatchStrData()
    :   tmpbdbCdata_(nullptr),
        bdbCdata_(nullptr),
        bdbCdescs_(nullptr),
        bdbCptrdescs_(nullptr),
        szbdbCdata_(0),
        szbdbCdescs_(0),
        nbdbCptrdescs_(0),
        maxdatasize_(0),
        maxdatalen_(0),
        maxnstrs_(0)
    {
        memset( bdbCpmbeg_, 0, pmv2DTotFlds * sizeof(void*));
        memset( bdbCpmend_, 0, pmv2DTotFlds * sizeof(void*));
        memset( bdbCpmendovhd_, 0, pmv2DTotFlds * sizeof(void*));
//         memset( szpm2dvf_, 0, pmv2DTotFlds * sizeof(size_t));
//         memset( szpm2dvfovhd_, 0, pmv2DTotFlds * sizeof(size_t));
    }

    /**
     * @brief 销毁 `PMBatchStrData`，释放其管理的结构数据读取与布局资源。
     * @par 参数
     * 无。
     * @return 无返回值；对象持有的资源在返回前完成释放。
     */
    ~PMBatchStrData() {}

    // *** METHODS ***
    /**
     * @brief 在结构数据读取与布局中分配 `AllocateSpace` 对应的数据。
     * @param chunkdatasize 控制当前步骤范围或规模的 `chunkdatasize`。
     * @param chunkdatalen 控制当前步骤范围或规模的 `chunkdatalen`。
     * @param chunknstrs 供该函数读取或更新的 `chunknstrs` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void AllocateSpace(size_t chunkdatasize, size_t chunkdatalen, size_t chunknstrs);

    //serialize/deserialize
    /**
     * @brief 在结构数据读取与布局中处理 `Serialize` 对应的数据。
     * @param string 供该函数读取或更新的 `string` 参数。
     * @param nagents 控制当前步骤范围或规模的 `nagents`。
     * @param eod 供该函数读取或更新的 `eod` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void Serialize(const std::string&, const size_t nagents, const size_t eod) const;
    /**
     * @brief 在结构数据读取与布局中处理 `Deserialize` 对应的数据。
     * @param string 供该函数读取或更新的 `string` 参数。
     * @param nagents 控制当前步骤范围或规模的 `nagents`。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    size_t Deserialize(const std::string&, const size_t nagents);

    /**
     * @brief 在结构数据读取与布局中处理 `Fallback` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void Fallback() {FallbackOvhdPtrs();}

    //sort by length; used when all chunk data has been compiled
    /**
     * @brief 在结构数据读取与布局中排序 `Sort` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void Sort();

    /**
     * @brief 在结构数据读取与布局中复制 `Copy` 对应的数据。
     * @param bdbdesc 描述参考结构的 `bdbdesc`。
     * @param bdbpmbeg 描述参考结构的 `bdbpmbeg`。
     * @param bdbpmend 描述参考结构的 `bdbpmend`。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void Copy(const char** bdbdesc, char* const * const bdbpmbeg, char* const * const bdbpmend);
    /**
     * @brief 在结构数据读取与布局中复制 `CopyFrom` 对应的数据。
     * @param bsd 供该函数读取或更新的 `bsd` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void CopyFrom(const PMBatchStrData& bsd);

    /**
     * @brief 在结构数据读取与布局中读取 `GetPMDataSize` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    size_t GetPMDataSize();
    /**
     * @brief 在结构数据读取与布局中读取 `GetPMDataSize1` 对应的数据。
     * @param structlen 控制当前步骤范围或规模的 `structlen`。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    static size_t GetPMDataSize1(size_t structlen);
    /**
     * @brief 在结构数据读取与布局中读取 `GetPMDataSizeUB` 对应的数据。
     * @param totallen 控制当前步骤范围或规模的 `totallen`。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    static size_t GetPMDataSizeUB(size_t totallen);

    /**
     * @brief 在结构数据读取与布局中处理 `ContainsData` 对应的数据。
     * @param bdbpmbeg 描述参考结构的 `bdbpmbeg`。
     * @param bdbpmend 描述参考结构的 `bdbpmend`。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    static bool ContainsData(char* const * const bdbpmbeg, char* const * const bdbpmend) {
        return 
            bdbpmbeg[pps2DLen] < bdbpmend[pps2DLen] &&
            bdbpmbeg[pps2DDist] < bdbpmend[pps2DDist] &&
            bdbpmbeg[pps2DType] < bdbpmend[pps2DType] &&
            bdbpmbeg[pmv2DCoords] < bdbpmend[pmv2DCoords];
    }

    /**
     * @brief 在结构数据读取与布局中处理 `ContainsData` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    bool ContainsData() const {//whether data is present
        return bdbCpmbeg_[pmv2DCoords] < bdbCpmend_[pmv2DCoords];
    }

    /**
     * @brief 在结构数据读取与布局中处理 `ContainsDataLast` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    bool ContainsDataLast() const {//whether data has been written for the last structure
        return bdbCpmend_[pmv2DCoords] < bdbCpmendovhd_[pmv2DCoords];
    }

    //change/pack structure pointers to include structures indexed in filterdata
    /**
     * @brief 在结构数据读取与布局中筛选 `FilterStructs` 对应的数据。
     * @param bdbdesc 描述参考结构的 `bdbdesc`。
     * @param bdbpmbeg 描述参考结构的 `bdbpmbeg`。
     * @param bdbpmend 描述参考结构的 `bdbpmend`。
     * @param filterdata 供该函数读取或更新的 `filterdata` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    static void FilterStructs(
        const char** bdbdesc, char** bdbpmbeg, char** bdbpmend,
        const unsigned int* filterdata);

    //length of the structure at the given position/index (NOTE:pointers assumed valid):
    /**
     * @brief 在结构数据读取与布局中读取 `GetLengthAt` 对应的数据。
     * @param bdbpmbeg 描述参考结构的 `bdbpmbeg`。
     * @param ndx 控制当前步骤范围或规模的 `ndx`。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    static size_t GetLengthAt(const char* const * const bdbpmbeg, int ndx) {
        return ((INTYPE*)(bdbpmbeg[pps2DLen]))[ndx];
    }

    //address of the structure at the given position/index (NOTE:pointers assumed valid):
    /**
     * @brief 在结构数据读取与布局中读取 `GetAddressAt` 对应的数据。
     * @param bdbpmbeg 描述参考结构的 `bdbpmbeg`。
     * @param ndx 控制当前步骤范围或规模的 `ndx`。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    static size_t GetAddressAt(const char* const * const bdbpmbeg, int ndx) {
        return ((LNTYPE*)(bdbpmbeg[pps2DDist]))[ndx];
    }

    //field of the structure at the given position/index (NOTE:pointers assumed valid):
    template<typename T, int F>
    /**
     * @brief 在结构数据读取与布局中读取 `GetFieldAt` 对应的数据。
     * @param bdbpmbeg 描述参考结构的 `bdbpmbeg`。
     * @param ndx 控制当前步骤范围或规模的 `ndx`。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    static T GetFieldAt(const char* const * const bdbpmbeg, int ndx) {
        return ((T*)(bdbpmbeg[F]))[ndx];
    }

    template<typename T, int F>
    /**
     * @brief 在结构数据读取与布局中读取 `GetFieldAt` 对应的数据。
     * @param ndx 控制当前步骤范围或规模的 `ndx`。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    T GetFieldAt(int ndx) const {return ((T*)(bdbCpmbeg_[F]))[ndx];}

    //set a structure field at the given position/index (NOTE:pointers assumed valid):
    template<typename T, int F>
    /**
     * @brief 在结构数据读取与布局中设置 `SetFieldAt` 对应的数据。
     * @param bdbpmbeg 描述参考结构的 `bdbpmbeg`。
     * @param ndx 控制当前步骤范围或规模的 `ndx`。
     * @param value 需要读取、写入或转换的值。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    static void SetFieldAt(char* const * const bdbpmbeg, int ndx, T value) {
        ((T*)(bdbpmbeg[F]))[ndx] = value;
    }

    template<typename T, int F>
    /**
     * @brief 在结构数据读取与布局中设置 `SetFieldAt` 对应的数据。
     * @param ndx 控制当前步骤范围或规模的 `ndx`。
     * @param value 需要读取、写入或转换的值。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void SetFieldAt(int ndx, T value) {((T*)(bdbCpmbeg_[F]))[ndx] = value;}

    //#structures written in the buffers (NOTE:pointers assumed valid):
    /**
     * @brief 在结构数据读取与布局中读取 `GetNoStructs` 对应的数据。
     * @param bdbpmbeg 描述参考结构的 `bdbpmbeg`。
     * @param bdbpmend 描述参考结构的 `bdbpmend`。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    static size_t GetNoStructs(char* const * const bdbpmbeg, char* const * const bdbpmend) {
        return (size_t)(bdbpmend[pps2DLen]-bdbpmbeg[pps2DLen]) / SZINTYPE;
    }

    /**
     * @brief 在结构数据读取与布局中读取 `GetNoStructsWritten` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    size_t GetNoStructsWritten() const {//#structures written in bdbCdata_
        return (size_t)(bdbCpmend_[pps2DLen]-bdbCpmbeg_[pps2DLen]) / SZINTYPE;
    }

    //#total positions written in the buffers (NOTE:pointers assumed valid):
    /**
     * @brief 在结构数据读取与布局中读取 `GetNoPosits` 对应的数据。
     * @param bdbpmbeg 描述参考结构的 `bdbpmbeg`。
     * @param bdbpmend 描述参考结构的 `bdbpmend`。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    static size_t GetNoPosits(char* const * const bdbpmbeg, char* const * const bdbpmend) {
        return (size_t)(bdbpmend[pmv2Drsd]-bdbpmbeg[pmv2Drsd]) / SZCHTYPE;
    }

    /**
     * @brief 在结构数据读取与布局中读取 `GetNoPositsWritten` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    size_t GetNoPositsWritten() const {//#positions written in bdbCdata_
        return (size_t)(bdbCpmend_[pmv2Drsd]-bdbCpmbeg_[pmv2Drsd]) / SZCHTYPE;
    }

    /**
     * @brief 在结构数据读取与布局中读取 `GetNoPositsOvhd` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    size_t GetNoPositsOvhd() const {//#positions in the structure being compiled
        return (size_t)(bdbCpmendovhd_[pmv2Drsd]-bdbCpmend_[pmv2Drsd]) / SZCHTYPE;
    }

    //{{these methods are for filling by residue one structure at a time
    /**
     * @brief 在结构数据读取与布局中处理 `AddOneResidue` 对应的数据。
     * @param maxstrlen 控制当前步骤范围或规模的 `maxstrlen`。
     * @param rsdcode 供该函数读取或更新的 `rsdcode` 参数。
     * @param resnum 供该函数读取或更新的 `resnum` 参数。
     * @param restype 供该函数读取或更新的 `restype` 参数。
     * @param inscode 供该函数读取或更新的 `inscode` 参数。
     * @param chain 供该函数读取或更新的 `chain` 参数。
     * @param chord 供该函数读取或更新的 `chord` 参数。
     * @param coords 供该函数读取或更新的 `coords` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    bool AddOneResidue(
        int maxstrlen,
        CHTYPE rsdcode, INTYPE resnum, int restype, char inscode, char chain, char chord, 
        FPTYPE coords[pmv2DNoElems]);

    /**
     * @brief 在结构数据读取与布局中处理 `FinalizeCrntStructure` 对应的数据。
     * @param filendx 控制当前步骤范围或规模的 `filendx`。
     * @param strndx 供该函数读取或更新的 `strndx` 参数。
     * @param globndx 供该函数读取或更新的 `globndx` 参数。
     * @param moltype 供该函数读取或更新的 `moltype` 参数。
     * @param description 供该函数读取或更新的 `description` 参数。
     * @param strchain 供该函数读取或更新的 `strchain` 参数。
     * @param strmodel 供该函数读取或更新的 `strmodel` 参数。
     * @param usechaininfo 供该函数读取或更新的 `usechaininfo` 参数。
     * @param usemodelinfo 供该函数读取或更新的 `usemodelinfo` 参数。
     * @param queryblocks 描述查询结构的 `queryblocks`。
     * @param querypmbegs 描述查询结构的 `querypmbegs`。
     * @param querypmends 描述查询结构的 `querypmends`。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    TPMBSDFinRetCode FinalizeCrntStructure(//finalize the structure being compiled
        size_t filendx, int strndx,
        int globndx, int moltype, const std::string& description, 
        const std::string& strchain, const std::string& strmodel, 
        bool usechaininfo, bool usemodelinfo,
        const int queryblocks,
        char* const * const * const querypmbegs,
        char* const * const * const querypmends);
    //}}

    /**
     * @brief 在结构数据读取与布局中处理 `FieldTypeValid` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    bool FieldTypeValid() const;//whether the field Type is valid across structures

    /**
     * @brief 在结构数据读取与布局中读取 `GetFileNdxAt` 对应的数据。
     * @param ndx 控制当前步骤范围或规模的 `ndx`。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    size_t GetFileNdxAt(int ndx) const {return tmpclustfilendxs_[ndx];}
    /**
     * @brief 在结构数据读取与布局中读取 `GetStructNdxAt` 对应的数据。
     * @param ndx 控制当前步骤范围或规模的 `ndx`。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    int GetStructNdxAt(int ndx) const {return tmpcluststrndxs_[ndx];}

    //get the description of the structure in the overhead buffer
    /**
     * @brief 在结构数据读取与布局中读取 `GetOvhdStrDescription` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    const char* GetOvhdStrDescription() const;

    //copy overhead to another batch object:
    /**
     * @brief 在结构数据读取与布局中复制 `CopyOvhdTo` 对应的数据。
     * @param PMBatchStrData 供该函数读取或更新的 `PMBatchStrData` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    TPMBSDFinRetCode CopyOvhdTo(PMBatchStrData&) const;

    //return true if the sequence similairty between the structure in Ovhd and 
    //any of the queries in given blocks is above the threshold
    /**
     * @brief 在结构数据读取与布局中处理 `SequenceSimilarityOvhd` 对应的数据。
     * @param queryblocks 描述查询结构的 `queryblocks`。
     * @param querypmbegs 描述查询结构的 `querypmbegs`。
     * @param querypmends 描述查询结构的 `querypmends`。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    bool SequenceSimilarityOvhd(
        const int queryblocks,
        char* const * const * const querypmbegs,
        char* const * const * const querypmends) const;

    template<int DIMD>
    /**
     * @brief 在结构数据读取与布局中检查 `CheckAlignmentScore` 对应的数据。
     * @param seqsimthrscore 当前步骤使用或写回的 `seqsimthrscore` 分数。
     * @param qrydst 描述查询结构的 `qrydst`。
     * @param qrylen 控制当前步骤范围或规模的 `qrylen`。
     * @param dbstrdst 描述参考结构的 `dbstrdst`。
     * @param dbstrlen 控制当前步骤范围或规模的 `dbstrlen`。
     * @param qrybegpos 描述查询结构的 `qrybegpos`。
     * @param rfnbegpos 描述参考结构的 `rfnbegpos`。
     * @param querypmbeg 查询结构打包字段的起始指针数组。
     * @param bdbCpmbeg 参考结构打包字段的起始指针数组。
     * @param rfnRE 描述参考结构的 `rfnRE`。
     * @param qryRE 描述查询结构的 `qryRE`。
     * @param scores 保存或读取对齐分数的缓冲区。
     * @param pxmins 供该函数读取或更新的 `pxmins` 参数。
     * @param tmp 供该函数读取或更新的 `tmp` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    static bool CheckAlignmentScore(
        const float seqsimthrscore,
        const int qrydst, const int qrylen,
        const int dbstrdst, const int dbstrlen,
        const int qrybegpos, const int rfnbegpos,
        const char* const * const __restrict querypmbeg,
        const char* const * const __restrict bdbCpmbeg,
        char rfnRE[DIMD], char qryRE[DIMD],
        float scores[DIMD], float pxmins[DIMD],
        float tmp[DIMD]);

    /**
     * @brief 在结构数据读取与布局中格式化输出 `Print` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void Print() const;
    /**
     * @brief 在结构数据读取与布局中格式化输出 `Print` 对应的数据。
     * @param bdbpmbeg 描述参考结构的 `bdbpmbeg`。
     * @param bdbpmend 描述参考结构的 `bdbpmend`。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    static void Print(char* const * const bdbpmbeg, char* const * const bdbpmend);

private:
    /**
     * @brief 在结构数据读取与布局中设置 `SetMaxDataLimits` 对应的数据。
     * @param maxdatasize 控制当前步骤范围或规模的 `maxdatasize`。
     * @param maxdatalen 控制当前步骤范围或规模的 `maxdatalen`。
     * @param maxnstrs 供该函数读取或更新的 `maxnstrs` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void SetMaxDataLimits(
        size_t maxdatasize, size_t maxdatalen, size_t maxnstrs)
    {
        maxdatasize_ = maxdatasize;
        maxdatalen_ = maxdatalen;
        maxnstrs_ = maxnstrs;
    }

    /**
     * @brief 在结构数据读取与布局中分配 `AllocateSpaceForData` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void AllocateSpaceForData();
    /**
     * @brief 在结构数据读取与布局中分配 `AllocateSpaceForDescriptions` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    bool AllocateSpaceForDescriptions();
    /**
     * @brief 在结构数据读取与布局中分配 `AllocateSpaceForDescPtrs` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void AllocateSpaceForDescPtrs();

    /**
     * @brief 在结构数据读取与布局中处理 `FallbackOvhdPtrs` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void FallbackOvhdPtrs() {//fall back overhead pointers
        for(int f = 0; f < pmv2DTotFlds; f++)
            bdbCpmendovhd_[f] = bdbCpmend_[f];
    }

    /**
     * @brief 在结构数据读取与布局中处理 `FormatDescription` 对应的数据。
     * @param description 供该函数读取或更新的 `description` 参数。
     * @param strchain 供该函数读取或更新的 `strchain` 参数。
     * @param strmodel 供该函数读取或更新的 `strmodel` 参数。
     * @param usechaininfo 供该函数读取或更新的 `usechaininfo` 参数。
     * @param usemodelinfo 供该函数读取或更新的 `usemodelinfo` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    std::string FormatDescription(
        std::string description,
        const std::string& strchain, const std::string& strmodel, 
        bool usechaininfo, bool usemodelinfo);

private:
    template<typename T, int field>
    /**
     * @brief 在结构数据读取与布局中排序 `sort_helper_ssfields_assign` 对应的数据。
     * @param nstts 控制当前步骤范围或规模的 `nstts`。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void sort_helper_ssfields_assign(size_t nstts);
    template<int field>
    /**
     * @brief 在结构数据读取与布局中排序 `sort_helper_psfields_assign` 对应的数据。
     * @param nstts 控制当前步骤范围或规模的 `nstts`。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void sort_helper_psfields_assign(size_t nstts);

private:
    std::unique_ptr<char,DRDataDeleter> tmpbdbCdata_;//temporary buffer for sorted data
    std::vector<INTYPE> tmpndxs_;//temporary index vector for structure lengths
    char tmpreadbuff_[PMBatchStrData_TMPBUFFSIZE];//tmp buffer for data read

public:
    std::vector<size_t> tmpclustfilendxs_;//temporary file index vector for clustering
    std::vector<int> tmpcluststrndxs_;//temporary structure-within-file index vector for clustering

public:
    // *** MEMBER VARIABLES ***
    std::unique_ptr<char,DRHostDataDeleter> bdbCdata_;
    char* bdbCpmbeg_[pmv2DTotFlds];//addresses of the beginnings of the fields in bdbCdata_
    char* bdbCpmend_[pmv2DTotFlds];//end addresses of the fields in bdbCdata_
    //end addresses of the fields in bdbCdata_
    // including one additional structure if any that does not fit into the memory limits:
    char* bdbCpmendovhd_[pmv2DTotFlds];
    //TODO: introduce bdbCpmendovhdbeg_ instead of using bdbCpmend_ for CPU version (ovhd copy)!

//     //sizes in bytes of the fields of complete structures written in bdbCdata_:
//     size_t szpm2dvf_[pmv2DTotFlds];
//     //sizes in bytes of the fields written in bdbCdata_, 
//     // including one additional structure if any that does not fit into the memory limits:
//     size_t szpm2dvfovhd_[pmv2DTotFlds];

    std::unique_ptr<char,DRDataDeleter> bdbCdescs_;//descriptions
    std::unique_ptr<char*[]> bdbCptrdescs_;//pointers to structure descriptions in bdbCdescs_

    size_t szbdbCdata_;//size allocated for bdbCdata_
    size_t szbdbCdescs_;//size allocated for bdbCdescs_
    size_t nbdbCptrdescs_;//number of slots allocated for bdbCptrdescs_

//     size_t szcrntsize_;//size currently used by bdbCdata_
//     size_t szcrntlen_;//total #positions (length) currently occupied by bdbCdata_
//     size_t szcrntstrs_;//#structures currently in bdbCdata_

    size_t maxdatasize_;//max size of data for bdbCdata_ to contain
    size_t maxdatalen_;//max total #positions (length) over all structures
    size_t maxnstrs_;//max #structures whose data bdbCdata_ can contain

    TSCounterVar cnt_;//thread-safe counter of how many agents access the data
};

// =========================================================================
// INLINES
// 
// AllocateSpace: allocate space for structure data and descriptions given the 
// limits of data chunk size, total number of positions (residues), and the 
// number of structures
//
/**
 * @brief 在结构数据读取与布局中分配 `PMBatchStrData::AllocateSpace` 对应的数据。
 * @param chunkdatasize 控制当前步骤范围或规模的 `chunkdatasize`。
 * @param chunkdatalen 控制当前步骤范围或规模的 `chunkdatalen`。
 * @param chunknstrs 供该函数读取或更新的 `chunknstrs` 参数。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
inline
void PMBatchStrData::AllocateSpace(
    size_t chunkdatasize, size_t chunkdatalen, size_t chunknstrs)
{
    //the following comes first:
    SetMaxDataLimits(chunkdatasize, chunkdatalen, chunknstrs);

    //allocate space for structure data plus space for the end addresses of 
    // structure descriptions once the limits have been set
    AllocateSpaceForData();

    AllocateSpaceForDescriptions();
    AllocateSpaceForDescPtrs();
}

// -------------------------------------------------------------------------
// GetPMDataSize: get the size of the structure model data written
/**
 * @brief 在结构数据读取与布局中读取 `PMBatchStrData::GetPMDataSize` 对应的数据。
 * @par 参数
 * 无。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline
size_t PMBatchStrData::GetPMDataSize()
{
    MYMSG("PMBatchStrData::GetPMDataSize",6);
    size_t size = 0;
    for(int n = 0; n < pmv2DTotFlds; n++)
        size += (size_t)(bdbCpmend_[n]-bdbCpmbeg_[n]);
    return size;
}

// -------------------------------------------------------------------------
// GetPMDataSize1: get the size of complete structure model data of one 
// structure 
/**
 * @brief 在结构数据读取与布局中读取 `PMBatchStrData::GetPMDataSize1` 对应的数据。
 * @param structlen 控制当前步骤范围或规模的 `structlen`。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline
size_t PMBatchStrData::GetPMDataSize1(size_t structlen)
{
    MYMSG("PMBatchStrData::GetPMDataSize1",6);
    size_t size = 0;
    if(structlen < 1) return size;
    int n = 0;
    for(; n < pps2DStrFlds; n++)
        size += TPM2DVectorFieldSize::szvfs_[n];
    for(; n < pmv2DTotFlds; n++)
        size += TPM2DVectorFieldSize::szvfs_[n] * structlen;
    return size;
}

// -------------------------------------------------------------------------
// GetPMDataSizeUB: get the max size of complete structure model data when 
// the total number of positions (residues) is totallen
/**
 * @brief 在结构数据读取与布局中读取 `PMBatchStrData::GetPMDataSizeUB` 对应的数据。
 * @param totallen 控制当前步骤范围或规模的 `totallen`。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline
size_t PMBatchStrData::GetPMDataSizeUB(size_t totallen)
{
    MYMSG("PMBatchStrData::GetPMDataSizeUB",6);
    size_t size = 0;
    for(int n = 0; n < pmv2DTotFlds; n++)
        //[pps2DLen]: when #structs=totallen
        //[pps2DDist]: when #structs=totallen
        size += TPM2DVectorFieldSize::szvfs_[n] * totallen;
    return size;
}

// -------------------------------------------------------------------------
// AddOneResidue: save data for one residue in the batch data object; 
// adjust accordingly pointers and running sizes;
// maxstrlen, max structure length (-1, use default max);
// rsdcode, one-letter residue code;
// resnum, residue serial number as appears in the structure file;
// inscode, residue insertion code;
// chain, chain id of the structure;
// chord, chain serial number (order) in the structure file;
// coords, residue coordinates
/**
 * @brief 在结构数据读取与布局中处理 `PMBatchStrData::AddOneResidue` 对应的数据。
 * @param maxstrlen 控制当前步骤范围或规模的 `maxstrlen`。
 * @param rsdcode 供该函数读取或更新的 `rsdcode` 参数。
 * @param resnum 供该函数读取或更新的 `resnum` 参数。
 * @param restype 供该函数读取或更新的 `restype` 参数。
 * @param inscode 供该函数读取或更新的 `inscode` 参数。
 * @param chain 供该函数读取或更新的 `chain` 参数。
 * @param chord 供该函数读取或更新的 `chord` 参数。
 * @param coords 供该函数读取或更新的 `coords` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline
bool PMBatchStrData::AddOneResidue(
    int maxstrlen,
    CHTYPE rsdcode, INTYPE resnum, int restype, char inscode, char chain, char chord,
    FPTYPE coords[pmv2DNoElems])
{
    //NOTE: structure-specific fields pps2DLen, pps2DType, and pps2DDist 
    // have to be updated on finalizing the data of the current structure

    //number of positions in this structure (not finalized) so far
    size_t nposits = GetNoPositsOvhd();

    if((0 < maxstrlen && (size_t)maxstrlen < nposits)/*max length violation*/||
       PMBSmaxonestructurelength < nposits/*buffer size exceeded*/||
       maxdatalen_ < nposits/*allocation too small*/)
    {
        //fallback, abandon the structure
        FallbackOvhdPtrs();
        return false;
    }

    //number of structures written so far:
    // INTYPE nstructs = GetNoStructsWritten();
    int n, f;

    for(n = 0, f = pmv2DCoords; n < pmv2DNoElems; n++, f++) {
        *(FPTYPE*)(bdbCpmendovhd_[f]) = coords[n];
        bdbCpmendovhd_[f] += TPM2DVectorFieldSize::szvfs_[f];
    }

    *(INTYPE*)(bdbCpmendovhd_[pmv2D_Ins_Ch_Ord]) = 
        (INTYPE)PM2D_MAKEINT_Ins_Ch_Ord(
            (unsigned int)inscode, (unsigned int)chain, (unsigned int)chord);
    WatermarkOnType(GetMoleculeType(restype), *(INTYPE*)(bdbCpmendovhd_[pmv2D_Ins_Ch_Ord]));
    bdbCpmendovhd_[pmv2D_Ins_Ch_Ord] += TPM2DVectorFieldSize::szvfs_[pmv2D_Ins_Ch_Ord];

    *(INTYPE*)(bdbCpmendovhd_[pmv2DResNumber]) = resnum;
    bdbCpmendovhd_[pmv2DResNumber] += TPM2DVectorFieldSize::szvfs_[pmv2DResNumber];

    *(CHTYPE*)(bdbCpmendovhd_[pmv2Drsd]) = rsdcode;
    bdbCpmendovhd_[pmv2Drsd] += TPM2DVectorFieldSize::szvfs_[pmv2Drsd];

    //NOTE: secondary structure assignments skipped; they will be 
    // calculated and filled in by an accelerator
    //*(CHTYPE*)(bdbCpmendovhd_[pmv2Dss]) = ...;
    bdbCpmendovhd_[pmv2Dss] += TPM2DVectorFieldSize::szvfs_[pmv2Dss];

    return true;
}

// -------------------------------------------------------------------------
// FinalizeCrntStructure: finalize the structure being compiled by 
// moving the data from the overhead buffer to the end of the data in the 
// chunk; moving data actually corresponds to changing the values of 
// pointers;
// return false (pmbsdfAbandoned) if the structure is too large to be 
// kept by the buffers associated with the chunk of data;
// filendx, file index for clustering;
// strndx, structure index within a file for clustering;
/**
 * @brief 在结构数据读取与布局中处理 `PMBatchStrData::FinalizeCrntStructure` 对应的数据。
 * @param filendx 控制当前步骤范围或规模的 `filendx`。
 * @param strndx 供该函数读取或更新的 `strndx` 参数。
 * @param globndx 供该函数读取或更新的 `globndx` 参数。
 * @param moltype 供该函数读取或更新的 `moltype` 参数。
 * @param description 供该函数读取或更新的 `description` 参数。
 * @param strchain 供该函数读取或更新的 `strchain` 参数。
 * @param strmodel 供该函数读取或更新的 `strmodel` 参数。
 * @param usechaininfo 供该函数读取或更新的 `usechaininfo` 参数。
 * @param usemodelinfo 供该函数读取或更新的 `usemodelinfo` 参数。
 * @param queryblocks 描述查询结构的 `queryblocks`。
 * @param querypmbegs 描述查询结构的 `querypmbegs`。
 * @param querypmends 描述查询结构的 `querypmends`。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline
PMBatchStrData::TPMBSDFinRetCode 
PMBatchStrData::FinalizeCrntStructure(
    size_t filendx, int strndx,
    int globndx, int moltype, const std::string& description, 
    const std::string& strchain, const std::string& strmodel, 
    bool usechaininfo, bool usemodelinfo,
    const int /* queryblocks */,
    char* const * const * const /* querypmbegs */,
    char* const * const * const /* querypmends */)
{
    MYMSG("PMBatchStrData::FinalizeCrntStructure",6);
    static const std::string preamb = "PMBatchStrData::FinalizeCrntStructure: ";
    static const size_t devminrlength = CLOptions::GetDEV_MINRLEN();
    //number of positions in this structure (not finalized)
    size_t nposits = GetNoPositsOvhd();//#positions to be written
    size_t npositswrt = GetNoPositsWritten();//#positions written
    size_t szstruct = PMBatchStrData::GetPMDataSize1(nposits);//size of data to be written
    size_t szstructswrt = GetPMDataSize();//size of data written
    size_t nstructswrt = GetNoStructsWritten();//number of structures written
    static constexpr size_t szalign = PMBSdatalignment * pmv2DTotFlds;//max size for data alignment

    if(nposits < 1) {
        FallbackOvhdPtrs();
        return pmbsdfEmpty;
    }
    if(nposits < devminrlength) {
        FallbackOvhdPtrs();
        return pmbsdfShort;
    }

    if(maxdatasize_ < szstruct + szalign || 
       maxdatalen_ < nposits)//allocation too small
    {
        //fallback, abandon the structure
        FallbackOvhdPtrs();
        return pmbsdfAbandoned;
    }

    // //verify mutual sequence similarity
    // //NOTE: verified at the computation level.
    // if(0 < queryblocks && querypmbegs && querypmends &&
    //    !SequenceSimilarityOvhd(queryblocks, querypmbegs, querypmends))
    // {
    //     //fallback, abandon the structure
    //     FallbackOvhdPtrs();
    //     return pmbsdfLowSimilarity;
    // }


    //{{ fill in accompanying data before returning with 
    // pmbsdfLimits (to be moved later) or pmbsdfWritten
    char** ptrdescs = bdbCptrdescs_.get();//assume a valid pointer

    std::string desc = 
        FormatDescription(description, strchain, strmodel, 
            usechaininfo, usemodelinfo);

    memcpy(ptrdescs[nstructswrt], desc.c_str(), desc.size());
    ptrdescs[nstructswrt][desc.size()] = 0;
    //}}


    //fill in structure-specific fields; use pps2DType as a type and ID simultaneously;
    *(INTYPE*)(bdbCpmendovhd_[pps2DLen]) = (INTYPE)nposits;
    *(INTYPE*)(bdbCpmendovhd_[pps2DType]) = (INTYPE)(globndx);
    *(LNTYPE*)(bdbCpmendovhd_[pps2DDist]) = (LNTYPE)npositswrt;
    //update the sign bit at the first structure position to indicate molecular type
    WatermarkOnType(GetMoleculeType(moltype), *(INTYPE*)(bdbCpmend_[pmv2D_Ins_Ch_Ord]));

    tmpclustfilendxs_.push_back(filendx);
    tmpcluststrndxs_.push_back(strndx);

    if(maxdatasize_ < szstructswrt + szstruct + szalign || 
       maxdatalen_ < npositswrt + nposits || 
       maxnstrs_ < nstructswrt + 1)
    {
        //this structure cannot be contained currently, leave it for the next round;
        //structure will be the 1st, distance = 0
        *(LNTYPE*)(bdbCpmendovhd_[pps2DDist]) = 0;
        return pmbsdfLimits;
    }

    if((size_t)INT_MAX <= npositswrt + nposits)
        throw MYRUNTIME_ERROR2( 
        preamb + "Overflow detected. Data chunk size must be reduced.", 
        CRITICAL);

    bdbCpmendovhd_[pps2DLen] += TPM2DVectorFieldSize::szvfs_[pps2DLen];
    bdbCpmendovhd_[pps2DType] += TPM2DVectorFieldSize::szvfs_[pps2DType];
    bdbCpmendovhd_[pps2DDist] += TPM2DVectorFieldSize::szvfs_[pps2DDist];

    //adjust pointers to the end of overhead data section
    for(int f = 0; f < pmv2DTotFlds; f++)
        bdbCpmend_[f] = bdbCpmendovhd_[f];

    return pmbsdfWritten;
}

// -------------------------------------------------------------------------
// FormatDescription: format structure description;
// return `move' string
/**
 * @brief 在结构数据读取与布局中处理 `PMBatchStrData::FormatDescription` 对应的数据。
 * @param description 供该函数读取或更新的 `description` 参数。
 * @param strchain 供该函数读取或更新的 `strchain` 参数。
 * @param strmodel 供该函数读取或更新的 `strmodel` 参数。
 * @param usechaininfo 供该函数读取或更新的 `usechaininfo` 参数。
 * @param usemodelinfo 供该函数读取或更新的 `usemodelinfo` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline
std::string PMBatchStrData::FormatDescription(
    std::string description,
    const std::string& strchain, const std::string& strmodel, 
    bool usechaininfo, bool usemodelinfo)
{
    if(usechaininfo && strchain.size() && strchain[0])
        description = description + " Chn:" + strchain;//chain
    if(usemodelinfo && strmodel.size() && strmodel[0])
        description = description + " (M:" + strmodel + ")";//model
    std::string desc = 
        (PMBSdatDEFDESCLEN <= description.size()) //4: "..." + 0
        ?   "..." + description.substr(description.size()-PMBSdatDEFDESCLEN+4)
        :   description;
    return desc;
}

// -------------------------------------------------------------------------
// FieldTypeValid: verify whether the Type field is valid across all 
// written structures;
// 
/**
 * @brief 在结构数据读取与布局中处理 `PMBatchStrData::FieldTypeValid` 对应的数据。
 * @par 参数
 * 无。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline
bool PMBatchStrData::FieldTypeValid() const
{
    int nstructswrt = (int)GetNoStructsWritten();
    for(int i = 0; i < nstructswrt; i++)
        if(GetFieldAt<INTYPE,pps2DType>(i) == INT_MAX)
            return false;
    return true;
}

// -------------------------------------------------------------------------
// CopyOvhdTo: move data (of one structure) from the overhead buffers to 
// another batch object;
// bsd, batch object to move data to;
// return false (pmbsdfAbandoned) if the structure is too large to be 
// kept (should not happen without a bug);
/**
 * @brief 在结构数据读取与布局中复制 `PMBatchStrData::CopyOvhdTo` 对应的数据。
 * @param bsd 供该函数读取或更新的 `bsd` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline
PMBatchStrData::TPMBSDFinRetCode
PMBatchStrData::CopyOvhdTo(PMBatchStrData& bsd) const
{
    MYMSG("PMBatchStrData::CopyOvhdTo",6);
    static const std::string preamb = "PMBatchStrData::CopyOvhdTo: ";
    //number of positions in this structure (not finalized)
    size_t nposits = GetNoPositsOvhd();//#positions to be written

    if(nposits < 1)
        return pmbsdfEmpty;

    size_t npositswrtbsd = bsd.GetNoPositsWritten();//#positions written
    size_t szstruct = PMBatchStrData::GetPMDataSize1(nposits);//size of data to be written
    size_t szstructswrtbsd = PMBatchStrData::GetPMDataSize1(npositswrtbsd);//size of data written
    size_t nstructswrtbsd = bsd.GetNoStructsWritten();//number of structures written
    static constexpr size_t szalign = PMBSdatalignment * pmv2DTotFlds;//max size for data alignment

    if(bsd.maxdatasize_ < szstruct+szalign || 
       bsd.maxdatalen_ < nposits)//allocation too small
    {
        //NOTE: this should not happen!
        return pmbsdfAbandoned;
    }

    if(bsd.maxdatasize_ < szstructswrtbsd + szstruct + szalign || 
       bsd.maxdatalen_ < npositswrtbsd + nposits || 
       bsd.maxnstrs_ < nstructswrtbsd + 1)
    {
        //NOTE: should not happen!
        return pmbsdfLimits;
    }

    if((size_t)INT_MAX <= npositswrtbsd + nposits)
        //NOTE: should not happen!
        throw MYRUNTIME_ERROR2( 
        preamb + "Overflow detected. Data chunk size must be reduced.", 
        CRITICAL);

    //{{actual copy
    int n = 0;
    for(; n < pps2DStrFlds; n++) {
        memcpy(bsd.bdbCpmend_[n], bdbCpmend_[n], TPM2DVectorFieldSize::szvfs_[n]);
        bsd.bdbCpmend_[n] += TPM2DVectorFieldSize::szvfs_[n];
        bsd.bdbCpmendovhd_[n] = bsd.bdbCpmend_[n];
    }
    for(; n < pmv2DTotFlds; n++) {
        memcpy(bsd.bdbCpmend_[n], bdbCpmend_[n], TPM2DVectorFieldSize::szvfs_[n] * nposits);
        bsd.bdbCpmend_[n] += TPM2DVectorFieldSize::szvfs_[n] * nposits;
        bsd.bdbCpmendovhd_[n] = bsd.bdbCpmend_[n];
    }
    //copy description
    char* const * ptrdescs = bdbCptrdescs_.get();//assume valid pointers
    char** bsdptrdescs = bsd.bdbCptrdescs_.get();
    size_t nstructswrt = GetNoStructsWritten();//number of structures written
    strcpy(bsdptrdescs[nstructswrtbsd], ptrdescs[nstructswrt]);
    //}}

    if(tmpclustfilendxs_.size() != tmpcluststrndxs_.size() ||
       tmpclustfilendxs_.size() != nstructswrt + 1)
        throw MYRUNTIME_ERROR2( 
        preamb + "Inconsistent file index size.", CRITICAL);

    if(tmpclustfilendxs_.size()) bsd.tmpclustfilendxs_.push_back(tmpclustfilendxs_.back());
    if(tmpcluststrndxs_.size()) bsd.tmpcluststrndxs_.push_back(tmpcluststrndxs_.back());

    return pmbsdfWritten;
}

// -------------------------------------------------------------------------
// GetOvhdStrDescription: get the description of the structure in the 
// overhead buffer
/**
 * @brief 在结构数据读取与布局中读取 `PMBatchStrData::GetOvhdStrDescription` 对应的数据。
 * @par 参数
 * 无。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline
const char* PMBatchStrData::GetOvhdStrDescription() const
{
    size_t nposits = GetNoPositsOvhd();//#positions to be written

    if(nposits < 1) return "";

    INTYPE nstructswrt = GetNoStructsWritten();
    char* const * ptrdescs = bdbCptrdescs_.get();

    return ptrdescs? ptrdescs[nstructswrt]: "";
}

// -------------------------------------------------------------------------

#endif//__PMBatchStrData_h__
