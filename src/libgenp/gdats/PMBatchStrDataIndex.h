/***************************************************************************
 *   Copyright (C) 2021-2023 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __PMBatchStrDataIndexIndex_h__
#define __PMBatchStrDataIndexIndex_h__

#include "libutil/mybase.h"

#include <string.h>

#include <vector>
#include <memory>

#include "libutil/CLOptions.h"
#include "PM2DVectorFields.h"
#include "DRDataDeleter.h"
#include "PMBatchStrData.h"

// -------------------------------------------------------------------------
// PMBatchStrDataIndex: Complete indexed batch structure data for parallel 
// processing
//
class PMBatchStrDataIndex {
public:
    /**
     * @brief 构造 `PMBatchStrDataIndex`，初始化其负责的结构数据读取与布局状态。
     * @par 参数
     * 无。
     * @return 无返回值；完成对象构造与初始状态设置。
     */
    PMBatchStrDataIndex()
    :   tmpbdbCdata_(nullptr),
        bdbCdata_(nullptr),
        szbdbCdata_(0)
    {
        memset(bdbCpmbeg_, 0, pmv2DTotIndexFlds * sizeof(void*));
        memset(bdbCpmend_, 0, pmv2DTotIndexFlds * sizeof(void*));
    }

    /**
     * @brief 销毁 `PMBatchStrDataIndex`，释放其管理的结构数据读取与布局资源。
     * @par 参数
     * 无。
     * @return 无返回值；对象持有的资源在返回前完成释放。
     */
    ~PMBatchStrDataIndex() {}

    // *** METHODS ***
    //index coordinates present in PMBatchStrData object
    /**
     * @brief 在结构数据读取与布局中构造 `MakeIndex` 对应的数据。
     * @param PMBatchStrData 供该函数读取或更新的 `PMBatchStrData` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void MakeIndex(const PMBatchStrData&);

    /**
     * @brief 在结构数据读取与布局中读取 `GetPMDataSize` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    size_t GetPMDataSize();
    /**
     * @brief 在结构数据读取与布局中读取 `GetPMDataSize` 对应的数据。
     * @param totallen 控制当前步骤范围或规模的 `totallen`。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    static size_t GetPMDataSize(size_t totallen);

    /**
     * @brief 在结构数据读取与布局中处理 `ContainsData` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    bool ContainsData() const {//whether data is present
        return bdbCpmbeg_[pmv2DNdxCoords] < bdbCpmend_[pmv2DNdxCoords] &&
            bdbCpmbeg_[pmv2DNdxLeft] < bdbCpmend_[pmv2DNdxLeft] &&
            bdbCpmbeg_[pmv2DNdxRight] < bdbCpmend_[pmv2DNdxRight];
    }

    //change/pack structure pointers to include structures indexed in filterdata
    /**
     * @brief 在结构数据读取与布局中筛选 `FilterStructs` 对应的数据。
     * @param bdbndxpmbeg 描述参考结构的 `bdbndxpmbeg`。
     * @param bdbndxpmend 描述参考结构的 `bdbndxpmend`。
     * @param bdbpmbeg 描述参考结构的 `bdbpmbeg`。
     * @param bdbpmend 描述参考结构的 `bdbpmend`。
     * @param filterdata 供该函数读取或更新的 `filterdata` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    static void FilterStructs(
        char** bdbndxpmbeg, char** bdbndxpmend,
        char* const* bdbpmbeg, char* const* bdbpmend,
        const unsigned int* filterdata);

    /**
     * @brief 在结构数据读取与布局中读取 `GetNoPositsWritten` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    size_t GetNoPositsWritten() const {//#positions written in bdbCdata_
        return (size_t)(bdbCpmend_[pmv2DNdxCoords]-bdbCpmbeg_[pmv2DNdxCoords]) / SZFPTYPE;
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

    /**
     * @brief 在结构数据读取与布局中搜索 `Search` 对应的数据。
     * @param PMBatchStrData 供该函数读取或更新的 `PMBatchStrData` 参数。
     * @param crds 供该函数读取或更新的 `crds` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void Search(const PMBatchStrData&, FPTYPE crds[]) const;

    /**
     * @brief 在结构数据读取与布局中格式化输出 `Print` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void Print() const;

private:
    constexpr static size_t ndims_ = 3;
    enum {ptrLeft_, ptrRight_, nPtrs_};

    /**
     * @brief 在结构数据读取与布局中分配 `AllocateSpace` 对应的数据。
     * @param param1 供该函数读取或更新的 `param1` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void AllocateSpace(size_t);

    /**
     * @brief 在结构数据读取与布局中处理 `ConstructKdtree` 对应的数据。
     * @param crds 供该函数读取或更新的 `crds` 参数。
     * @param ptrs 供该函数读取或更新的 `ptrs` 参数。
     * @param begin 供该函数读取或更新的 `begin` 参数。
     * @param end 供该函数读取或更新的 `end` 参数。
     * @param dimndx 供该函数读取或更新的 `dimndx` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    int ConstructKdtree(
        FPTYPE* crds[], INTYPE* ptrs[],
        int begin, int end, int dimndx);

    /**
     * @brief 在结构数据读取与布局中处理 `NNRecursive` 对应的数据。
     * @param address 供该函数读取或更新的 `address` 参数。
     * @param root 供该函数读取或更新的 `root` 参数。
     * @param crds 供该函数读取或更新的 `crds` 参数。
     * @param dimndx 供该函数读取或更新的 `dimndx` 参数。
     * @param nvisited 控制当前步骤范围或规模的 `nvisited`。
     * @param nestndx 控制当前步骤范围或规模的 `nestndx`。
     * @param nestdst 控制当前步骤范围或规模的 `nestdst`。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void NNRecursive(
        int address, int root, FPTYPE crds[ndims_], int dimndx,
        int& nvisited, int& nestndx, float& nestdst) const;

    /**
     * @brief 在结构数据读取与布局中处理 `NNIterative` 对应的数据。
     * @param address 供该函数读取或更新的 `address` 参数。
     * @param root 供该函数读取或更新的 `root` 参数。
     * @param crds 供该函数读取或更新的 `crds` 参数。
     * @param dimndx 供该函数读取或更新的 `dimndx` 参数。
     * @param nvisited 控制当前步骤范围或规模的 `nvisited`。
     * @param nestndx 控制当前步骤范围或规模的 `nestndx`。
     * @param nestdst 控制当前步骤范围或规模的 `nestdst`。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void NNIterative(
        int address, int root, FPTYPE crds[ndims_], int dimndx,
        int& nvisited, int& nestndx, float& nestdst) const;

    /**
     * @brief 在结构数据读取与布局中处理 `NNNaive` 对应的数据。
     * @param address 供该函数读取或更新的 `address` 参数。
     * @param len 控制当前步骤范围或规模的 `len`。
     * @param crds 供该函数读取或更新的 `crds` 参数。
     * @param nvisited 控制当前步骤范围或规模的 `nvisited`。
     * @param nestndx 控制当前步骤范围或规模的 `nestndx`。
     * @param nestdst 控制当前步骤范围或规模的 `nestdst`。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void NNNaive(
        int address, int len, FPTYPE crds[ndims_],
        int& nvisited, int& nestndx, float& nestdst) const;

    //tree branch (left/right) at a given position (node assumed to be valid):
    /**
     * @brief 在结构数据读取与布局中读取 `GetOrgndxAt` 对应的数据。
     * @param address 供该函数读取或更新的 `address` 参数。
     * @param node 控制当前步骤范围或规模的 `node`。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    INTYPE GetOrgndxAt(int address, int node) const {
        return ((INTYPE*)(bdbCpmbeg_[pmv2DNdxOrgndx]))[address + node];
    }

    //tree branch (left/right) at a given position (node assumed to be valid):
    /**
     * @brief 在结构数据读取与布局中读取 `GetBranchAt` 对应的数据。
     * @param address 供该函数读取或更新的 `address` 参数。
     * @param node 控制当前步骤范围或规模的 `node`。
     * @param branch 供该函数读取或更新的 `branch` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    INTYPE GetBranchAt(int address, int node, int branch) const {
        return ((INTYPE*)
            (bdbCpmbeg_[(branch == ptrLeft_)? pmv2DNdxLeft: pmv2DNdxRight]))
                [address + node];
    }

    //coordinate value at a given position (node assumed to be valid):
    /**
     * @brief 在结构数据读取与布局中读取 `GetCoordinateAt` 对应的数据。
     * @param address 供该函数读取或更新的 `address` 参数。
     * @param node 控制当前步骤范围或规模的 `node`。
     * @param dim 供该函数读取或更新的 `dim` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    FPTYPE GetCoordinateAt(int address, int node, int dim) const {
        return ((FPTYPE*)(bdbCpmbeg_[pmv2DNdxCoords+dim]))[address + node];
    }

    //squared distance to the point (node assumed to be valid):
    /**
     * @brief 在结构数据读取与布局中读取 `GetDistance2` 对应的数据。
     * @param address 供该函数读取或更新的 `address` 参数。
     * @param node 控制当前步骤范围或规模的 `node`。
     * @param crds 供该函数读取或更新的 `crds` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    float GetDistance2(int address, int node, FPTYPE crds[]) const {
        return
            SQRD(GetCoordinateAt(address, node, 0) - crds[0]) +
            SQRD(GetCoordinateAt(address, node, 1) - crds[1]) +
            SQRD(GetCoordinateAt(address, node, 2) - crds[2]);
    }

private:
    template<typename T, int field>
    /**
     * @brief 在结构数据读取与布局中处理 `index_helper_psfields_assign` 对应的数据。
     * @param nposs 控制当前步骤范围或规模的 `nposs`。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void index_helper_psfields_assign(size_t nposs);

    template <typename T, int field>
    /**
     * @brief 在结构数据读取与布局中处理 `index_helper_psfields_adjend` 对应的数据。
     * @param nposs 控制当前步骤范围或规模的 `nposs`。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void index_helper_psfields_adjend(size_t nposs);

private:
    std::unique_ptr<char,DRDataDeleter> tmpbdbCdata_;//temporary buffer for indexed data
    std::vector<INTYPE> tmpndxs_;//temporary index vector

public:
    // *** MEMBER VARIABLES ***
    std::unique_ptr<char,DRHostDataDeleter> bdbCdata_;
    char* bdbCpmbeg_[pmv2DTotIndexFlds];//addresses of the beginnings of the fields in bdbCdata_
    char* bdbCpmend_[pmv2DTotIndexFlds];//end addresses of the fields in bdbCdata_

    size_t szbdbCdata_;//size allocated for bdbCdata_
};

// =========================================================================
// INLINES
// -------------------------------------------------------------------------
// GetPMDataSize: get the size of the indexed structure data
/**
 * @brief 在结构数据读取与布局中读取 `PMBatchStrDataIndex::GetPMDataSize` 对应的数据。
 * @par 参数
 * 无。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline
size_t PMBatchStrDataIndex::GetPMDataSize()
{
    MYMSG("PMBatchStrDataIndex::GetPMDataSize",6);
    size_t size = 0;
    for(int n = 0; n < pmv2DTotIndexFlds; n++)
        size += (size_t)(bdbCpmend_[n]-bdbCpmbeg_[n]);
    return size;
}

// -------------------------------------------------------------------------
// GetPMDataSize: get the size of complete indexed structure model data;
// totallen, the total number of positions (residues/atoms);
/**
 * @brief 在结构数据读取与布局中读取 `PMBatchStrDataIndex::GetPMDataSize` 对应的数据。
 * @param totallen 控制当前步骤范围或规模的 `totallen`。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline
size_t PMBatchStrDataIndex::GetPMDataSize(size_t totallen)
{
    MYMSG("PMBatchStrDataIndex::GetPMDataSize [static]",6);
    size_t size = 0;
    for(int n = 0; n < pmv2DTotIndexFlds; n++)
        size += TPM2DIndexFieldSize::szvfs_[n] * totallen;
    return size;
}

#endif//__PMBatchStrDataIndexIndex_h__
