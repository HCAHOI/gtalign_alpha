/***************************************************************************
 *   Copyright (C) 2021-2023 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "libutil/mybase.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <string>

#include "libutil/alpha.h"
#include "libgenp/gdats/PM2DVectorFields.h"
#include "libmympbase/mplayout/CuMemoryBase.h"

#include "libmycu/cucom/cucommon.h"
#include "libmycu/cuproc/Devices.h"
#include "cuconstant.cuh"
#include "CuDeviceMemory.cuh"

// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
// constructor
//
/**
 * @brief 构造 `CuDeviceMemory`，初始化其负责的CUDA 内存布局状态。
 * @param dprop 供该函数读取或更新的 `dprop` 参数。
 * @param deviceallocsize 控制当前步骤范围或规模的 `deviceallocsize`。
 * @param nareas 控制当前步骤范围或规模的 `nareas`。
 * @return 无返回值；完成对象构造与初始状态设置。
 */
CuDeviceMemory::CuDeviceMemory(
    DeviceProperties dprop, 
    size_t deviceallocsize,
    int nareas)
:
    CuMemoryBase(deviceallocsize, nareas),
    deviceProp_(dprop),
    d_heap_(NULL)
{
    MYMSG("CuDeviceMemory::CuDeviceMemory", 4);
    Initialize();
}

// -------------------------------------------------------------------------
// destructor
//
/**
 * @brief 销毁 `CuDeviceMemory`，释放其管理的CUDA 内存布局资源。
 * @par 参数
 * 无。
 * @return 无返回值；对象持有的资源在返回前完成释放。
 */
CuDeviceMemory::~CuDeviceMemory()
{
    MYMSG("CuDeviceMemory::~CuDeviceMemory", 4);
}

// -------------------------------------------------------------------------
// CacheCompleteData: transfer complete required data to device
// 
/**
 * @brief 在CUDA 内存布局中处理 `CuDeviceMemory::CacheCompleteData` 对应的数据。
 * @par 参数
 * 无。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
void CuDeviceMemory::CacheCompleteData()
{
    MYMSG("CuDeviceMemory::CacheCompleteData", 4);
    MYCUDACHECK(
        cudaMemcpyToSymbol(
            dc_Gonnet_scores_, GONNET_SCORES.data_, NEA * NEA * sizeof(float))
    );
    MYCUDACHECKLAST;
}





// =========================================================================
//device addresses of vectors representing structure data and constant data:

__constant__ char* dc_pm2dvfields_[sztot_pm2dvfields];
__constant__ float dc_Gonnet_scores_[NEA * NEA];

// =========================================================================
// -------------------------------------------------------------------------
// AllocateHeap: allocate device memory
/**
 * @brief 在CUDA 内存布局中分配 `CuDeviceMemory::AllocateHeap` 对应的数据。
 * @par 参数
 * 无。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
inline
void CuDeviceMemory::AllocateHeap()
{
    MYMSG("CuDeviceMemory::AllocateHeap", 6);

    MYCUDACHECK(cudaSetDevice(deviceProp_.devid_));
    MYCUDACHECKLAST;

    MYCUDACHECK(cudaMalloc((void**)&d_heap_, GetAllocSize()));
    MYCUDACHECKLAST;
}

// -------------------------------------------------------------------------
// FreeDevicePtr: free device pointer
/**
 * @brief 在CUDA 内存布局中释放 `CuDeviceMemory::FreeDevicePtr` 对应的数据。
 * @param d_ptr 供该函数读取或更新的 `d_ptr` 参数。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
inline
void CuDeviceMemory::FreeDevicePtr( char*& d_ptr )
{
    MYMSG("CuDeviceMemory::FreeDevicePtr", 6);
    if( d_ptr ) {
        MYCUDACHECK(cudaFree(d_ptr));
        d_ptr = NULL;
        MYCUDACHECKLAST;
    }
}

// -------------------------------------------------------------------------
// CopyCPMDataToDevice: copy packed structure data to device;
// NFIELDS, template parameter, #fields corresponding to their memory addresses in cmem;
// bdbCpmbeg, host pointers for beginning addresses of (query/target) db structures;
// bdbCpmend, host pointers for terminal addresses of (query/target) db structures;
// dev_pckdpm, device pointer to the beginning the data to be copied;
// szmaxsize, max memory size allowed to be occupied by the (arguments) data;
// cmbegndx, beginning index of device memory addresses in constant memory;
// NOTE: memory is allocated for device pointer dev_pckdpm 
//
template <size_t NFIELDS>
/**
 * @brief 在CUDA 内存布局中复制 `CuDeviceMemory::CopyCPMDataToDevice` 对应的数据。
 * @param bdbCpmbeg 参考结构打包字段的起始指针数组。
 * @param bdbCpmend 参考结构打包字段的结束指针数组。
 * @param dev_pckdpm 供该函数读取或更新的 `dev_pckdpm` 参数。
 * @param szmaxsize 控制当前步骤范围或规模的 `szmaxsize`。
 * @param cmbegndx 供该函数读取或更新的 `cmbegndx` 参数。
 * @param sinfo 供该函数读取或更新的 `sinfo` 参数。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
void CuDeviceMemory::CopyCPMDataToDevice(
    char** bdbCpmbeg,
    char** bdbCpmend,
    char* dev_pckdpm,
    const size_t szmaxsize,
    const size_t cmbegndx,
    const char* sinfo)
{
    MYMSG("CuDeviceMemory::CopyCPMDataToDevice", 4);
    static const std::string preamb = "CuDeviceMemory::CopyCPMDataToDevice: ";
    const size_t cszalnment = GetMemAlignment();
    myruntime_error mre;

    if( bdbCpmbeg == NULL || bdbCpmend == NULL )
        throw MYRUNTIME_ERROR(preamb + "Null data structures.");

    char msgbuf[BUF_MAX];
    size_t i, szfldna, szfld, sztot = 0UL;
    //device pointers to the fields of structure data packed to one buffer
    char* d_pckdpmdatflds[NFIELDS];

    try {
        for(i = 0; i < NFIELDS; i++ ) {
            if(bdbCpmbeg[i] == NULL || bdbCpmend[i] == NULL || 
               bdbCpmend[i] <= bdbCpmbeg[i]) {
                sprintf( msgbuf, "Invalid pointers to data field %zu: %p - %p",
                        i, bdbCpmbeg[i], bdbCpmend[i]);
                throw MYRUNTIME_ERROR(preamb + msgbuf);
            }

            szfldna = (size_t)(bdbCpmend[i] - bdbCpmbeg[i]);
            szfld = ALIGN_UP(szfldna, cszalnment);

            //NOTE: increment before sztot gets updated
            d_pckdpmdatflds[i] = dev_pckdpm + sztot;

            sztot += szfld;

            if(szmaxsize < sztot) {
                sprintf(msgbuf,
                        "%sInsufficient allocated device memory: "
                        "%zuB requested vs %zuB allowed.",
                        preamb.c_str(), sztot, szmaxsize);
                throw MYRUNTIME_ERROR(preamb + msgbuf);
            }

            //NOTE: secondary structure not calculated in host: no copy;
            if(i != pmv2Dss)
                //synchronous copy; have to wait for data to be present on device
                MYCUDACHECK(cudaMemcpy(d_pckdpmdatflds[i], bdbCpmbeg[i], szfldna,
                            cudaMemcpyHostToDevice));
            MYCUDACHECKLAST;
        }

        //save device addresses to constant memory
        MYCUDACHECK(cudaMemcpyToSymbol(
            dc_pm2dvfields_, d_pckdpmdatflds, NFIELDS * sizeof(char*), 
                        cmbegndx * sizeof(char*)/*offset*/));
        MYCUDACHECKLAST;

        MYMSGBEGl(3)
            sprintf(msgbuf, "%sChunk of %s transferred:  "
                "total data size= %zu", preamb.c_str(), sinfo, sztot);
            MYMSG(msgbuf, 3);
        MYMSGENDl

    } catch( myexception const& ex ) {
        mre = ex;
    }

    if( mre.isset())
        throw mre;
}

// -------------------------------------------------------------------------
// TransferQueryPMDataToDevice: transfer a chunk of query db structure 
// data to device
// 
/**
 * @brief 在CUDA 内存布局中处理 `CuDeviceMemory::TransferQueryPMDataToDevice` 对应的数据。
 * @param querypmbeg 查询结构打包字段的起始指针数组。
 * @param querypmend 查询结构打包字段的结束指针数组。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
void CuDeviceMemory::TransferQueryPMDataToDevice(
    char** querypmbeg, char** querypmend)
{
    MYMSG("CuDeviceMemory::TransferQueryPMDataToDevice", 4);
    static const std::string preamb = "CuDeviceMemory::TransferQueryPMDataToDevice: ";

    if( sz_heapsections_[0][ddsEndOfQrsChunk] <= sz_heapsections_[0][ddsBegOfQrsChunk])
        throw MYRUNTIME_ERROR(preamb + "Unallocated memory for queries.");

    //beginning of relative data
    char* d_qrschunkdata = d_heap_ + sz_heapsections_[0][ddsBegOfQrsChunk];

    //maxsize for data
    size_t szmaxallowed = 
        sz_heapsections_[0][ddsEndOfQrsChunk] - sz_heapsections_[0][ddsBegOfQrsChunk];

    CopyCPMDataToDevice<pmv2DTotFlds>(
        querypmbeg, querypmend,
        d_qrschunkdata,
        szmaxallowed,
        ndx_qrs_dc_pm2dvfields_,
        "query structures");
}

// -------------------------------------------------------------------------
// TransferQueryPMIndexToDevice: transfer indexed query db structure 
// data to device
// 
/**
 * @brief 在CUDA 内存布局中处理 `CuDeviceMemory::TransferQueryPMIndexToDevice` 对应的数据。
 * @param queryndxpmbeg 描述查询结构的 `queryndxpmbeg`。
 * @param queryndxpmend 描述查询结构的 `queryndxpmend`。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
void CuDeviceMemory::TransferQueryPMIndexToDevice(
    char** queryndxpmbeg, char** queryndxpmend)
{
    MYMSG("CuDeviceMemory::TransferQueryPMIndexToDevice", 4);
    static const std::string preamb = "CuDeviceMemory::TransferQueryPMIndexToDevice: ";

    if( sz_heapsections_[0][ddsEndOfQrsIndex] <= sz_heapsections_[0][ddsBegOfQrsIndex])
        throw MYRUNTIME_ERROR(preamb + "Unallocated memory for query index.");

    //beginning of relative data
    char* d_qrsindexdata = d_heap_ + sz_heapsections_[0][ddsBegOfQrsIndex];

    //maxsize for data
    size_t szmaxallowed = 
        sz_heapsections_[0][ddsEndOfQrsIndex] - sz_heapsections_[0][ddsBegOfQrsIndex];

    CopyCPMDataToDevice<pmv2DTotIndexFlds>(
        queryndxpmbeg, queryndxpmend,
        d_qrsindexdata,
        szmaxallowed,
        ndx_qrs_dc_pm2dvndxfds_,
        "query index");
}

// -------------------------------------------------------------------------
// TransferCPMDataToDevice: transfer a chunk of target db structure data to 
// device
// 
/**
 * @brief 在CUDA 内存布局中处理 `CuDeviceMemory::TransferCPMDataToDevice` 对应的数据。
 * @param bdbCpmbeg 参考结构打包字段的起始指针数组。
 * @param bdbCpmend 参考结构打包字段的结束指针数组。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
void CuDeviceMemory::TransferCPMDataToDevice(
    char** bdbCpmbeg,
    char** bdbCpmend)
{
    MYMSG("CuDeviceMemory::TransferCPMDataToDevice", 4);
    static const std::string preamb = "CuDeviceMemory::TransferCPMDataToDevice: ";

    if( sz_heapsections_[0][ddsEndOfDbsChunk] <= sz_heapsections_[0][ddsBegOfDbsChunk])
        throw MYRUNTIME_ERROR(preamb + "Unallocated memory for target db data.");

    //beginning of relative data
    char* d_dbschunkdata = d_heap_ + sz_heapsections_[0][ddsBegOfDbsChunk];

    //maxsize for data
    size_t szmaxallowed = 
        sz_heapsections_[0][ddsEndOfDbsChunk] - sz_heapsections_[0][ddsBegOfDbsChunk];

    CopyCPMDataToDevice<pmv2DTotFlds>(
        bdbCpmbeg, bdbCpmend,
        d_dbschunkdata,
        szmaxallowed,
        ndx_dbs_dc_pm2dvfields_,
        "Db structures");
}

// -------------------------------------------------------------------------
// TransferCPMIndexToDevice: transfer indexed target db structure data to device
// 
/**
 * @brief 在CUDA 内存布局中处理 `CuDeviceMemory::TransferCPMIndexToDevice` 对应的数据。
 * @param bdbCndxpmbeg 描述参考结构的 `bdbCndxpmbeg`。
 * @param bdbCndxpmend 描述参考结构的 `bdbCndxpmend`。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
void CuDeviceMemory::TransferCPMIndexToDevice(
    char** bdbCndxpmbeg, char** bdbCndxpmend)
{
    MYMSG("CuDeviceMemory::TransferCPMIndexToDevice", 4);
    static const std::string preamb = "CuDeviceMemory::TransferCPMIndexToDevice: ";

    if( sz_heapsections_[0][ddsEndOfDbsIndex] <= sz_heapsections_[0][ddsBegOfDbsIndex])
        throw MYRUNTIME_ERROR(preamb + "Unallocated memory for target db index.");

    //beginning of relative data
    char* d_dbsindexdata = d_heap_ + sz_heapsections_[0][ddsBegOfDbsIndex];

    //maxsize for data
    size_t szmaxallowed = 
        sz_heapsections_[0][ddsEndOfDbsIndex] - sz_heapsections_[0][ddsBegOfDbsIndex];

    CopyCPMDataToDevice<pmv2DTotIndexFlds>(
        bdbCndxpmbeg, bdbCndxpmend,
        d_dbsindexdata,
        szmaxallowed,
        ndx_dbs_dc_pm2dvndxfds_,
        "Db index");
}

