/***************************************************************************
 *   Copyright (C) 2021-2023 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __Devices_h__
#define __Devices_h__

#include <stdio.h>

#include <string>
#include <memory>
#include <vector>

#include "libutil/mybase.h"
#include "tsafety/TSCounterVar.h"

#define MAXNDEVS 32

class Devices;

extern Devices DEVPROPs;

// -------------------------------------------------------------------------
//
struct DeviceProperties {
    enum {
        DEVMEMORYRESERVE = 256 * ONEM
    };
    /**
     * @brief 构造 `DeviceProperties`，初始化其负责的CUDA 设备管理状态。
     * @par 参数
     * 无。
     * @return 无返回值；完成对象构造与初始状态设置。
     */
    DeviceProperties():
        shdcnt_(new TSCounterVar),
        devid_(-1),
        ccmajor_(0),ccminor_(0),
        totmem_(0),reqmem_(0),
        textureAlignment_(0),maxTexture1DLinear_(0),
        deviceOverlap_(0),asyncEngineCount_(0),
        computeMode_(0)
    {
        memset(maxGridSize_, 0, 3 * sizeof(int));
    };
    /**
     * @brief 构造 `DeviceProperties`，初始化其负责的CUDA 设备管理状态。
     * @param dprop 供该函数读取或更新的 `dprop` 参数。
     * @return 无返回值；完成对象构造与初始状态设置。
     */
    DeviceProperties( const DeviceProperties& dprop ):
        shdcnt_(dprop.shdcnt_),
        devid_(dprop.devid_),
        ccmajor_(dprop.ccmajor_),
        ccminor_(dprop.ccminor_),
        totmem_(dprop.totmem_),
        reqmem_(dprop.reqmem_),
        textureAlignment_(dprop.textureAlignment_),
        maxTexture1DLinear_(dprop.maxTexture1DLinear_),
        deviceOverlap_(dprop.deviceOverlap_),
        asyncEngineCount_(dprop.asyncEngineCount_),
        computeMode_(dprop.computeMode_),
        name_(dprop.name_)
    {
        memcpy(maxGridSize_, dprop.maxGridSize_, 3 * sizeof(int));
    };
    /**
     * @brief 在CUDA 设备管理中处理 `DevidValid` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    bool DevidValid() const {return devid_ < 0? false: true;}
    /**
     * @brief 在CUDA 设备管理中处理 `GridMaxXdim` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    int GridMaxXdim() const {return maxGridSize_[0];}
    /**
     * @brief 在CUDA 设备管理中处理 `GridMaxYdim` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    int GridMaxYdim() const {return maxGridSize_[1];}
    //shared device counter for synchronization between worker threads 
    // communicating with the same device:
    std::shared_ptr<TSCounterVar> shdcnt_;
    int devid_;//device id
    int ccmajor_;//major compute capability
    int ccminor_;//minor compute capability
    size_t totmem_;//total global memory
    size_t reqmem_;//requested amount of global memory
    size_t textureAlignment_;//alignment size for textures
    int maxTexture1DLinear_;//1D texture size
    int maxGridSize_[3];//maximumm grid size in each dimension
    int deviceOverlap_;//concurrent memory transfer and kernel execution
    int asyncEngineCount_;//number of asynchronous memory transfer engines
    int computeMode_;//enum:
    //cudaComputeModeDefault = 0
    //cudaComputeModeExclusive = 1
    //cudaComputeModeProhibited = 2
    //cudaComputeModeExclusiveProcess = 3
    std::string name_;
};

// _________________________________________________________________________
// Class Devices
//
// Implementation of queries for device properties
//
class Devices
{
public:
    /**
     * @brief 构造 `Devices`，初始化其负责的CUDA 设备管理状态。
     * @param maxdvs 供该函数读取或更新的 `maxdvs` 参数。
     * @return 无返回值；完成对象构造与初始状态设置。
     */
    Devices( int maxdvs = MAXNDEVS ): 
        maxdevs_(maxdvs),
        maxmem_(-1L)
    {};
    /**
     * @brief 销毁 `Devices`，释放其管理的CUDA 设备管理资源。
     * @par 参数
     * 无。
     * @return 无返回值；对象持有的资源在返回前完成释放。
     */
    ~Devices() {};

    /**
     * @brief 在CUDA 设备管理中设置 `SetMaxMemoryAmount` 对应的数据。
     * @param mem_in_mb 供当前步骤读取或更新的 `mem_in_mb` 缓冲区。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void SetMaxMemoryAmount( ssize_t mem_in_mb ) { maxmem_ = mem_in_mb * ONEM; };

    /**
     * @brief 在CUDA 设备管理中格式化输出 `PrintDevices` 对应的数据。
     * @param FILE 供该函数读取或更新的 `FILE` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void PrintDevices( FILE* );
    /**
     * @brief 在CUDA 设备管理中处理 `PrettyPrintUsedDevices` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void PrettyPrintUsedDevices();

    /**
     * @brief 在CUDA 设备管理中处理 `RegisterDevices` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void RegisterDevices();

    /**
     * @brief 在CUDA 设备管理中读取 `GetNDevices` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    int GetNDevices() const;
    /**
     * @brief 在CUDA 设备管理中读取 `GetDevicePropertiesAt` 对应的数据。
     * @param n 控制当前步骤范围或规模的 `n`。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    const DeviceProperties* GetDevicePropertiesAt(size_t n) const;

    /**
     * @brief 在CUDA 设备管理中读取 `GetDevIdWithMinRequestedMem` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    int GetDevIdWithMinRequestedMem() const;

protected:
    /**
     * @brief 在CUDA 设备管理中读取 `ReadDevices` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void ReadDevices();
    /**
     * @brief 在CUDA 设备管理中排序 `SortDevices` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void SortDevices();
    /**
     * @brief 在CUDA 设备管理中处理 `RegisterDeviceProperties` 对应的数据。
     * @param devid 供该函数读取或更新的 `devid` 参数。
     * @param maxmem 供当前步骤读取或更新的 `maxmem` 缓冲区。
     * @param checkduplicates 控制该处理分支是否启用的 `checkduplicates` 标志。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    bool RegisterDeviceProperties( int devid, ssize_t maxmem, bool checkduplicates );
    /**
     * @brief 在CUDA 设备管理中处理 `PruneRegisteredDevices` 对应的数据。
     * @param ndevs 控制当前步骤范围或规模的 `ndevs`。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void PruneRegisteredDevices( int ndevs );

private:
    int maxdevs_;//maximum number of devices
    ssize_t maxmem_;//maximum memory amount for all devices
    std::vector<DeviceProperties> devs_;//list of devices
};

////////////////////////////////////////////////////////////////////////////
// INLINES
//
// -------------------------------------------------------------------------
//
/**
 * @brief 在CUDA 设备管理中读取 `Devices::GetNDevices` 对应的数据。
 * @par 参数
 * 无。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline
int Devices::GetNDevices() const
{
    MYMSG( "Devices::GetNDevices", 7 );
    return (int)devs_.size();
}

// -------------------------------------------------------------------------
//
/**
 * @brief 在CUDA 设备管理中读取 `Devices::GetDevicePropertiesAt` 对应的数据。
 * @param n 控制当前步骤范围或规模的 `n`。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline
const DeviceProperties* Devices::GetDevicePropertiesAt(size_t n) const
{
    MYMSG( "Devices::GetDevicePropertiesAt", 6 );
    if( devs_.size() <= n || devs_.data() == NULL )
        return NULL;
    return devs_.data() + n;
}

#endif//__Devices_h__
