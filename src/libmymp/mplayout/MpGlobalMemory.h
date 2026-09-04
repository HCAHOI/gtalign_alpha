/***************************************************************************
 *   Copyright (C) 2021-2023 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __MpGlobalMemory_h__
#define __MpGlobalMemory_h__

#include "libutil/mybase.h"

#include <stdio.h>
#include <stdlib.h>
#include <string>

#include "libutil/CLOptions.h"
#include "libmympbase/mplayout/CuMemoryBase.h"


////////////////////////////////////////////////////////////////////////////
// CLASS MpGlobalMemory
// global memory arrangement for structure search and alignment computation
//
class MpGlobalMemory: public CuMemoryBase
{
public:
    /**
     * @brief 构造 `MpGlobalMemory`，初始化其负责的CPU 对齐流水线状态。
     * @param deviceallocsize 控制当前步骤范围或规模的 `deviceallocsize`。
     * @param nareas 控制当前步骤范围或规模的 `nareas`。
     * @return 无返回值；完成对象构造与初始状态设置。
     */
    MpGlobalMemory(
        size_t deviceallocsize,
        int nareas );

    /**
     * @brief 销毁 `MpGlobalMemory`，释放其管理的CPU 对齐流水线资源。
     * @par 参数
     * 无。
     * @return 无返回值；对象持有的资源在返回前完成释放。
     */
    virtual ~MpGlobalMemory();

    /**
     * @brief 在CPU 对齐流水线中读取 `GetDeviceName` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    const std::string& GetDeviceName() const { return devname_;}

    /**
     * @brief 在CPU 对齐流水线中读取 `GetHeap` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    virtual char* GetHeap() const {return g_heap_;}

    /**
     * @brief 在CPU 对齐流水线中读取 `GetMemAlignment` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    virtual size_t GetMemAlignment() const {
        // return CuMemoryBase::GetMemAlignment();
        size_t cszalnment = 512UL;
        return cszalnment;
    }

protected:

    /**
     * @brief 在CPU 对齐流水线中分配 `AllocateHeap` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    virtual void AllocateHeap();
    /**
     * @brief 在CPU 对齐流水线中释放 `DeallocateHeap` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    virtual void DeallocateHeap();

private:
    std::string devname_;
    char* g_heap_;//global heap containing all data written, generated, and read
};

// -------------------------------------------------------------------------
// INLINES ...
//

#endif//__MpGlobalMemory_h__
