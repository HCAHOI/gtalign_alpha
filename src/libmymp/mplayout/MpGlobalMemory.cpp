/***************************************************************************
 *   Copyright (C) 2021-2023 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "libutil/mybase.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cstdlib>

#include "libmympbase/mplayout/CuMemoryBase.h"
#include "MpGlobalMemory.h"

// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
// constructor
//
/**
 * @brief 构造 `MpGlobalMemory`，初始化其负责的CPU 对齐流水线状态。
 * @param deviceallocsize 控制当前步骤范围或规模的 `deviceallocsize`。
 * @param nareas 控制当前步骤范围或规模的 `nareas`。
 * @return 无返回值；完成对象构造与初始状态设置。
 */
MpGlobalMemory::MpGlobalMemory(
    size_t deviceallocsize,
    int nareas)
:
    CuMemoryBase(deviceallocsize, nareas),
    devname_("n/a"),
    g_heap_(NULL)
{
    MYMSG("MpGlobalMemory::MpGlobalMemory", 4);
    Initialize();
}

// -------------------------------------------------------------------------
// destructor
//
/**
 * @brief 销毁 `MpGlobalMemory`，释放其管理的CPU 对齐流水线资源。
 * @par 参数
 * 无。
 * @return 无返回值；对象持有的资源在返回前完成释放。
 */
MpGlobalMemory::~MpGlobalMemory()
{
    MYMSG("MpGlobalMemory::~MpGlobalMemory", 4);
}





// =========================================================================
// AllocateHeap: allocate device memory
/**
 * @brief 在CPU 对齐流水线中分配 `MpGlobalMemory::AllocateHeap` 对应的数据。
 * @par 参数
 * 无。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
inline
void MpGlobalMemory::AllocateHeap()
{
    MYMSG("MpGlobalMemory::AllocateHeap", 6);
    g_heap_ = (char*)my_aligned_alloc(GetMemAlignment(), GetAllocSize());
    if(g_heap_ == NULL)
        throw MYRUNTIME_ERROR(
        "MpGlobalMemory::AllocateHeap: Not enough memory.");
}

// -------------------------------------------------------------------------
// FreeDevicePtr: free device pointer
/**
 * @brief 在CPU 对齐流水线中释放 `MpGlobalMemory::DeallocateHeap` 对应的数据。
 * @par 参数
 * 无。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
inline
void MpGlobalMemory::DeallocateHeap()
{
    MYMSG("MpGlobalMemory::DeallocateHeap", 6);
    if(g_heap_) {
        my_aligned_free(g_heap_);
        g_heap_ = NULL;
    }
}
