/***************************************************************************
 *   Copyright (C) 2021-2023 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __DRDataDeleter_h__
#define __DRDataDeleter_h__

#include "libutil/mybase.h"

#include <memory>

#if defined(GPUINUSE) && 0
#   include <cuda_runtime_api.h>
#endif

#include "libutil/CLOptions.h"

// -------------------------------------------------------------------------
//
struct DRDataDeleter {
    /**
     * @brief 在结构数据读取与布局中处理 `operator()` 对应的数据。
     * @param p 供该函数读取或更新的 `p` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void operator()(void* p) const {
        if(p)
            std::free(p);
    };
};

struct DRHostDataDeleter {
    /**
     * @brief 在结构数据读取与布局中处理 `operator()` 对应的数据。
     * @param p 供该函数读取或更新的 `p` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void operator()(void* p) const {
        if(p) {
#if defined(GPUINUSE) && 0
            if(CLOptions::GetIO_UNPINNED() == 0)
                cudaFreeHost(p);
            else
#endif
                // std::free(p);
                my_aligned_free(p);
        }
    };
};

// -------------------------------------------------------------------------

#endif//__DRDataDeleter_h__
