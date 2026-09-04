/***************************************************************************
 *   Copyright (C) 2021-2023 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __myatomic_cuh__
#define __myatomic_cuh__

// -------------------------------------------------------------------------
// atomicMinFloat: atomic min function for float; based on CUDA programming 
// guide (B14.Atomic Functions) and Bonsai github
// 
__device__ __forceinline__
/**
 * @brief 在CUDA 通用原语中处理 `atomicMinFloat` 对应的数据。
 * @param address 供该函数读取或更新的 `address` 参数。
 * @param val 供该函数读取或更新的 `val` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
float atomicMinFloat(float* address, float val)
{
    int* address_as_int = (int*)address;
    int old = __float_as_int(*address);

    while(val < __int_as_float(old))
    {
        int assumed = old;
        //NOTE: use integer comparison to avoid hang in case of NaN (since NaN != NaN)
        if((old = atomicCAS(address_as_int, assumed, __float_as_int(val))) == assumed)
            break;
    }

    return __int_as_float(old);
}

#endif//__myatomic_cuh__
