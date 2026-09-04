/***************************************************************************
 *   Copyright (C) 2021-2023 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __alpha_h__
#define __alpha_h__

#include "platform.h"

#if defined (GPUINUSE) || defined(OS_MS_WINDOWS)
#define ALPHA_OMPDECLARE
#else 
#define ALPHA_OMPDECLARE _Pragma("omp declare simd notinbranch ")
#endif

// cardinality of the English alphabet
#define NEA 26

//{{SS states
enum SSSTATES {
    SS_C,
    SS_E,
    SS_H,
    SS_NSTATES
};
extern const char* gSSAlphabet;//{"CEH"}
extern const char* gSSlcAlphabet;//{" eh"}
//}}

//NOTE: map functions for 3-letter residue names!
/**
 * @brief 在通用工具中处理 `ResName2Code` 对应的数据。
 * @param param1 供该函数读取或更新的 `param1` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
char ResName2Code(const char*);
/**
 * @brief 在通用工具中处理 `ResCode2Name` 对应的数据。
 * @param param1 供该函数读取或更新的 `param1` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
const char* ResCode2Name(char);

//for testing
/**
 * @brief 在通用工具中处理 `testmyresnamehash` 对应的数据。
 * @par 参数
 * 无。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
void testmyresnamehash();

// GONNET residue scores
struct _GONNET_SCORES_
{
    /**
     * @brief 构造 `_GONNET_SCORES_`，初始化其负责的通用工具状态。
     * @par 参数
     * 无。
     * @return 无返回值；完成对象构造与初始状态设置。
     */
    _GONNET_SCORES_();
    /**
     * @brief 在通用工具中读取 `get` 对应的数据。
     * @param res1 供该函数读取或更新的 `res1` 参数。
     * @param res2 供该函数读取或更新的 `res2` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    float get(char res1, char res2);
    float data_[NEA * NEA];
};

extern _GONNET_SCORES_  GONNET_SCORES;

// -------------------------------------------------------------------------
// _GONNET_SCORES_::get: get a pairwise residue score obtained from the
// Gonnet frequencies
//
ALPHA_OMPDECLARE
/**
 * @brief 在通用工具中读取 `_GONNET_SCORES_::get` 对应的数据。
 * @param res1 供该函数读取或更新的 `res1` 参数。
 * @param res2 供该函数读取或更新的 `res2` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline
float _GONNET_SCORES_::get(char res1, char res2)
{
    if('A' <= res1 && res1 <= 'Z' && 'A' <= res2 && res2 <= 'Z')
        return data_[(res1 - 'A') * NEA + (res2 - 'A')];
    return 0.0f;
}

// -------------------------------------------------------------------------

#endif//__alpha_h__
