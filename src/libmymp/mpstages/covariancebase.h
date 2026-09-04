/***************************************************************************
 *   Copyright (C) 2021-2023 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __covariancebase_h__
#define __covariancebase_h__

#include "libutil/cnsts.h"
#include "libutil/CLOptions.h"
#include "libmycu/cucom/cudef.h"
#include "libmycu/cucom/cutemplates.h"

// -------------------------------------------------------------------------
// GetFragStepSize_frg_deep: get the the step size which corresponds to 
// the depth of superposition exploration based on fragments dependent upon 
// structure lengths; version of more extensive parallelization;
// NOTE: deep depth;
// length, structure length;
//
/**
 * @brief 在CPU 刚体拟合与评分中读取 `GetFragStepSize_frg_deep` 对应的数据。
 * @param length 控制当前步骤范围或规模的 `length`。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
__HDINLINE__
int GetFragStepSize_frg_deep(int length)
{
    if(length > 150) return 15;
    int leno3 = myhdmax(1, (int)((float)length * oneTHIRDf));
    if(leno3 < 15) return leno3;
    return 15;
}
// NOTE: high depth;
/**
 * @brief 在CPU 刚体拟合与评分中读取 `GetFragStepSize_frg_high` 对应的数据。
 * @param length 控制当前步骤范围或规模的 `length`。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
__HDINLINE__
int GetFragStepSize_frg_high(int length)
{
    if(length > 150) return 23;
    int leno3 = myhdmax(1, (int)((float)length * oneTHIRDf));
    if(leno3 < 15) return leno3;
    return 15;
}
// NOTE: medium depth;
/**
 * @brief 在CPU 刚体拟合与评分中读取 `GetFragStepSize_frg_medium` 对应的数据。
 * @param length 控制当前步骤范围或规模的 `length`。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
__HDINLINE__
int GetFragStepSize_frg_medium(int length)
{
    if(length > 250) return 45;
    if(length > 200) return 35;
    if(length > 150) return 25;
    int leno3 = myhdmax(1, (int)((float)length * oneTHIRDf));
    if(leno3 < 15) return leno3;
    return 15;
}
// NOTE: shallow depth;
/**
 * @brief 在CPU 刚体拟合与评分中读取 `GetFragStepSize_frg_shallow_factor` 对应的数据。
 * @par 参数
 * 无。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
__HDINLINE__
constexpr int GetFragStepSize_frg_shallow_factor() {return 2;}
//
/**
 * @brief 在CPU 刚体拟合与评分中读取 `GetFragStepSize_frg_shallow` 对应的数据。
 * @param length 控制当前步骤范围或规模的 `length`。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
__HDINLINE__
int GetFragStepSize_frg_shallow(int length)
{
    return
        GetFragStepSize_frg_medium(length) *
        GetFragStepSize_frg_shallow_factor();
}

// -------------------------------------------------------------------------
// GetQryRfnStepsize2: return the stepsize for query and reference when 
// calculating transformation matrices from variable-length fragments;
//
/**
 * @brief 在CPU 刚体拟合与评分中读取 `GetQryRfnStepsize2` 对应的数据。
 * @param depth 供该函数读取或更新的 `depth` 参数。
 * @param qrylen 控制当前步骤范围或规模的 `qrylen`。
 * @param dbstrlen 控制当前步骤范围或规模的 `dbstrlen`。
 * @param qrystepsz 描述查询结构的 `qrystepsz`。
 * @param rfnstepsz 描述参考结构的 `rfnstepsz`。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
__HDINLINE__
void GetQryRfnStepsize2(
    const int depth,
    int qrylen, int dbstrlen,
    int* qrystepsz, int* rfnstepsz)
{
    *qrystepsz = GetFragStepSize_frg_shallow(qrylen);
    *rfnstepsz = GetFragStepSize_frg_shallow(dbstrlen);
    if(depth == CLOptions::csdDeep) {
        *qrystepsz = GetFragStepSize_frg_deep(qrylen);
        *rfnstepsz = GetFragStepSize_frg_deep(dbstrlen);
    } else if(depth == CLOptions::csdHigh) {
        *qrystepsz = GetFragStepSize_frg_high(qrylen);
        *rfnstepsz = GetFragStepSize_frg_high(dbstrlen);
    } else if(depth == CLOptions::csdMedium) {
        *qrystepsz = GetFragStepSize_frg_medium(qrylen);
        *rfnstepsz = GetFragStepSize_frg_medium(dbstrlen);
    }
}

// -------------------------------------------------------------------------
// GetNAlnPoss_frg: return the maximum number of alignment positions 
// (fragment size) given the lengths and start positions of the query and 
// reference structures for length-dependent fragments;
// version of more extensive parallelization;
// arg3 is fragndx, fragment index determining the fragment size dependent 
// upon lengths;
/**
 * @brief 在CPU 刚体拟合与评分中读取 `GetNAlnPoss_frg` 对应的数据。
 * @param qrylen 控制当前步骤范围或规模的 `qrylen`。
 * @param dbstrlen 控制当前步骤范围或规模的 `dbstrlen`。
 * @param qrypos 描述查询结构的 `qrypos`。
 * @param rfnpos 描述参考结构的 `rfnpos`。
 * @param arg1 供该函数读取或更新的 `arg1` 参数。
 * @param arg2 供该函数读取或更新的 `arg2` 参数。
 * @param arg3 供该函数读取或更新的 `arg3` 参数。
 * @param seedapproachstruct 供该函数读取或更新的 `seedapproachstruct` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
__HDINLINE__
int GetNAlnPoss_frg(
    int qrylen, int dbstrlen,
    int /*qrypos*/, int /*rfnpos*/,
    int /*arg1*/, int /*arg2*/, int arg3,
    const int seedapproachstruct = 0)
{
    const int minlen = myhdmin(qrylen, dbstrlen);
    const int minleno2 = (minlen >> 1);

    if(seedapproachstruct) return myhdmin(seedapproachstruct/*64/128*/, minleno2);

    if(arg3 == 0) {
        int leno3 = myhdmax(1, (int)((float)minlen * oneTHIRDf));
        return myhdmin(20, leno3);
    }

    // arg3 == 1
    return myhdmin(100, minleno2);
}

// -------------------------------------------------------------------------

#endif//__covariancebase_h__
