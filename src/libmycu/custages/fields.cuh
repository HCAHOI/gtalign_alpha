/***************************************************************************
 *   Copyright (C) 2021-2026 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __fields_cuh__
#define __fields_cuh__

#include "libutil/alpha.h"
#include "libgenp/gdats/PM2DVectorFields.h"
#include "libmycu/culayout/cuconstant.cuh"

// =========================================================================
// template parameter values: query and reference structures
#define FLDS_STRUCTS_QRIES 0
#define FLDS_STRUCTS_REFNS 1

// -------------------------------------------------------------------------
// GetDbStrField: get a field value of a given reference (db) structure;
template<typename T, int field>
__device__ __forceinline__
/**
 * @brief 在CUDA 刚体拟合与评分中读取 `GetDbStrField` 对应的数据。
 * @param dbstrndx 描述参考结构的 `dbstrndx`。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
T GetDbStrField(int dbstrndx)
{
    return ((T*)(dc_pm2dvfields_[ndx_dbs_dc_pm2dvfields_+field]))[dbstrndx];
}

// SetDbStrField: set a field value for a given reference (db) structure;
template<typename T, int field>
__device__ __forceinline__
/**
 * @brief 在CUDA 刚体拟合与评分中设置 `SetDbStrField` 对应的数据。
 * @param dbstrndx 描述参考结构的 `dbstrndx`。
 * @param value 需要读取、写入或转换的值。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
void SetDbStrField(int dbstrndx, T value)
{
    ((T*)(dc_pm2dvfields_[ndx_dbs_dc_pm2dvfields_+field]))[dbstrndx] = value;
}



// -------------------------------------------------------------------------
// GetQueryStrField: get a field value of a given query structure;
template<typename T, int field>
__device__ __forceinline__
/**
 * @brief 在CUDA 刚体拟合与评分中读取 `GetQueryStrField` 对应的数据。
 * @param qrystrndx 描述查询结构的 `qrystrndx`。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
T GetQueryStrField(int qrystrndx)
{
    return ((T*)(dc_pm2dvfields_[ndx_qrs_dc_pm2dvfields_+field]))[qrystrndx];
}



// -------------------------------------------------------------------------
// GetDbStrLength: get the length of the given reference (db) structure;
__device__ __forceinline__
/**
 * @brief 在CUDA 刚体拟合与评分中读取 `GetDbStrLength` 对应的数据。
 * @param dbstrndx 描述参考结构的 `dbstrndx`。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
INTYPE GetDbStrLength(int dbstrndx)
{
    return ((INTYPE*)(dc_pm2dvfields_[ndx_dbs_dc_pm2dvfields_+pps2DLen]))[dbstrndx];
}

// GetDbStrDst: get the distance (address) of the given reference structure;
__device__ __forceinline__
/**
 * @brief 在CUDA 刚体拟合与评分中读取 `GetDbStrDst` 对应的数据。
 * @param dbstrndx 描述参考结构的 `dbstrndx`。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
LNTYPE GetDbStrDst(int dbstrndx)
{
    return ((LNTYPE*)(dc_pm2dvfields_[ndx_dbs_dc_pm2dvfields_+pps2DDist]))[dbstrndx];
}

// GetQueryLength: get the length of the given query structure;
__device__ __forceinline__
/**
 * @brief 在CUDA 刚体拟合与评分中读取 `GetQueryLength` 对应的数据。
 * @param qrystrndx 描述查询结构的 `qrystrndx`。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
INTYPE GetQueryLength(int qrystrndx)
{
    return ((INTYPE*)(dc_pm2dvfields_[ndx_qrs_dc_pm2dvfields_+pps2DLen]))[qrystrndx];
}

// GetQueryDst: get the distance (address) of the given query structure;
__device__ __forceinline__
/**
 * @brief 在CUDA 刚体拟合与评分中读取 `GetQueryDst` 对应的数据。
 * @param qrystrndx 描述查询结构的 `qrystrndx`。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
LNTYPE GetQueryDst(int qrystrndx)
{
    return ((LNTYPE*)(dc_pm2dvfields_[ndx_qrs_dc_pm2dvfields_+pps2DDist]))[qrystrndx];
}

// -------------------------------------------------------------------------
// GetQryRfnField: get the value of field fld of the given query or 
// reference structure;
//
template<int STRUCTS>
__device__ __forceinline__
/**
 * @brief 在CUDA 刚体拟合与评分中读取 `GetQryRfnField` 对应的数据。
 * @param strndx 供该函数读取或更新的 `strndx` 参数。
 * @param fld 供该函数读取或更新的 `fld` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
INTYPE GetQryRfnField(int strndx, int fld)
{
    return
    (STRUCTS == FLDS_STRUCTS_QRIES)
    ? ((INTYPE*)(dc_pm2dvfields_[ndx_qrs_dc_pm2dvfields_ + fld]))[strndx]
    : ((INTYPE*)(dc_pm2dvfields_[ndx_dbs_dc_pm2dvfields_ + fld]))[strndx];
}


// -------------------------------------------------------------------------
// GetDbStrLenDst: get the length and distance from the beginning of 
// data of the given reference (db) structure;
// NOTE: pps2DLen and pps2DDist assumed to be adjacent: see PM2DVectorFields.h!
// NOTE: pps2DLen and pps2DDist are written in the first slots of cache;
//
__device__ __forceinline__
/**
 * @brief 在CUDA 刚体拟合与评分中读取 `GetDbStrLenDst` 对应的数据。
 * @param dbstrndx 描述参考结构的 `dbstrndx`。
 * @param cache 供该函数读取或更新的 `cache` 参数。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
void GetDbStrLenDst(int dbstrndx, int* __restrict__ cache)
{
    cache[threadIdx.x] = 
        ((INTYPE*)(dc_pm2dvfields_[
            ndx_dbs_dc_pm2dvfields_+pps2DLen+threadIdx.x]))[dbstrndx];
}

// GetQueryLenDst: get the length and distance for the given query structure;
// the same notes hold;
//
__device__ __forceinline__
/**
 * @brief 在CUDA 刚体拟合与评分中读取 `GetQueryLenDst` 对应的数据。
 * @param qrystrndx 描述查询结构的 `qrystrndx`。
 * @param cache 供该函数读取或更新的 `cache` 参数。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
void GetQueryLenDst(int qrystrndx, int* __restrict__ cache)
{
    cache[threadIdx.x] = 
        ((INTYPE*)(dc_pm2dvfields_[
            ndx_qrs_dc_pm2dvfields_+pps2DLen+threadIdx.x]))[qrystrndx];
}

// -------------------------------------------------------------------------
// GetDbStrSS: get the secondary structure assignment at the given reference 
// structure position;
__device__ __forceinline__
/**
 * @brief 在CUDA 刚体拟合与评分中读取 `GetDbStrRsd` 对应的数据。
 * @param pos 供该函数读取或更新的 `pos` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
CHTYPE GetDbStrRsd(int pos)
{
    return ((CHTYPE*)(dc_pm2dvfields_[ndx_dbs_dc_pm2dvfields_ + pmv2Drsd]))[pos];
}

// GetQuerySS: get the secondary structure assignment at the given query 
// structure position;
__device__ __forceinline__
/**
 * @brief 在CUDA 刚体拟合与评分中读取 `GetQueryRsd` 对应的数据。
 * @param pos 供该函数读取或更新的 `pos` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
CHTYPE GetQueryRsd(int pos)
{
    return ((CHTYPE*)(dc_pm2dvfields_[ndx_qrs_dc_pm2dvfields_ + pmv2Drsd]))[pos];
}

// -------------------------------------------------------------------------
// GetDbStrSS: get the secondary structure assignment at the given reference 
// structure position;
__device__ __forceinline__
/**
 * @brief 在CUDA 刚体拟合与评分中读取 `GetDbStrSS` 对应的数据。
 * @param pos 供该函数读取或更新的 `pos` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
CHTYPE GetDbStrSS(int pos)
{
    return ((CHTYPE*)(dc_pm2dvfields_[ndx_dbs_dc_pm2dvfields_ + pmv2Dss]))[pos];
}

// GetQuerySS: get the secondary structure assignment at the given query 
// structure position;
__device__ __forceinline__
/**
 * @brief 在CUDA 刚体拟合与评分中读取 `GetQuerySS` 对应的数据。
 * @param pos 供该函数读取或更新的 `pos` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
CHTYPE GetQuerySS(int pos)
{
    return ((CHTYPE*)(dc_pm2dvfields_[ndx_qrs_dc_pm2dvfields_ + pmv2Dss]))[pos];
}

// GetQryRfnSS: get the value of secondary structure of the given query or 
// reference structure at the given pos;
//
template<int STRUCTS>
__device__ __forceinline__
/**
 * @brief 在CUDA 刚体拟合与评分中读取 `GetQryRfnSS` 对应的数据。
 * @param pos 供该函数读取或更新的 `pos` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
CHTYPE GetQryRfnSS(int pos)
{
    return (STRUCTS == FLDS_STRUCTS_QRIES)? GetQuerySS(pos): GetDbStrSS(pos);
}

// -------------------------------------------------------------------------
// GetDbStrCoord: get coordinate at the given position of reference 
// structure;
//
template<int CRD>
__device__ __forceinline__
/**
 * @brief 在CUDA 刚体拟合与评分中读取 `GetDbStrCoord` 对应的数据。
 * @param pos 供该函数读取或更新的 `pos` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
FPTYPE GetDbStrCoord(int pos)
{
    return ((FPTYPE*)(dc_pm2dvfields_[
        ndx_dbs_dc_pm2dvfields_ + pmv2DCoords + CRD]))[pos];
}
__device__ __forceinline__
FPTYPE GetDbStrCoord(int crd, int pos)
{
    return ((FPTYPE*)(dc_pm2dvfields_[
        ndx_dbs_dc_pm2dvfields_ + pmv2DCoords + crd]))[pos];
}

// GetQueryCoord: get coordinate at the given position of a query;
//
template<int CRD>
__device__ __forceinline__
/**
 * @brief 在CUDA 刚体拟合与评分中读取 `GetQueryCoord` 对应的数据。
 * @param pos 供该函数读取或更新的 `pos` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
FPTYPE GetQueryCoord(int pos)
{
    return ((FPTYPE*)(dc_pm2dvfields_[
        ndx_qrs_dc_pm2dvfields_ + pmv2DCoords + CRD]))[pos];
}
__device__ __forceinline__
FPTYPE GetQueryCoord(int crd, int pos)
{
    return ((FPTYPE*)(dc_pm2dvfields_[
        ndx_qrs_dc_pm2dvfields_ + pmv2DCoords + crd]))[pos];
}

__device__ __forceinline__
/**
 * @brief 在CUDA 刚体拟合与评分中读取 `GetQryRfnCoord` 对应的数据。
 * @param fldsection 供该函数读取或更新的 `fldsection` 参数。
 * @param crd 供该函数读取或更新的 `crd` 参数。
 * @param pos 供该函数读取或更新的 `pos` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
FPTYPE GetQryRfnCoord(int fldsection, int crd, int pos)
{
    return ((FPTYPE*)(dc_pm2dvfields_[fldsection + pmv2DCoords + crd]))[pos];
}

// -------------------------------------------------------------------------
// SetDbStrSecStr: set secondary structure assignment code at the given 
// reference (db) structure position;
__device__ __forceinline__
/**
 * @brief 在CUDA 刚体拟合与评分中设置 `SetDbStrSecStr` 对应的数据。
 * @param pos 供该函数读取或更新的 `pos` 参数。
 * @param ss 供该函数读取或更新的 `ss` 参数。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
void SetDbStrSecStr(int pos, char ss)
{
    ((CHTYPE*)(dc_pm2dvfields_[ndx_dbs_dc_pm2dvfields_ + pmv2Dss]))[pos] = ss;
}

// SetQuerySecStr: set secondary structure assignment code at the given 
// query structure position;
__device__ __forceinline__
/**
 * @brief 在CUDA 刚体拟合与评分中设置 `SetQuerySecStr` 对应的数据。
 * @param pos 供该函数读取或更新的 `pos` 参数。
 * @param ss 供该函数读取或更新的 `ss` 参数。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
void SetQuerySecStr(int pos, char ss)
{
    ((CHTYPE*)(dc_pm2dvfields_[ndx_qrs_dc_pm2dvfields_ + pmv2Dss]))[pos] = ss;
}



// --- INDEX ---------------------------------------------------------------
// -------------------------------------------------------------------------
// GetIndxdDbStrField: get a field value at the given position of indexed
// reference (db) structure;
template<typename T, int field>
__device__ __forceinline__
/**
 * @brief 在CUDA 刚体拟合与评分中读取 `GetIndxdDbStrField` 对应的数据。
 * @param pos 供该函数读取或更新的 `pos` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
T GetIndxdDbStrField(int pos)
{
    return ((T*)(dc_pm2dvfields_[ndx_dbs_dc_pm2dvndxfds_ + field]))[pos];
}

// SsetIndxdDbStrField: set a field value at the given position of indexed
// reference (db) structure;
template<typename T, int field>
__device__ __forceinline__
/**
 * @brief 在CUDA 刚体拟合与评分中设置 `SetIndxdDbStrField` 对应的数据。
 * @param pos 供该函数读取或更新的 `pos` 参数。
 * @param value 需要读取、写入或转换的值。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
void SetIndxdDbStrField(int pos, T value)
{
    ((T*)(dc_pm2dvfields_[ndx_dbs_dc_pm2dvndxfds_ + field]))[pos] = value;
}



// -------------------------------------------------------------------------
// GetIndxdDbstrCoord: get coordinate at the given position of indexed 
// reference structure;
template<int CRD>
__device__ __forceinline__
/**
 * @brief 在CUDA 刚体拟合与评分中读取 `GetIndxdDbStrCoord` 对应的数据。
 * @param pos 供该函数读取或更新的 `pos` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
FPTYPE GetIndxdDbStrCoord(int pos)
{
    return ((FPTYPE*)(dc_pm2dvfields_[ndx_dbs_dc_pm2dvndxfds_ + pmv2DNdxCoords + CRD]))[pos];
}
// -------------------------------------------------------------------------
// GetIndxdQueryCoord: get coordinate at the given position of indexed query
// structure;
//
template<int CRD>
__device__ __forceinline__
/**
 * @brief 在CUDA 刚体拟合与评分中读取 `GetIndxdQueryCoord` 对应的数据。
 * @param pos 供该函数读取或更新的 `pos` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
FPTYPE GetIndxdQueryCoord(int pos)
{
    return ((FPTYPE*)(dc_pm2dvfields_[ndx_qrs_dc_pm2dvndxfds_ + pmv2DNdxCoords + CRD]))[pos];
}
__device__ __forceinline__
FPTYPE GetIndxdQueryCoord(int crd, int pos)
{
    return ((FPTYPE*)(dc_pm2dvfields_[ndx_qrs_dc_pm2dvndxfds_ + pmv2DNdxCoords + crd]))[pos];
}

// -------------------------------------------------------------------------
// GetIndxdDbStrOrgndx: get the original index at the given position of 
// indexed reference structure;
//
__device__ __forceinline__
/**
 * @brief 在CUDA 刚体拟合与评分中读取 `GetIndxdDbStrOrgndx` 对应的数据。
 * @param pos 供该函数读取或更新的 `pos` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
INTYPE GetIndxdDbStrOrgndx(int pos)
{
    return ((INTYPE*)(dc_pm2dvfields_[ndx_dbs_dc_pm2dvndxfds_ + pmv2DNdxOrgndx]))[pos];
}

// GetIndxdQueryOrgndx: get the original index at the given position of 
// indexed query structure;
//
__device__ __forceinline__
/**
 * @brief 在CUDA 刚体拟合与评分中读取 `GetIndxdQueryOrgndx` 对应的数据。
 * @param pos 供该函数读取或更新的 `pos` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
INTYPE GetIndxdQueryOrgndx(int pos)
{
    return ((INTYPE*)(dc_pm2dvfields_[ndx_qrs_dc_pm2dvndxfds_ + pmv2DNdxOrgndx]))[pos];
}

// -------------------------------------------------------------------------
// GetIndxdDbStrBranch: get branch (left/right) at the given position of 
// indexed reference structure;
//
__device__ __forceinline__
/**
 * @brief 在CUDA 刚体拟合与评分中读取 `GetIndxdDbStrBranchLeft` 对应的数据。
 * @param pos 供该函数读取或更新的 `pos` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
INTYPE GetIndxdDbStrBranchLeft(int pos)
{
    return ((INTYPE*)(dc_pm2dvfields_[ndx_dbs_dc_pm2dvndxfds_ + pmv2DNdxLeft]))[pos];
}
__device__ __forceinline__
/**
 * @brief 在CUDA 刚体拟合与评分中读取 `GetIndxdDbStrBranchRight` 对应的数据。
 * @param pos 供该函数读取或更新的 `pos` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
INTYPE GetIndxdDbStrBranchRight(int pos)
{
    return ((INTYPE*)(dc_pm2dvfields_[ndx_dbs_dc_pm2dvndxfds_ + pmv2DNdxRight]))[pos];
}

// -------------------------------------------------------------------------
// GetIndxdQueryBranch: get branch (left/right) at the given position of 
// indexed query structure;
// NOTE: BRANCH should be pmv2DNdxLeft or pmv2DNdxRight;
//
template<int BRANCH>
__device__ __forceinline__
/**
 * @brief 在CUDA 刚体拟合与评分中读取 `GetIndxdQueryBranch` 对应的数据。
 * @param pos 供该函数读取或更新的 `pos` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
INTYPE GetIndxdQueryBranch(int pos)
{
    return ((INTYPE*)(dc_pm2dvfields_[ndx_qrs_dc_pm2dvndxfds_ + BRANCH]))[pos];
}
__device__ __forceinline__
INTYPE GetIndxdQueryBranch(int branch, int pos)
{
    return ((INTYPE*)(dc_pm2dvfields_[ndx_qrs_dc_pm2dvndxfds_ + branch]))[pos];
}

__device__ __forceinline__
/**
 * @brief 在CUDA 刚体拟合与评分中读取 `GetIndxdQueryBranchLeft` 对应的数据。
 * @param pos 供该函数读取或更新的 `pos` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
INTYPE GetIndxdQueryBranchLeft(int pos)
{
    return ((INTYPE*)(dc_pm2dvfields_[ndx_qrs_dc_pm2dvndxfds_ + pmv2DNdxLeft]))[pos];
}
__device__ __forceinline__
/**
 * @brief 在CUDA 刚体拟合与评分中读取 `GetIndxdQueryBranchRight` 对应的数据。
 * @param pos 供该函数读取或更新的 `pos` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
INTYPE GetIndxdQueryBranchRight(int pos)
{
    return ((INTYPE*)(dc_pm2dvfields_[ndx_qrs_dc_pm2dvndxfds_ + pmv2DNdxRight]))[pos];
}

// -------------------------------------------------------------------------
// GetGonnetScore: get a pairwise residue score obtained from the Gonnet
// frequencies
//
__device__ __forceinline__
/**
 * @brief 在CUDA 刚体拟合与评分中读取 `GetGonnetScore` 对应的数据。
 * @param res1 供该函数读取或更新的 `res1` 参数。
 * @param res2 供该函数读取或更新的 `res2` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
float GetGonnetScore(char res1, char res2)
{
    if('A' <= res1 && res1 <= 'Z' && 'A' <= res2 && res2 <= 'Z')
        return dc_Gonnet_scores_[(res1 - 'A') * NEA + (res2 - 'A')];
    return 0.0f;
}

// -------------------------------------------------------------------------

#endif//__fields_cuh__
