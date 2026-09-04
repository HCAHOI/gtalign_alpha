/***************************************************************************
 *   Copyright (C) 2021-2023 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include <stdio.h>
#include <string>

#include "platform.h"
#include "mylimits.h"
#include "myexception.h"

// -------------------------------------------------------------------------
// CLASS myruntime_error
//
// pretty_format: place exception information on a string object
//
/**
 * @brief 在通用工具中处理 `myruntime_error::pretty_format` 对应的数据。
 * @param preamb 供该函数读取或更新的 `preamb` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
std::string myruntime_error::pretty_format( std::string preamb ) const throw()
{
    char buf[BUF_MAX];
    bool emp = _errdsc.empty();
    if( emp )
        preamb += myexception::what();
    else
        preamb += _errdsc;
    int bfile = 0, bline = 0, bfunc = 0;
    if( file())
        bfile = 1;
    if( line())
        bline = 2;
    if( function())
        bfunc = 4;
    if( bfile | bline | bfunc ) {
        preamb += NL;
        preamb += "    (";
        if( file()) {
            preamb += "File: ";
            preamb += file();
            if( bline | bfunc ) {
                preamb += "; ";
            }
        }
        if( line()) {
            preamb += "Line: ";
            sprintf( buf, "%u", line());
            preamb += buf;
            if( bfunc ) {
                preamb += "; ";
            }
        }
        if( function()) {
            preamb += "Function: ";
            preamb += function();
        }
        preamb += ")";
    }
    return preamb;
}

