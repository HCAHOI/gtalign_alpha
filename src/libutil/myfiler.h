/***************************************************************************
 *   Copyright (C) 2021-2023 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __myfiler_h__
#define __myfiler_h__

#include <stdio.h>
#include <string>

#include <sys/types.h>
#include <sys/stat.h>
#include "platform.h"

#ifdef OS_MS_WINDOWS
#	include <Windows.h>
#endif

struct TCharStream {
    /**
     * @brief 构造 `TCharStream`，初始化其负责的通用工具状态。
     * @par 参数
     * 无。
     * @return 无返回值；完成对象构造与初始状态设置。
     */
    TCharStream()
    :
#ifdef OS_MS_WINDOWS
        hMapFile_(NULL),
#endif
        data_(NULL), datlen_(0), curpos_(0),
        pagenr_(0), pagesize_(0), pageoff_(0)
    {}
    /**
     * @brief 构造 `TCharStream`，初始化其负责的通用工具状态。
     * @param hmapfile 供该函数读取或更新的 `hmapfile` 参数。
     * @param data 供该函数读取或更新的 `data` 参数。
     * @param datlen 控制当前步骤范围或规模的 `datlen`。
     * @param curpos 供该函数读取或更新的 `curpos` 参数。
     * @param pagenr 供该函数读取或更新的 `pagenr` 参数。
     * @param pagesize 控制当前步骤范围或规模的 `pagesize`。
     * @param pageoff 供该函数读取或更新的 `pageoff` 参数。
     * @return 无返回值；完成对象构造与初始状态设置。
     */
    TCharStream( 
#ifdef OS_MS_WINDOWS
        HANDLE hmapfile,
#endif
        char* data, size_t datlen, size_t curpos,
        size_t pagenr, size_t pagesize, size_t pageoff)
    :
#ifdef OS_MS_WINDOWS
        hMapFile_(hmapfile),
#endif
        data_(data), datlen_(datlen), curpos_(curpos),
        pagenr_(pagenr), pagesize_(pagesize), pageoff_(pageoff)
    {}
    /**
     * @brief 在通用工具中处理 `operator=` 对应的数据。
     * @param chstr 供该函数读取或更新的 `chstr` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    TCharStream& operator=(const TCharStream& chstr) {
        datlen_ = chstr.datlen_;
        curpos_ = chstr.curpos_;
        pagenr_ = chstr.pagenr_;
        pageoff_ = chstr.pageoff_;
        return *this;
    }
    /**
     * @brief 在通用工具中重置 `ResetCounters` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void ResetCounters() {
        datlen_ = 0;
        curpos_ = 0;
        pagenr_ = 0;
        pageoff_ = 0;
    }
    /**
     * @brief 在通用工具中处理 `incpos` 对应的数据。
     * @param by 供该函数读取或更新的 `by` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void incpos(size_t by) {curpos_ += by;}
    /**
     * @brief 在通用工具中处理 `incposnl` 对应的数据。
     * @param by 供该函数读取或更新的 `by` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void incposnl(size_t by) {
        curpos_ += by;
        const char* pbeg = data_ + curpos_;
        const char* pend = data_ + datlen_;
        const char* p = pbeg;
        for(; p < pend && *p != '\n' && *p != '\r'; p++);
        for(; p < pend && (*p == '\n' || *p == '\r'); p++);
        curpos_ += (size_t)(p-pbeg);
    }
    /**
     * @brief 在通用工具中处理 `eof` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    bool eof() const {return datlen_ <= curpos_;}
    //
#ifdef OS_MS_WINDOWS
    HANDLE hMapFile_;//handle of the file mapping object
#endif
    char* data_;//data
    size_t datlen_;//data size
    size_t curpos_;//next position in the current data to read
    size_t pagenr_;//number of the page where the data begins
    size_t pagesize_;//page size in bytes
    size_t pageoff_;//offset of page pagenr_
};

// file/directory routines
/**
 * @brief 在通用工具中处理 `file_exists` 对应的数据。
 * @param param1 供该函数读取或更新的 `param1` 参数。
 * @param mode 供该函数读取或更新的 `mode` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
bool file_exists( const char*,
#ifdef OS_MS_WINDOWS
    unsigned short mode = _S_IFREG
#else
    mode_t mode = S_IFREG
#endif
);
/**
 * @brief 在通用工具中处理 `directory_exists` 对应的数据。
 * @param pathname 供该函数读取或更新的 `pathname` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline
bool directory_exists( const char* pathname ) {
    return file_exists( pathname, 
#ifdef OS_MS_WINDOWS
    _S_IFDIR
#else
    S_IFDIR
#endif
    );
}
/**
 * @brief 在通用工具中处理 `file_size` 对应的数据。
 * @param param1 供该函数读取或更新的 `param1` 参数。
 * @param param2 供该函数读取或更新的 `param2` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
int file_size(const char*, size_t*);
/**
 * @brief 在通用工具中处理 `mymkdir` 对应的数据。
 * @param pathname 供该函数读取或更新的 `pathname` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
int mymkdir(const char* pathname);


//functor for reading from a stream into the string buffer
struct CopyToBuffer {
    /**
     * @brief 在通用工具中处理 `operator()` 对应的数据。
     * @param srcptr 提供输入数据的 `srcptr`。
     * @param srclen 控制当前步骤范围或规模的 `srclen`。
     * @param buf 供当前步骤读取或更新的 `buf` 缓冲区。
     * @param param4 供该函数读取或更新的 `param4` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void operator()(char* srcptr, size_t srclen, std::string* buf, size_t*) {
        *buf = std::string(srcptr, srclen);
    }
};
//functor for reading from a stream and setting the address where the data begins
struct GetPtr {
    /**
     * @brief 在通用工具中处理 `operator()` 对应的数据。
     * @param srcptr 提供输入数据的 `srcptr`。
     * @param srclen 控制当前步骤范围或规模的 `srclen`。
     * @param outptr 接收当前步骤输出的 `outptr`。
     * @param outlen 控制当前步骤范围或规模的 `outlen`。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void operator()(char* srcptr, size_t srclen, char** outptr, size_t* outlen) {
        *outptr = srcptr;
        *outlen = srclen;
    }
};


/**
 * @brief 在通用工具中处理 `skip_comments` 对应的数据。
 * @param fp 供该函数读取或更新的 `fp` 参数。
 * @param buffer 供当前步骤读取或更新的 `buffer` 缓冲区。
 * @param bufsize 控制当前步骤范围或规模的 `bufsize`。
 * @param readlen 控制当前步骤范围或规模的 `readlen`。
 * @param cc 供该函数读取或更新的 `cc` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
int skip_comments( FILE* fp, char* buffer, size_t bufsize, size_t* readlen, char cc = '#');
/**
 * @brief 在通用工具中处理 `skip_comments` 对应的数据。
 * @param fp 供该函数读取或更新的 `fp` 参数。
 * @param buffer 供当前步骤读取或更新的 `buffer` 缓冲区。
 * @param cc 供该函数读取或更新的 `cc` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
int skip_comments( FILE* fp, std::string& buffer, char cc = '#');
/**
 * @brief 在通用工具中处理 `skip_comments` 对应的数据。
 * @param source 供该函数读取或更新的 `source` 参数。
 * @param buffer 供当前步骤读取或更新的 `buffer` 缓冲区。
 * @param cc 供该函数读取或更新的 `cc` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
int skip_comments( TCharStream* source, std::string& buffer, char cc  = '#');
/**
 * @brief 在通用工具中处理 `skip_comments` 对应的数据。
 * @param source 供该函数读取或更新的 `source` 参数。
 * @param ptr 供该函数读取或更新的 `ptr` 参数。
 * @param param3 供该函数读取或更新的 `param3` 参数。
 * @param len 控制当前步骤范围或规模的 `len`。
 * @param cc 供该函数读取或更新的 `cc` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
int skip_comments( TCharStream* source, char*& ptr, size_t, size_t* len, char cc  = '#');
/**
 * @brief 在通用工具中读取 `read_double` 对应的数据。
 * @param readfrom 供该函数读取或更新的 `readfrom` 参数。
 * @param readlen 控制当前步骤范围或规模的 `readlen`。
 * @param membuf 供当前步骤读取或更新的 `membuf` 缓冲区。
 * @param rbytes 供该函数读取或更新的 `rbytes` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
int read_double( const char* readfrom, size_t readlen, double* membuf, size_t* rbytes );
/**
 * @brief 在通用工具中读取 `read_float` 对应的数据。
 * @param readfrom 供该函数读取或更新的 `readfrom` 参数。
 * @param readlen 控制当前步骤范围或规模的 `readlen`。
 * @param membuf 供当前步骤读取或更新的 `membuf` 缓冲区。
 * @param rbytes 供该函数读取或更新的 `rbytes` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
int read_float( const char* readfrom, size_t readlen, float* membuf, size_t* rbytes );
/**
 * @brief 在通用工具中读取 `read_integer` 对应的数据。
 * @param readfrom 供该函数读取或更新的 `readfrom` 参数。
 * @param readlen 控制当前步骤范围或规模的 `readlen`。
 * @param membuf 供当前步骤读取或更新的 `membuf` 缓冲区。
 * @param rbytes 供该函数读取或更新的 `rbytes` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
int read_integer( const char* readfrom, size_t readlen, int* membuf, size_t* rbytes );
/**
 * @brief 在通用工具中读取 `read_llinteger` 对应的数据。
 * @param readfrom 供该函数读取或更新的 `readfrom` 参数。
 * @param readlen 控制当前步骤范围或规模的 `readlen`。
 * @param membuf 供当前步骤读取或更新的 `membuf` 缓冲区。
 * @param rbytes 供该函数读取或更新的 `rbytes` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
int read_llinteger( const char* readfrom, size_t readlen, long long int* membuf, size_t* rbytes );
/**
 * @brief 在通用工具中读取 `read_symbol` 对应的数据。
 * @param readfrom 供该函数读取或更新的 `readfrom` 参数。
 * @param readlen 控制当前步骤范围或规模的 `readlen`。
 * @param membuf 供当前步骤读取或更新的 `membuf` 缓冲区。
 * @param rbytes 供该函数读取或更新的 `rbytes` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
int read_symbol( const char* readfrom, size_t readlen, char* membuf, size_t* rbytes );

/**
 * @brief 在通用工具中处理 `feof` 对应的数据。
 * @param chs 供该函数读取或更新的 `chs` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline int feof(TCharStream* chs ) {
    if( chs == NULL || chs->data_ == NULL || 
        chs->datlen_ < 1 || chs->datlen_ <= chs->curpos_ )
        return 1;
    return 0;
}

#define ERR_SC_MACC ( 111 )
#define ERR_SC_READ ( 113 )
#define ERR_RD_MACC ( 115 )
#define ERR_RD_NOVL ( 117 )
#define ERR_RD_INVL ( 119 )
#define ERR_RI_MACC ( 131 )
#define ERR_RI_NOVL ( 133 )
#define ERR_RI_INVL ( 135 )
#define ERR_RL_MACC ( 137 )
#define ERR_RL_NOVL ( 139 )
#define ERR_RL_INVL ( 151 )
#define ERR_RS_MACC ( 153 )
#define ERR_RS_INVL ( 155 )

//error messages
/**
 * @brief 在通用工具中处理 `TranslateReadError` 对应的数据。
 * @param code 供该函数读取或更新的 `code` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline const char* TranslateReadError( int code )
{
    switch( code ) {
        case 0: return "OK";
        case ERR_SC_MACC: return "skip_comments: Memory access error";
        case ERR_SC_READ: return "skip_comments: Reading error";
        case ERR_RD_MACC: return "read_double: Memory access error";
        case ERR_RD_NOVL: return "read_double: No double value read";
        case ERR_RD_INVL: return "read_double: Invalid double value";
        case ERR_RI_MACC: return "read_integer: Memory access error";
        case ERR_RI_NOVL: return "read_integer: No integer read";
        case ERR_RI_INVL: return "read_integer: Invalid integer value";
        case ERR_RL_MACC: return "read_llinteger: Memory access error";
        case ERR_RL_NOVL: return "read_llinteger: No integer read";
        case ERR_RL_INVL: return "read_llinteger: Invalid integer value";
        case ERR_RS_MACC: return "read_symbol: Memory access error";
        case ERR_RS_INVL: return "read_symbol: Not a single character";
    }
    return "Unknown";
}

#endif//__myfiler_h__
