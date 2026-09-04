/***************************************************************************
 *   Copyright (C) 2021-2023 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __mygetopt_h__
#define __mygetopt_h__

#include <string>
#include <vector>

enum TArgFlag {
    my_required_argument,
    my_optional_argument,
    my_no_argument,
    my_n_targflags
};

// predefined return values
enum TRetValues {
    my_gorv_term = -1,//terminated processing of command line
    my_gorv_value = '*',//value, not option
    my_gorv_noarg = ':',//argument missing for option
    my_gorv_notfound = '?',//option not found
    my_gorv_illformed = '!',//ill-formed option
};

struct myoption {
    const char* sLongname_;//option's name
    TArgFlag    eFlag_;//option's property
    int         nRet_;//return value
};

// _________________________________________________________________________
// CLASS MyGetopt
// for parsing command line
//
class MyGetopt {
public:
    /**
     * @brief 构造 `MyGetopt`，初始化其负责的通用工具状态。
     * @param myoption 供该函数读取或更新的 `myoption` 参数。
     * @param argv 命令行参数字符串数组。
     * @param argc 命令行参数个数。
     * @return 无返回值；完成对象构造与初始状态设置。
     */
    MyGetopt( const myoption*, const char *argv[], int argc );
    /**
     * @brief 销毁 `MyGetopt`，释放其管理的通用工具资源。
     * @par 参数
     * 无。
     * @return 无返回值；对象持有的资源在返回前完成释放。
     */
    ~MyGetopt();

    /**
     * @brief 在通用工具中读取 `GetNextOption` 对应的数据。
     * @param argument 供该函数读取或更新的 `argument` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    int GetNextOption( std::string* argument );

protected:
    /**
     * @brief 在通用工具中初始化 `Init` 对应的数据。
     * @param myoption 供该函数读取或更新的 `myoption` 参数。
     * @param argv 命令行参数字符串数组。
     * @param argc 命令行参数个数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void Init( const myoption*, const char *argv[], int argc );
    /**
     * @brief 在通用工具中重置 `Reset` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void Reset() { ncnt_ = 0; stopped_ = false; };

private:
    int ncnt_;//word counter in options; max value is argc
    bool stopped_;//flag of stopping processing command line
    const char** argv_;//command line
    int argc_;//number of words in command line
    std::vector<myoption> srtoptions_;//sorted options
};

// --- INLINES -------------------------------------------------------------
// comparison operators for two option names
//
/**
 * @brief 在通用工具中处理 `operator<` 对应的数据。
 * @param left 供该函数读取或更新的 `left` 参数。
 * @param right 供该函数读取或更新的 `right` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline
bool operator<(const myoption& left, const myoption& right )
{
    if( left.sLongname_ == NULL || right.sLongname_ == NULL )
        return false;
    return strcmp(left.sLongname_, right.sLongname_) < 0;
}

/**
 * @brief 在通用工具中处理 `operator==` 对应的数据。
 * @param left 供该函数读取或更新的 `left` 参数。
 * @param right 供该函数读取或更新的 `right` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline
bool operator==(const myoption& left, const myoption& right )
{
    if( left.sLongname_ == NULL || right.sLongname_ == NULL )
        return false;
    return strcmp(left.sLongname_, right.sLongname_) == 0;
}

// -------------------------------------------------------------------------

#endif//__mygetopt_h__
