/***************************************************************************
 *   Copyright (C) 2021-2023 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __myexception_h__
#define __myexception_h__

#include <string>
#include "debug.h"

#define MYRUNTIME_ERROR( MSG )  myruntime_error( MSG, __EXCPOINT__ )
#define MYRUNTIME_ERROR2( MSG, CLS )  myruntime_error( MSG, __EXCPOINT__, CLS )

#define TRY try {
#define CATCH_ERROR_RETURN(STATEMENTS) \
        STATEMENTS; \
    } catch( myruntime_error const& ex ) { \
        STATEMENTS; \
        error( ex.pretty_format().c_str()); \
        return EXIT_FAILURE; \
    } catch( myexception const& ex ) { \
        STATEMENTS; \
        error( ex.what()); \
        return EXIT_FAILURE; \
    } catch( ... ) { \
        STATEMENTS; \
        error("Unknown exception caught."); \
        return EXIT_FAILURE; \
    }

enum {
    NOCLASS,
    SCALING,
    CRITICAL
};

// _________________________________________________________________________
// CLASS myexception
// for exception handling
//
class myexception
{
public:
    /**
     * @brief 构造 `myexception`，初始化其负责的通用工具状态。
     * @par 参数
     * 无。
     * @return 无返回值；完成对象构造与初始状态设置。
     */
    myexception() throw() {}
    /**
     * @brief 构造 `myexception`，初始化其负责的通用工具状态。
     * @param ex 供该函数读取或更新的 `ex` 参数。
     * @return 无返回值；完成对象构造与初始状态设置。
     */
    myexception( const myexception& ex ) throw() { operator=(ex);};
    /**
     * @brief 销毁 `myexception`，释放其管理的通用工具资源。
     * @par 参数
     * 无。
     * @return 无返回值；对象持有的资源在返回前完成释放。
     */
    virtual ~myexception() throw() {};
    /**
     * @brief 在通用工具中处理 `operator=` 对应的数据。
     * @param myexception 供该函数读取或更新的 `myexception` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    virtual myexception& operator=( const myexception& ) throw() { return *this;};
    /**
     * @brief 在通用工具中处理 `what` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    virtual const char* what() const throw();
    /**
     * @brief 在通用工具中处理 `eclass` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    virtual int         eclass() const throw();
};

// _________________________________________________________________________
// CLASS myruntime_error
// runtime error exception
//
class myruntime_error: public myexception
{
public:
    /**
     * @brief 构造 `myruntime_error`，初始化其负责的通用工具状态。
     * @par 参数
     * 无。
     * @return 无返回值；完成对象构造与初始状态设置。
     */
    myruntime_error() throw();
    /**
     * @brief 构造 `myruntime_error`，初始化其负责的通用工具状态。
     * @param ex 供该函数读取或更新的 `ex` 参数。
     * @return 无返回值；完成对象构造与初始状态设置。
     */
    myruntime_error( const myexception& ex ) throw(){ operator=(ex);};
    /**
     * @brief 构造 `myruntime_error`，初始化其负责的通用工具状态。
     * @param mre 供该函数读取或更新的 `mre` 参数。
     * @return 无返回值；完成对象构造与初始状态设置。
     */
    myruntime_error( const myruntime_error& mre ) throw(): myexception(mre) { operator=(mre);};
    /**
     * @brief 构造 `myruntime_error`，初始化其负责的通用工具状态。
     * @param arg 供该函数读取或更新的 `arg` 参数。
     * @param ecl 供该函数读取或更新的 `ecl` 参数。
     * @return 无返回值；完成对象构造与初始状态设置。
     */
    explicit myruntime_error( const std::string& arg, int ecl = NOCLASS ) throw();
    /**
     * @brief 构造 `myruntime_error`，初始化其负责的通用工具状态。
     * @param arg 供该函数读取或更新的 `arg` 参数。
     * @param fl 供该函数读取或更新的 `fl` 参数。
     * @param ln 供该函数读取或更新的 `ln` 参数。
     * @param func 供该函数读取或更新的 `func` 参数。
     * @param ecl 供该函数读取或更新的 `ecl` 参数。
     * @return 无返回值；完成对象构造与初始状态设置。
     */
    myruntime_error( const std::string& arg, 
            const char* fl, unsigned int ln, const char* func, int ecl = NOCLASS ) throw();
    /**
     * @brief 构造 `myruntime_error`，初始化其负责的通用工具状态。
     * @param fl 供该函数读取或更新的 `fl` 参数。
     * @param ln 供该函数读取或更新的 `ln` 参数。
     * @param func 供该函数读取或更新的 `func` 参数。
     * @param ecl 供该函数读取或更新的 `ecl` 参数。
     * @return 无返回值；完成对象构造与初始状态设置。
     */
    myruntime_error( const char* fl, unsigned int ln, const char* func, int ecl = NOCLASS ) throw();
    /**
     * @brief 销毁 `myruntime_error`，释放其管理的通用工具资源。
     * @par 参数
     * 无。
     * @return 无返回值；对象持有的资源在返回前完成释放。
     */
    virtual ~myruntime_error() throw();

    /**
     * @brief 在通用工具中处理 `operator=` 对应的数据。
     * @param myexception 供该函数读取或更新的 `myexception` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    virtual myexception& operator=( const myexception& ) throw();
    /**
     * @brief 在通用工具中处理 `operator=` 对应的数据。
     * @param myruntime_error 供该函数读取或更新的 `myruntime_error` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    virtual myruntime_error& operator=( const myruntime_error& ) throw();

    /**
     * @brief 在通用工具中处理 `isset` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    virtual bool isset() const throw();

    /**
     * @brief 在通用工具中处理 `what` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    virtual const char*     what() const throw();       //cause of the error
    /**
     * @brief 在通用工具中处理 `eclass` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    virtual int             eclass() const throw();     //error class
    /**
     * @brief 在通用工具中处理 `file` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    virtual const char*     file() const throw();       //filename
    /**
     * @brief 在通用工具中处理 `line` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    virtual unsigned int    line() const throw();       //line number
    /**
     * @brief 在通用工具中处理 `function` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    virtual const char*     function() const throw();   //function name

    /**
     * @brief 在通用工具中处理 `pretty_format` 对应的数据。
     * @param string 供该函数读取或更新的 `string` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    virtual std::string     pretty_format( std::string = std::string()) const throw();

private:
    std::string     _errdsc;    //error description
    int             _class;     //error class
    const char*     _file;      //filename
    unsigned int    _line;      //line number
    const char*     _function;  //function name
};

// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
// CLASS myexception implementation
//
// what: error description
//
/**
 * @brief 在通用工具中处理 `myexception::what` 对应的数据。
 * @par 参数
 * 无。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline const char* myexception::what() const throw()
{
    return "{myexception}";
}

// eclass: error class
//
/**
 * @brief 在通用工具中处理 `myexception::eclass` 对应的数据。
 * @par 参数
 * 无。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline int myexception::eclass() const throw()
{
    return NOCLASS;
}

// -------------------------------------------------------------------------
// CLASS myruntime_error
//
// Constructors
//
/**
 * @brief 构造 `myruntime_error`，初始化其负责的通用工具状态。
 * @par 参数
 * 无。
 * @return 无返回值；完成对象构造与初始状态设置。
 */
inline myruntime_error::myruntime_error() throw()
:   myexception(),
    _errdsc(),
    _class( NOCLASS ),
    _file( NULL ),
    _line(0),
    _function( NULL )
{
}

/**
 * @brief 构造 `myruntime_error`，初始化其负责的通用工具状态。
 * @param arg 供该函数读取或更新的 `arg` 参数。
 * @param ecldesc 供该函数读取或更新的 `ecldesc` 参数。
 * @return 无返回值；完成对象构造与初始状态设置。
 */
inline myruntime_error::myruntime_error( const std::string& arg, int ecldesc ) throw()
:   myexception(),
    _errdsc( arg ),
    _class( ecldesc ),
    _file( NULL ),
    _line(0),
    _function( NULL )
{
}

/**
 * @brief 构造 `myruntime_error`，初始化其负责的通用工具状态。
 * @param arg 供该函数读取或更新的 `arg` 参数。
 * @param fl 供该函数读取或更新的 `fl` 参数。
 * @param ln 供该函数读取或更新的 `ln` 参数。
 * @param func 供该函数读取或更新的 `func` 参数。
 * @param ecldesc 供该函数读取或更新的 `ecldesc` 参数。
 * @return 无返回值；完成对象构造与初始状态设置。
 */
inline myruntime_error::myruntime_error( const std::string& arg, 
    const char* fl, unsigned int ln, const char* func, int ecldesc ) throw()
:   myexception(),
    _errdsc( arg ),
    _class( ecldesc ),
    _file( fl ),
    _line( ln ),
    _function( func )
{
}

/**
 * @brief 构造 `myruntime_error`，初始化其负责的通用工具状态。
 * @param fl 供该函数读取或更新的 `fl` 参数。
 * @param ln 供该函数读取或更新的 `ln` 参数。
 * @param func 供该函数读取或更新的 `func` 参数。
 * @param ecldesc 供该函数读取或更新的 `ecldesc` 参数。
 * @return 无返回值；完成对象构造与初始状态设置。
 */
inline myruntime_error::myruntime_error( 
    const char* fl, unsigned int ln, const char* func, int ecldesc ) throw()
:   myexception(),
    _errdsc(),
    _class( ecldesc ),
    _file( fl ),
    _line( ln ),
    _function( func )
{
}

// Destructor
//
/**
 * @brief 销毁 `myruntime_error`，释放其管理的通用工具资源。
 * @par 参数
 * 无。
 * @return 无返回值；对象持有的资源在返回前完成释放。
 */
inline myruntime_error::~myruntime_error() throw()
{
}

// -------------------------------------------------------------------------
// operator=: assignment operators
//
/**
 * @brief 在通用工具中处理 `myruntime_error::operator=` 对应的数据。
 * @param me 供该函数读取或更新的 `me` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline myexception& myruntime_error::operator=( const myexception& me ) throw()
{
    const myruntime_error* pme = dynamic_cast<const myruntime_error*>(&me);
    if( pme ) {
        operator=(*pme);
    }
    else {
        _errdsc = myexception::what();
        _class = myexception::eclass();
        _file = NULL;
        _line = 0;
        _function = NULL;
    }
    return *this;
}

/**
 * @brief 在通用工具中处理 `myruntime_error::operator=` 对应的数据。
 * @param mre 供该函数读取或更新的 `mre` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline myruntime_error& myruntime_error::operator=( const myruntime_error& mre ) throw()
{
    _errdsc = mre._errdsc;
    _class = mre.eclass();
    _file = mre.file();
    _line = mre.line();
    _function = mre.function();
    return *this;
}

// -------------------------------------------------------------------------
// isset: check whether the error is set
//
/**
 * @brief 在通用工具中处理 `myruntime_error::isset` 对应的数据。
 * @par 参数
 * 无。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline bool myruntime_error::isset() const throw()
{
    return !_errdsc.empty() || _class != NOCLASS || 
          _file || _line || _function;
}

// -------------------------------------------------------------------------
// what: error description
//
/**
 * @brief 在通用工具中处理 `myruntime_error::what` 对应的数据。
 * @par 参数
 * 无。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline const char* myruntime_error::what() const throw()
{
    return _errdsc.c_str();
}
// -------------------------------------------------------------------------
// eclass: error class 
//
/**
 * @brief 在通用工具中处理 `myruntime_error::eclass` 对应的数据。
 * @par 参数
 * 无。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline int myruntime_error::eclass() const throw()
{
    return _class;
}
// -------------------------------------------------------------------------
// file: filename
//
/**
 * @brief 在通用工具中处理 `myruntime_error::file` 对应的数据。
 * @par 参数
 * 无。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline const char* myruntime_error::file() const throw()
{
    return _file;
}
// -------------------------------------------------------------------------
// line: line number
//
/**
 * @brief 在通用工具中处理 `myruntime_error::line` 对应的数据。
 * @par 参数
 * 无。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline unsigned int myruntime_error::line() const throw()
{
    return _line;
}
// -------------------------------------------------------------------------
// function: function name
//
/**
 * @brief 在通用工具中处理 `myruntime_error::function` 对应的数据。
 * @par 参数
 * 无。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline const char* myruntime_error::function() const throw()
{
    return _function;
}


#endif//__myexception_h__
