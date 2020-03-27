#pragma once

#include <stdexcept>

class StackOverflowException : public std::exception
{
private:
    std::string msg;

public:
    StackOverflowException(const std::string& message = "") : msg(message)
    {
    }

    const char * what() const throw()
    {
        return msg.c_str();
    }
};

class OutOfMemoryException : public std::exception
{
private:
    std::string msg;

public:
    OutOfMemoryException(const std::string& message = "") : msg(message)
    {
    }

    const char * what() const throw()
    {
        return msg.c_str();
    }
};

class ThreadAbortException : public std::exception
{
private:
    std::string msg;

public:
    ThreadAbortException(const std::string& message = "") : msg(message)
    {
    }

    const char * what() const throw()
    {
        return msg.c_str();
    }
};

class DivideByZeroException : public std::exception
{
private:
    std::string msg;

public:
    DivideByZeroException(const std::string& message = "") : msg(message)
    {
    }

    const char * what() const throw()
    {
        return msg.c_str();
    }
};
