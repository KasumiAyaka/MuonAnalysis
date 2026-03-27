#ifndef CUF_KEYWORD_ARGS_H
#define CUF_KEYWORD_ARGS_H

#include <limits>
#include <tuple>
#include <string>

namespace l2c
{

#ifdef _MSC_VER
class Exception : public std::exception
{
public:
	//Exception(int e) noexcept : std::exception(ToHexString(e).c_str()) {}
	Exception(const char* c) noexcept : std::exception(c) {}
	Exception(const std::string& message)  noexcept : std::exception(message.c_str()) {}
	Exception(std::string&& message)  noexcept : std::exception(std::move(message).c_str()) {}
	virtual ~Exception() = default;

	virtual std::string GetErrorMessage() const { return what(); }
};
#else
class Exception : public std::exception
{
public:
	//Exception(int e)  noexcept : mMessage(ToHexString(e)) {}
	Exception(const char* c) noexcept : mMessage(c) {}
	Exception(const std::string& message) noexcept : mMessage(message) {}
	Exception(std::string&& message) noexcept : mMessage(std::move(message)) {}
	virtual ~Exception() = default;

	virtual const char* what() const { return mMessage.c_str(); }
	std::string GetErrorMessage() const { return mMessage; }

protected:
	std::string mMessage;
};
#endif

class InvalidArg : public Exception
{
public:
	using Exception::Exception;
	virtual std::string GetErrorMessage() const override { return std::string("INVALID ARG : ") + what(); }
};

template <class Derived, template <class ...> class Base>
class IsBasedOn_T
{
private:
	typedef char  Yes;
	typedef struct { char c[2]; } No;

	template <class ...U>
	static constexpr Yes check(const Base<U...>&);
	static constexpr No  check(...);

	static const Derived& d;
public:
	static constexpr bool value = sizeof(check(d)) == sizeof(Yes);
};
template <bool ...BoolSequence>
struct AndOperationSeq;
template <bool Head, bool ...Body>
struct AndOperationSeq<Head, Body...>
{
	static constexpr bool value = Head && AndOperationSeq<Body...>::value;
};
template <bool Head>
struct AndOperationSeq<Head>
{
	static constexpr bool value = Head;
};
template <>
struct AndOperationSeq<>
{
	static constexpr bool value = false;
};

class AnyTypeKeyword {};

namespace detail
{

template <class Name_, class Type_, class Tag_>
struct KeywordValue
{
	using Name = Name_;
	using Type = Type_;
	using Tag = Tag_;

	KeywordValue(Type v) : mValue(std::forward<Type>(v)) {}
	Type Get() { return std::forward<Type>(mValue); }
	template <class T>
	constexpr bool Is() const { return std::is_same<T, Type>::value; }
private:
	Type mValue;
};

template <class Name, class Type, class Tag>
struct KeywordName;
template <class Name_, class Type_, class Tag_>
struct KeywordName
{
	using Name = Name_;
	using Type = Type_;
	using Tag = Tag_;
	using Value = KeywordValue<Name, Type, Tag>;

	constexpr KeywordName() {}
	constexpr Value operator=(Type v) const { return Value(std::forward<Type>(v)); }
};
template <class Name_, class Tag_>
struct KeywordName<Name_, AnyTypeKeyword, Tag_>
{
	using Name = Name_;
	using Tag = Tag_;

	constexpr KeywordName() {}
	template <class Type>
	constexpr KeywordValue<Name, Type, Tag> operator=(Type&& v) const
	{
		return KeywordValue<Name, Type, Tag>;
		(std::forward<Type>(v));
	}
};
template <class Name_, class Tag_>
struct KeywordName<Name_, bool, Tag_>
{
	using Name = Name_;
	using Type = bool;
	using Tag = Tag_;
	using Value = KeywordValue<Name, bool, Tag>;

	static constexpr bool Get() { return true; }
	constexpr Value operator=(bool v) const { return Value(v); }
};

template <class Keyword, class Arg1,
	std::enable_if_t<!std::is_constructible<typename Keyword::Type, Arg1>::value, std::nullptr_t> = nullptr>
	std::remove_reference_t<typename Keyword::Type> GetDefault(Keyword /*keyword*/, Arg1&& /*arg1*/)
{
	throw InvalidArg("Default value does not exist.");
}
template <class Keyword, class Arg1,
	std::enable_if_t<std::is_constructible<typename Keyword::Type, Arg1>::value, std::nullptr_t> = nullptr>
	std::remove_reference_t<typename Keyword::Type> GetDefault(Keyword /*keyword*/, Arg1&& arg1)
{
	return std::forward<Arg1>(arg1);
}
template <class Keyword, class Arg1, class ...Args,
	bool B = (sizeof...(Args) != 0),
	std::enable_if_t<B, std::nullptr_t> = nullptr>
	std::remove_reference_t<typename Keyword::Type> GetDefault(Keyword keyword, Arg1&&, Args&& ...args)
{
	return GetDefault(std::forward<Keyword>(keyword), std::forward<Args>(args)...);
}
template <class Keyword>
std::remove_reference_t<typename Keyword::Type> GetDefault(Keyword /*keyword*/)
{
	throw InvalidArg("Default value does not exist.");
}

template <std::size_t KeyIndex>
struct GetKeywordArg_impl
{
	template <class Keyword, class ...Args>
	static typename Keyword::Type f(Keyword, Args&& ...args)
	{
		return std::get<KeyIndex>(std::forward_as_tuple(std::forward<Args>(args)...)).Get();
	}
};
template <>
struct GetKeywordArg_impl<std::numeric_limits<std::size_t>::max()>
{
	template <class Keyword, class ...Args>
	static std::remove_reference_t<typename Keyword::Type> f(Keyword k, Args&& ...args)
	{
		return detail::GetDefault(k, std::forward<Args>(args)...);
	}
};

}

template <class Type>
struct IsKeyword_impl
{
	static const bool value = IsBasedOn_T<std::remove_cv_t<Type>, detail::KeywordName>::value ||
		IsBasedOn_T<std::remove_cv_t<Type>, detail::KeywordValue>::value;
};
template <class Type>
struct IsKeyword
{
	static const bool value = IsKeyword_impl<Type>::value;
};
template <class ...Types>
struct AreAllKeywords
{
	static const bool value = AndOperationSeq<IsKeyword<Types>::value...>::value;
};
template <>
struct AreAllKeywords<>
{
	static const bool value = true;
};

#define L2C_DEFINE_KEYWORD_OPTION(NAME)\
constexpr auto NAME = l2c::detail::KeywordName<struct _##NAME, bool, void>();

#define L2C_DEFINE_TAGGED_KEYWORD_OPTION(NAME, TAG)\
constexpr auto NAME = l2c::detail::KeywordName<struct _##NAME, bool, TAG>();

#define L2C_DEFINE_KEYWORD_OPTION_WITH_VALUE(NAME, TYPE)\
constexpr auto NAME = l2c::detail::KeywordName<struct _##NAME, TYPE, void>();

#define L2C_DEFINE_TAGGED_KEYWORD_OPTION_WITH_VALUE(NAME, TYPE, TAG)\
constexpr auto NAME = l2c::detail::KeywordName<struct _##NAME, TYPE, TAG>();

namespace detail
{

template <class Option, class Tag, bool B = IsKeyword<Option>::value>
struct IsTaggedWith
{
	static constexpr bool value = false;
};
template <class Option, class Tag>
struct IsTaggedWith<Option, Tag, true>
{
	static constexpr bool value = std::is_base_of<typename Option::Tag, Tag>::value;
};

template <class Keyword, class Arg>
struct IsSameName
{
	static constexpr bool value = false;
};
template <class Keyword, class Name, class Type, class Tag>
struct IsSameName<Keyword, KeywordName<Name, Type, Tag>>
{
	static constexpr bool value = std::is_same<typename Keyword::Name, Name>::value;
};
template <class Keyword, class Name, class Type, class Tag>
struct IsSameName<Keyword, KeywordValue<Name, Type, Tag>>
{
	static constexpr bool value = std::is_same<typename Keyword::Name, Name>::value;
};

template <std::size_t N, class Keyword, class ...Args>
struct FindKeyword_impl;
template <std::size_t N, class Keyword>
struct FindKeyword_impl<N, Keyword>
{
	static constexpr std::size_t Index = std::numeric_limits<std::size_t>::max();
};
template <std::size_t N, class Keyword, class Arg, class ...Args>
struct FindKeyword_impl<N, Keyword, Arg, Args...>
{
	static constexpr std::size_t Index =
		IsSameName<Keyword, Arg>::value ? N : FindKeyword_impl<N + 1, Keyword, Args...>::Index;
};

template <class Keyword, class ...Args>
struct FindKeyword
{
	static constexpr std::size_t Index = FindKeyword_impl<0, Keyword, Args...>::Index;
	static constexpr bool value = Index != std::numeric_limits<std::size_t>::max();
};

}

#define CUF_TAGGED_ARGS_ENABLER(OPTIONS, TAG)\
bool CUF_TAG_IS_BASE_OF_ARGTAG = (sizeof...(OPTIONS) == 0 ||\
								  AndOperationSeq<adapt::detail::IsTaggedWith<std::remove_reference_t<OPTIONS>, TAG>::value...>::value),\
std::enable_if_t<CUF_TAG_IS_BASE_OF_ARGTAG, std::nullptr_t> = nullptr

template <class Keyword, class ...Args>
const bool KeywordExists(Keyword, Args&& ...)
{
	return detail::FindKeyword<Keyword, std::remove_cv_t<Args>...>::value;
}
template <class Keyword, class ...Args>
const bool KeywordExists(Keyword, std::tuple<Args...>)
{
	return detail::FindKeyword<Keyword, std::remove_cv_t<Args>...>::value;
}
template <class Keyword, class ...Args>
constexpr bool KeywordExists()
{
	return detail::FindKeyword<std::remove_cv_t<Keyword>, std::remove_cv_t<Args>...>::value;
}

template <class Keyword, class ...Args>
decltype(auto) GetKeywordArg(Keyword k, Args&& ...args)
{
	return detail::GetKeywordArg_impl<
		detail::FindKeyword<Keyword, std::remove_cv_t<Args>...>::Index>::
		f(k, std::forward<Args>(args)...);
}

}

#endif