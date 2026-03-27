#ifndef VXX_KEYWORD_ARGS_H
#define VXX_KEYWORD_ARGS_H

#include "Template.h"
#include <tuple>

#define VXX_TIE_ARGS(...) __VA_ARGS__

namespace vxx
{

//関数の引数にPythonっぽいキーワード可変長引数を作るための補助クラス。
//func(param1 = x, param2 = y);みたいな呼び出しをしたい.

namespace detail
{

template <class Name, class Type>
struct KeywordValue
{
	using _Type = Type;

	KeywordValue(_Type v) : mValue(std::forward<_Type>(v)) {}
	_Type Get() { return std::forward<_Type>(mValue); }
private:
	_Type mValue;
};
template <class Name, class Type>
struct KeywordName;
template <class Name, class Type>
struct KeywordName
{
	using _Type = Type;
	using _Value = KeywordValue<Name, Type>;

	constexpr KeywordName() {}
	constexpr _Value operator=(_Type v) const { return _Value(std::forward<_Type>(v)); }
};
template <class Name>
struct KeywordName<Name, void>
{
	using _Type = void;
	using _Value = KeywordName<Name, void>;
};

template <class Keyword>
std::remove_reference_t<typename Keyword::_Type> GetDefault(Keyword keyword)
{
	static_assert(0, "default value does not exist.");
}
template <class Keyword, class Arg1>
std::remove_reference_t<typename Keyword::_Type> GetDefault(Keyword keyword, Arg1&& arg1)
{
	static_assert(IsBasedOn<Arg1, typename Keyword::_Type>::value, "default value does not exist.");
	return std::forward<Arg1>(arg1);
}
template <class Keyword, class Arg1, class ...Args,
		  bool B = (sizeof...(Args) != 0),
		  std::enable_if_t<B, std::nullptr_t> = nullptr>
std::remove_reference_t<typename Keyword::_Type> GetDefault(Keyword keyword, Arg1&&, Args&& ...args)
{
	return GetDefault(std::forward<Keyword>(keyword), std::forward<Args>(args)...);
}

template <class keyword>
struct GetKeywordArgInfo
{
	//using Name = void;
	using Value = void;
	using Arg = void;
};
template <class Name_>
struct GetKeywordArgInfo<KeywordName<Name_, void>>
{
	//using Name = KeywordName<Name_, void>;
	using Value = void;
	using Arg = KeywordName<Name_, void>;
};
template <class Name_, class Type_>
struct GetKeywordArgInfo<KeywordName<Name_, Type_>>
{
	//using Name = KeywordName<Name_, Type_>;
	using Value = typename KeywordName<Name_, Type_>::_Value;
	using Arg = Value;
};

template <std::size_t KeyIndex>
struct GetKeywordArg_impl
{
	//キーワード引数が与えられている場合。
	template <class Keyword, class ...Args>
	static typename Keyword::_Type f(Keyword, Args&& ...args)
	{
		return std::get<KeyIndex>(std::forward_as_tuple(std::forward<Args>(args)...)).Get();
	}
};
template <>
struct GetKeywordArg_impl<std::numeric_limits<std::size_t>::max()>
{
	//キーワード引数が与えられていない場合。indexはstd::size_tの最大値になる。
	//デフォルト値が引数の最後に与えられている場合はそれを返し、
	//なければstatic_assertでコンパイルエラーにする。
	template <class Keyword, class ...Args>
	static std::remove_reference_t<typename Keyword::_Type> f(Keyword k, Args&& ...args)
	{
		return detail::GetDefault(k, std::forward<Args>(args)...);
	}
};

}

template <class Type>
struct IsKeyword_impl
{
	static const bool value = IsBasedOn_T<std::remove_reference_t<std::remove_const_t<Type>>, detail::KeywordName>::value ||
		IsBasedOn_T<std::remove_reference_t<std::remove_const_t<Type>>, detail::KeywordValue>::value;
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

#define VXX_DEFINE_KEYWORD_OPTION(NAME)\
constexpr auto NAME = vxx::detail::KeywordName<struct _##NAME, void>();
#define VXX_DEFINE_KEYWORD_OPTIONAL_ARG(NAME, TYPE)\
constexpr auto NAME = vxx::detail::KeywordName<struct _##NAME, TYPE>();

template <class Keyword, class ...Args>
constexpr bool KeywordExists(Keyword, Args&& ...args)
{
	return Find<typename detail::GetKeywordArgInfo<Keyword>::Arg, std::remove_const_t<std::remove_reference_t<Args>>...>::value;
}

template <class Keyword, class ...Args>
typename decltype(auto) GetKeywordArg(Keyword k, Args&& ...args)
{
	//キーワード引数が与えられている場合に呼ばれる。
	//該当するキーワードから値を取り出して返す。
	return detail::GetKeywordArg_impl<
		Find<typename detail::GetKeywordArgInfo<Keyword>::Value, std::remove_const_t<std::remove_reference_t<Args>>...>::Index>::
		f(k, std::forward<Args>(args)...);
}

}

#endif