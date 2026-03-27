#define NOMINMAX

#include <vector>
#include <set>
#include <cassert>
#include "KeywordArgs_l2c.h"

#if defined L2CX_DLL_BUILD
#define L2CX_DLL_EXPORT __declspec(dllexport)
#else
#define L2CX_DLL_EXPORT __declspec(dllimport)
#endif

namespace l2c
{
struct Group;
}

namespace netscan
{
template<typename T>
class SynchronizedQueue;
struct Group;
void trans_group_async(int32_t npl,
					   SynchronizedQueue< std::pair<std::pair<int64_t, int64_t>, l2c::Group> >& qbuf,
					   SynchronizedQueue<std::pair<int64_t, int64_t>>& nchain_sum,
					   const Group& gr, int64_t gid, int64_t gnchain_outlim);
}

namespace l2c
{

struct Linklet
{
	Linklet(int32_t aPos1, int32_t aPos2, int64_t aId1, int64_t aId2)
		: pos1(aPos1), pos2(aPos2), id1(aId1), id2(aId2) {}
	int32_t pos1, pos2;
	int64_t id1, id2;
};

struct BaseTrackID
{
	BaseTrackID(int64_t id) : bit(id) {}
	int32_t GetPL() const {
		return bit >> 48;
	}
	int64_t GetRawID() const {
		return bit & 0x0000ffffffffffff;
	}
	operator int64_t() const { return bit; }

	bool IsEmpty() const { return bit == -1; }

private:

	int64_t bit;
};

struct Chain
{
	Chain(int64_t aId, int32_t aNPL, int32_t aSPL, int32_t aEPL, int64_t aNSeg, const BaseTrackID* aChainPtr)
		: id(aId), npl(aNPL), spl(aSPL), epl(aEPL), nseg(aNSeg), ch(aChainPtr)
	{}

	int64_t GetID() const { return id; }
	int16_t GetStartPL() const { return spl; }
	int16_t GetEndPL() const { return epl; }
	int64_t GetNSeg() const { return nseg; }

	const BaseTrackID& GetBaseTrack(size_t pl) const
	{
		assert(pl <= npl);
		return ch[pl];
	}

private:

	int64_t id;
	int32_t npl;
	int32_t spl, epl;
	int64_t nseg;
	const BaseTrackID* ch;
};

struct Group
{
	static constexpr size_t Offset = 3;

	friend void netscan::trans_group_async(int32_t npl,
										   netscan::SynchronizedQueue< std::pair<std::pair<int64_t, int64_t>, l2c::Group> >& qbuf,
										   netscan::SynchronizedQueue<std::pair<int64_t, int64_t>>& nchain_sum,
										   const netscan::Group& gr, int64_t gid, int64_t gnchain_outlim);

	Group(int64_t aID, int64_t aNChain, int32_t aNPL, int32_t aSPL, int32_t aEPL)
		: id(aID), nch(aNChain), npl(aNPL), spl(aSPL), epl(aEPL)
	{}
	Group(int64_t aID, int32_t aNChain, int32_t aNPL)
		: id(aID), nch(aNChain), npl(aNPL), spl(0), epl(0)
	{}

	Group(const Group&) = default;
	Group(Group&& g) noexcept
		: id(g.id), nch(g.nch), npl(g.npl), spl(g.spl), epl(g.epl),
		ch(std::move(g.ch))
	{
		bt = std::move(g.bt);
		ll = std::move(g.ll);
	}
	Group& operator=(const Group&) = default;
	Group& operator=(Group&&) noexcept = default;

	std::set<BaseTrackID>& GetBaseTracks() { return bt; }
	const std::set<BaseTrackID>& GetBaseTracks() const { return bt; }
	std::set<std::pair<BaseTrackID, BaseTrackID>>& GetLinklets() { return ll; }
	const std::set<std::pair<BaseTrackID, BaseTrackID>>& GetLinklets() const { return ll; }

	int64_t GetID() const { return id; }

	bool IsOverUpperLim() const { return ch.empty(); }
	int64_t GetNumOfChains() const { return nch; }

	int32_t GetStartPL() const { return spl; }
	int32_t GetEndPL() const { return epl; }

	Chain GetChain(size_t n) const
	{
		assert(!IsOverUpperLim());
		size_t index = n * (npl + Offset);
		const BaseTrackID* c = &(ch[index]);
		int64_t id = *c; ++c;
		int64_t x = *c; ++c;
		int32_t spl = x >> 32;
		int32_t epl = 0x00000000ffffffff & x;
		int64_t nseg = *c; ++c;
		return Chain(id, npl, spl, epl, nseg, c);
	}

private:

	void AddBaseTrack(int64_t id) { bt.insert(id); };
	void AddLinklet(std::pair<int64_t, int64_t> id) { ll.insert(id); };
	BaseTrackID* InsertChain(int64_t cid, int16_t cspl, int16_t cepl, int32_t cnseg)
	{
		auto it = ch.insert(ch.end(), npl + Offset, -1);
		int64_t* p = (int64_t*)(&(*it));
		*p = cid; ++p;
		*p = (int64_t(cspl) << 32) + cepl; ++p;
		*p = cnseg; ++p;
		return (BaseTrackID*)p;
	}

	int64_t id;
	int64_t nch;
	int32_t npl;
	int32_t spl, epl;
	std::set<BaseTrackID> bt;
	std::set<std::pair<BaseTrackID, BaseTrackID>> ll;
	std::vector<BaseTrackID> ch;
};

struct Cdat
{
	Cdat(int64_t n) : npl(n) {}
	Cdat(int64_t n, std::vector<Group>&& gr) : npl(n), gr(std::move(gr)) {}

	size_t GetNumOfGroups() const { return gr.size(); }
	Group& GetGroup(size_t index) { return gr[index]; }
	const Group& GetGroup(size_t index) const { return gr[index]; }
	std::vector<Group>& GetGroups() { return gr; }
	const std::vector<Group>& GetGroups() const { return gr; }

private:
	std::vector<Group> gr;
	int64_t npl;
};

namespace opt
{

L2C_DEFINE_KEYWORD_OPTION_WITH_VALUE(usepos, const std::vector<int>&)
L2C_DEFINE_KEYWORD_OPTION_WITH_VALUE(upperlim, int64_t)
L2C_DEFINE_KEYWORD_OPTION(output_isolated_linklet)

struct Options
{
	std::vector<int32_t> mUsePos;
	int64_t mUpperLim;
	bool mOutputIsolatedLinklet;
};

}

Cdat L2CX_DLL_EXPORT MakeCdat(const std::vector<Linklet>& ltlist, opt::Options ops);

template <class ...Ops>
Cdat MakeCdat(const std::vector<Linklet>& ltlist, Ops&& ...ops)
{
	opt::Options opspack;
	constexpr bool b = KeywordExists<decltype(opt::usepos), Ops...>();
	opspack.mUsePos = GetKeywordArg(opt::usepos, std::forward<Ops>(ops)...);
	opspack.mUpperLim = GetKeywordArg(opt::upperlim, std::forward<Ops>(ops)..., -1);
	opspack.mOutputIsolatedLinklet = GetKeywordArg(opt::output_isolated_linklet, std::forward<Ops>(ops)..., false);

	return MakeCdat(ltlist, opspack);

}
Cdat L2CX_DLL_EXPORT MakeCdat_no_out(const std::vector<Linklet>& ltlist, opt::Options ops);

template <class ...Ops>
Cdat MakeCdat_no_out(const std::vector<Linklet>& ltlist, Ops&& ...ops)
{
	opt::Options opspack;
	constexpr bool b = KeywordExists<decltype(opt::usepos), Ops...>();
	opspack.mUsePos = GetKeywordArg(opt::usepos, std::forward<Ops>(ops)...);
	opspack.mUpperLim = GetKeywordArg(opt::upperlim, std::forward<Ops>(ops)..., -1);
	opspack.mOutputIsolatedLinklet = GetKeywordArg(opt::output_isolated_linklet, std::forward<Ops>(ops)..., false);

	return MakeCdat_no_out(ltlist, opspack);

}

}