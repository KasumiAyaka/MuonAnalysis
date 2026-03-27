//muonのchainに対して、linkletから再度chainを再生成する
//1group 1chainが前提
#define _CRT_SECURE_NO_WARNINGS
#define DEFAULT_CHAIN_UPPERLIM 1000000000000000

#pragma comment(lib, "VxxReader.lib")
#include "VxxReader.h"
#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.h>
#pragma comment(lib, "lib_l2cx.lib")
#include <LibL2c-x.h>
//これを202行目に追記すると動く
/*template <class ...Ops>
Cdat MakeCdat_no_out(const std::vector<Linklet>& ltlist, Ops&& ...ops)
{
	opt::Options opspack;
	constexpr bool b = KeywordExists<decltype(opt::usepos), Ops...>();
	opspack.mUsePos = GetKeywordArg(opt::usepos, std::forward<Ops>(ops)...);
	opspack.mUpperLim = GetKeywordArg(opt::upperlim, std::forward<Ops>(ops)..., -1);
	opspack.mOutputIsolatedLinklet = GetKeywordArg(opt::output_isolated_linklet, std::forward<Ops>(ops)..., false);

	return MakeCdat_no_out(ltlist, opspack);

}
*/

#include <filesystem>
#include <set>
#include <netscan_data_types_ui.h>
#include <algorithm>
#include <iostream>
#include <unordered_map>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <stack>
#include <unordered_set>

#include <Chain_convolution_2.h>
#include <bipartite_graph_enumeration.h>
#include <Cycle_enumerate.h>

#include <omp.h>
using namespace l2c;

class output_format {
public:
	int groupid, trackid;
	int num_comfirmed_path, num_cut_path, num_select_path;
	std::vector<std::pair<int, std::tuple<int, int, int, int>>> comfirmed_path;
	std::vector<std::pair<int, std::tuple<int, int, int, int>>> cut_path;
	std::vector<std::pair<int, std::tuple<int, int, int, int>>> select_path;
};
//bool operator<(const Point &left, const Point &right){return left.y < right.y;}
// 自分で作った型をunorderdコンテナに入れたいときは、operator== の他に
// 型と同じ名前空間で hash_value 関数を定義
namespace mfile0 {
	bool operator<(const mfile0::M_Base& left, const mfile0::M_Base& right)
	{
		if (left.pos == right.pos)return left.rawid < right.rawid;
		return left.pos < right.pos;
	}
	bool operator==(const mfile0::M_Base& right, const mfile0::M_Base& left)
	{
		return right.pos == left.pos && right.rawid == left.rawid;
	}
	size_t hash_value(const mfile0::M_Base& b)
	{
		// 複数の値のハッシュ値を組み合わせてハッシュ値を計算するには、
		// boost::hash_combine を使います。
		size_t h = 0;
		boost::hash_combine(h, b.pos);
		boost::hash_combine(h, b.rawid);
		return h;
	}

}
//bool operator<(const std::pair<int, int>&right, const std::pair<int, int>& left) {
//	if (right.first == left.first)return right.second < left.second;
//	return right.first < left.first;
//}
bool Base_set_compare(std::set<mfile0::M_Base>& b0, std::set<mfile0::M_Base>& b1) {
	if (b0.size() != b1.size())return false;
	auto itr_b0 = b0.begin();
	auto itr_b1 = b1.begin();
	for (int i = 0; i < b0.size(); i++) {
		if (*itr_b0 == *itr_b1) {
			itr_b0++;
			itr_b1++;
			continue;
		}
		return false;

	}
	return true;
}
std::vector< output_format> read_chain_file(std::string filename, std::multimap<int, mfile0::M_Base>& base_map);
std::map<int, corrmap0::Corrmap> read_corrmap_abs(std::string file_in_ECC);
void read_basetrack(std::string file_in_ECC, std::multimap<int, mfile0::M_Base>& base_map);
void read_basetrack_id(std::string file_in_ECC, std::multimap<int, mfile0::M_Base>& base_map);
void read_basetrack2(std::string file_in_base, std::multimap<int, mfile0::M_Base>& base_map);

void apply_corrmap(std::multimap<int, mfile0::M_Base>& base_map, std::map<int, corrmap0::Corrmap>& corr_abs, std::map<int, double>& z_map);
void apply_corrmap_local(std::multimap<int, mfile0::M_Base>& base_map, std::map<int, std::vector<corrmap_3d::align_param2>>& corrmap_dd);

std::vector<std::vector<std::pair<mfile0::M_Base, mfile0::M_Base>>> divide_connected(output_format& g, std::map<std::pair<int, int>, mfile0::M_Base>& base);
l2c::Cdat l2c_x(std::set<std::pair<int32_t, int64_t>>& btset, std::vector<Linklet>& ltlist, std::vector<int32_t>& usepos, int64_t upperlim = DEFAULT_CHAIN_UPPERLIM, bool output = false);
std::vector<std::pair<mfile0::M_Base, mfile0::M_Base>> select_path(std::vector<std::pair<mfile0::M_Base, mfile0::M_Base>>& path, std::map<std::pair<int, int>, mfile0::M_Base>& base_map_id);
std::vector<std::vector<std::pair<mfile0::M_Base, mfile0::M_Base>>> divide_connected(output_format& g, std::map<std::pair<int, int>, mfile0::M_Base>& base);
std::vector<std::vector<std::pair<mfile0::M_Base, mfile0::M_Base>>> divide_connected(std::vector<std::pair<int, std::tuple<int, int, int, int>>>& path_list, std::map<std::pair<int, int>, mfile0::M_Base>& base);
std::vector<std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>> divide_connected(std::vector<std::pair<int, std::tuple<int, int, int, int>>>& path_list);

std::vector<std::vector<std::pair<mfile0::M_Base, mfile0::M_Base>>> select_path_v2(std::vector<std::pair<mfile0::M_Base, mfile0::M_Base>>& path);
std::map<std::pair<mfile0::M_Base, mfile0::M_Base>, Chain_path> path_compress_next(std::set<mfile0::M_Base>& all_vertex, boost::unordered_multimap<mfile0::M_Base, mfile0::M_Base>& path_next, boost::unordered_multimap<mfile0::M_Base, mfile0::M_Base>& path_prev);
std::map<std::pair<mfile0::M_Base, mfile0::M_Base>, Chain_path> path_compress_prev(std::set<mfile0::M_Base>& all_vertex, boost::unordered_multimap<mfile0::M_Base, mfile0::M_Base>& path_next, boost::unordered_multimap<mfile0::M_Base, mfile0::M_Base>& path_prev);

std::set<mfile0::M_Base> select_start(std::set<mfile0::M_Base>& all_vertex, boost::unordered_multimap<mfile0::M_Base, mfile0::M_Base>& path_next, boost::unordered_multimap<mfile0::M_Base, mfile0::M_Base>& path_prev);
std::vector < std::pair<std::set<mfile0::M_Base>, std::set<mfile0::M_Base>>> Divide_bipartite_graph(output_format& g, std::vector<std::pair<mfile0::M_Base, mfile0::M_Base>>& path);

void output_linklet(std::string filename, std::vector<output_format>& out);
void output_linklet_compress(std::string filename, std::vector<output_format>& out);
bool Judge_bipartite_graph(int target, boost::unordered_multimap <int, int>& path_prev, boost::unordered_multimap <int, int>& path_next, std::set<int>& up, std::set<int>& down);
void result_hangnail_cycle(output_format& g, std::map<std::pair<int, int>, mfile0::M_Base>& base, int& path_id);
std::map<int, int> add_depth_information(std::set<int>vertex, std::set<int>up, std::set<int>down, boost::unordered_multimap <int, int>& path_next);

void result_bipartite_graph_0(
	std::vector < std::pair<std::set<mfile0::M_Base>, std::set<mfile0::M_Base>>>& bipartite_graph_v,
	std::list<std::pair<int, std::tuple<int, int, int, int>>>& single_path,
	std::list<std::pair<int, std::tuple<int, int, int, int>>>& branch_path,
	std::map<std::pair<int, int>, mfile0::M_Base>& base_map_id, int& path_id);
void path_organize(std::vector<std::pair<mfile0::M_Base, mfile0::M_Base>>& best_path,
	boost::unordered_multimap <mfile0::M_Base, mfile0::M_Base>& path_next,
	boost::unordered_multimap <mfile0::M_Base, mfile0::M_Base>& path_prev,
	std::list<std::pair<int, std::tuple<int, int, int, int>>>& single_path,
	std::list<std::pair<int, std::tuple<int, int, int, int>>>& branch_path,
	int& path_id
);
std::vector<std::pair<mfile0::M_Base, mfile0::M_Base>> Select_best_combination(std::vector<std::vector<std::pair<mfile0::M_Base, mfile0::M_Base>>>& enumerat_all_path_base, std::map<mfile0::M_Base, Chain_path> chain_map);
double Calc_combination_value(std::vector<std::pair<Chain_path, Chain_path>>& chains);
void result_no_bipartite_graph(output_format& group, std::vector<std::pair<mfile0::M_Base, mfile0::M_Base>>& path, int& path_id);
void result_cycle_cross(output_format& g, std::map<std::pair<int, int>, mfile0::M_Base>& base, int& path_id);

void make_to_next_chain(std::map<mfile0::M_Base, Chain_path>& chain_map, std::set<mfile0::M_Base>& base, boost::unordered_multimap <mfile0::M_Base, mfile0::M_Base>& path_next);
void make_to_prev_chain(std::map<mfile0::M_Base, Chain_path>& chain_map, std::set<mfile0::M_Base>& base, boost::unordered_multimap <mfile0::M_Base, mfile0::M_Base>& path_prev);
bool check_insert_path(std::pair< int, int >& path, boost::unordered_multimap <int, int>& path_next);
void path_id_reroll(output_format& g);
void enumerate_path_dfs(std::vector<mfile0::M_Base>& path, std::vector<std::vector<mfile0::M_Base>>& all_hist, mfile0::M_Base now, mfile0::M_Base& goal, std::set<mfile0::M_Base>& seen, std::multimap<mfile0::M_Base, mfile0::M_Base>& all_path);
void result_no_bipartite_graph(output_format& group, std::vector< std::pair<std::pair<int, int>, std::pair<int, int>>>& path, int& path_id);
int use_thread(double ratio, bool output);

int main(int argc, char** argv) {
	if (argc != 5) {
		fprintf(stderr, "usage:file-in-chain file-in-ECC file-out mode\n");
		fprintf(stderr, "mode=0:result complete bipartite graph \n");
		fprintf(stderr, "mode=1:result bipartite graph \n");
		fprintf(stderr, "mode=2:result one segment cycle\n");
		exit(1);
	}

	std::string file_in_chain = argv[1];
	std::string file_in_ECC = argv[2];
	std::string file_out = argv[3];
	int mode = std::stoi(argv[4]);

	std::multimap<int, mfile0::M_Base>base_map;
	std::vector< output_format> group = read_chain_file(file_in_chain, base_map);
	printf("group size = %d\n", group.size());

	//corrmap absを使う場合
	//std::map<int, corrmap0::Corrmap> corrmap = read_corrmap_abs(file_in_ECC);
	////gap nominal read
	//std::stringstream structure_path;
	//structure_path << file_in_ECC << "\\st\\st.dat";
	//printf("%s\n", structure_path.str().c_str());
	//chamber1::Chamber chamber;
	//chamber1::read_structure(structure_path.str(), chamber);
	//std::map<int, double> z_map = chamber1::base_z_convert(chamber);

	//新alignmentを使う場合
	std::string file_in_corrmap = file_in_ECC + "\\Area0\\0\\align\\fine\\local\\corrmap-local-abs.lst";

	std::map<int, std::vector<corrmap_3d::align_param>>corrmap = corrmap_3d::read_ali_param_abs(file_in_corrmap, 1);
	std::map<int, std::vector<corrmap_3d::align_param2>>corrmap_dd = corrmap_3d::DelaunayDivide_map(corrmap);

	//basetrack情報のpick up
	read_basetrack(file_in_ECC, base_map);
	//read_basetrack_id(file_in_ECC, base_map);
	//read_basetrack2(file_in_base, base_map);
	printf("read fin\n");

	//corrmapの適用
	//corrmap absを使う場合
	//apply_corrmap(base_map, corrmap, z_map);
	apply_corrmap_local(base_map, corrmap_dd);

	std::map<std::pair<int, int>, mfile0::M_Base>base_map_id;
	for (auto itr = base_map.begin(); itr != base_map.end(); itr++) {
		base_map_id.insert(std::make_pair(std::make_pair(itr->first, itr->second.rawid), itr->second));
	}
	printf("apply corr fin\n");

	//完全2部グラフの解消
	for (int i = 0; i < group.size(); i++) {
		printf("Result complete bipartite graph :group %d/%d %d-%d\r", i, group.size(), group[i].groupid, group[i].trackid);
		int path_id = 0;
		for (auto& p : group[i].comfirmed_path) {
			path_id = std::max(path_id, p.first);
		}
		for (auto& p : group[i].cut_path) {
			path_id = std::max(path_id, p.first);
		}
		for (auto& p : group[i].select_path) {
			path_id = std::max(path_id, p.first);
		}
		path_id += 1;
		//idを使ってる-->高効率
		std::vector<std::vector<std::pair<mfile0::M_Base, mfile0::M_Base>>> connect_path = divide_connected(group[i], base_map_id);
		for (int j = 0; j < connect_path.size(); j++) {
			//printf("path %d/%d\n", j,connect_path.size());
			auto result = select_path_v2(connect_path[j]);
			for (int k = 0; k < result.size(); k++) {
				for (int l = 0; l < result[k].size(); l++) {
					group[i].comfirmed_path.push_back(std::make_pair(path_id,
						std::make_tuple(
							result[k][l].first.pos,
							result[k][l].second.pos,
							result[k][l].first.rawid,
							result[k][l].second.rawid
						)));
				}
				path_id++;
			}
		}
		group[i].select_path.clear();

		group[i].num_comfirmed_path = group[i].comfirmed_path.size();
		group[i].num_cut_path = group[i].cut_path.size();
		group[i].num_select_path = group[i].select_path.size();
	}
	printf("\n");
	mode -= 1;
	if (mode < 0) {
		output_linklet(file_out, group);
		return 0;
	}

	//2部グラフでない場所を2部グラフにする
	//2部グラフでない頂点+pathを書き出す(全頂点で実行)
	//pathを足す
	//非2部グラフの解消
	for (int i = 0; i < group.size(); i++) {
		printf("Result no bipartite graph :group %d-%d\r", i, group.size(), group[i].groupid, group[i].trackid);
		//////////////////////
		int path_id = 0;
		for (auto& p : group[i].comfirmed_path) {
			path_id = std::max(path_id, p.first);
		}
		for (auto& p : group[i].cut_path) {
			path_id = std::max(path_id, p.first);
		}
		for (auto& p : group[i].select_path) {
			path_id = std::max(path_id, p.first);
		}
		path_id += 1;

		//連結成分に分解
		std::vector<std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>> connect_path = divide_connected(group[i].cut_path);
		for (int j = 0; j < connect_path.size(); j++) {
			//printf("connect part %d/%d\r", j, connect_path.size());

			result_no_bipartite_graph(group[i], connect_path[j], path_id);

		}
		path_id_reroll(group[i]);

	}
	printf("\n");

	/*
	for (int i = 0; i < group.size(); i++) {
		printf("Result no bipartite graph : group %d/%d\r", i, group.size());

		//printf("result no bipartite_graph group %d\n", i);
		//printf("group %d\n", i);
		int path_id = 0;
		for (auto &p : group[i].comfirmed_path) {
			path_id = std::max(path_id, p.first);
		}
		for (auto &p : group[i].cut_path) {
			path_id = std::max(path_id, p.first);
		}
		for (auto &p : group[i].select_path) {
			path_id = std::max(path_id, p.first);
		}
		path_id += 1;
		//idを使わない-->普遍的
		std::vector<std::vector<std::pair<mfile0::M_Base, mfile0::M_Base>>> connect_path = divide_connected(group[i].cut_path, base_map_id);
		for (int j = 0; j < connect_path.size(); j++) {
			//printf("connect part %d/%d\n", j, connect_path.size());

			result_no_bipartite_graph(group[i], connect_path[j], path_id);

		}
	}
	printf("\n");
	*/
	//解決できそうな場所を随時解消
	for (int i = 0; i < group.size(); i++) {
		printf("Result bipartite graph : group %d-%d\r", i, group.size(), group[i].groupid, group[i].trackid);
		int path_id = 0;
		for (auto& p : group[i].comfirmed_path) {
			path_id = std::max(path_id, p.first);
		}
		for (auto& p : group[i].cut_path) {
			path_id = std::max(path_id, p.first);
		}
		for (auto& p : group[i].select_path) {
			path_id = std::max(path_id, p.first);
		}
		path_id += 1;
		//1-n-1閉路の解決
		result_hangnail_cycle(group[i], base_map_id, path_id);

		//idを使わない-->普遍的
		std::vector<std::vector<std::pair<mfile0::M_Base, mfile0::M_Base>>> connect_path = divide_connected(group[i].cut_path, base_map_id);
		for (int j = 0; j < connect_path.size(); j++) {
			//printf("connect part %d/%d\n", j, connect_path.size());


			//pathを2部グラフに分解
			auto bipartite_graph_v = Divide_bipartite_graph(group[i], connect_path[j]);
			//2部グラフを一個ずつ解決
			std::list<std::pair<int, std::tuple<int, int, int, int>>> single_path(group[i].comfirmed_path.begin(), group[i].comfirmed_path.end());
			std::list<std::pair<int, std::tuple<int, int, int, int>>> branch_path(group[i].cut_path.begin(), group[i].cut_path.end());

			result_bipartite_graph_0(bipartite_graph_v, single_path, branch_path, base_map_id, path_id);

			std::vector<std::pair<int, std::tuple<int, int, int, int>>> single_path_v(single_path.begin(), single_path.end());
			std::vector<std::pair<int, std::tuple<int, int, int, int>>> branch_path_v(branch_path.begin(), branch_path.end());
			group[i].comfirmed_path = single_path_v;
			group[i].cut_path = branch_path_v;
			//1-n-1閉路の解決
			result_hangnail_cycle(group[i], base_map_id, path_id);
		}
		group[i].select_path.clear();

		group[i].num_comfirmed_path = group[i].comfirmed_path.size();
		group[i].num_cut_path = group[i].cut_path.size();
		group[i].num_select_path = group[i].select_path.size();

	}
	printf("\n");
	mode -= 1;
	if (mode < 0) {
		output_linklet(file_out, group);
		return 0;
	}

	for (int i = 0; i < group.size(); i++) {
		path_id_reroll(group[i]);
		printf("Result cycle : group %d-%d\r", i, group.size(), group[i].groupid, group[i].trackid);
		int path_id = 0;
		for (auto& p : group[i].comfirmed_path) {
			path_id = std::max(path_id, p.first);
		}
		for (auto& p : group[i].cut_path) {
			path_id = std::max(path_id, p.first);
		}
		for (auto& p : group[i].select_path) {
			path_id = std::max(path_id, p.first);
		}
		path_id += 1;
		result_cycle_cross(group[i], base_map_id, path_id);
		path_id_reroll(group[i]);
	}
	printf("\n");
	output_linklet(file_out, group);

}

l2c::Cdat l2c_x(std::set<std::pair<int32_t, int64_t>>& btset, std::vector<Linklet>& ltlist, std::vector<int32_t>& usepos, int64_t upperlim, bool output) {
	try
	{
		size_t possize = usepos.size();

		//l2c関数本体を呼び出す。
		//出力されるCdatは基本的にl2c-xと同等のものになっているはずだが、若干の違いがある。
		//1. upperlimを超過したGroupは、Chainは持たず、全BaseTrackの情報だけを持った状態で出力される。
		//2. upperlim超過のGroupも含め、属す全BaseTrackの情報を持っている。
		//またplate枚数は最大64枚。変更するにはLibL2c-xのリビルドを要する。
		//デフォルトでoutput_isolated_linkletはfalse、upperlimは無効になっている。useposは与えないとエラーになる。
		if (output) {
			l2c::Cdat cdat = l2c::MakeCdat(ltlist, opt::usepos = usepos, opt::upperlim = upperlim, opt::output_isolated_linklet = true);
			return cdat;
		}
		l2c::Cdat cdat = l2c::MakeCdat_no_out(ltlist, opt::usepos = usepos, opt::upperlim = upperlim, opt::output_isolated_linklet = true);
		return cdat;

	}
	catch (const std::exception& e)
	{
		fprintf(stderr, "%s\n", e.what());
	}

}

std::vector< output_format> read_chain_file(std::string filename, std::multimap<int, mfile0::M_Base>& base_map) {
	std::vector< output_format> ret;
	std::map<int, std::map<int, mfile0::M_Base>>base_map_pl;
	std::ifstream ifs(filename);
	output_format out;
	std::tuple<int, int, int, int> path;
	int id;
	int pl, rawid;
	while (ifs >> out.groupid >> out.trackid >> out.num_comfirmed_path >> out.num_cut_path >> out.num_select_path) {
		for (int i = 0; i < out.num_comfirmed_path; i++) {
			ifs >> id >> std::get<0>(path) >> std::get<1>(path) >> std::get<2>(path) >> std::get<3>(path);
			out.comfirmed_path.push_back(std::make_pair(id, path));

			pl = std::get<0>(path) / 10;
			rawid = std::get<2>(path);
			auto res = base_map_pl.find(pl);
			mfile0::M_Base b;
			b.pos = pl * 10;
			b.rawid = rawid;
			if (res == base_map_pl.end()) {
				std::map<int, mfile0::M_Base> base_map_tmp;
				base_map_tmp.insert(std::make_pair(rawid, b));
				base_map_pl.insert(std::make_pair(pl, base_map_tmp));
			}
			else {
				res->second.insert(std::make_pair(rawid, b));
			}

			pl = std::get<1>(path) / 10;
			rawid = std::get<3>(path);
			res = base_map_pl.find(pl);
			b.pos = pl * 10;
			b.rawid = rawid;
			if (res == base_map_pl.end()) {
				std::map<int, mfile0::M_Base> base_map_tmp;
				base_map_tmp.insert(std::make_pair(rawid, b));
				base_map_pl.insert(std::make_pair(pl, base_map_tmp));
			}
			else {
				res->second.insert(std::make_pair(rawid, b));
			}

		}
		for (int i = 0; i < out.num_cut_path; i++) {
			ifs >> id >> std::get<0>(path) >> std::get<1>(path) >> std::get<2>(path) >> std::get<3>(path);
			out.cut_path.push_back(std::make_pair(id, path));

			pl = std::get<0>(path) / 10;
			rawid = std::get<2>(path);
			auto res = base_map_pl.find(pl);
			mfile0::M_Base b;
			b.pos = pl * 10;
			b.rawid = rawid;
			if (res == base_map_pl.end()) {
				std::map<int, mfile0::M_Base> base_map_tmp;
				base_map_tmp.insert(std::make_pair(rawid, b));
				base_map_pl.insert(std::make_pair(pl, base_map_tmp));
			}
			else {
				res->second.insert(std::make_pair(rawid, b));
			}

			pl = std::get<1>(path) / 10;
			rawid = std::get<3>(path);
			res = base_map_pl.find(pl);
			b.pos = pl * 10;
			b.rawid = rawid;
			if (res == base_map_pl.end()) {
				std::map<int, mfile0::M_Base> base_map_tmp;
				base_map_tmp.insert(std::make_pair(rawid, b));
				base_map_pl.insert(std::make_pair(pl, base_map_tmp));
			}
			else {
				res->second.insert(std::make_pair(rawid, b));
			}

		}
		for (int i = 0; i < out.num_select_path; i++) {
			ifs >> id >> std::get<0>(path) >> std::get<1>(path) >> std::get<2>(path) >> std::get<3>(path);
			out.select_path.push_back(std::make_pair(id, path));

			pl = std::get<0>(path) / 10;
			rawid = std::get<2>(path);
			auto res = base_map_pl.find(pl);
			mfile0::M_Base b;
			b.pos = pl * 10;
			b.rawid = rawid;
			if (res == base_map_pl.end()) {
				std::map<int, mfile0::M_Base> base_map_tmp;
				base_map_tmp.insert(std::make_pair(rawid, b));
				base_map_pl.insert(std::make_pair(pl, base_map_tmp));
			}
			else {
				res->second.insert(std::make_pair(rawid, b));
			}

			pl = std::get<1>(path) / 10;
			rawid = std::get<3>(path);
			res = base_map_pl.find(pl);
			b.pos = pl * 10;
			b.rawid = rawid;
			if (res == base_map_pl.end()) {
				std::map<int, mfile0::M_Base> base_map_tmp;
				base_map_tmp.insert(std::make_pair(rawid, b));
				base_map_pl.insert(std::make_pair(pl, base_map_tmp));
			}
			else {
				res->second.insert(std::make_pair(rawid, b));
			}

		}

		ret.push_back(out);
		out.comfirmed_path.clear();
		out.cut_path.clear();
		out.select_path.clear();
	}

	for (auto itr = base_map_pl.begin(); itr != base_map_pl.end(); itr++) {
		for (auto itr2 = itr->second.begin(); itr2 != itr->second.end(); itr2++) {
			base_map.insert(std::make_pair(itr->first, itr2->second));
			//printf("%d %d\n", itr->first, itr2->first);
		}
	}
	return ret;
}

std::map<int, corrmap0::Corrmap> read_corrmap_abs(std::string file_in_ECC) {
	std::stringstream file_in_corr;
	file_in_corr << file_in_ECC << "\\Area0\\0\\align\\corrmap-abs.lst";

	std::vector<corrmap0::Corrmap> corr;
	corrmap0::read_cormap(file_in_corr.str(), corr);

	std::map<int, corrmap0::Corrmap> ret;
	for (auto itr = corr.begin(); itr != corr.end(); itr++) {
		ret.insert(std::make_pair(itr->pos[0] / 10, *itr));
	}
	return ret;
}

void read_basetrack(std::string file_in_ECC, std::multimap<int, mfile0::M_Base>& base_map) {
	std::vector<int> read_pl;
	for (auto itr = base_map.begin(); itr != base_map.end(); itr++) {
		int count = base_map.count(itr->first);
		read_pl.push_back(itr->first);
		itr = std::next(itr, count - 1);
	}

#pragma omp parallel for num_threads(10) schedule(dynamic,1)
	for (int i = 0; i < read_pl.size(); i++) {

		std::stringstream file_in_base;
		file_in_base << file_in_ECC << "\\Area0\\PL"
			<< std::setw(3) << std::setfill('0') << read_pl[i] << "\\b"
			<< std::setw(3) << std::setfill('0') << read_pl[i] << ".sel.cor.vxx";
		std::vector<vxx::base_track_t >base;
		vxx::BvxxReader br;
		base = br.ReadAll(file_in_base.str(), read_pl[i], 0);
		std::map<int, vxx::base_track_t> base_rawmap;
		for (auto itr = base.begin(); itr != base.end(); itr++) {
			base_rawmap.insert(std::make_pair(itr->rawid, *itr));
		}
		printf("count = %d\n", base_map.count(read_pl[i]));
		auto range = base_map.equal_range(read_pl[i]);
		for (auto res = range.first; res != range.second; res++) {
			auto b = base_rawmap.find(res->second.rawid);
			if (b == base_rawmap.end()) {
				fprintf(stderr, "exception PL%03d rawid=%d not found\n", read_pl[i], res->second.rawid);
			}
			res->second.ax = b->second.ax;
			res->second.ay = b->second.ay;
			res->second.x = b->second.x;
			res->second.y = b->second.y;
			res->second.ph = b->second.m[0].ph + b->second.m[1].ph;
			res->second.flg_d[0] = 0;
			res->second.flg_d[1] = 0;
			res->second.flg_i[0] = 0;
			res->second.flg_i[1] = 0;
			res->second.flg_i[2] = 0;
			res->second.flg_i[3] = 0;
		}
	}


}

void read_basetrack_id(std::string file_in_ECC, std::multimap<int, mfile0::M_Base>& base_map) {


	std::vector< mfile0::M_Base*>base_p;
	std::vector<int> read_pl;
	for (auto itr = base_map.begin(); itr != base_map.end(); itr++) {
		base_p.push_back(&(itr->second));
	}

#pragma omp parallel for  num_threads(10) 
	for (int i = 0; i < base_p.size(); i++) {
		int pl = base_p[i]->pos / 10;
		int rawid = base_p[i]->rawid;
		std::stringstream file_in_base;
		file_in_base << file_in_ECC << "\\Area0\\PL"
			<< std::setw(3) << std::setfill('0') << pl << "\\b"
			<< std::setw(3) << std::setfill('0') << pl << ".sel.cor.vxx";
		std::vector<vxx::base_track_t >base;
		vxx::BvxxReader br;
		std::array<int, 2> index = { rawid,rawid + 1 };//1234<=rawid<=5678であるようなものだけを読む。
		base = br.ReadAll(file_in_base.str(), pl, 0, vxx::opt::index = index);

		base_p[i]->ax = base[0].ax;
		base_p[i]->ay = base[0].ay;
		base_p[i]->x = base[0].x;
		base_p[i]->y = base[0].y;
		base_p[i]->ph = base[0].m[0].ph + base[0].m[1].ph;
		base_p[i]->flg_d[0] = 0;
		base_p[i]->flg_d[1] = 0;
		base_p[i]->flg_i[0] = 0;
		base_p[i]->flg_i[1] = 0;
		base_p[i]->flg_i[2] = 0;
		base_p[i]->flg_i[3] = 0;
	}



}

void read_basetrack2(std::string file_in_base, std::multimap<int, mfile0::M_Base>& base_map) {


	std::map<std::pair<int, int>, mfile0::M_Base*>base_map_id;
	for (auto itr = base_map.begin(); itr != base_map.end(); itr++) {
		base_map_id.insert(std::make_pair(std::make_pair(itr->first, itr->second.rawid), &(itr->second)));
	}
	//std::ofstream ofs("U:\\Linklet\\use_base2.bin", std::ios::binary);

	std::ifstream ifs(file_in_base, std::ios::binary);
	//filesize取得
	ifs.seekg(0, std::ios::end);
	int64_t eofpos = ifs.tellg();
	ifs.clear();
	ifs.seekg(0, std::ios::beg);
	int64_t begpos = ifs.tellg();
	int64_t nowpos = ifs.tellg();
	int64_t size2 = eofpos - begpos;
	int64_t GB = size2 / (1000 * 1000 * 1000);
	int64_t MB = (size2 - GB * 1000 * 1000 * 1000) / (1000 * 1000);
	int64_t KB = (size2 - GB * 1000 * 1000 * 1000 - MB * 1000 * 1000) / (1000);
	if (GB > 0) {
		std::cout << "FILE size :" << GB << "." << MB << " [GB]" << std::endl;
	}
	else {
		std::cout << "FILE size :" << MB << "." << KB << " [MB]" << std::endl;
	}
	int64_t count = 0;
	netscan::base_track_t base;
	while (ifs.read((char*)&base, sizeof(netscan::base_track_t))) {
		if (count % 10000 == 0) {
			nowpos = ifs.tellg();
			auto size1 = nowpos - begpos;
			std::cerr << std::right << std::fixed << "\r now reading ..." << std::setw(4) << std::setprecision(1) << size1 * 100. / size2 << "%";
		}
		count++;

		auto b = base_map_id.find(std::make_pair(base.pl, base.rawid));
		if (b == base_map_id.end())continue;
		b->second->ax = base.ax;
		b->second->ay = base.ay;
		b->second->x = base.x;
		b->second->y = base.y;
		b->second->z = base.z;
		b->second->ph = base.m[0].ph + base.m[1].ph;
		b->second->flg_d[0] = 0;
		b->second->flg_d[1] = 0;
		b->second->flg_i[0] = 0;
		b->second->flg_i[1] = 0;
		b->second->flg_i[2] = 0;
		b->second->flg_i[3] = 0;
		//ofs.write((char*)& base, sizeof(netscan::base_track_t));

	}
	auto size1 = eofpos - begpos;
	std::cerr << "\r now reading ..." << std::setw(4) << std::setprecision(1) << size1 * 100. / size2 << "%" << std::endl;;
	if (count == 0) {
		fprintf(stderr, "%s no linklet!\n", file_in_base.c_str());
		exit(1);
	}

}
void apply_corrmap(std::multimap<int, mfile0::M_Base>& base_map, std::map<int, corrmap0::Corrmap>& corr_abs, std::map<int, double>& z_map) {


	double tmp_x, tmp_y;

	for (auto itr = base_map.begin(); itr != base_map.end(); itr++) {
		int count = base_map.count(itr->first);
		auto range = base_map.equal_range(itr->first);
		auto z = z_map.find(itr->first);
		auto param = corr_abs.find(itr->first);
		if (z == z_map.end()) {
			printf("exception PL%03d z not found\n", itr->first);
		}
		if (param == corr_abs.end()) {
			printf("exception PL%03d corr abs not found\n", itr->first);
		}

		for (auto res = range.first; res != range.second; res++) {
			tmp_x = res->second.x;
			tmp_y = res->second.y;
			res->second.x = param->second.position[0] * tmp_x + param->second.position[1] * tmp_y + param->second.position[4];
			res->second.y = param->second.position[2] * tmp_x + param->second.position[3] * tmp_y + param->second.position[5];
			tmp_x = res->second.ax;
			tmp_y = res->second.ay;
			res->second.ax = param->second.angle[0] * tmp_x + param->second.angle[1] * tmp_y + param->second.angle[4];
			res->second.ay = param->second.angle[2] * tmp_x + param->second.angle[3] * tmp_y + param->second.angle[5];

			res->second.z = z->second + param->second.dz;

		}
		itr = std::next(itr, count - 1);
	}
}

void apply_corrmap_local(std::multimap<int, mfile0::M_Base>& base_map, std::map<int, std::vector<corrmap_3d::align_param2>>& corr) {

	//ここを書く
	std::set<int> all_pl_set;
	for (auto itr = base_map.begin(); itr != base_map.end(); itr++) {
		all_pl_set.insert(itr->first);
	}
	std::vector<int> all_pl;
	for (auto itr = all_pl_set.begin(); itr != all_pl_set.end(); itr++) {
		all_pl.push_back(*itr);
	}

	int all = all_pl.size(), count = 0;
#pragma omp parallel for num_threads(use_thread(0.4,true)) schedule(dynamic,1)
	for (int i = 0; i < all_pl.size(); i++) {
		//for (auto itr = base_map_single.begin(); itr != base_map_single.end(); itr++) {
			//count = base_map_single.count(itr->first);
		int pl = all_pl[i];
#pragma omp critical
		{
			printf("PL%03d basetrack tans %3d/%3d\n", pl, count, all);
			count++;
			if (corr.count(pl) == 0) {
				fprintf(stderr, "PL%03d corrmap not found\n", pl);
				exit(1);
			}
			if (base_map.count(pl) == 0) {
				fprintf(stderr, "PL%03d basetrack not found\n", pl);
				exit(1);
			}
		}

		std::vector<corrmap_3d::align_param2> param = corr.at(pl);
		auto range = base_map.equal_range(pl);
		std::vector< mfile0::M_Base*> base_trans;
		base_trans.reserve(base_map.count(pl));
		for (auto itr = range.first; itr != range.second; itr++) {
			base_trans.push_back(&(itr->second));
		}
		std::vector <std::pair<mfile0::M_Base*, corrmap_3d::align_param2*>> base_trans_map = corrmap_3d::track_affineparam_correspondence(base_trans, param);
		trans_base_all(base_trans_map);

	}
}
std::set<std::pair<int, int>> edge_point(std::set<std::tuple<int, int, int, int>>& path) {
	std::set<std::pair<int, int>> base0, base1;
	for (auto itr = path.begin(); itr != path.end(); itr++) {
		base0.insert(std::make_pair(std::get<0>(*itr), std::get<2>(*itr)));
		base1.insert(std::make_pair(std::get<1>(*itr), std::get<3>(*itr)));
	}
	std::set<std::pair<int, int>> ret;
	for (auto itr = base0.begin(); itr != base0.end(); itr++) {
		if (base1.count(*itr) == 1)continue;
		ret.insert(*itr);
	}
	return ret;
}

std::vector<std::vector<std::pair<mfile0::M_Base, mfile0::M_Base>>> divide_connected(output_format& g, std::map<std::pair<int, int>, mfile0::M_Base>& base) {
	std::multimap<int, std::pair<mfile0::M_Base, mfile0::M_Base>>connected_id;
	for (auto itr = g.select_path.begin(); itr != g.select_path.end(); itr++) {
		std::pair<int, int> id0, id1;
		id0 = std::make_pair(std::get<0>(itr->second) / 10, std::get<2>(itr->second));
		id1 = std::make_pair(std::get<1>(itr->second) / 10, std::get<3>(itr->second));
		auto base0 = base.find(id0);
		auto base1 = base.find(id1);
		if (base0 == base.end() || base1 == base.end()) {
			printf("basetrack not found\n");
		}
		connected_id.insert(std::make_pair(itr->first, std::make_pair(base0->second, base1->second)));
	}
	std::vector<std::vector<std::pair<mfile0::M_Base, mfile0::M_Base>>> ret;
	for (auto itr = connected_id.begin(); itr != connected_id.end(); itr++) {
		int count = connected_id.count(itr->first);
		std::vector<std::pair<mfile0::M_Base, mfile0::M_Base>> ret_v;
		auto range = connected_id.equal_range(itr->first);
		for (auto res = range.first; res != range.second; res++) {
			ret_v.push_back(res->second);
		}
		ret.push_back(ret_v);

		itr = std::next(itr, count - 1);
	}
	return ret;
}

std::vector<std::vector<std::pair<mfile0::M_Base, mfile0::M_Base>>> divide_connected(std::vector<std::pair<int, std::tuple<int, int, int, int>>>& path_list, std::map<std::pair<int, int>, mfile0::M_Base>& base) {
	//idを使わない
	std::multimap<mfile0::M_Base, mfile0::M_Base>connected_list;
	std::set<mfile0::M_Base>all_edge;
	for (auto itr = path_list.begin(); itr != path_list.end(); itr++) {
		std::pair<int, int> id0, id1;
		id0 = std::make_pair(std::get<0>(itr->second) / 10, std::get<2>(itr->second));
		id1 = std::make_pair(std::get<1>(itr->second) / 10, std::get<3>(itr->second));
		auto base0 = base.find(id0);
		auto base1 = base.find(id1);
		if (base0 == base.end() || base1 == base.end()) {
			printf("basetrack not found\n");
		}
		connected_list.insert(std::make_pair(base0->second, base1->second));
		connected_list.insert(std::make_pair(base1->second, base0->second));
		all_edge.insert(base0->second);
		all_edge.insert(base1->second);
	}
	std::vector<std::vector<std::pair<mfile0::M_Base, mfile0::M_Base>>> ret;
	std::set<mfile0::M_Base>finished_edge;
	bool loop_flg = true;
	while (true) {
		bool loop_flg = false;
		for (auto itr = all_edge.begin(); itr != all_edge.end(); itr++) {
			if (finished_edge.count(*itr) == 1) continue;
			loop_flg = true;

			std::set<std::pair<mfile0::M_Base, mfile0::M_Base>> connect_part;
			std::set < mfile0::M_Base> search_list;
			std::set < mfile0::M_Base> add_search_list;
			search_list.insert(*itr);
			while (true) {
				for (auto itr2 = search_list.begin(); itr2 != search_list.end(); itr2++) {
					if (connected_list.count(*itr2) == 0) {
						finished_edge.insert(*itr2);
						continue;
					}
					auto range = connected_list.equal_range(*itr2);
					for (auto res = range.first; res != range.second; res++) {
						auto edge0 = *itr2;
						auto edge1 = res->second;
						if (edge0.pos > edge1.pos)std::swap(edge0, edge1);
						connect_part.insert(std::make_pair(edge0, edge1));
						if (finished_edge.count(res->second) == 0) {
							add_search_list.insert(res->second);
						}
					}
					finished_edge.insert(*itr2);
				}
				search_list = add_search_list;
				if (add_search_list.size() == 0)break;
				add_search_list.clear();
			}
			std::vector<std::pair<mfile0::M_Base, mfile0::M_Base>> connect_part_v(connect_part.begin(), connect_part.end());
			ret.push_back(connect_part_v);
		}
		if (!loop_flg)break;
	}

	//for (auto &vv : ret) {
	//	for (auto &v : vv) {
	//		Print_path(v.first, v.second);
	//		printf("\n");
	//	}
	//	printf("\n");
	//}

	return ret;

}

bool chain_matching(std::vector<mfile0::M_Base>& c0, std::vector<mfile0::M_Base>& c1) {
	std::set<std::pair<int, int>> c0_id;
	for (auto itr = c0.begin(); itr != c0.end(); itr++) {
		c0_id.insert(std::make_pair(itr->pos, itr->rawid));
	}
	for (auto itr = c1.begin(); itr != c1.end(); itr++) {
		if (c0_id.count(std::make_pair(itr->pos, itr->rawid)) == 1)return true;
	}
	return false;
}
std::vector<std::pair<mfile0::M_Base, mfile0::M_Base>> select_path(std::vector<std::pair<mfile0::M_Base, mfile0::M_Base>>& path, std::map<std::pair<int, int>, mfile0::M_Base>& base_map_id) {

	std::set<std::pair<int32_t, int64_t>> calc_path_btset;
	std::set<std::tuple<int, int, int, int>> set_calc_path_ltlist;
	std::vector<Linklet> calc_path_ltlist;
	std::set<int32_t> set_calc_path_usepos;
	std::vector<int32_t> calc_path_usepos;


	for (auto itr = path.begin(); itr != path.end(); itr++) {
		set_calc_path_ltlist.insert(std::make_tuple(itr->first.pos, itr->second.pos, itr->first.rawid, itr->second.rawid));
		set_calc_path_usepos.insert(itr->first.pos);
		set_calc_path_usepos.insert(itr->second.pos);
		calc_path_btset.insert(std::make_pair(itr->first.pos, itr->first.rawid));
		calc_path_btset.insert(std::make_pair(itr->second.pos, itr->second.rawid));
	}
	calc_path_usepos.reserve(set_calc_path_usepos.size());
	for (auto itr = set_calc_path_usepos.begin(); itr != set_calc_path_usepos.end(); itr++) {
		calc_path_usepos.push_back(*itr);
	}
	calc_path_ltlist.reserve(set_calc_path_ltlist.size());
	for (auto itr = set_calc_path_ltlist.begin(); itr != set_calc_path_ltlist.end(); itr++) {
		calc_path_ltlist.emplace_back(std::get<0>(*itr), std::get<1>(*itr), std::get<2>(*itr), std::get<3>(*itr));
	}
	std::set<std::pair<int, int>> edge_base = edge_point(set_calc_path_ltlist);

	l2c::Cdat cdat_all_path = l2c_x(calc_path_btset, calc_path_ltlist, calc_path_usepos);
	size_t grsize = cdat_all_path.GetNumOfGroups();
	size_t possize = calc_path_usepos.size();
	/*
	std::vector<std::vector<std::vector<mfile0::M_Base>>> chain_combination;
	for (size_t grid = 0; grid < grsize; ++grid)
	{
		const Group& gr = cdat_all_path.GetGroup(grid);
		size_t chsize = gr.GetNumOfChains();
		if (grid != gr.GetID()) throw std::exception("grid != gr.GetID().");//gridはgr.GetID()の戻り値と基本的に等しいはずである。
		int32_t spl = gr.GetStartPL();
		int32_t epl = gr.GetEndPL();
		if (gr.IsOverUpperLim())
		{
			//upperlimを超過している場合、chainの情報はない。
			//ただしchainの本数はGetNumOfChainsで正しく取得できる。
			//また属すBaseTrackの情報は保持しているので、GetBaseTracksで全BaseTrackを取得できる。
			fprintf(stdout, "over upperlim. nchain:%-9lld\n", chsize);
			exit(1);
		}
		else {
			for (size_t ich = 0; ich < chsize; ++ich)
			{
				std::vector<mfile0::M_Base>chain;
				Chain ch = gr.GetChain(ich);//仕様上、Chainオブジェクトは一時変数として戻ってくる。
				for (size_t pl = 0; pl < possize; ++pl)
				{
					BaseTrackID bt = ch.GetBaseTrack(pl);
					if (bt.IsEmpty()) continue;//IsEmptyがtrueのときは、このBaseTrackは歯抜け。
					int32_t btpl = bt.GetPL();
					int64_t btid = bt.GetRawID();
					if (btpl != pl) throw std::exception("BaseTrackID::PL != pl.");
					if (calc_path_btset.find(std::make_pair(calc_path_usepos[btpl], btid)) == calc_path_btset.end())
					{
						//エラーチェック。
						//Chain内の全てのBaseTrackはbtset（ファイルから読み込まれた全BaseTrack）に必ず含まれているはず。
						//含まれていなければエラー。
						throw std::exception("BaseTrack is not found in btset.");
					}
					auto res = base_map_id.find(std::make_pair(calc_path_usepos[btpl] / 10, btid));
					if (res == base_map_id.end()) {
						printf("PL%03d %d not found\n", calc_path_usepos[btpl] / 10, btid);
					}
					chain.push_back(res->second);
				}
				bool flg = true;
				bool input_flg = false;
				for (int i = 0; i < chain_combination.size(); i++) {
					if (chain_combination[i].size() == edge_base.size())continue;
					flg = true;
					for (int j = 0; j < chain_combination[i].size(); j++) {
						//chain_combination[i][j]とchainの比較
						if (!flg)break;
						if (!chain_matching(chain_combination[i][j], chain)) {
							flg = false;
						}
					}
					if (!flg) {
						chain_combination[i].push_back(chain);
						input_flg = true;
						break;
					}
				}
				if (!input_flg) {
					std::vector<std::vector<mfile0::M_Base>> chain_v;
					chain_v.push_back(chain);
					chain_combination.push_back(chain_v);
				}
			}
		}
	}

	for (int i = 0; i < chain_combination.size(); i++) {
		for (int j = 0; j < chain_combination[i].size(); j++) {
			for (auto itr = chain_combination[i][j].begin(); itr != chain_combination[i][j].end(); itr++) {
				printf("(%d %d) ", itr->pos, itr->rawid);
			}
			printf("\n");
		}
		printf("\n");
	}
	*/
	std::vector<std::pair<mfile0::M_Base, mfile0::M_Base>> ret;
	return ret;
}



//使用関数群
double minimum_distance(std::set<mfile0::M_Base>& base) {
	if (base.size() < 2)return -1;
	double dist = -1;
	double dist_tmp;
	for (auto itr = base.begin(); itr != base.end(); itr++) {
		for (auto itr2 = std::next(itr, 1); itr2 != base.end(); itr2++) {
			//printf("%d %.1lf %d %.1lf\n", itr->pos, itr->z, itr2->pos, itr2->z);
			double dx = itr->x - itr2->x + (itr->ax + itr2->ax) / 2 * (itr->z - itr2->z);
			double dy = itr->y - itr2->y + (itr->ay + itr2->ay) / 2 * (itr->z - itr2->z);
			dx = itr->x - itr2->x;
			dy = itr->y - itr2->y;
			dist_tmp = sqrt(pow(dx, 2) + pow(dy, 2));
			if (dist<0 || dist>dist_tmp) {
				dist = dist_tmp;
			}
		}
	}
	return dist;
}
std::vector<Chain_path> select_start(std::set<mfile0::M_Base>& all_vertex, boost::unordered_multimap<mfile0::M_Base, mfile0::M_Base>& path_next, boost::unordered_multimap<mfile0::M_Base, mfile0::M_Base>& path_prev, bool& to_next) {

	std::set<mfile0::M_Base> up_base, down_base;
	for (auto itr = all_vertex.begin(); itr != all_vertex.end(); itr++) {
		if (path_next.count(*itr) == 0)down_base.insert(*itr);
		if (path_prev.count(*itr) == 0)up_base.insert(*itr);
	}
	//デバッグ用
	//if (minimum_distance(up_base) > minimum_distance(down_base)) {
	if (minimum_distance(up_base) < minimum_distance(down_base)) {
		to_next = false;
		std::vector<Chain_path> ret;
		for (auto itr = down_base.begin(); itr != down_base.end(); itr++) {
			auto search_down_base = *itr;
			if (path_prev.count(search_down_base) == 0) {
				fprintf(stderr, "exception function [select_start]\n");
				exit(1);
			}
			Chain_path path_tmp;
			//printf("start (%d,%d) %d\n", search_down_base.pos, search_down_base.rawid, path_prev.count(search_down_base));
			while (path_prev.count(search_down_base) == 1) {
				//printf("(%d,%d) %d\n", search_down_base.pos, search_down_base.rawid, path_prev.count(search_down_base));
				path_tmp.Add_Basetrack(search_down_base);
				auto res = path_prev.find(search_down_base);
				search_down_base = res->second;
			}
			path_tmp.Add_Basetrack(search_down_base);
			ret.push_back(path_tmp);
		}
		return ret;
	}

	to_next = true;
	std::vector<Chain_path> ret;
	for (auto itr = up_base.begin(); itr != up_base.end(); itr++) {
		auto search_up_base = *itr;
		if (path_next.count(search_up_base) == 0) {
			fprintf(stderr, "exception function [select_start]\n");
			exit(1);
		}
		Chain_path path_tmp;
		while (path_next.count(search_up_base) == 1) {
			path_tmp.Add_Basetrack(search_up_base);
			auto res = path_next.find(search_up_base);
			search_up_base = res->second;
		}
		path_tmp.Add_Basetrack(search_up_base);
		ret.push_back(path_tmp);
	}
	return ret;
}
std::map<std::pair<mfile0::M_Base, mfile0::M_Base>, Chain_path> path_compress_next(std::set<mfile0::M_Base>& all_vertex, boost::unordered_multimap<mfile0::M_Base, mfile0::M_Base>& path_next, boost::unordered_multimap<mfile0::M_Base, mfile0::M_Base>& path_prev)
{
	std::map<std::pair<mfile0::M_Base, mfile0::M_Base>, Chain_path>  ret;

	std::multimap<int, Chain_path> chain_path_vec;

	std::set<mfile0::M_Base> finished_vertex;
	std::set<std::pair<mfile0::M_Base, mfile0::M_Base>>finished_path;
	for (auto itr = all_vertex.begin(); itr != all_vertex.end(); itr++) {
		//上から下に。 上が一本の場合は上の時にやる
		if (path_next.count(*itr) == 1 && path_prev.count(*itr) == 1)continue;
		if (path_next.count(*itr) == 0)continue;
		auto range = path_next.equal_range(*itr);
		for (auto res = range.first; res != range.second; res++) {
			mfile0::M_Base now = *itr;
			std::set<mfile0::M_Base> path;
			//下まで続ける
			path.insert(now);
			now = res->second;
			path.insert(now);
			while (path_next.count(now) == 1) {
				now = path_next.find(now)->second;
				path.insert(now);
			}
			auto res_fin = finished_path.insert(std::make_pair(*path.begin(), *path.rbegin()));
			if (res_fin.second) {
				auto input_path = Chain_path(path);
				chain_path_vec.insert(std::make_pair(input_path.Get_all_base().size(), input_path));
			}
		}
	}

	//chain_path_vecの圧縮-->ret.insert(std::make_pair(edge, input_path));
	for (auto itr = chain_path_vec.rbegin(); itr != chain_path_vec.rend(); itr++) {
		auto input_all_base = itr->second.Get_all_base();
		bool flg = false;
		for (auto itr2 = ret.begin(); itr2 != ret.end(); itr2++) {
			auto all_base = itr2->second.Get_all_base();
			std::set<mfile0::M_Base> all_base_set;
			for (auto itr3 = all_base.begin(); itr3 != all_base.end(); itr3++) {
				all_base_set.insert(*itr3);
			}

			flg = true;
			for (auto itr3 = input_all_base.begin(); itr3 != input_all_base.end(); itr3++) {
				if (all_base_set.count(*itr3) == 0)flg = false;
			}
			if (flg)break;
		}
		if (!flg) {
			auto edge = itr->second.Get_path_edge();
			ret.insert(std::make_pair(edge, itr->second));
		}

	}

	path_next.clear();
	path_prev.clear();
	all_vertex.clear();
	for (auto itr = ret.begin(); itr != ret.end(); itr++) {
		auto edge = itr->second.Get_path_edge();
		path_next.insert(edge);
		path_prev.insert(std::make_pair(edge.second, edge.first));
		all_vertex.insert(edge.first);
		all_vertex.insert(edge.second);
	}
	//printf("path next\n");
	//for (auto itr = path_next.begin(); itr != path_next.end(); itr++) {
	//	printf("%d,%d -- %d,%d\n", itr->first.pos, itr->first.rawid, itr->second.pos, itr->second.rawid);
	//}
	//printf("path prev\n");
	//for (auto itr = path_prev.begin(); itr != path_prev.end(); itr++) {
	//	printf("%d,%d -- %d,%d\n", itr->first.pos, itr->first.rawid, itr->second.pos, itr->second.rawid);
	//}

	return ret;

}
std::map<std::pair<mfile0::M_Base, mfile0::M_Base>, Chain_path> path_compress_prev(std::set<mfile0::M_Base>& all_vertex, boost::unordered_multimap<mfile0::M_Base, mfile0::M_Base>& path_next, boost::unordered_multimap<mfile0::M_Base, mfile0::M_Base>& path_prev)
{
	std::map<std::pair<mfile0::M_Base, mfile0::M_Base>, Chain_path>  ret;

	std::multimap<int, Chain_path> chain_path_vec;

	std::set<mfile0::M_Base> finished_vertex;
	std::set<std::pair<mfile0::M_Base, mfile0::M_Base>>finished_path;
	for (auto itr = all_vertex.begin(); itr != all_vertex.end(); itr++) {
		//下から上に。 下が一本の場合は下の時にやる
		if (path_next.count(*itr) == 1 && path_prev.count(*itr) == 1)continue;
		if (path_prev.count(*itr) == 0)continue;
		auto range = path_prev.equal_range(*itr);
		for (auto res = range.first; res != range.second; res++) {
			mfile0::M_Base now = *itr;
			std::set<mfile0::M_Base> path;
			//下まで続ける
			path.insert(now);
			now = res->second;
			path.insert(now);
			while (path_prev.count(now) == 1) {
				now = path_prev.find(now)->second;
				path.insert(now);
			}
			auto res_fin = finished_path.insert(std::make_pair(*path.begin(), *path.rbegin()));
			if (res_fin.second) {
				auto input_path = Chain_path(path);
				chain_path_vec.insert(std::make_pair(input_path.Get_all_base().size(), input_path));
			}
		}
	}

	//chain_path_vecの圧縮-->ret.insert(std::make_pair(edge, input_path));
	for (auto itr = chain_path_vec.rbegin(); itr != chain_path_vec.rend(); itr++) {
		auto input_all_base = itr->second.Get_all_base();
		bool flg = false;
		for (auto itr2 = ret.begin(); itr2 != ret.end(); itr2++) {
			auto all_base = itr2->second.Get_all_base();
			std::set<mfile0::M_Base> all_base_set;
			for (auto itr3 = all_base.begin(); itr3 != all_base.end(); itr3++) {
				all_base_set.insert(*itr3);
			}

			flg = true;
			for (auto itr3 = input_all_base.begin(); itr3 != input_all_base.end(); itr3++) {
				if (all_base_set.count(*itr3) == 0)flg = false;
			}
			if (flg)break;
		}
		if (!flg) {
			auto edge = itr->second.Get_path_edge();
			ret.insert(std::make_pair(edge, itr->second));
		}

	}

	path_next.clear();
	path_prev.clear();
	all_vertex.clear();
	for (auto itr = ret.begin(); itr != ret.end(); itr++) {
		auto edge = itr->second.Get_path_edge();
		path_next.insert(edge);
		path_prev.insert(std::make_pair(edge.second, edge.first));
		all_vertex.insert(edge.first);
		all_vertex.insert(edge.second);
	}
	//printf("path next\n");
	//for (auto itr = path_next.begin(); itr != path_next.end(); itr++) {
	//	printf("%d,%d -- %d,%d\n", itr->first.pos, itr->first.rawid, itr->second.pos, itr->second.rawid);
	//}
	//printf("path prev\n");
	//for (auto itr = path_prev.begin(); itr != path_prev.end(); itr++) {
	//	printf("%d,%d -- %d,%d\n", itr->first.pos, itr->first.rawid, itr->second.pos, itr->second.rawid);
	//}

	return ret;

}

bool judge_loop_fin(std::set<Chain_path>& selected_path, std::vector < std::pair<Chain_path, Chain_path>>& path_pair) {
	std::vector < std::pair<Chain_path, Chain_path>> ret;
	for (int i = 0; i < path_pair.size(); i++) {
		//片方のpathが選択済みであればスルー
		if (selected_path.count(path_pair[i].first) == 1)continue;
		if (selected_path.count(path_pair[i].second) == 1)continue;
		ret.push_back(path_pair[i]);
	}
	path_pair = ret;
	return path_pair.size() != 0;

}
bool result_cross_path_next(std::vector<Chain_path>& path, std::map<std::pair<mfile0::M_Base, mfile0::M_Base>, Chain_path>& path_map, boost::unordered_multimap<mfile0::M_Base, mfile0::M_Base>& path_next, boost::unordered_multimap<mfile0::M_Base, mfile0::M_Base>& path_prev) {
	std::vector<Chain_path> ret;
	std::map<mfile0::M_Base, std::set<mfile0::M_Base>> next_branch;
	for (auto itr = path.begin(); itr != path.end(); itr++) {
		auto edge = itr->Get_path_edge();
		std::set<mfile0::M_Base> br;
		//終わり
		if (path_next.count(edge.second) == 0)continue;
		auto range = path_next.equal_range(edge.second);
		for (auto res = range.first; res != range.second; res++) {
			br.insert(res->second);
		}
		next_branch.insert(std::make_pair(edge.second, br));
	}
	if (next_branch.size() == 0)return false;
	else if (next_branch.size() == 1) {
		if (next_branch.begin()->second.size() == 1) {
			mfile0::M_Base start = next_branch.begin()->first;
			mfile0::M_Base end = *next_branch.begin()->second.begin();
			for (int i = 0; i < path.size(); i++) {
				auto edge = path[i].Get_path_edge();
				if (edge.second == start) {
					auto res = path_map.find(std::make_pair(start, end));
					if (res == path_map.end()) {
						printf("error function[result_cross_path_next]\n");
						exit(1);
					}
					auto add_base = res->second.Get_all_base();
					for (auto itr = add_base.begin(); itr != add_base.end(); itr++) {
						path[i].Add_Basetrack(*itr);
					}
					return true;
				}
			}
			printf("error function[result_cross_path_next]\n");
			exit(1);
		}
		else {
			printf("error function[result_cross_path_next]\n");
			exit(1);
		}
	}

	//printf("next branch\n");
	//for (auto itr = next_branch.begin(); itr != next_branch.end(); itr++) {
	//	printf("(%d,%d)\n", itr->first.pos, itr->first.rawid);
	//	for (auto itr2 = itr->second.begin(); itr2 != itr->second.end(); itr2++) {
	//		printf("\t(%d,%d)\n", itr2->pos, itr2->rawid);
	//	}
	//}

	std::vector<std::pair<std::set<mfile0::M_Base>, std::set<mfile0::M_Base>>> cross_path;
	for (auto itr = next_branch.begin(); itr != next_branch.end(); itr++) {
		for (auto itr2 = std::next(itr, 1); itr2 != next_branch.end(); itr2++) {
			//set の一致確認
			if (!Base_set_compare(itr->second, itr2->second))continue;
			auto fin_base_set = itr->second;
			auto start_base0 = itr->first;
			auto start_base1 = itr2->first;
			bool flg = false;
			for (int i = 0; i < cross_path.size(); i++) {
				if (!Base_set_compare(fin_base_set, cross_path[i].second))continue;
				cross_path[i].first.insert(start_base0);
				cross_path[i].first.insert(start_base1);
				flg = true;
				break;
			}
			if (!flg) {
				std::pair < std::set<mfile0::M_Base>, std::set<mfile0::M_Base>> cross_path_pair;
				cross_path_pair.first.insert(start_base0);
				cross_path_pair.first.insert(start_base1);
				cross_path_pair.second = fin_base_set;
				cross_path.push_back(cross_path_pair);
			}
		}
	}
	for (int i = 0; i < cross_path.size(); i++) {
		//i番目の完全2部グラフ部分の解決(基本的にはi=0のみ)
		//printf("cross path %d\n", i);
		if (false) {
			std::vector <mfile0::M_Base> edge_now(cross_path[i].first.begin(), cross_path[i].first.end());
			std::vector <mfile0::M_Base> edge_next(cross_path[i].second.begin(), cross_path[i].second.end());
			for (int j = 0; j < edge_now.size(); j++) {
				printf("(%d,%d)  (%d,%d)\n", edge_now[j].pos, edge_now[j].rawid, edge_next[j].pos, edge_next[j].rawid);
			}
		}
		//Kn,nになっていない-->まだ解決していない場所がある
		if (cross_path[i].first.size() != cross_path[i].second.size())continue;
		//マッチングの列挙
		int path_num = cross_path[i].first.size();
		//path_numが多い時は最良経路から選択する
		while (path_num > 5) {
			//未実装
			break;
		}
		int combination_num = 1;
		for (int j = 1; j <= path_num; j++) {
			combination_num *= j;
		}
		std::vector<std::vector<std::pair<Chain_path, Chain_path>>> all_matching;
		auto path0 = path[0];
		auto path1 = path[0];
		for (int j = 0; j < combination_num; j++) {
			std::vector<std::pair<Chain_path, Chain_path>> match_tmp;
			for (int k = 0; k < path_num; k++) {
				match_tmp.push_back(std::make_pair(path0, path1));
			}
			all_matching.push_back(match_tmp);
		}
		std::vector <mfile0::M_Base> edge_now(cross_path[i].first.begin(), cross_path[i].first.end());
		std::vector <mfile0::M_Base> edge_next(cross_path[i].second.begin(), cross_path[i].second.end());
		int count_permutation = 0;
		for (int j = 0; j < edge_now.size(); j++) {
			for (int k = 0; k < path.size(); k++) {
				if (path[k].Get_path_edge().second == edge_now[j]) {
					path0 = path[k];
					break;
				}
			}
			for (int k = 0; k < all_matching.size(); k++) {
				all_matching[k][j].first = path0;
			}
		}
		sort(edge_next.begin(), edge_next.end());
		count_permutation = 0;
		do {
			for (int j = 0; j < edge_next.size(); j++) {
				//printf("%d %d\n", count_permutation, j);
				auto res = path_map.find(std::make_pair(edge_now[j], edge_next[j]));
				if (res == path_map.end()) {
					printf("error function[result_cross_path_next]\n");
					exit(1);
				}
				auto path1 = res->second;
				all_matching[count_permutation][j].second = path1;
			}
			count_permutation++;
		} while (next_permutation(edge_next.begin(), edge_next.end()));

		//for (int j = 0; j < all_matching.size(); j++) {
		//	printf("combination %d\n", j);
		//	for (int k = 0; k < all_matching[j].size(); k++) {
		//		auto edge0 = all_matching[j][k].first.Get_path_edge();
		//		auto edge1 = all_matching[j][k].second.Get_path_edge();

		//		printf("%d,%d: (%d,%d)-(%d,%d) ==  (%d,%d)-(%d,%d)\n", j, k,
		//			edge0.first.pos, edge0.first.rawid, edge0.second.pos, edge0.second.rawid,
		//			edge1.first.pos, edge1.first.rawid, edge1.second.pos, edge1.second.rawid
		//		);
		//	}
		//}

		//basetrack --> 経路への変換 + 重みの計算
		std::map < std::pair<Chain_path, Chain_path>, double> path_pair;
		for (auto itr_now = cross_path[i].first.begin(); itr_now != cross_path[i].first.end(); itr_now++) {
			for (auto itr_next = cross_path[i].second.begin(); itr_next != cross_path[i].second.end(); itr_next++) {
				auto res = path_map.find(std::make_pair(*itr_now, *itr_next));
				if (res == path_map.end()) {
					printf("error function[result_cross_path_next]\n");
					exit(1);
				}
				for (int j = 0; j < path.size(); j++) {
					if (path[j].Get_path_edge().second == *itr_now) {
						//ここで角度差の計算までする
						auto path0 = path[j];
						auto path1 = res->second;
						auto edge0 = path0.Get_path_edge();
						auto edge1 = path1.Get_path_edge();
						path0.Line_Fit(edge0.second.pos / 10 - 5, edge0.second.pos / 10);
						path1.Line_Fit(edge1.first.pos / 10, edge1.first.pos / 10 + 5);
						double angle_diff = sqrt(pow(path0.Get_line_ax() - path1.Get_line_ax(), 2) + pow(path0.Get_line_ay() - path1.Get_line_ay(), 2));
						path_pair.insert(std::make_pair(std::make_pair(path0, path1), angle_diff));
					}
				}
			}
		}

		std::vector<std::pair<Chain_path, Chain_path>> best_path_combination;
		double val_min = -1;
		for (int j = 0; j < all_matching.size(); j++) {
			//printf("combination %d\n", j);
			double val_tmp = 0;
			for (int k = 0; k < all_matching[j].size(); k++) {
				auto edge0 = all_matching[j][k].second.Get_path_edge();

				//auto path0 = all_matching[j][k].first;
				//auto path1 = all_matching[j][k].second;
				//edge0 = path0.Get_path_edge();
				//auto edge1 = path1.Get_path_edge();
				//printf("(%d,%d)-(%d-%d) == (%d,%d)-(%d-%d)\n",
				//	edge0.first.pos, edge0.first.rawid, edge0.second.pos, edge0.second.rawid,
				//	edge1.first.pos, edge1.first.rawid, edge1.second.pos, edge1.second.rawid);

				//auto edge_out = all_matching[j][k].second.Get_path_edge();
				//printf("%d %d %d %d\n", edge_out.first.pos, edge_out.first.rawid, edge_out.second.pos, edge_out.second.rawid);
				auto res = path_pair.find(all_matching[j][k]);
				if (res == path_pair.end()) {
					printf("error function[result_cross_path_next]\n");
					exit(1);
				}
				val_tmp += res->second;
				//printf("(%4d,%10d) - (%4d,%10d) : %.4lf\n", edge0.first.pos, edge0.first.rawid, edge0.second.pos, edge0.second.rawid, res->second);
			}
			if (val_min < 0 || val_tmp < val_min) {
				val_min = val_tmp;
				best_path_combination = all_matching[j];
			}
		}

		for (int j = 0; j < best_path_combination.size(); j++) {
			for (auto& b : best_path_combination[j].second.Get_all_base()) {
				best_path_combination[j].first.Add_Basetrack(b);
			}
			ret.push_back(best_path_combination[j].first);
		}
	}

	//最終結果
	//for (int i = 0; i < path.size(); i++) {
	//	auto edge = path[i].Get_path_edge();
	//	printf("(%d,%d) - (%d,%d)\n", edge.first.pos, edge.first.rawid, edge.second.pos, edge.second.rawid);
	//}
	//printf("\n");
	//for (int i = 0; i < ret.size(); i++) {
	//	auto edge = ret[i].Get_path_edge();
	//	printf("(%d,%d) - (%d,%d)\n", edge.first.pos, edge.first.rawid, edge.second.pos, edge.second.rawid);
	//}
	//printf("\n");
	//printf("fin\n");

	for (auto itr = path.begin(); itr != path.end(); itr++) {
		auto edge0 = itr->Get_path_edge();
		bool flg = true;
		for (auto itr1 = ret.begin(); itr1 != ret.end(); itr1++) {
			auto edge1 = itr1->Get_path_edge();
			if (edge0.first == edge1.first)flg = false;
		}
		if (flg)ret.push_back(*itr);
	}

	path = ret;
	return true;
}
bool result_cross_path_prev(std::vector<Chain_path>& path, std::map<std::pair<mfile0::M_Base, mfile0::M_Base>, Chain_path>& path_map, boost::unordered_multimap<mfile0::M_Base, mfile0::M_Base>& path_next, boost::unordered_multimap<mfile0::M_Base, mfile0::M_Base>& path_prev) {

	std::vector<Chain_path> ret;
	std::map<mfile0::M_Base, std::set<mfile0::M_Base>> prev_branch;
	for (auto itr = path.begin(); itr != path.end(); itr++) {
		auto edge = itr->Get_path_edge();
		std::set<mfile0::M_Base> br;
		//終わり
		if (path_prev.count(edge.first) == 0)continue;
		auto range = path_prev.equal_range(edge.first);
		for (auto res = range.first; res != range.second; res++) {
			br.insert(res->second);
		}
		prev_branch.insert(std::make_pair(edge.first, br));
	}
	if (prev_branch.size() == 0)return false;
	else if (prev_branch.size() == 1) {
		if (prev_branch.begin()->second.size() == 1) {
			//start pl < end pl
			mfile0::M_Base end = prev_branch.begin()->first;
			mfile0::M_Base start = *prev_branch.begin()->second.begin();
			for (int i = 0; i < path.size(); i++) {
				auto edge = path[i].Get_path_edge();
				if (edge.first == end) {
					auto res = path_map.find(std::make_pair(start, end));
					if (res == path_map.end()) {
						printf("error function[result_cross_path_prev]\n");
						exit(1);
					}
					auto add_base = res->second.Get_all_base();
					for (auto itr = add_base.begin(); itr != add_base.end(); itr++) {
						path[i].Add_Basetrack(*itr);
					}
					return true;
				}
			}
			printf("error function[result_cross_path_prev]\n");
			exit(1);
		}
		else {
			printf("error function[result_cross_path_prev]\n");
			exit(1);
		}
	}
	////ここまで
	//for (auto itr = prev_branch.begin(); itr != prev_branch.end(); itr++) {
	//	for (auto itr2 = itr->second.begin(); itr2 != itr->second.end(); itr2++) {
	//		Print_path(itr->first, *itr2);
	//		printf("\n");
	//	}
	//}
	//printf("-----------\n");


	std::vector<std::pair<std::set<mfile0::M_Base>, std::set<mfile0::M_Base>>> cross_path;
	for (auto itr = prev_branch.begin(); itr != prev_branch.end(); itr++) {
		for (auto itr2 = std::next(itr, 1); itr2 != prev_branch.end(); itr2++) {
			//set の一致確認
			if (!Base_set_compare(itr->second, itr2->second))continue;
			auto fin_base_set = itr->second;
			auto start_base0 = itr->first;
			auto start_base1 = itr2->first;
			bool flg = false;
			for (int i = 0; i < cross_path.size(); i++) {
				if (!Base_set_compare(fin_base_set, cross_path[i].second))continue;
				cross_path[i].first.insert(start_base0);
				cross_path[i].first.insert(start_base1);
				flg = true;
				break;
			}
			if (!flg) {
				std::pair < std::set<mfile0::M_Base>, std::set<mfile0::M_Base>> cross_path_pair;
				cross_path_pair.first.insert(start_base0);
				cross_path_pair.first.insert(start_base1);
				cross_path_pair.second = fin_base_set;
				cross_path.push_back(cross_path_pair);
			}
		}
	}
	for (int i = 0; i < cross_path.size(); i++) {
		//i番目の完全2部グラフ部分の解決(基本的にはi=0のみ)
		//printf("cross path %d\n", i);
		if (false) {
			std::vector <mfile0::M_Base> edge_now(cross_path[i].first.begin(), cross_path[i].first.end());
			std::vector <mfile0::M_Base> edge_next(cross_path[i].second.begin(), cross_path[i].second.end());
			printf("edge now=%d , edge next=%d\n", edge_now.size(), edge_next.size());
			for (int j = 0; j < edge_now.size(); j++) {
				printf("(%d,%d)  (%d,%d)\n", edge_now[j].pos, edge_now[j].rawid, edge_next[j].pos, edge_next[j].rawid);
			}
		}
		//Kn,nになっていない-->まだ解決していない場所がある
		if (cross_path[i].first.size() != cross_path[i].second.size())continue;

		//マッチングの列挙
		int path_num = cross_path[i].first.size();
		//path_numが多い時は最良経路から選択する
		while (path_num > 5) {
			//未実装
			break;
		}
		int combination_num = 1;
		for (int j = 1; j <= path_num; j++) {
			combination_num *= j;
		}
		std::vector<std::vector<std::pair<Chain_path, Chain_path>>> all_matching;
		auto path0 = path[0];
		auto path1 = path[0];
		for (int j = 0; j < combination_num; j++) {
			std::vector<std::pair<Chain_path, Chain_path>> match_tmp;
			for (int k = 0; k < path_num; k++) {
				match_tmp.push_back(std::make_pair(path0, path1));
			}
			all_matching.push_back(match_tmp);
		}
		std::vector <mfile0::M_Base> edge_now(cross_path[i].first.begin(), cross_path[i].first.end());
		std::vector <mfile0::M_Base> edge_prev(cross_path[i].second.begin(), cross_path[i].second.end());
		int count_permutation = 0;
		for (int j = 0; j < edge_now.size(); j++) {
			for (int k = 0; k < path.size(); k++) {
				if (path[k].Get_path_edge().first == edge_now[j]) {
					path0 = path[k];
					break;
				}
			}
			for (int k = 0; k < all_matching.size(); k++) {
				all_matching[k][j].first = path0;
			}
		}
		sort(edge_prev.begin(), edge_prev.end());
		count_permutation = 0;
		do {
			for (int j = 0; j < edge_prev.size(); j++) {
				//printf("%d %d\n", count_permutation, j);
				auto res = path_map.find(std::make_pair(edge_prev[j], edge_now[j]));
				if (res == path_map.end()) {
					printf("error function[result_cross_path_prev]\n");
					exit(1);
				}
				auto path1 = res->second;
				all_matching[count_permutation][j].second = path1;
			}
			count_permutation++;
		} while (next_permutation(edge_prev.begin(), edge_prev.end()));

		//for (int j = 0; j < all_matching.size(); j++) {
		//	printf("combination %d\n", j);
		//	for (int k = 0; k < all_matching[j].size(); k++) {
		//		auto edge1 = all_matching[j][k].first.Get_path_edge();
		//		auto edge0 = all_matching[j][k].second.Get_path_edge();

		//		printf("%d,%d: (%d,%d)-(%d,%d) ==  (%d,%d)-(%d,%d)\n", j, k,
		//			edge0.first.pos, edge0.first.rawid, edge0.second.pos, edge0.second.rawid,
		//			edge1.first.pos, edge1.first.rawid, edge1.second.pos, edge1.second.rawid
		//		);
		//	}
		//}

		//basetrack --> 経路への変換 + 重みの計算
		std::map < std::pair<Chain_path, Chain_path>, double> path_pair;
		for (auto itr_now = cross_path[i].first.begin(); itr_now != cross_path[i].first.end(); itr_now++) {
			for (auto itr_prev = cross_path[i].second.begin(); itr_prev != cross_path[i].second.end(); itr_prev++) {
				auto res = path_map.find(std::make_pair(*itr_prev, *itr_now));
				if (res == path_map.end()) {
					printf("error function[result_cross_path_prev]\n");
					exit(1);
				}
				for (int j = 0; j < path.size(); j++) {
					if (path[j].Get_path_edge().first == *itr_now) {
						//ここで角度差の計算までする
						auto path1 = path[j];
						auto path0 = res->second;
						auto edge0 = path0.Get_path_edge();
						auto edge1 = path1.Get_path_edge();
						path0.Line_Fit(edge0.second.pos / 10 - 5, edge0.second.pos / 10);
						path1.Line_Fit(edge1.first.pos / 10, edge1.first.pos / 10 + 5);
						double angle_diff = sqrt(pow(path0.Get_line_ax() - path1.Get_line_ax(), 2) + pow(path0.Get_line_ay() - path1.Get_line_ay(), 2));
						path_pair.insert(std::make_pair(std::make_pair(path1, path0), angle_diff));
					}
				}
			}
		}

		std::vector<std::pair<Chain_path, Chain_path>> best_path_combination;
		double val_min = -1;
		for (int j = 0; j < all_matching.size(); j++) {
			//printf("combination %d\n", j);
			double val_tmp = 0;
			for (int k = 0; k < all_matching[j].size(); k++) {
				auto edge0 = all_matching[j][k].second.Get_path_edge();

				//auto path0 = all_matching[j][k].first;
				//auto path1 = all_matching[j][k].second;
				//edge0 = path0.Get_path_edge();
				//auto edge1 = path1.Get_path_edge();
				//printf("(%d,%d)-(%d-%d) == (%d,%d)-(%d-%d)\n",
				//	edge0.first.pos, edge0.first.rawid, edge0.second.pos, edge0.second.rawid,
				//	edge1.first.pos, edge1.first.rawid, edge1.second.pos, edge1.second.rawid);

				//auto edge_out = all_matching[j][k].second.Get_path_edge();
				//printf("%d %d %d %d\n", edge_out.first.pos, edge_out.first.rawid, edge_out.second.pos, edge_out.second.rawid);
				auto res = path_pair.find(all_matching[j][k]);
				if (res == path_pair.end()) {
					printf("error function[result_cross_path_prev]\n");
					exit(1);
				}
				val_tmp += res->second;
				//printf("(%4d,%10d) - (%4d,%10d) : %.4lf\n", edge0.first.pos, edge0.first.rawid, edge0.second.pos, edge0.second.rawid, res->second);
			}
			if (val_min < 0 || val_tmp < val_min) {
				val_min = val_tmp;
				best_path_combination = all_matching[j];
			}
		}

		for (int j = 0; j < best_path_combination.size(); j++) {
			for (auto& b : best_path_combination[j].second.Get_all_base()) {
				best_path_combination[j].first.Add_Basetrack(b);
			}
			ret.push_back(best_path_combination[j].first);
		}
	}

	//最終結果
	//for (int i = 0; i < path.size(); i++) {
	//	auto edge = path[i].Get_path_edge();
	//	printf("(%d,%d) - (%d,%d)\n", edge.first.pos, edge.first.rawid, edge.second.pos, edge.second.rawid);
	//}
	//printf("\n");
	//for (int i = 0; i < ret.size(); i++) {
	//	auto edge = ret[i].Get_path_edge();
	//	printf("(%d,%d) - (%d,%d)\n", edge.first.pos, edge.first.rawid, edge.second.pos, edge.second.rawid);
	//}
	//printf("\n");
	for (auto itr = path.begin(); itr != path.end(); itr++) {
		auto edge0 = itr->Get_path_edge();
		bool flg = true;
		for (auto itr1 = ret.begin(); itr1 != ret.end(); itr1++) {
			auto edge1 = itr1->Get_path_edge();
			if (edge0.second == edge1.second)flg = false;
		}
		if (flg)ret.push_back(*itr);
	}


	path = ret;
	return true;
}

//本体
std::vector<std::vector<std::pair<mfile0::M_Base, mfile0::M_Base>>> select_path_v2(std::vector<std::pair<mfile0::M_Base, mfile0::M_Base>>& path) {

	boost::unordered_multimap<mfile0::M_Base, mfile0::M_Base> path_next, path_prev;
	std::set<mfile0::M_Base> all_vertex;

	for (auto itr = path.begin(); itr != path.end(); itr++) {
		path_next.insert(std::make_pair(itr->first, itr->second));
		path_prev.insert(std::make_pair(itr->second, itr->first));
		all_vertex.insert(itr->first);
		all_vertex.insert(itr->second);
	}
	//一本路の圧縮
	//true --> PLが大きいほうへ進む
	bool dirction_flg;
	//std::set<mfile0::M_Base> start_edge= select_start(all_vertex, path_next, path_prev, dirction_flg);
	std::vector<Chain_path> start_path = select_start(all_vertex, path_next, path_prev, dirction_flg);
	//printf("start path\n");
	//for (int i = 0; i < start_path.size(); i++) {
	//	auto edge = start_path[i].Get_path_edge();
	//	printf("path %d :%d,%d -- %d,%d\n", i, edge.first.pos, edge.first.rawid, edge.second.pos, edge.second.rawid);
	//}

	if (dirction_flg) {
		//printf("next\n");
		std::map<std::pair<mfile0::M_Base, mfile0::M_Base>, Chain_path>  path_map = path_compress_next(all_vertex, path_next, path_prev);
		//for (auto itr = path_map.begin(); itr != path_map.end(); itr++) {
		//	printf("(%d,%d)-(%d-%d)\n", itr->first.first.pos, itr->first.first.rawid, itr->first.second.pos, itr->first.second.rawid);
		//}
		while (result_cross_path_next(start_path, path_map, path_next, path_prev));
	}
	else {
		//printf("prev\n");
		std::map<std::pair<mfile0::M_Base, mfile0::M_Base>, Chain_path>  path_map = path_compress_prev(all_vertex, path_next, path_prev);
		//for (auto itr = path_map.begin(); itr != path_map.end(); itr++) {
		//	printf("(%d,%d)-(%d-%d)\n", itr->first.first.pos, itr->first.first.rawid, itr->first.second.pos, itr->first.second.rawid);
		//}
		while (result_cross_path_prev(start_path, path_map, path_next, path_prev));
	}


	std::vector<std::vector<std::pair<mfile0::M_Base, mfile0::M_Base>>> ret;
	for (int i = 0; i < start_path.size(); i++) {
		std::vector<std::pair<mfile0::M_Base, mfile0::M_Base>> one_path;
		auto save_base = start_path[i].Get_all_base();
		for (int j = 0; j < save_base.size(); j++) {
			if (j + 1 == save_base.size())continue;
			one_path.push_back(std::make_pair(save_base[j], save_base[j + 1]));
		}
		ret.push_back(one_path);
	}
	return ret;
}

//pathを2部グラフに分解
std::vector < std::pair<std::set<mfile0::M_Base>, std::set<mfile0::M_Base>>> Divide_bipartite_graph(output_format& g, std::vector<std::pair<mfile0::M_Base, mfile0::M_Base>>& path) {
	std::vector < std::pair<std::set<mfile0::M_Base>, std::set<mfile0::M_Base>>>bipartite_graph_v;
	boost::unordered_multimap <int, int> path_next;
	boost::unordered_multimap <int, int> path_prev;
	std::map<int, mfile0::M_Base>index_base;
	std::map< mfile0::M_Base, int>base_index;
	int count = 0;
	for (auto itr = path.begin(); itr != path.end(); itr++) {
		auto base0 = itr->first;
		auto base1 = itr->second;
		auto res0 = base_index.find(base0);
		auto res1 = base_index.find(base1);
		if (res0 == base_index.end()) {
			index_base.insert(std::make_pair(count, base0));
			base_index.insert(std::make_pair(base0, count));
			res0 = base_index.find(base0);
			count++;
		}
		if (res1 == base_index.end()) {
			index_base.insert(std::make_pair(count, base1));
			base_index.insert(std::make_pair(base1, count));
			res1 = base_index.find(base1);
			count++;
		}
		path_next.insert(std::make_pair(res0->second, res1->second));
		path_prev.insert(std::make_pair(res1->second, res0->second));
	}
	std::set<int> finished;
	for (auto itr = index_base.begin(); itr != index_base.end(); itr++) {
		int target = itr->first;
		if (finished.count(target) == 1)continue;
		std::set<int>up, down;
		//2部グラフだった場合
		if (Judge_bipartite_graph(target, path_prev, path_next, up, down)) {
			//下端の場合
			if (down.size() == 0) {
				finished.insert(target);
				continue;
			}
			std::pair<std::set<mfile0::M_Base>, std::set<mfile0::M_Base>> bipartite_graph;
			for (auto itr = up.begin(); itr != up.end(); itr++) {
				auto res = index_base.find(*itr);
				if (res == index_base.end()) {
					fprintf(stderr, "error function [Divide_bipartite_graph]\n");
					exit(1);
				}
				bipartite_graph.first.insert(res->second);
				finished.insert(*itr);
			}
			for (auto itr = down.begin(); itr != down.end(); itr++) {
				auto res = index_base.find(*itr);
				if (res == index_base.end()) {
					fprintf(stderr, "error function [Divide_bipartite_graph]\n");
					exit(1);
				}
				bipartite_graph.second.insert(res->second);
				//下流側は入れない
				//finished.insert(*itr);
			}
			bipartite_graph_v.push_back(bipartite_graph);
		}
		else {
			//2部グラフでない場合
			auto res = index_base.find(target);
			if (res == index_base.end()) {
				fprintf(stderr, "error function [Divide_bipartite_graph]\n");
				exit(1);
			}

			printf("%d %d no bipartite graph\n", res->second.pos, res->second.rawid);
			finished.insert(target);
		}


	}


	return bipartite_graph_v;
}


//角度差等使わなくていいもののみ
//matchingの列挙
void result_bipartite_graph_0(
	std::vector < std::pair<std::set<mfile0::M_Base>, std::set<mfile0::M_Base>>>& bipartite_graph_v,
	std::list<std::pair<int, std::tuple<int, int, int, int>>>& single_path,
	std::list<std::pair<int, std::tuple<int, int, int, int>>>& branch_path,
	std::map<std::pair<int, int>, mfile0::M_Base>& base_map_id, int& path_id) {

	//bipartite_graph_v　2部グラフの頂点
	//single_path&branch_path 2部グラフを作るpath

	//pathの記憶
	boost::unordered_multimap <mfile0::M_Base, mfile0::M_Base> path_next;
	boost::unordered_multimap <mfile0::M_Base, mfile0::M_Base> path_prev;
	for (auto itr = single_path.begin(); itr != single_path.end(); itr++) {
		int pos[2], rawid[2];
		pos[0] = std::get<0>(itr->second);
		pos[1] = std::get<1>(itr->second);
		rawid[0] = std::get<2>(itr->second);
		rawid[1] = std::get<3>(itr->second);
		auto base0 = base_map_id.find(std::make_pair(pos[0] / 10, rawid[0]));
		auto base1 = base_map_id.find(std::make_pair(pos[1] / 10, rawid[1]));
		if (base0 == base_map_id.end()) {
			fprintf(stderr, "error PL%03d %d not found 0\n", pos[0] / 10, rawid[0]);
			fprintf(stderr, "function [result_bipartite_graph_0]\n");
			exit(1);
		}
		if (base1 == base_map_id.end()) {
			fprintf(stderr, "error PL%03d %d not found 1\n", pos[1] / 10, rawid[1]);
			fprintf(stderr, "function [result_bipartite_graph_0]\n");
			exit(1);
		}
		path_next.insert(std::make_pair(base0->second, base1->second));
		path_prev.insert(std::make_pair(base1->second, base0->second));
	}

	for (auto& bg : bipartite_graph_v) {
		std::map<int, mfile0::M_Base> index_to_vertex;
		std::map< mfile0::M_Base, int> vertex_to_index;
		std::set<int>up, down;
		int count = 0;
		//頂点にiDを振り分ける
		for (auto itr = bg.first.begin(); itr != bg.first.end(); itr++) {
			mfile0::M_Base base0 = *itr;
			auto res = vertex_to_index.insert(std::make_pair(*itr, count));
			if (!res.second)continue;
			else {
				index_to_vertex.insert(std::make_pair(count, *itr));
				up.insert(count);
				count++;
			}
		}
		for (auto itr = bg.second.begin(); itr != bg.second.end(); itr++) {
			mfile0::M_Base base0 = *itr;
			auto res = vertex_to_index.insert(std::make_pair(*itr, count));
			if (!res.second)continue;
			else {
				index_to_vertex.insert(std::make_pair(count, *itr));
				down.insert(count);
				count++;
			}
		}

		std::vector<std::pair<int, int>> all_path;
		//pathを記憶する
		for (auto itr = branch_path.begin(); itr != branch_path.end(); itr++) {
			int pos[2], rawid[2];
			pos[0] = std::get<0>(itr->second);
			pos[1] = std::get<1>(itr->second);
			rawid[0] = std::get<2>(itr->second);
			rawid[1] = std::get<3>(itr->second);
			auto base0 = base_map_id.find(std::make_pair(pos[0] / 10, rawid[0]));
			auto base1 = base_map_id.find(std::make_pair(pos[1] / 10, rawid[1]));
			if (base0 == base_map_id.end()) {
				fprintf(stderr, "error PL%03d %d not found\n", pos[0] / 10, rawid[0]);
				fprintf(stderr, "function [result_bipartite_graph_0]\n");
				exit(1);
			}
			if (base1 == base_map_id.end()) {
				fprintf(stderr, "error PL%03d %d not found\n", pos[1] / 10, rawid[1]);
				fprintf(stderr, "function [result_bipartite_graph_0]\n");
				exit(1);
			}

			auto index0 = vertex_to_index.find(base0->second);
			auto index1 = vertex_to_index.find(base1->second);
			if (index0 == vertex_to_index.end() || index1 == vertex_to_index.end())continue;
			if (up.count(index0->second) == 0)continue;
			if (down.count(index1->second) == 0)continue;

			all_path.push_back(std::make_pair(index0->second, index1->second));
		}
		if (false) {
			printf("index vertex\n");
			for (auto itr = index_to_vertex.begin(); itr != index_to_vertex.end(); itr++) {
				printf("%d:(%d,%d)\n", itr->first, itr->second.pos, itr->second.rawid);
			}
			printf("all path\n");
			for (auto itr = all_path.begin(); itr != all_path.end(); itr++) {
				printf("%d-%d\n", itr->first, itr->second);
			}
		}


		//とれる組み合わせの列挙
		std::vector<std::vector<std::pair<int, int>>> enumerat_all_path = Enumeration(all_path);
		//basetrackとchainの対応
		std::map<mfile0::M_Base, Chain_path> chain_map;
		make_to_next_chain(chain_map, bg.second, path_next);
		make_to_prev_chain(chain_map, bg.first, path_prev);
		//列挙した組み合わせ-->basetrackに変換
		std::vector<std::vector<std::pair<mfile0::M_Base, mfile0::M_Base>>>enumerat_all_path_base;
		for (int i = 0; i < enumerat_all_path.size(); i++) {
			std::vector<std::pair<mfile0::M_Base, mfile0::M_Base>> path_combination;
			for (int j = 0; j < enumerat_all_path[i].size(); j++) {
				auto base0 = index_to_vertex.find(enumerat_all_path[i][j].first);
				auto base1 = index_to_vertex.find(enumerat_all_path[i][j].second);
				if (base0 == index_to_vertex.end()) {
					fprintf(stderr, "function [result_bipartite_graph_0]\n");
					exit(1);
				}
				if (base1 == index_to_vertex.end()) {
					fprintf(stderr, "function [result_bipartite_graph_0]\n");
					exit(1);
				}
				path_combination.push_back(std::make_pair(base0->second, base1->second));
			}
			enumerat_all_path_base.push_back(path_combination);
		}

		//select best combination
		std::vector<std::pair<mfile0::M_Base, mfile0::M_Base>> best_path = Select_best_combination(enumerat_all_path_base, chain_map);

		path_organize(best_path, path_next, path_prev, single_path, branch_path, path_id);

	}
	return;
}
void path_organize(std::vector<std::pair<mfile0::M_Base, mfile0::M_Base>>& best_path,
	boost::unordered_multimap <mfile0::M_Base, mfile0::M_Base>& path_next,
	boost::unordered_multimap <mfile0::M_Base, mfile0::M_Base>& path_prev,
	std::list<std::pair<int, std::tuple<int, int, int, int>>>& single_path,
	std::list<std::pair<int, std::tuple<int, int, int, int>>>& branch_path,
	int& path_id
) {
	//printf("best path\n");
	//for (int i = 0; i < best_path.size(); i++) {

	//	Print_path(best_path[i].first, best_path[i].second);
	//	printf("\n");
	//}


	//追加
	//signle path
	//path_next
	//path_prev
	std::multimap<mfile0::M_Base, mfile0::M_Base> best_path_next;
	std::multimap<mfile0::M_Base, mfile0::M_Base> best_path_prev;
	for (auto itr = best_path.begin(); itr != best_path.end(); itr++) {
		best_path_next.insert(std::make_pair(itr->first, itr->second));
		best_path_prev.insert(std::make_pair(itr->second, itr->first));
	}
	std::set<std::pair<mfile0::M_Base, mfile0::M_Base>> new_single_path;
	for (auto itr = best_path_next.begin(); itr != best_path_next.end(); itr++) {
		if (best_path_next.count(itr->first) != 1)continue;
		if (best_path_prev.count(itr->second) != 1)continue;
		if (new_single_path.count(*itr) != 0)continue;
		new_single_path.insert(*itr);
		path_next.insert(std::make_pair(itr->first, itr->second));
		path_prev.insert(std::make_pair(itr->second, itr->first));
		single_path.push_back(std::make_pair(path_id, std::make_tuple(
			itr->first.pos, itr->second.pos, itr->first.rawid, itr->second.rawid)));
		path_id++;
	}


	//削除
	//branch path
	std::multimap<std::pair<int, int>, std::pair<int, int>> best_path_next_id;
	std::multimap<std::pair<int, int>, std::pair<int, int>> best_path_prev_id;
	for (auto itr = best_path.begin(); itr != best_path.end(); itr++) {
		best_path_next_id.insert(std::make_pair(std::make_pair(itr->first.pos, itr->first.rawid), std::make_pair(itr->second.pos, itr->second.rawid)));
		best_path_prev_id.insert(std::make_pair(std::make_pair(itr->second.pos, itr->second.rawid), std::make_pair(itr->first.pos, itr->first.rawid)));
	}

	for (auto itr = branch_path.begin(); itr != branch_path.end();) {
		std::pair<int, int>b0, b1;
		b0.first = std::get<0>(itr->second);
		b0.second = std::get<2>(itr->second);
		b1.first = std::get<1>(itr->second);
		b1.second = std::get<3>(itr->second);

		if (best_path_next_id.count(b0) == 0) {
			itr++;
			continue;
		}
		//signle pathへ行く
		if (best_path_next_id.count(b0) == 1 && best_path_prev_id.count(b1) == 1) {
			itr = branch_path.erase(itr);
			continue;
		}
		auto range = best_path_next_id.equal_range(b0);
		bool flg = false;
		for (auto res = range.first; res != range.second; res++) {
			if (res->first == b0 && res->second == b1)flg = true;
		}
		if (!flg) {
			itr = branch_path.erase(itr);
		}
		else {
			itr++;
		}

	}
	return;


}
std::vector<std::pair<mfile0::M_Base, mfile0::M_Base>> Select_best_combination(std::vector<std::vector<std::pair<mfile0::M_Base, mfile0::M_Base>>>& enumerat_all_path_base, std::map<mfile0::M_Base, Chain_path> chain_map) {
	double best_val = -1;
	double val = 0;
	std::vector<std::pair<mfile0::M_Base, mfile0::M_Base>> best_path;
	for (int i = 0; i < enumerat_all_path_base.size(); i++) {
		std::vector<std::pair<Chain_path, Chain_path>>chains;
		for (int j = 0; j < enumerat_all_path_base[i].size(); j++) {
			if (chain_map.count(enumerat_all_path_base[i][j].first) == 0 ||
				chain_map.count(enumerat_all_path_base[i][j].second) == 0) {
				fprintf(stderr, "error function[Select_best_combination]\n");
				exit(1);
			}
			auto c0 = chain_map.find(enumerat_all_path_base[i][j].first)->second;
			auto c1 = chain_map.find(enumerat_all_path_base[i][j].second)->second;
			c1.Add_Basetrack(enumerat_all_path_base[i][j].first);
			chains.push_back(std::make_pair(c0, c1));
		}
		//printf("Path %d\n", i);
		val = Calc_combination_value(chains);
		if (best_val < 0 || val < best_val) {
			best_val = val;
			best_path = enumerat_all_path_base[i];
		}
	}
	return best_path;
}
double Calc_combination_value(std::vector<std::pair<Chain_path, Chain_path>>& chains) {
	double val = 0;
	for (int i = 0; i < chains.size(); i++) {
		int pl0 = chains[i].first.Get_path_edge().second.pos / 10;
		int pl1 = chains[i].second.Get_path_edge().first.pos / 10;
		//Print_path(chains[i].first.Get_path_edge().first, chains[i].first.Get_path_edge().second);
		//printf("--");
		//Print_path(chains[i].second.Get_path_edge().first, chains[i].second.Get_path_edge().second);
		//printf("\n");
		chains[i].first.Line_Fit(pl0 - 5, pl0);
		chains[i].second.Line_Fit(pl1, pl1 + 5);
		val += pow(chains[i].first.Get_line_ax() - chains[i].second.Get_line_ax(), 2);
		val += pow(chains[i].first.Get_line_ay() - chains[i].second.Get_line_ay(), 2);
	}
	return val;

}
//分岐-->下流へつながる1本のchain
void make_to_next_chain(std::map<mfile0::M_Base, Chain_path>& chain_map, std::set<mfile0::M_Base>& base, boost::unordered_multimap <mfile0::M_Base, mfile0::M_Base>& path_next) {
	for (auto itr = base.begin(); itr != base.end(); itr++) {
		mfile0::M_Base start = *itr;
		Chain_path c;
		c.Add_Basetrack(start);
		while (true) {
			if (path_next.count(start) != 1)break;
			auto res = path_next.find(start);
			start = res->second;
			c.Add_Basetrack(start);
		}
		chain_map.insert(std::make_pair(*itr, c));
	}
}
//上流-->分岐へつながる1本のchain
void make_to_prev_chain(std::map<mfile0::M_Base, Chain_path>& chain_map, std::set<mfile0::M_Base>& base, boost::unordered_multimap <mfile0::M_Base, mfile0::M_Base>& path_prev) {
	for (auto itr = base.begin(); itr != base.end(); itr++) {
		mfile0::M_Base start = *itr;
		Chain_path c;
		c.Add_Basetrack(start);
		while (true) {
			if (path_prev.count(start) != 1)break;
			auto res = path_prev.find(start);
			start = res->second;
			c.Add_Basetrack(start);
		}
		chain_map.insert(std::make_pair(*itr, c));
	}
}

//ささくれ閉路の解決
void result_hangnail_cycle(output_format& g, std::map<std::pair<int, int>, mfile0::M_Base>& base, int& path_id) {
	//idを使わない
	auto path_list = g.cut_path;
	std::multimap<mfile0::M_Base, mfile0::M_Base>path_next;
	std::multimap<mfile0::M_Base, mfile0::M_Base>path_prev;
	std::set<mfile0::M_Base>all_edge;
	for (auto itr = path_list.begin(); itr != path_list.end(); itr++) {
		std::pair<int, int> id0, id1;
		id0 = std::make_pair(std::get<0>(itr->second) / 10, std::get<2>(itr->second));
		id1 = std::make_pair(std::get<1>(itr->second) / 10, std::get<3>(itr->second));
		auto base0 = base.find(id0);
		auto base1 = base.find(id1);
		if (base0 == base.end() || base1 == base.end()) {
			printf("basetrack not found\n");
		}
		path_next.insert(std::make_pair(base0->second, base1->second));
		path_prev.insert(std::make_pair(base1->second, base0->second));
		all_edge.insert(base0->second);
		all_edge.insert(base1->second);
	}
	std::vector<std::vector<std::pair<mfile0::M_Base, mfile0::M_Base>>> ret;
	std::set<mfile0::M_Base>finished_edge;
	for (auto itr = all_edge.begin(); itr != all_edge.end(); itr++) {
		if (finished_edge.count(*itr) == 1) continue;
		if (path_next.count(*itr) < 2) {
			finished_edge.insert(*itr);
			continue;
		}
		bool flg = true;
		auto range = path_next.equal_range(*itr);
		std::set<mfile0::M_Base> center;
		std::set<mfile0::M_Base> end_point;
		for (auto res = range.first; res != range.second; res++) {
			if (path_next.count(res->second) != 1) {
				flg = false;
				break;
			}
			if (path_prev.count(res->second) != 1) {
				flg = false;
				break;
			}
			center.insert(res->second);
			auto res2 = path_next.find(res->second);
			end_point.insert(res2->second);
			if (end_point.size() > 1) {
				flg = false;
				break;
			}
		}
		//if (g.groupid == 1538) {
		//	printf("start:%d,%d\n", itr->pos, itr->rawid);
		//	for (auto itr2 = center.begin(); itr2 != center.end(); itr2++) {
		//		printf("%d,%d\n", itr2->pos, itr2->rawid);
		//	}
		//	for (auto itr2 = end_point.begin(); itr2 != end_point.end(); itr2++) {
		//		printf("end:%d,%d\n", itr2->pos, itr2->rawid);
		//	}
		//	printf("flg=%d\n", flg);
		//}
		//1-n-1閉路
		if (flg) {
			double dist = -1;
			mfile0::M_Base select_base;
			auto end_p = *(end_point.begin());
			for (auto itr2 = center.begin(); itr2 != center.end(); itr2++) {
				double dx = itr->x - (itr->x - end_p.x) / (itr->z - end_p.z) * (itr->z - itr2->z);
				double dy = itr->y - (itr->y - end_p.y) / (itr->z - end_p.z) * (itr->z - itr2->z);
				dx = dx - itr2->x;
				dy = dy - itr2->y;
				if (dist<0 || dist>dx * dx + dy * dy) {
					dist = dx * dx + dy * dy;
					select_base = *itr2;
				}
			}
			for (auto itr2 = center.begin(); itr2 != center.end(); itr2++) {
				auto path0 = std::make_tuple(itr->pos, itr2->pos, (int)itr->rawid, (int)itr2->rawid);
				auto path1 = std::make_tuple(itr2->pos, end_point.begin()->pos, (int)itr2->rawid, (int)end_point.begin()->rawid);

				if (*itr2 == select_base) {
					g.comfirmed_path.push_back(std::make_pair(path_id, path0));
					g.comfirmed_path.push_back(std::make_pair(path_id, path1));
					path_id++;
				}

				for (auto itr3 = g.cut_path.begin(); itr3 != g.cut_path.end(); ) {
					if (itr3->second == path0) {
						itr3 = g.cut_path.erase(itr3);
					}
					else {
						itr3++;
					}
				}
				for (auto itr3 = g.cut_path.begin(); itr3 != g.cut_path.end(); ) {
					if (itr3->second == path1) {
						itr3 = g.cut_path.erase(itr3);
					}
					else {
						itr3++;
					}
				}
			}
			for (auto itr2 = center.begin(); itr2 != center.end(); itr2++) {
				finished_edge.insert(*itr2);
			}
		}
		finished_edge.insert(*itr);
	}

	g.num_comfirmed_path = g.comfirmed_path.size();
	g.num_cut_path = g.cut_path.size();
	g.num_select_path = g.select_path.size();

}
//ささくれ閉路の解決
//--o-o-o-o-o--
//  |---o---|
bool check_simple_cycle(std::vector<std::vector<mfile0::M_Base>>& all_hist, std::multimap<mfile0::M_Base, mfile0::M_Base>& all_path) {
	std::set<mfile0::M_Base>next_point;
	mfile0::M_Base start_point = all_hist[0][0];
	for (int i = 0; i < all_hist.size(); i++) {
		if (all_hist[i].size() < 2)continue;
		next_point.insert(all_hist[i][1]);
	}
	//解決済み
	if (all_path.count(start_point) == 1)return false;
	auto range = all_path.equal_range(start_point);
	for (auto res = range.first; res != range.second; res++) {
		if (next_point.count(res->second) == 0)return false;
	}
	return true;
}
void result_cycle_cross(output_format& g, std::map<std::pair<int, int>, mfile0::M_Base>& base, int& path_id) {
	//idを使わない
	auto path_list = g.cut_path;
	std::multimap<mfile0::M_Base, mfile0::M_Base>path_next;
	std::multimap<mfile0::M_Base, mfile0::M_Base>path_prev;
	std::set<mfile0::M_Base>all_edge;
	for (auto itr = g.cut_path.begin(); itr != g.cut_path.end(); itr++) {
		std::pair<int, int> id0, id1;
		id0 = std::make_pair(std::get<0>(itr->second) / 10, std::get<2>(itr->second));
		id1 = std::make_pair(std::get<1>(itr->second) / 10, std::get<3>(itr->second));
		auto base0 = base.find(id0);
		auto base1 = base.find(id1);
		if (base0 == base.end() || base1 == base.end()) {
			printf("basetrack not found\n");
		}
		path_next.insert(std::make_pair(base0->second, base1->second));
		path_prev.insert(std::make_pair(base1->second, base0->second));
		all_edge.insert(base0->second);
		all_edge.insert(base1->second);
	}
	for (auto itr = g.comfirmed_path.begin(); itr != g.comfirmed_path.end(); itr++) {
		std::pair<int, int> id0, id1;
		id0 = std::make_pair(std::get<0>(itr->second) / 10, std::get<2>(itr->second));
		id1 = std::make_pair(std::get<1>(itr->second) / 10, std::get<3>(itr->second));
		auto base0 = base.find(id0);
		auto base1 = base.find(id1);
		if (base0 == base.end() || base1 == base.end()) {
			printf("basetrack not found\n");
		}
		path_next.insert(std::make_pair(base0->second, base1->second));
		path_prev.insert(std::make_pair(base1->second, base0->second));
		all_edge.insert(base0->second);
		all_edge.insert(base1->second);
	}

	std::set<mfile0::M_Base>start_edge;

	std::set<mfile0::M_Base>end_edge;

	//閉路の解決
	for (auto itr = all_edge.begin(); itr != all_edge.end(); itr++) {
		if (path_next.count(*itr) > 1) {
			start_edge.insert(*itr);
		}
		if (path_prev.count(*itr) > 1) {
			end_edge.insert(*itr);
		}
	}
	std::vector<mfile0::M_Base>path;
	std::set<mfile0::M_Base>seen;
	for (auto itr = start_edge.begin(); itr != start_edge.end(); itr++) {
		//if (g.groupid == 2931) {
		//	printf("start %d,%d\n", itr->pos, itr->rawid);
		//}
		mfile0::M_Base now = *itr;
		for (auto itr2 = end_edge.begin(); itr2 != end_edge.end(); itr2++) {
			mfile0::M_Base goal = *itr2;
			//if (g.groupid == 2931) {
			//	printf("goal %d,%d\n", itr2->pos, itr2->rawid);
			//}
			std::vector<std::vector<mfile0::M_Base>>all_hist;
			path.clear();
			seen.clear();
			enumerate_path_dfs(path, all_hist, now, goal, seen, path_next);
			//if (g.groupid == 2931) {
			//	printf("\n");
			//	printf("path:(%d,%d)-(%d,%d)\n", itr->pos, itr->rawid, itr2->pos, itr2->rawid);
			//	for (int i = 0; i < all_hist.size(); i++) {
			//		for (int j = 0; j < all_hist[i].size(); j++) {
			//			if (j + 1 != all_hist[i].size()) {
			//				printf("(%d,%d)-->", all_hist[i][j].pos, all_hist[i][j].rawid);
			//			}
			//			else {
			//				printf("(%d,%d)\n", all_hist[i][j].pos, all_hist[i][j].rawid);
			//			}
			//		}
			//	}
			//}
			if (all_hist.size() < 2)continue;
			//startから1本線が出ている(2本線ではない)-->閉路がつぶされている-->使わない
			if (!check_simple_cycle(all_hist, path_next))continue;
			for (int i = 0; i < all_hist.size(); i++) {
				if (all_hist[i].size() == 3) {
					//このpathは消す
					std::pair<mfile0::M_Base, mfile0::M_Base>del_path0 = std::make_pair(all_hist[i][0], all_hist[i][1]);
					std::pair<mfile0::M_Base, mfile0::M_Base>del_path1 = std::make_pair(all_hist[i][1], all_hist[i][2]);
					//std::multimap<mfile0::M_Base, mfile0::M_Base>path_next;
					//std::multimap<mfile0::M_Base, mfile0::M_Base>path_prev;
					//if (g.groupid == 2931) {
					//	printf("del_path0:(%d,%d)-(%d,%d)\n", del_path0.first.pos, del_path0.second.pos, del_path0.first.rawid, del_path0.second.rawid);
					//	printf("del_path1:(%d,%d)-(%d,%d)\n", del_path1.first.pos, del_path1.second.pos, del_path1.first.rawid, del_path1.second.rawid);
					//}
					for (auto itr3 = path_next.begin(); itr3 != path_next.end(); ) {
						//if (g.groupid == 2931) {
						//	printf("path next:(%d,%d)-(%d,%d) flg0 %d,flg1 %d\n", itr->first.pos, itr->second.pos, itr->first.rawid, itr->second.rawid, del_path0.first == itr->first&&del_path0.second == itr->second, del_path1.first == itr->first&&del_path1.second == itr->second);
						//}
						if ((del_path0.first == itr3->first && del_path0.second == itr3->second) ||
							(del_path1.first == itr3->first && del_path1.second == itr3->second)) {
							itr3 = path_next.erase(itr3);
						}
						else {
							itr3++;
						}
					}
					for (auto itr3 = path_prev.begin(); itr3 != path_prev.end(); ) {
						//if (g.groupid == 2931) {
						//	printf("path prev:(%d,%d)-(%d,%d) flg0 %d,flg1 %d\n", itr->first.pos, itr->second.pos, itr->first.rawid, itr->second.rawid, del_path0.second == itr->first&&del_path0.first == itr->second, del_path1.second == itr->first&&del_path1.first == itr->second);
						//}

						if ((del_path0.second == itr3->first && del_path0.first == itr3->second) ||
							(del_path1.second == itr3->first && del_path1.first == itr3->second)) {
							itr3 = path_prev.erase(itr3);
						}
						else {
							itr3++;
						}
					}
					std::tuple<int, int, int, int> del_path_tuple0 = std::make_tuple(del_path0.first.pos, del_path0.second.pos, (int)del_path0.first.rawid, (int)del_path0.second.rawid);
					std::tuple<int, int, int, int> del_path_tuple1 = std::make_tuple(del_path1.first.pos, del_path1.second.pos, (int)del_path1.first.rawid, (int)del_path1.second.rawid);
					//if (g.groupid == 2931) {
					//	printf("del_path0:(%d,%d)-(%d,%d)\n", std::get<0>(del_path_tuple0), std::get<1>(del_path_tuple0), std::get<2>(del_path_tuple0), std::get<3>(del_path_tuple0));
					//	printf("del_path1:(%d,%d)-(%d,%d)\n", std::get<0>(del_path_tuple1), std::get<1>(del_path_tuple1), std::get<2>(del_path_tuple1), std::get<3>(del_path_tuple1));
					//}
					for (auto itr3 = g.cut_path.begin(); itr3 != g.cut_path.end();) {
						//if (g.groupid == 2931) {
						//	printf("path:(%d,%d)-(%d,%d) flg %d\n", std::get<0>(itr->second), std::get<1>(itr->second), std::get<2>(itr->second), std::get<3>(itr->second), itr->second == del_path_tuple0 || itr->second == del_path_tuple1);
						//}
						if (itr3->second == del_path_tuple0 || itr3->second == del_path_tuple1) {
							itr3 = g.cut_path.erase(itr3);
						}
						else {
							itr3++;
						}
					}
				}
			}
		}
	}
}
void enumerate_path_dfs(std::vector<mfile0::M_Base>& path, std::vector<std::vector<mfile0::M_Base>>& all_hist, mfile0::M_Base now, mfile0::M_Base& goal, std::set<mfile0::M_Base>& seen, std::multimap<mfile0::M_Base, mfile0::M_Base>& all_path) {
	path.push_back(now);
	if (now == goal) {
		//printf("hit\n");
		//Print_hist(path);
		all_hist.push_back(path);
		path.pop_back();
		return;
	}
	if (all_path.count(now) == 0) {
		seen.insert(now);
		path.pop_back();
		return;
	}
	auto range = all_path.equal_range(now);
	for (auto res = range.first; res != range.second; res++) {
		if (seen.count(res->second) == 1)continue;
		if (goal.pos < res->second.pos)continue;
		enumerate_path_dfs(path, all_hist, res->second, goal, seen, all_path);
	}
	seen.insert(now);
	path.pop_back();

}
bool Judge_bipartite_graph(int target, boost::unordered_multimap <int, int>& path_prev, boost::unordered_multimap <int, int>& path_next, std::set<int>& up, std::set<int>& down) {
	up.insert(target);
	int sum = 0, sum_p = -1;
	while (sum != sum_p) {
		sum_p = sum;
		for (auto itr = up.begin(); itr != up.end(); itr++) {
			if (path_next.count(*itr) == 0)continue;
			auto range = path_next.equal_range(*itr);
			for (auto res = range.first; res != range.second; res++) {
				down.insert(res->second);
			}
		}
		for (auto itr = down.begin(); itr != down.end(); itr++) {
			if (path_prev.count(*itr) == 0)continue;
			auto range = path_prev.equal_range(*itr);
			for (auto res = range.first; res != range.second; res++) {
				up.insert(res->second);
			}
		}
		sum = up.size() + down.size();
	}

	//2部グラフの判定
	bool flg = true;
	for (auto itr = up.begin(); itr != up.end(); itr++) {
		if (down.count(*itr) == 1)flg = false;
	}
	for (auto itr = down.begin(); itr != down.end(); itr++) {
		if (up.count(*itr) == 1)flg = false;
	}

	return flg;
}

void result_no_bipartite_graph(output_format& group, std::vector<std::pair<mfile0::M_Base, mfile0::M_Base>>& path, int& path_id) {
	boost::unordered_multimap <int, int> path_next;
	boost::unordered_multimap <int, int> path_prev;
	std::map<int, mfile0::M_Base>index_base;
	std::map< mfile0::M_Base, int>base_index;
	int count = 0;
	for (auto itr = path.begin(); itr != path.end(); itr++) {
		auto base0 = itr->first;
		auto base1 = itr->second;
		auto res0 = base_index.find(base0);
		auto res1 = base_index.find(base1);
		if (res0 == base_index.end()) {
			index_base.insert(std::make_pair(count, base0));
			base_index.insert(std::make_pair(base0, count));
			res0 = base_index.find(base0);
			count++;
		}
		if (res1 == base_index.end()) {
			index_base.insert(std::make_pair(count, base1));
			base_index.insert(std::make_pair(base1, count));
			res1 = base_index.find(base1);
			count++;
		}
		path_next.insert(std::make_pair(res0->second, res1->second));
		path_prev.insert(std::make_pair(res1->second, res0->second));
	}
	std::set<int> finished;
	for (auto itr = index_base.begin(); itr != index_base.end(); ) {
		int target = itr->first;
		if (finished.count(target) == 1) {
			itr++;
			continue;
		}
		std::set<int>up, down;
		//2部グラフだった場合
		if (Judge_bipartite_graph(target, path_prev, path_next, up, down)) {
			for (auto itr2 = up.begin(); itr2 != up.end(); itr2++) {
				finished.insert(*itr2);
			}
			for (auto itr2 = down.begin(); itr2 != down.end(); itr2++) {
				finished.insert(*itr2);
			}
			itr++;
			continue;
		}
		//2部グラフでない場合
		std::set<int>all_vertex;
		for (auto itr2 = up.begin(); itr2 != up.end(); itr2++) {
			all_vertex.insert(*itr2);
		}
		for (auto itr2 = down.begin(); itr2 != down.end(); itr2++) {
			all_vertex.insert(*itr2);
		}
		//深さ0-1でつながっている部分を抽出
		//深さ0-2でつながっている部分を抽出-->cut
		//深さ1-2でつなげる
		std::map<int, int> all_vertex_depth = add_depth_information(all_vertex, up, down, path_next);
		for (auto itr2 = all_vertex_depth.begin(); itr2 != all_vertex_depth.end(); itr2++) {
			if (itr2->second != 0) 				continue;

			boost::unordered_multimap <int, int> path_01, path_02;
			auto range = path_next.equal_range(itr2->first);
			for (auto res = range.first; res != range.second; res++) {
				auto res_depth = all_vertex_depth.find(res->second);
				if (res_depth == all_vertex_depth.end())continue;
				if (res_depth->second == 1)path_01.insert(std::make_pair(itr2->first, res->second));
				if (res_depth->second == 2)path_02.insert(std::make_pair(itr2->first, res->second));
			}
			if (path_01.size() == 0 || path_02.size() == 0)continue;
			for (auto itr_02 = path_02.begin(); itr_02 != path_02.end(); itr_02++) {
				auto res0 = index_base.find(itr_02->first);
				auto res2 = index_base.find(itr_02->second);
				if (res0 == index_base.end()) {
					fprintf(stderr, "error function[result_no_bipartite_graph]\n");
					exit(1);
				}
				if (res2 == index_base.end()) {
					fprintf(stderr, "error function[result_no_bipartite_graph]\n");
					exit(1);
				}
				std::pair<mfile0::M_Base, mfile0::M_Base> base02 = std::make_pair(res0->second, res2->second);
				for (auto itr_01 = path_01.begin(); itr_01 != path_01.end(); itr_01++) {
					auto res1 = index_base.find(itr_01->second);
					if (res1 == index_base.end()) {
						fprintf(stderr, "error function[result_no_bipartite_graph]\n");
						exit(1);
					}
					std::pair<mfile0::M_Base, mfile0::M_Base> base12 = std::make_pair(res1->second, res2->second);
					std::pair< int, int > add_path = std::make_pair(res1->first, res2->first);
					if (!check_insert_path(add_path, path_next))continue;
					//printf("add path:");
					//Print_path(base12.first, base12.second);
					//printf("\n");
					//1->2のpathを追加
					//0->2のpathを削除
					//path,group
					path.push_back(base12);
					group.cut_path.push_back(std::make_pair(path_id, std::make_tuple(
						base12.first.pos, base12.second.pos, base12.first.rawid, base12.second.rawid)));
					path_id++;
					//boost::unordered_multimap <int, int> path_next;
					//boost::unordered_multimap <int, int> path_prev;
					path_next.insert(std::make_pair(itr_01->second, itr_02->second));
					path_prev.insert(std::make_pair(itr_02->second, itr_01->second));
				}
				//printf("del path:");
				//Print_path(base02.first, base02.second);
				//printf("\n");
				for (auto itr3 = path.begin(); itr3 != path.end();) {
					if (*itr3 == base02) itr3 = path.erase(itr3);
					else itr3++;
				}
				auto del_path = std::make_tuple(base02.first.pos, base02.second.pos, base02.first.rawid, base02.second.rawid);
				for (auto itr3 = group.cut_path.begin(); itr3 != group.cut_path.end();) {
					if (itr3->second == del_path) itr3 = group.cut_path.erase(itr3);
					else itr3++;
				}
				for (auto itr3 = path_next.begin(); itr3 != path_next.end();) {
					if (itr3->first == itr_02->first && itr3->second == itr_02->second) itr3 = path_next.erase(itr3);
					else itr3++;
				}
				for (auto itr3 = path_prev.begin(); itr3 != path_prev.end();) {
					if (itr3->first == itr_02->second && itr3->second == itr_02->first) itr3 = path_prev.erase(itr3);
					else itr3++;
				}
			}
		}
		//処理後イテレータは最初に戻す
		itr = index_base.begin();
		finished.clear();
	}

}
bool check_insert_path(std::pair< int, int >& path, boost::unordered_multimap <int, int>& path_next) {
	if (path_next.count(path.first) == 0)return true;
	auto range = path_next.equal_range(path.first);
	for (auto res = range.first; res != range.second; res++) {
		if (res->first == path.first && res->second == path.second)return false;
	}
	return true;
}
std::map<int, int> add_depth_information(std::set<int>vertex, std::set<int>up, std::set<int>down, boost::unordered_multimap <int, int>& path_next) {
	int count = 0;
	std::map<int, int> ret;
	std::set<int> next;
	std::set<int> now;

	std::set<int>start;
	for (auto itr = up.begin(); itr != up.end(); itr++) {
		if (down.count(*itr) == 0) {
			ret.insert(std::make_pair(*itr, count));
			next.insert(*itr);
		}
	}
	while (next.size() != 0) {
		now = next;
		next.clear();
		count++;
		for (auto itr = now.begin(); itr != now.end(); itr++) {
			if (path_next.count(*itr) == 0)continue;
			auto range = path_next.equal_range(*itr);
			for (auto res = range.first; res != range.second; res++) {
				if (vertex.count(res->second) == 0)continue;
				auto res_in = ret.insert(std::make_pair(res->second, count));
				if (!res_in.second) {
					res_in.first->second = count;
				}
				next.insert(res->second);
			}
		}
	}
	return ret;
}

void output_linklet(std::string filename, std::vector<output_format>& out) {
	std::ofstream ofs(filename);
	double weight = 1;
	int link_num = 0;
	int all = out.size();
	int count = 0;
	for (auto itr = out.begin(); itr != out.end(); itr++) {

		if (count % 1000 == 0) {
			printf("\r write path %d/%d(%4.1lf%%)", count, all, count * 100. / all);
		}
		count++;


		ofs << std::fixed << std::right
			<< std::setw(5) << std::setprecision(0) << itr->groupid << " "
			<< std::setw(5) << std::setprecision(0) << itr->trackid << " "
			<< std::setw(3) << std::setprecision(0) << itr->num_comfirmed_path << " "
			<< std::setw(3) << std::setprecision(0) << itr->num_cut_path << " "
			<< std::setw(3) << std::setprecision(0) << itr->num_select_path << std::endl;
		for (auto itr2 = itr->comfirmed_path.begin(); itr2 != itr->comfirmed_path.end(); itr2++) {
			ofs << std::fixed << std::right
				<< std::setw(3) << std::setprecision(0) << itr2->first << " "
				<< std::setw(5) << std::setprecision(0) << std::get<0>(itr2->second) << " "
				<< std::setw(5) << std::setprecision(0) << std::get<1>(itr2->second) << " "
				<< std::setw(10) << std::setprecision(0) << std::get<2>(itr2->second) << " "
				<< std::setw(10) << std::setprecision(0) << std::get<3>(itr2->second) << std::endl;
		}

		for (auto itr2 = itr->cut_path.begin(); itr2 != itr->cut_path.end(); itr2++) {
			ofs << std::fixed << std::right
				<< std::setw(3) << std::setprecision(0) << itr2->first << " "
				<< std::setw(5) << std::setprecision(0) << std::get<0>(itr2->second) << " "
				<< std::setw(5) << std::setprecision(0) << std::get<1>(itr2->second) << " "
				<< std::setw(10) << std::setprecision(0) << std::get<2>(itr2->second) << " "
				<< std::setw(10) << std::setprecision(0) << std::get<3>(itr2->second) << std::endl;
		}

		for (auto itr2 = itr->select_path.begin(); itr2 != itr->select_path.end(); itr2++) {
			ofs << std::fixed << std::right
				<< std::setw(3) << std::setprecision(0) << itr2->first << " "
				<< std::setw(5) << std::setprecision(0) << std::get<0>(itr2->second) << " "
				<< std::setw(5) << std::setprecision(0) << std::get<1>(itr2->second) << " "
				<< std::setw(10) << std::setprecision(0) << std::get<2>(itr2->second) << " "
				<< std::setw(10) << std::setprecision(0) << std::get<3>(itr2->second) << std::endl;
		}
	}
	printf("\r write path %d/%d(%4.1lf%%)\n", count, all, count * 100. / all);


}


std::vector<std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>> divide_connected(std::vector<std::pair<int, std::tuple<int, int, int, int>>>& path_list) {
	//idを使わない
	std::multimap<std::pair<int, int>, std::pair<int, int>>connected_list;
	std::set< std::pair<int, int>>all_edge;
	for (auto itr = path_list.begin(); itr != path_list.end(); itr++) {
		std::pair<int, int> edge0, edge1;
		edge0 = std::make_pair(std::get<0>(itr->second), std::get<2>(itr->second));
		edge1 = std::make_pair(std::get<1>(itr->second), std::get<3>(itr->second));

		connected_list.insert(std::make_pair(edge0, edge1));
		connected_list.insert(std::make_pair(edge1, edge0));
		all_edge.insert(edge0);
		all_edge.insert(edge1);
	}
	std::vector<std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>> ret;
	std::set<std::pair<int, int>>finished_edge;
	bool loop_flg = true;
	while (true) {
		bool loop_flg = false;
		for (auto itr = all_edge.begin(); itr != all_edge.end(); itr++) {
			if (finished_edge.count(*itr) == 1) continue;
			loop_flg = true;

			std::set<std::pair<std::pair<int, int>, std::pair<int, int>>> connect_part;
			std::set < std::pair<int, int>> search_list;
			std::set < std::pair<int, int>> add_search_list;
			search_list.insert(*itr);
			while (true) {
				for (auto itr2 = search_list.begin(); itr2 != search_list.end(); itr2++) {
					if (connected_list.count(*itr2) == 0) {
						finished_edge.insert(*itr2);
						continue;
					}
					auto range = connected_list.equal_range(*itr2);
					for (auto res = range.first; res != range.second; res++) {
						auto edge0 = *itr2;
						auto edge1 = res->second;
						if (edge0.first > edge1.first)std::swap(edge0, edge1);
						connect_part.insert(std::make_pair(edge0, edge1));
						if (finished_edge.count(res->second) == 0) {
							add_search_list.insert(res->second);
						}
					}
					finished_edge.insert(*itr2);
				}
				search_list = add_search_list;
				if (add_search_list.size() == 0)break;
				add_search_list.clear();
			}
			std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> connect_part_v(connect_part.begin(), connect_part.end());
			ret.push_back(connect_part_v);
		}
		if (!loop_flg)break;
	}

	//for (auto &vv : ret) {
	//	for (auto &v : vv) {
	//		Print_path(v.first, v.second);
	//		printf("\n");
	//	}
	//	printf("\n");
	//}

	return ret;

}
//use
bool Judge_bipartite_graph_one(std::pair<int, int> target, boost::unordered_multimap <std::pair<int, int>, std::pair<int, int>>& path_prev, boost::unordered_multimap <std::pair<int, int>, std::pair<int, int>>& path_next, std::set<std::pair<int, int>>& up, std::set<std::pair<int, int>>& down) {
	up.insert(target);
	int sum = 0, sum_p = -1;
	while (sum != sum_p) {
		sum_p = sum;
		for (auto itr = up.begin(); itr != up.end(); itr++) {
			if (path_next.count(*itr) == 0)continue;
			auto range = path_next.equal_range(*itr);
			for (auto res = range.first; res != range.second; res++) {
				down.insert(res->second);
			}
		}
		for (auto itr = down.begin(); itr != down.end(); itr++) {
			if (path_prev.count(*itr) == 0)continue;
			auto range = path_prev.equal_range(*itr);
			for (auto res = range.first; res != range.second; res++) {
				up.insert(res->second);
			}
		}
		sum = up.size() + down.size();
	}

	////2部グラフの判定
	bool flg = true;
	for (auto itr = up.begin(); itr != up.end(); itr++) {
		if (down.count(*itr) == 1)flg = false;
	}
	for (auto itr = down.begin(); itr != down.end(); itr++) {
		if (up.count(*itr) == 1)flg = false;
	}
	return flg;
}
//use
bool Judge_bipartite_graph(int target, std::multimap <int, int>& path_prev, std::multimap<int, int>& path_next, std::set<int>& up, std::set<int>& down) {
	up.insert(target);
	int sum = 0, sum_p = -1;
	while (sum != sum_p) {
		sum_p = sum;
		for (auto itr = up.begin(); itr != up.end(); itr++) {
			if (path_next.count(*itr) == 0)continue;
			auto range = path_next.equal_range(*itr);
			for (auto res = range.first; res != range.second; res++) {
				down.insert(res->second);
			}
		}
		for (auto itr = down.begin(); itr != down.end(); itr++) {
			if (path_prev.count(*itr) == 0)continue;
			auto range = path_prev.equal_range(*itr);
			for (auto res = range.first; res != range.second; res++) {
				up.insert(res->second);
			}
		}
		sum = up.size() + down.size();
	}

	//2部グラフの判定
	bool flg = true;
	for (auto itr = up.begin(); itr != up.end(); itr++) {
		if (down.count(*itr) == 1)flg = false;
	}
	for (auto itr = down.begin(); itr != down.end(); itr++) {
		if (up.count(*itr) == 1)flg = false;
	}

	return flg;
}
//use
void decide_path_del_add(std::multimap<int, int>& path_next_id, std::multimap<int, int>& path_prev_id, std::set<std::pair<int, int>>& del_path, std::set < std::pair<int, int>>& add_path, std::map<int, std::pair<int, int>>& index_to_pair, std::pair<int, int>& path_01, std::pair<int, int>& path_02, std::pair<int, int>& path_12) {
	//02を削除
	//12を追加
	//12のposが同一-->何もしない
	//printf("path 01:(%d,%d)-(%d,%d)\n",
	//	index_to_pair.find(path_01.first)->second.first, index_to_pair.find(path_01.first)->second.second,
	//	index_to_pair.find(path_01.second)->second.first, index_to_pair.find(path_01.second)->second.second);
	//printf("path 12:(%d,%d)-(%d,%d)\n",
	//	index_to_pair.find(path_12.first)->second.first, index_to_pair.find(path_12.first)->second.second,
	//	index_to_pair.find(path_12.second)->second.first, index_to_pair.find(path_12.second)->second.second);
	//printf("path 02:(%d,%d)-(%d,%d)\n",
	//	index_to_pair.find(path_02.first)->second.first, index_to_pair.find(path_02.first)->second.second,
	//	index_to_pair.find(path_02.second)->second.first, index_to_pair.find(path_02.second)->second.second);

	auto vertex0 = index_to_pair.find(path_01.first);
	auto vertex1 = index_to_pair.find(path_12.first);
	auto vertex2 = index_to_pair.find(path_12.second);
	if (vertex0 == index_to_pair.end() || vertex1 == index_to_pair.end() || vertex2 == index_to_pair.end()) {
		fprintf(stderr, "error decide_path_del_add\n");
		exit(1);
	}
	//上下関係が成り立っているもののみ
	if (vertex0->second.first >= vertex1->second.first)return;
	if (vertex1->second.first >= vertex2->second.first)return;
	if (vertex0->second.first >= vertex2->second.first)return;

	bool add_flg1 = true;
	bool add_flg2 = true;
	for (auto itr3 = path_next_id.begin(); itr3 != path_next_id.end(); itr3++) {
		if (itr3->first == path_12.first && itr3->second == path_12.second)add_flg1 = false;
	}
	if (add_flg1) {
		path_next_id.insert(std::make_pair(path_12.first, path_12.second));
		path_prev_id.insert(std::make_pair(path_12.second, path_12.first));
		add_path.insert(path_12);
	}
	for (auto itr3 = path_next_id.begin(); itr3 != path_next_id.end(); itr3++) {
		if (itr3->first == path_01.first && itr3->second == path_01.second)add_flg2 = false;
	}
	if (add_flg2) {
		path_next_id.insert(std::make_pair(path_01.first, path_01.second));
		path_prev_id.insert(std::make_pair(path_01.second, path_01.first));
		add_path.insert(path_01);
	}

	del_path.insert(path_02);
	for (auto itr3 = path_next_id.begin(); itr3 != path_next_id.end(); ) {
		if (itr3->first == path_02.first && itr3->second == path_02.second) {
			itr3 = path_next_id.erase(itr3);
		}
		else {
			itr3++;
		}
	}
	for (auto itr3 = path_prev_id.begin(); itr3 != path_prev_id.end(); ) {
		if (itr3->first == path_02.second && itr3->second == path_02.first) {
			itr3 = path_prev_id.erase(itr3);
		}
		else {
			itr3++;
		}
	}

	//printf("del path:%d-%d\n", path_02.first, path_02.second);
	//printf("add path:%d-%d:%d\n", path_01.first, path_01.second, add_flg2);
	//printf("add path:%d-%d:%d\n", path_12.first, path_12.second, add_flg1);

}
//use
void decide_path_del_add(std::multimap<int, int>& path_next_id, std::multimap<int, int>& path_prev_id, std::set<std::pair<int, int>>& del_path, std::set < std::pair<int, int>>& add_path, std::map<int, std::pair<int, int>>& index_to_pair) {
	bool flg = true;
	while (flg) {
		flg = false;
		//printf("input id\n");
		//for (auto itr = path_next_id.begin(); itr != path_next_id.end(); itr++) {
		//	printf("%d->%d\n", itr->first, itr->second);
		//}
		std::vector<std::vector<std::pair<int, int>>> all_cycle = cycle_enumerate(path_next_id);
		//first-secondは任意
		std::multimap<int, std::vector<std::pair<int, int>>>cycle_map;
		for (int i = 0; i < all_cycle.size(); i++) {
			if (all_cycle[i].size() % 2 == 1) {
				cycle_map.insert(std::make_pair(all_cycle[i].size(), all_cycle[i]));
			}
			//for (int j = 0; j < all_cycle[i].size(); j++) {
			//	if (j + 1 == all_cycle[i].size()) {
			//		printf("[%d-->%d]\n", all_cycle[i][j].first, all_cycle[i][j].second);
			//	}
			//	else {
			//		printf("[%d-->%d]::", all_cycle[i][j].first, all_cycle[i][j].second);
			//	}
			//}
		}
		for (auto itr = cycle_map.begin(); itr != cycle_map.end(); itr++) {
			//解決できるサイクルを解決
			//最初へ戻る
			std::multimap <int, int> path_next;
			std::multimap <int, int> path_prev;
			std::set<int>up, down, fin, all;
			//path,vertexの設定
			for (auto itr2 = itr->second.begin(); itr2 != itr->second.end(); itr2++) {
				all.insert(itr2->first);
				all.insert(itr2->second);
				if (path_next_id.count(itr2->first) != 0) {
					auto range = path_next_id.equal_range(itr2->first);
					for (auto res = range.first; res != range.second; res++) {
						if (res->first == itr2->first && res->second == itr2->second) {
							path_next.insert(std::make_pair(itr2->first, itr2->second));
							path_prev.insert(std::make_pair(itr2->second, itr2->first));
							break;
						}
					}
				}
				if (path_next_id.count(itr2->second) != 0) {
					auto range = path_next_id.equal_range(itr2->second);
					for (auto res = range.first; res != range.second; res++) {
						if (res->first == itr2->second && res->second == itr2->first) {
							path_next.insert(std::make_pair(itr2->second, itr2->first));
							path_prev.insert(std::make_pair(itr2->first, itr2->second));
							break;
						}
					}
				}
			}
			//print
			//printf("hoge\n");
			//printf("vertex:%d\n", all.size());
			//printf("path_next:%d\n", path_next.size());
			//printf("path_prev:%d\n", path_prev.size());
			//for (auto itr2 = path_next.begin(); itr2 != path_next.end(); itr2++) {
			//	if (std::next(itr2, 1) != path_next.end()) {
			//		printf("%d - %d::", itr2->first, itr2->second);
			//	}
			//	else {
			//		printf("%d - %d\n", itr2->first, itr2->second);
			//	}
			//}
			//for (auto itr2 = path_prev.begin(); itr2 != path_prev.end(); itr2++) {
			//	if (std::next(itr2, 1) != path_prev.end()) {
			//		printf("%d - %d::", itr2->first, itr2->second);
			//	}
			//	else {
			//		printf("%d - %d\n", itr2->first, itr2->second);
			//	}
			//}


			for (auto itr2 = all.begin(); itr2 != all.end(); itr2++) {
				if (fin.count(*itr2) == 1)continue;
				up.clear();
				down.clear();
				if (Judge_bipartite_graph(*itr2, path_prev, path_next, up, down)) {
					for (auto itr = up.begin(); itr != up.end(); itr++) {
						fin.insert(*itr);
					}
					continue;
				}
				int merge_count = 0;
				int target;
				for (auto itr3 = up.begin(); itr3 != up.end(); itr3++) {
					if (down.count(*itr3) == 1) {
						merge_count += 1;
						target = *itr3;
					}
				}
				//printf("merge count=%d\n", merge_count);
				if (merge_count == 1) {
					//targetを起点に深さ付け
					//一周回ったときに+1 or -1でどっちにつながぐか決まる

					//解消
					std::pair<int, int>path_next01, path_next02, path_next12;
					std::pair<int, int>path_prev01, path_prev02, path_prev12;
					auto range_prev = path_prev.equal_range(target);
					auto range_next = path_next.equal_range(target);
					for (auto res = range_prev.first; res != range_prev.second; res++) {
						if (up.count(res->second) == 1) {
							path_next01.first = res->second;
							path_next01.second = res->first;
						}
					}
					for (auto res = range_next.first; res != range_next.second; res++) {
						if (down.count(res->second) == 1) {
							path_prev12.first = res->first;
							path_prev12.second = res->second;
						}
					}

					range_next = path_next.equal_range(path_next01.first);
					range_prev = path_prev.equal_range(path_prev12.second);
					for (auto res = range_next.first; res != range_next.second; res++) {
						if (down.count(res->second) == 1 && path_next01.second != res->second) {
							path_next02.first = res->first;
							path_next02.second = res->second;
						}
					}
					for (auto res = range_prev.first; res != range_prev.second; res++) {
						if (up.count(res->second) == 1 && path_prev12.first != res->second) {
							path_prev02.first = res->second;
							path_prev02.second = res->first;
						}
					}
					path_next12.first = path_next01.second;
					path_next12.second = path_next02.second;
					path_prev01.first = path_prev02.first;
					path_prev01.second = path_prev12.first;
					//pathの操作
					decide_path_del_add(path_next_id, path_prev_id, del_path, add_path, index_to_pair, path_prev01, path_prev02, path_prev12);
					decide_path_del_add(path_next_id, path_prev_id, del_path, add_path, index_to_pair, path_next01, path_next02, path_next12);
					//system("pause");
					flg = true;
					break;
				}

			}
			if (flg)break;
		}
	}

	//printf("loop fin\n");
	//for (auto itr = add_path.begin(); itr != add_path.end(); itr++) {
	//	printf("add:%d %d\n", itr->first, itr->second);
	//}
	//for (auto itr = del_path.begin(); itr != del_path.end(); itr++) {
	//	printf("del:%d %d\n", itr->first, itr->second);
	//}


}
void result_no_bipartite_graph(output_format& group, std::vector< std::pair<std::pair<int, int>, std::pair<int, int>>>& path, int& path_id) {
	boost::unordered_multimap <std::pair<int, int>, std::pair<int, int>> path_next;
	boost::unordered_multimap <std::pair<int, int>, std::pair<int, int>> path_prev;
	std::set<std::pair<int, int>> all_vertex;
	for (auto itr = path.begin(); itr != path.end(); itr++) {
		auto base0 = itr->first;
		auto base1 = itr->second;
		path_next.insert(std::make_pair(base0, base1));
		path_prev.insert(std::make_pair(base1, base0));

		all_vertex.insert(base0);
		all_vertex.insert(base1);
	}
	std::set<std::pair<int, int>> finished;
	for (auto itr = all_vertex.begin(); itr != all_vertex.end(); ) {
		std::pair<int, int> target = *itr;
		if (finished.count(target) == 1) {
			itr++;
			continue;
		}
		std::set<std::pair<int, int>>up, down;
		//2部グラフだった場合
		if (Judge_bipartite_graph_one(target, path_prev, path_next, up, down)) {
			for (auto itr2 = up.begin(); itr2 != up.end(); itr2++) {
				finished.insert(*itr2);
			}
			itr++;
			continue;
		}
		//2部グラフでない場合
		std::multimap<int, int>path_next_id;
		std::multimap<int, int>path_prev_id;
		std::map<int, std::pair<int, int>>index_to_pair;
		std::map< std::pair<int, int>, int>pair_to_index;
		int count = 0;
		for (auto itr2 = up.begin(); itr2 != up.end(); itr2++) {
			if (path_next.count(*itr2) == 0)continue;
			auto range = path_next.equal_range(*itr2);
			for (auto res = range.first; res != range.second; res++) {
				if (up.count(res->first) != 1 || down.count(res->second) != 1)continue;
				auto res_id0 = pair_to_index.insert(std::make_pair(res->first, count));
				if (res_id0.second) {
					index_to_pair.insert(std::make_pair(count, res->first));
					count++;
				}
				auto res_id1 = pair_to_index.insert(std::make_pair(res->second, count));
				if (res_id1.second) {
					index_to_pair.insert(std::make_pair(count, res->second));
					count++;
				}
				path_next_id.insert(std::make_pair(res_id0.first->second, res_id1.first->second));
				path_prev_id.insert(std::make_pair(res_id1.first->second, res_id0.first->second));
			}
		}
		//for (auto itr2 = index_to_pair.begin(); itr2 != index_to_pair.end(); itr2++) {
		//	printf("%d:%d,%d\n", itr2->first, itr2->second.first, itr2->second.second);
		//}
		//for (auto itr2 = path_next_id.begin(); itr2 != path_next_id.end(); itr2++) {
		//	printf("path:%d->%d\n", itr2->first, itr2->second);
		//}
		std::set<std::pair<int, int>>del_path, add_path;
		decide_path_del_add(path_next_id, path_prev_id, del_path, add_path, index_to_pair);
		//add path
		//del path
		//の処理
		//printf("path add %d\n", add_path.size());
		//printf("path del %d\n", del_path.size());
		//printf("add_path\n");
		for (auto itr2 = add_path.begin(); itr2 != add_path.end(); itr2++) {
			bool insert_flg = true;
			auto res0 = index_to_pair.find(itr2->first)->second;
			auto res1 = index_to_pair.find(itr2->second)->second;
			//printf("(%d,%d)-(%d,%d)\n", res0.first, res0.second, res1.first, res1.second);
			if (path_next.count(res0) == 0) insert_flg = true;
			else {
				auto range = path_next.equal_range(res0);
				for (auto res = range.first; res != range.second; res++) {
					if (res->second == res1)insert_flg = false;
				}
			}
			if (insert_flg) {
				std::pair<std::pair<int, int>, std::pair<int, int>>base12 = std::make_pair(res0, res1);
				path.push_back(base12);
				group.cut_path.push_back(std::make_pair(path_id, std::make_tuple(
					base12.first.first, base12.second.first, base12.first.second, base12.second.second)));
				path_id++;
				//boost::unordered_multimap <int, int> path_next;
				//boost::unordered_multimap <int, int> path_prev;
				path_next.insert(std::make_pair(res0, res1));
				path_prev.insert(std::make_pair(res1, res0));
			}
		}
		for (auto itr2 = del_path.begin(); itr2 != del_path.end(); itr2++) {
			auto res0 = index_to_pair.find(itr2->first)->second;
			auto res1 = index_to_pair.find(itr2->second)->second;
			std::pair<std::pair<int, int>, std::pair<int, int>>base02 = std::make_pair(res0, res1);

			for (auto itr3 = path.begin(); itr3 != path.end();) {
				if (*itr3 == base02) itr3 = path.erase(itr3);
				else itr3++;
			}
			auto del_path = std::make_tuple(base02.first.first, base02.second.first, base02.first.second, base02.second.second);
			for (auto itr3 = group.cut_path.begin(); itr3 != group.cut_path.end();) {
				if (itr3->second == del_path) itr3 = group.cut_path.erase(itr3);
				else itr3++;
			}
			for (auto itr3 = path_next.begin(); itr3 != path_next.end();) {
				if (itr3->first == res0 && itr3->second == res1) itr3 = path_next.erase(itr3);
				else itr3++;
			}
			for (auto itr3 = path_prev.begin(); itr3 != path_prev.end();) {
				if (itr3->first == res1 && itr3->second == res0) itr3 = path_prev.erase(itr3);
				else itr3++;
			}
		}
		//処理後イテレータは最初に戻す
		itr = all_vertex.begin();
		finished.clear();
		//system("pause");
	}

}


void path_id_reroll(output_format& g) {
	boost::unordered_multimap <std::pair<int, int>, std::pair<int, int>> path_next;
	boost::unordered_multimap<std::pair<int, int>, std::pair<int, int>> path_prev;
	std::set<std::pair<int, int>>all_vertex;
	std::pair<int, int>ver0, ver1;

	for (auto itr = g.comfirmed_path.begin(); itr != g.comfirmed_path.end(); itr++) {
		ver0.first = std::get<0>(itr->second);
		ver0.second = std::get<2>(itr->second);
		ver1.first = std::get<1>(itr->second);
		ver1.second = std::get<3>(itr->second);
		path_next.insert(std::make_pair(ver0, ver1));
		path_prev.insert(std::make_pair(ver1, ver0));
		all_vertex.insert(ver0);
		all_vertex.insert(ver1);
	}
	for (auto itr = g.cut_path.begin(); itr != g.cut_path.end(); itr++) {
		ver0.first = std::get<0>(itr->second);
		ver0.second = std::get<2>(itr->second);
		ver1.first = std::get<1>(itr->second);
		ver1.second = std::get<3>(itr->second);

		path_next.insert(std::make_pair(ver0, ver1));
		path_prev.insert(std::make_pair(ver1, ver0));
		all_vertex.insert(ver0);
		all_vertex.insert(ver1);
	}
	int path_id = 0;
	g.cut_path.clear();
	g.comfirmed_path.clear();

	std::set<std::pair<int, int>>start;
	for (auto itr = all_vertex.begin(); itr != all_vertex.end(); itr++) {
		if (path_prev.count(*itr) != 1) {
			start.insert(*itr);
			continue;
		}
		auto res = path_prev.find(*itr);
		if (path_next.count(res->second) > 1) {
			start.insert(*itr);
			continue;
		}
	}
	int count = 0;
	std::set < std::pair<std::pair<int, int>, std::pair<int, int>>>finished;

	for (auto itr = start.begin(); itr != start.end(); itr++) {
		auto now = *itr;
		count = 0;
		while (path_next.count(now) == 1) {
			auto res = path_next.find(now);
			if (path_prev.count(res->second) > 1)break;
			g.comfirmed_path.push_back(std::make_pair(path_id, std::make_tuple(
				res->first.first, res->second.first, res->first.second, res->second.second
			)));
			finished.insert(std::make_pair(res->first, res->second));
			now = res->second;
			count++;
		}
		if (count > 0) path_id++;
	}

	for (auto itr = path_next.begin(); itr != path_next.end(); itr++) {
		if (finished.count(*itr) == 1)continue;
		g.cut_path.push_back(std::make_pair(path_id, std::make_tuple(
			itr->first.first, itr->second.first, itr->first.second, itr->second.second
		)));
		path_id++;
	}
	g.num_comfirmed_path = g.comfirmed_path.size();
	g.num_cut_path = g.cut_path.size();
	g.select_path.clear();
	g.num_select_path = 0;

}

int use_thread(double ratio, bool output) {
	if (ratio >= 1) {
		ratio = 1;
	}
	else if (ratio < 0.1) {
		ratio = 0.1;
	}
	int num_all_thread = omp_get_max_threads();
	if (output) {
		printf("max thread = %d\n", num_all_thread);
		printf("ratio      = %3.2lf\n", ratio);
		printf("using... %d thread\n", int(num_all_thread * ratio));
	}
	return (int)(num_all_thread * ratio);
}
