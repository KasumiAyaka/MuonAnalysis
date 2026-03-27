//muonのchainに対して、linkletから再度chainを再生成する
//1group 1chainが前提
#define _CRT_SECURE_NO_WARNINGS
//大体2^50
#define DEFAULT_CHAIN_UPPERLIM 1000000000000000

//#pragma comment(lib, "VxxReader.lib")
//#include "VxxReader.h"
#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>
#pragma comment(lib, "LibL2c-x.lib")
#include <LibL2c-x.h>

static decltype(std::chrono::system_clock::now()) abort_time;
static bool time_limit_exist;
std::chrono::seconds max_time_limit(120);

class TimeupError : public std::runtime_error {
public:
	TimeupError(std::string s) : runtime_error(s) {}
};


#include "I:\NINJA\E71a\work\kasumi\ECC\MuonAnalysis\Chain_convolution_2.cpp\bipartite_graph_enumeration.h"
#include "I:\NINJA\E71a\work\kasumi\ECC\MuonAnalysis\Chain_convolution_2.cpp\Cycle_enumerate.h"

//#include "../Chain_convolution_1/bipartite_graph_enumeration.h"
//#include "../Chain_convolution_1/Cycle_enumerate.h"

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
#include <chrono>


using namespace l2c;
// 自分で作った型
struct Segment {
	int16_t pos;
	int32_t rawid;
	bool operator==(const Segment& rhs) const
	{
		return pos == rhs.pos && rawid == rhs.rawid;
	}
	bool operator<(const Segment& rhs) const {
		if (pos == rhs.pos) {
			return rawid < rhs.rawid;
		}
		return pos < rhs.pos;
	}
};
class Line {
public:
	std::vector<int> points;
	int prev_edge;
	int next_edge;

};

// 自分で作った型をunorderdコンテナに入れたいときは、operator== の他に
// 型と同じ名前空間で hash_value 関数を定義

size_t hash_value(const Segment& d)
{
	// 複数の値のハッシュ値を組み合わせてハッシュ値を計算するには、
	// boost::hash_combine を使います。
	size_t h = 0;
	boost::hash_combine(h, d.pos);
	boost::hash_combine(h, d.rawid);
	return h;
}

class Chain_baselist {
public:
	int64_t groupid, trackid;
	std::pair<int32_t, int64_t> target_track;
	std::set<std::pair<int32_t, int64_t>> btset;
	std::set<std::tuple<int, int, int, int>>ltlist;
	std::vector<int32_t> usepos;
	int chain_num;
	std::vector<Linklet> make_link_list();
	void set_usepos();
};
class Chain_baselist_compress : public Chain_baselist
{
public:
	std::set<std::pair<int32_t, int64_t>> comp_btset;
	std::multimap<std::tuple<int, int, int, int>, std::vector<std::pair<int32_t, int64_t>>>comp_ltlist;
	std::set<std::tuple<int, int, int, int>>cut_ltlist;
	std::vector<Linklet> make_comp_link_list();
	std::vector<Linklet> make_comp_link_list(std::set<std::tuple<int, int, int, int>>& cut_list);
	void set_comp_usepos();

};

class output_format {
public:
	int groupid, trackid;
	int num_comfirmed_path, num_cut_path, num_select_path;
	std::vector<std::pair<int, std::tuple<int, int, int, int>>> comfirmed_path;
	std::vector<std::pair<int, std::tuple<int, int, int, int>>> cut_path;
	std::vector<std::pair<int, std::tuple<int, int, int, int>>> select_path;
};

bool sort_M_Base(const mfile0::M_Base& left, const mfile0::M_Base& right) {

	if (left.pos == right.pos) {
		return left.rawid < right.rawid;
	}
	return left.pos < right.pos;
}
//型変換
std::pair<std::pair<int, int>, std::pair<int, int>> ltlist_to_pair(std::tuple<int, int, int, int> ltlist) {
	std::pair<std::pair<int, int>, std::pair<int, int>> ret;
	ret.first.first = std::get<0>(ltlist);
	ret.second.first = std::get<1>(ltlist);
	ret.first.second = std::get<2>(ltlist);
	ret.second.second = std::get<3>(ltlist);
	return ret;
}
std::tuple<int, int, int, int> pair_to_ltlist(std::pair<int, int>pair0, std::pair<int, int>pair1) {
	return std::make_tuple(
		pair0.first, pair1.first, pair0.second, pair1.second
	);
}
void Print_path(std::pair<int, int>edge0, std::pair<int, int>edge1) {
	printf("(%d,%d)-(%d,%d)", edge0.first, edge0.second, edge1.first, edge1.second);
}



std::vector < Chain_baselist > read_linklet_list2(std::string filename);
l2c::Cdat l2c_x(std::set<std::pair<int32_t, int64_t>>& btset, std::vector<Linklet>& ltlist, std::vector<int32_t>& usepos, int64_t upperlim = DEFAULT_CHAIN_UPPERLIM, bool output = false);
void path_id_reroll(output_format& g);
output_format change_format(Chain_baselist& c);

std::vector<std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>> divide_connected(std::vector<std::pair<int, std::tuple<int, int, int, int>>>& path_list);
void result_no_bipartite_graph(output_format& group, std::vector< std::pair<std::pair<int, int>, std::pair<int, int>>>& path, int& path_id);

Chain_baselist result_uniqiuebase_cycle(output_format& group, int& path_id);
Chain_baselist remove_branch_hangnail(Chain_baselist b);
Chain_baselist link_convolution(Chain_baselist b, l2c::Cdat& cdat);

std::vector < std::pair<std::set<std::pair<int, int>>, std::set<std::pair<int, int>>>> Divide_bipartite_graph(output_format& g, std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>& path);
void result_bipartite_graph_0(
	std::vector < std::pair<std::set<std::pair<int, int>>, std::set<std::pair<int, int>>>>& bipartite_graph_v,
	std::list<std::pair<int, std::tuple<int, int, int, int>>>& branch_path,
	std::set< std::tuple<int, int, int, int>>& del_path);


void path_id_reroll(output_format& g, l2c::Cdat& cdat, std::vector<int32_t>& usepos);
void gragh_cut(output_format& group, std::set<std::pair<int32_t, int64_t>>& btset, std::vector<Linklet>& ltlist, std::vector<int32_t>& usepos);





void output_linklet_list(std::string filename, std::vector<Chain_baselist>& chain_list);
void output_linklet_list_comp1(std::string filename, std::vector<Chain_baselist_compress>& chain_list);
void output_linklet_list_comp1_cut(std::string filename, std::vector<Chain_baselist_compress>& chain_list);
void output_linklet(std::string filename, std::vector<output_format>& out);

bool judeg_overflow(l2c::Cdat& cdat);
Chain_baselist link_convolution_closed_path(Chain_baselist b, l2c::Cdat& cdat);
Chain_baselist link_convolution_closed_path2(Chain_baselist b);
std::vector<Segment> cycle_pickup(mfile0::M_Chain& c0, mfile0::M_Chain& c1);
bool judeg_search(mfile0::M_Chain& c0, mfile0::M_Chain& c1);
void insert_all_cycles(std::vector<std::vector<Segment>>& all_cycles, std::vector<Segment>& input_seg);
bool judge_cycle_uniquebase(std::vector<Segment>& seg);
void link_convolution_closed_path3(Chain_baselist b);
void DFS_all_path(const std::vector < std::vector<int >>& G, int v, int p, std::vector<bool>& seen, std::stack<int>& hist, std::vector<bool>& finished, std::vector<std::vector<int>>& path);
void detect_branch_hangnail(boost::unordered_multimap<Segment, Segment>& link, boost::unordered_multimap<Segment, Segment>& link_inv, std::set<Segment>& seg_edge, std::set<Segment>& remove_vertex, std::set<std::tuple<int, int, int, int>>& remove_path);
Chain_baselist_compress path_compress(Chain_baselist b);
Chain_baselist_compress gragh_cut(Chain_baselist_compress& b);
bool Judge_complete_bipartite_graph(boost::unordered_multimap <int, int>& path_prev, boost::unordered_multimap <int, int>& path_next, std::set<int>& up, std::set<int>& down);
bool Judge_bipartite_graph(int target, boost::unordered_multimap<int, int>& path_prev, boost::unordered_multimap<int, int>& path_next, std::set<int>& up, std::set<int>& down);
bool Judge_bipartite_graph(int target, std::multimap <int, int>& path_prev, std::multimap<int, int>& path_next, std::set<int>& up, std::set<int>& down);

l2c::Cdat l2c_x_one_chain(std::set<std::pair<int32_t, int64_t>>& btset, std::vector<Linklet>& ltlist, std::vector<int32_t>& usepos);
Chain_baselist_compress select_cross_path(Chain_baselist_compress& b, l2c::Cdat& cdat);
output_format change_format(Chain_baselist_compress& b, l2c::Cdat& cdat);
output_format cut_path_organize(output_format& out);
void result_no_bipartite_graph_debug(output_format& group, std::vector< std::pair<std::pair<int, int>, std::pair<int, int>>>& path, int& path_id);
int count_linklet_max(std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>& connected_part);
std::vector<std::pair<int, std::tuple<int, int, int, int>>> cut_path_remove(std::vector<std::pair<int, std::tuple<int, int, int, int>>>& path_list);


void cleanup_aborted(std::vector<output_format>& out_form, std::set<int>& aborted_indexes)
{
	int n = out_form.size();
	for (int i = n - 1; i >= 0; i--) {
		if (aborted_indexes.find(i) != aborted_indexes.end()) {
			out_form.erase(out_form.begin() + i);

			printf("skipped: %d %d-%d\n", i, out_form[i].groupid, out_form[i].trackid);
		}
	}

	aborted_indexes.clear();

	return;
}

int main(int argc, char** argv) {
	time_limit_exist = false;

	if (argc != 4) {
		fprintf(stderr, "usage:prg file_in_link file_out_link mode\n");
		fprintf(stderr, "mode=0:linklet convollution\n");
		fprintf(stderr, "mode=1:no-bipartite graph remove\n");
		fprintf(stderr, "mode=2:unique base cycle remove\n");
		fprintf(stderr, "mode=3:1seg hangnail remove\n");
		fprintf(stderr, "mode=4:bipartite graph(no tangled) remove\n");
		fprintf(stderr, "mode=5:Complete bipartite graph decomposition\n");

		exit(1);
	}

	std::set<int> aborted_indexes;

	std::string file_in_link_list = argv[1];
	std::string file_out_link_list = argv[2];
	int mode = std::stoi(argv[3]);
	std::vector<Chain_baselist> chain_list = read_linklet_list2(file_in_link_list);

	std::vector<Chain_baselist_compress> chain_list_2;
	std::vector<output_format> out_form;
	int count = 0;
	//通常のlinklet畳み込み
	int interval = 1;
	for (int i = 0; i < chain_list.size(); i++) {
		if (i % interval == 0 || i == chain_list.size() - 1) {
			printf("%d/%d %d-%d\r", i, chain_list.size(), chain_list[i].groupid, chain_list[i].trackid);
		}

		//**ここは省略可能?**//
		chain_list[i].set_usepos();
		std::vector<Linklet> link = chain_list[i].make_link_list();
		l2c::Cdat cdat = l2c_x(chain_list[i].btset, link, chain_list[i].usepos, 1, false);
		chain_list[i] = link_convolution(chain_list[i], cdat);

		output_format group = change_format(chain_list[i]);
		out_form.push_back(group);
	}
	printf("\n");
	mode -= 1;
	if (mode < 0) {
		output_linklet(file_out_link_list, out_form);
		return 0;
	}

	count = 0;
	//非2部グラフの解消
	for (int i = 0; i < out_form.size(); i++) {
		if (count % interval == 0 || count + 1 == out_form.size()) {
			printf("%d/%d %d-%d\r", i, out_form.size(), out_form[i].groupid, out_form[i].trackid);
		}
		count++;
		//////////////////////
		int path_id = 0;
		for (auto& p : out_form[i].comfirmed_path) {
			path_id = std::max(path_id, p.first);
		}
		for (auto& p : out_form[i].cut_path) {
			path_id = std::max(path_id, p.first);
		}

		// 2025/05/05時点の仕様では、この段階でsekect_pathには何も入っていないと思われる。
		// 元々2引数をとるchange_formatの中でselect_pathを生成しているが、
		// 現在のコードでは上記関数はそもそも使われていない。
		for (auto& p : out_form[i].select_path) {
			path_id = std::max(path_id, p.first);
		}
		path_id += 1;

		out_form[i].cut_path = cut_path_remove(out_form[i].cut_path);

		//連結成分に分解
		std::vector<std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>> connect_path = divide_connected(out_form[i].cut_path);
		for (int j = 0; j < connect_path.size(); j++) {
			//printf("connect part %d/%d\n", j, connect_path.size());
			//printf("link num = %d\n", connect_path[j].size());
			int linklet_num_max = count_linklet_max(connect_path[j]);

			if (connect_path[j].size() > 133 * 2 || linklet_num_max > 50) {
				//printf("\nconnect part %d/%d\n", j, connect_path.size());
				//printf("link num = %d\n", connect_path[j].size());


				continue;
			}
			abort_time = std::chrono::system_clock::now() + max_time_limit;
			time_limit_exist = true;
			try {
				result_no_bipartite_graph(out_form[i], connect_path[j], path_id);
			}
			catch (TimeupError& e) {
				aborted_indexes.insert(i);
				std::cerr << "\nTime Up - skipped." << std::endl;
			}
			time_limit_exist = false;
		}
		path_id_reroll(out_form[i]);

	}
	printf("\n");

	cleanup_aborted(out_form, aborted_indexes);

	mode -= 1;
	if (mode < 0) {
		output_linklet(file_out_link_list, out_form);
		return 0;
	}

	//閉路の解消
	count = 0;
	std::vector<Chain_baselist> chain_list2;
	for (int i = 0; i < out_form.size(); i++) {
		if (count % interval == 0 || count + 1 == out_form.size()) {
			printf("%d/%d %d-%d\r", i, out_form.size(), out_form[i].groupid, out_form[i].trackid);

		}
		count++;
		//if (chain_list[i].groupid != 359700000)continue;

		//////////////////////
		int path_id = 0;
		for (auto& p : out_form[i].comfirmed_path) {
			path_id = std::max(path_id, p.first);
		}
		for (auto& p : out_form[i].cut_path) {
			path_id = std::max(path_id, p.first);
		}
		for (auto& p : out_form[i].select_path) {
			path_id = std::max(path_id, p.first);
		}
		path_id += 1;
		if (i != 2573) {
			chain_list2.push_back(result_uniqiuebase_cycle(out_form[i], path_id));
		}
		else {
			chain_list2.push_back(result_uniqiuebase_cycle(out_form[i], path_id));
		}
		//連結成分に分解
	}
	printf("\n");
	mode -= 1;
	if (mode < 0) {
		for (int i = 0; i < chain_list2.size(); i++) {
			out_form[i] = change_format(chain_list2[i]);
		}
		output_linklet(file_out_link_list, out_form);
		return 0;
	}

	//ささくれの除去
	count = 0;
	for (int i = 0; i < chain_list2.size(); i++) {
		if (count % interval == 0 || count + 1 == chain_list2.size()) {
			printf("res hangnail %d/%d %d-%d\r", i, chain_list2.size(), chain_list2[i].groupid, chain_list2[i].trackid);

		}
		count++;
		chain_list2[i] = remove_branch_hangnail(chain_list2[i]);
		out_form[i] = change_format(chain_list2[i]);
	}
	printf("\n");
	mode -= 1;
	if (mode < 0) {
		output_linklet(file_out_link_list, out_form);
		return 0;
	}

	//非2部グラフの解消
	count = 0;
	for (int i = 0; i < out_form.size(); i++) {
		if (count % interval == 0 || count + 1 == out_form.size()) {
			printf("res no_bipartite graph %d/%d %d-%d\r", i, out_form.size(), out_form[i].groupid, out_form[i].trackid);
		}
		count++;

		//////////////////////
		int path_id = 0;
		for (auto& p : out_form[i].comfirmed_path) {
			path_id = std::max(path_id, p.first);
		}
		for (auto& p : out_form[i].cut_path) {
			path_id = std::max(path_id, p.first);
		}
		for (auto& p : out_form[i].select_path) {
			path_id = std::max(path_id, p.first);
		}
		path_id += 1;

		//連結成分に分解
		std::vector<std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>> connect_path = divide_connected(out_form[i].cut_path);
		for (int j = 0; j < connect_path.size(); j++) {
			//if (out_form[i].groupid == 11520&&out_form[i].trackid==144) {
			//	printf("\nconnect part %d/%d\n", j, connect_path.size());
			//	printf("link num = %d\n", connect_path[j].size());
			//	for (int k = 0; k < connect_path[j].size(); k++) {
			//		printf("%d %d - %d %d\n", connect_path[j][k].first.first, connect_path[j][k].first.second, connect_path[j][k].second.first, connect_path[j][k].second.second);
			//	}
			//}

			abort_time = std::chrono::system_clock::now() + max_time_limit;
			time_limit_exist = true;
			try {
				result_no_bipartite_graph(out_form[i], connect_path[j], path_id);
			}
			catch (TimeupError& e) {
				aborted_indexes.insert(i);
				std::cerr << "\nTime Up - skipped." << std::endl;
			}
			time_limit_exist = false;

		}
		path_id_reroll(out_form[i]);

	}
	printf("\n");

	cleanup_aborted(out_form, aborted_indexes);

	//選び方1通りの2部グラフの解消
	count = 0;
	for (int i = 0; i < out_form.size(); i++) {
		if (count % interval == 0 || count + 1 == out_form.size()) {
			printf("%d/%d %d-%d\r", i, out_form.size(), out_form[i].groupid, out_form[i].trackid);
		}
		count++;

		//////////////////////
		int path_id = 0;
		for (auto& p : out_form[i].comfirmed_path) {
			path_id = std::max(path_id, p.first);
		}
		for (auto& p : out_form[i].cut_path) {
			path_id = std::max(path_id, p.first);
		}
		for (auto& p : out_form[i].select_path) {
			path_id = std::max(path_id, p.first);
		}
		path_id += 1;
		//連結成分に分解
		std::vector<std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>> connect_path = divide_connected(out_form[i].cut_path);
		std::set< std::tuple<int, int, int, int>>del_path;
		std::list<std::pair<int, std::tuple<int, int, int, int>>> single_path(out_form[i].comfirmed_path.begin(), out_form[i].comfirmed_path.end());
		std::list<std::pair<int, std::tuple<int, int, int, int>>> branch_path(out_form[i].cut_path.begin(), out_form[i].cut_path.end());

		abort_time = std::chrono::system_clock::now() + max_time_limit;
		time_limit_exist = true;
		try {

			for (int j = 0; j < connect_path.size(); j++) {
				auto bipartite_graph_v = Divide_bipartite_graph(out_form[i], connect_path[j]);
				//2部グラフを一個ずつ解決

				result_bipartite_graph_0(bipartite_graph_v, branch_path, del_path);

			}
			//不要なpathを消去
			for (auto itr = branch_path.begin(); itr != branch_path.end();) {
				if (del_path.count(itr->second) == 1) {
					itr = branch_path.erase(itr);
				}
				else {
					itr++;
				}
			}
			std::vector<std::pair<int, std::tuple<int, int, int, int>>> single_path_v(single_path.begin(), single_path.end());
			std::vector<std::pair<int, std::tuple<int, int, int, int>>> branch_path_v(branch_path.begin(), branch_path.end());
			out_form[i].comfirmed_path = single_path_v;
			out_form[i].cut_path = branch_path_v;
			path_id_reroll(out_form[i]);
		}
		catch (TimeupError& e) {
			aborted_indexes.insert(i);
			std::cerr << "\nTime Up - skipped." << std::endl;
		}
		time_limit_exist = false;

	}
	printf("\n");

	cleanup_aborted(out_form, aborted_indexes);

	mode -= 1;
	if (mode < 0) {
		output_linklet(file_out_link_list, out_form);
		return 0;
	}


	//閉路の解消
	count = 0;
	for (int i = 0; i < out_form.size(); i++) {
		if (count % interval == 0 || count + 1 == out_form.size()) {
			printf("%d/%d %d-%d\r", i, out_form.size(), out_form[i].groupid, out_form[i].trackid);
		}
		count++;

		//////////////////////
		int path_id = 0;
		for (auto& p : out_form[i].comfirmed_path) {
			path_id = std::max(path_id, p.first);
		}
		for (auto& p : out_form[i].cut_path) {
			path_id = std::max(path_id, p.first);
		}
		for (auto& p : out_form[i].select_path) {
			path_id = std::max(path_id, p.first);
		}
		path_id += 1;
		auto chain_tmp = result_uniqiuebase_cycle(out_form[i], path_id);
		out_form[i] = change_format(chain_tmp);

		//連結成分に分解
	}

	//完全2部グラフ分解
	count = 0;
	for (int i = 0; i < out_form.size(); i++) {
		if (count % interval == 0 || count + 1 == out_form.size()) {
			printf("%d/%d %d-%d\r", i, out_form.size(), out_form[i].groupid, out_form[i].trackid);
		}
		count++;
		//完全2部グラフ以外の切断
		std::set<std::pair<int32_t, int64_t>> btset;
		std::vector<Linklet> ltlist;
		std::vector<int32_t> usepos;
		gragh_cut(out_form[i], btset, ltlist, usepos);
		auto cdat = l2c_x(btset, ltlist, usepos, 1, false);
		path_id_reroll(out_form[i], cdat, usepos);
	}
	printf("\n");

	output_linklet(file_out_link_list, out_form);



	return 0;
	//printf("%d\n", chain_list[i].groupid);
	//chain_list[i].set_usepos();
	//std::vector<Linklet> link = chain_list[i].make_link_list();
	//if (chain_list[i].ltlist.size() == 0) {

	//}
	//else {
	//	//linkletの畳み込み
	//	l2c::Cdat cdat = l2c_x(chain_list[i].btset, link, chain_list[i].usepos, 1, false);
	//	chain_list[i] = link_convolution(chain_list[i], cdat);
	//	//閉路の畳み込み
	//	auto res = link_convolution_closed_path2(chain_list[i]);

	//	link = res.make_link_list();
	//	cdat = l2c_x(res.btset, link, res.usepos, 1, false);
	//	res = link_convolution(res, cdat);
	//	//ささくれの除去
	//	res = remove_branch_hangnail(res);
	//	//1本路の圧縮
	//	Chain_baselist_compress res2 = path_compress(res);
	//	//完全2部グラフ以外の切断
	//	res2 = gragh_cut(res2);
	//	link = res2.make_comp_link_list(res2.cut_ltlist);
	//	res2.set_comp_usepos();

	//	//連結成分の抽出
	//	cdat = l2c_x(res2.btset, link, res2.usepos, 1, false);
	//	//format変更
	//	auto out_tmp = change_format(res2, cdat);
	//	//圧縮の解除、切断の連結
	//	out_tmp = cut_path_organize(out_tmp);

	//	out_form.push_back(out_tmp);

	//	chain_list_2.push_back(res2);

	//}



	//output_linklet_list_comp1_cut(file_out_link_list, chain_list_2);
	//output_linklet(file_out_link_list, out_form);

}

std::vector < Chain_baselist > read_linklet_list2(std::string filename) {
	std::vector < Chain_baselist > ret;
	std::ifstream ifs(filename);
	int gid, tid, muon_num, t_pl, t_rawid;
	int64_t link_num;
	int pl, raw;
	double weight;
	std::tuple<int, int, int, int> link;
	int count = 0;

	while (ifs >> gid >> tid >> t_pl >> t_rawid >> link_num) {
		if (count % 1000 == 0) {
			printf("\r read group %d gid=%d link=%lld", count, gid, link_num);
		}
		count++;
		Chain_baselist c;
		c.groupid = gid;
		c.trackid = tid;
		c.target_track.first = t_pl;
		c.target_track.second = t_rawid;
		for (int64_t i = 0; i < link_num; i++) {
			ifs >> std::get<0>(link) >> std::get<1>(link) >> std::get<2>(link) >> std::get<3>(link);
			c.btset.insert(std::make_pair(std::get<0>(link), std::get<2>(link)));
			c.btset.insert(std::make_pair(std::get<1>(link), std::get<3>(link)));
			c.ltlist.insert(link);
		}
		ret.push_back(c);
	}
	printf("\r read group %d gid=%d link=%lld\n", count, gid, link_num);

	return ret;

}

// Chainのlinklet一覧を元に、1対1で結びつくlinkletとそうでないものの情報を取り出す。
// confirmed_pathは1対1。cut_pathは少なくとも一方が多。
output_format change_format(Chain_baselist& c) {
	output_format ret;
	ret.groupid = c.groupid;
	ret.trackid = c.trackid;
	boost::unordered_multimap<std::pair<int, int>, std::pair<int, int>> path_next;
	boost::unordered_multimap<std::pair<int, int>, std::pair<int, int>> path_prev;

	for (auto itr = c.ltlist.begin(); itr != c.ltlist.end(); itr++) {
		// hayakawa memo
		// lilistには2本のbasetrackの情報がstd::tupleに{ pl0, pl1, rawid0, rawid1 }という形式で格納されている。
		// これをltlist_to_pairによって、{ { pl0, rawid0 }, { pl1, rawid1 } }という形式に変換している。
		// なお、pl0 < pl1の関係にあると思われる。
		auto path = ltlist_to_pair(*itr);
		path_next.insert(std::make_pair(path.first, path.second));
		path_prev.insert(std::make_pair(path.second, path.second));
	}

	int path_id = 0;
	// path_nextはbt0をキーとするmultimapなので、走査するときもbt0によって整列されている。
	for (auto itr = path_next.begin(); itr != path_next.end(); itr++) {
		if (path_next.count(itr->first) == 1 && path_prev.count(itr->second) == 1) {
			// hayakawa memo
			// bt0とbt1から成るlinkletが1対1で結びつけられている場合。
			// つまり、bt0とbt1はともに、互い以外にlinkletを持たない。
			ret.comfirmed_path.push_back(std::make_pair(path_id, pair_to_ltlist(itr->first, itr->second)));
		}
		else {
			// hayakawa memo
			// bt0とbt1とで、互いの他に何かしら繋がるbasetrackが存在する場合。
			ret.cut_path.push_back(std::make_pair(path_id, pair_to_ltlist(itr->first, itr->second)));
		}
		path_id++;
	}

	path_id_reroll(ret);
	return ret;
}

std::vector<std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>> divide_connected(std::vector<std::pair<int, std::tuple<int, int, int, int>>>& path_list) {
	//idを使わない
	std::multimap<std::pair<int, int>, std::pair<int, int>>connected_list;//Basetrack間の双方向接続、つまり上下流双方向のLinklet
	std::set< std::pair<int, int>>all_edge;//全エッジ、という変数名だが、格納しているのはBasetrackのpl, rawidペア。なので全ノードという方が正しいような。
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


std::vector<std::pair<int, std::tuple<int, int, int, int>>> cut_path_remove(std::vector<std::pair<int, std::tuple<int, int, int, int>>>& path_list) {
	std::vector<std::pair<int, std::tuple<int, int, int, int>>> ret;
	std::set<std::tuple<int, int, int, int>>cut_path_list;

	//連結成分に分解
	std::vector<std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>> connect_path = divide_connected(path_list);

	for (int j = 0; j < connect_path.size(); j++) {
		double pl_link_average;
		int linknum = connect_path[j].size(), pl_min, pl_max;
		for (int i = 0; i < connect_path[j].size(); i++) {

			if (i == 0) {
				pl_min = connect_path[j][i].first.first / 10;
				pl_max = connect_path[j][i].second.first / 10;
			}
			pl_min = std::min(pl_min, connect_path[j][i].first.first / 10);
			pl_max = std::max(pl_max, connect_path[j][i].second.first / 10);
		}
		double pl_link_ave = linknum * 1. / (pl_max - pl_min + 1);
		int linklet_num_max = count_linklet_max(connect_path[j]);

		//printf("PL%03d-%03d %d : ave %.2lf\n", pl_min, pl_max, linknum, linknum*1. / (pl_max - pl_min + 1));
		//printf("connect part %d/%d\n", j, connect_path.size());
		//printf("link num = %d\n", connect_path[j].size());
		//printf("link num max= %d\n", linklet_num_max);

		if (pl_link_ave > 5 && linknum > 100 && linklet_num_max > 5) {
			//printf("\nconnect part %d/%d\n", j, connect_path.size());
			//printf("link num = %d\n", connect_path[j].size());
			std::multimap<std::pair<int, int>, std::pair<int, int>>path_next, path_prev;
			std::multimap<std::pair<int, int>, std::pair<int, int>>path_next_2x, path_prev_2x;
			for (auto& path : connect_path[j]) {
				if (path.second.first - path.first.first <= 20) {
					path_next.insert(std::make_pair(path.first, path.second));
					path_prev.insert(std::make_pair(path.second, path.first));
				}
				else {
					path_next_2x.insert(std::make_pair(path.first, path.second));
					path_prev_2x.insert(std::make_pair(path.second, path.first));
				}
			}
			for (auto itr = path_next_2x.begin(); itr != path_next_2x.end(); itr++) {
				//それしかpathがない場合は保存
				if (path_next.count(itr->first) == 0)continue;
				if (path_prev.count(itr->second) == 0)continue;
				cut_path_list.insert(std::make_tuple(itr->first.first, itr->second.first, itr->first.second, itr->second.second));
			}
			//printf("cut path num %d\n", cut_path_list.size());

		}
	}

	for (auto& p : path_list) {
		if (cut_path_list.count(p.second) == 0) {
			ret.push_back(p);
		}
	}
	//printf("path:%d-->%d\n", path_list.size(), ret.size());

	return ret;

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
void enumerate_path_dfs(std::vector<int>& path, std::vector<std::vector<int>>& all_hist, int now, int& goal, std::set<int>& seen, boost::unordered_multimap<int, int>& all_path) {
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
		enumerate_path_dfs(path, all_hist, res->second, goal, seen, all_path);
	}
	seen.insert(now);
	path.pop_back();

}
void enumerate_path_dfs(std::vector<std::pair<int, int>>& path, std::vector<std::vector<std::pair<int, int>>>& all_hist, std::pair<int, int> now, std::pair<int, int> start, std::pair<int, int> goal, std::set<std::pair<int, int>>& seen, boost::unordered_multimap<std::pair<int, int>, std::pair<int, int>>& all_path) {
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
		if (res->second.first < start.first)continue;
		if (res->second.first > goal.first)continue;


		enumerate_path_dfs(path, all_hist, res->second, start, goal, seen, all_path);
	}
	seen.insert(now);
	path.pop_back();

}
bool check_simple_cycle(std::vector<int>& path0, std::vector<int>& path1) {
	std::set<int> hist;
	for (int i = 1; i < path0.size() - 1; i++) {
		hist.insert(path0[i]);
	}
	for (int i = 1; i < path1.size() - 1; i++) {
		if (hist.count(path1[i]) == 1)return false;
	}
	return true;
}
bool judge_start(std::pair<int, int>start, std::set<std::pair<int, int>>& goal, boost::unordered_multimap <std::pair<int, int>, std::pair<int, int>>& path_next) {
	auto range = path_next.equal_range(start);
	std::set<int> pos_end;
	for (auto res = range.first; res != range.second; res++) {
		pos_end.insert(res->second.first);
	}
	if (pos_end.size() < 2)return false;
	int min_pos = *pos_end.begin();
	for (auto res = range.first; res != range.second; res++) {
		if (min_pos == res->second.first)continue;
		goal.insert(res->second);
	}
	return true;
}
//use
bool decide_path_del_add(std::multimap<int, int>& path_next_id, std::multimap<int, int>& path_prev_id, std::set<std::pair<int, int>>& del_path, std::set < std::pair<int, int>>& add_path, std::map<int, std::pair<int, int>>& index_to_pair, std::pair<int, int>& path_01, std::pair<int, int>& path_02, std::pair<int, int>& path_12) {
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
	if (vertex0->second.first >= vertex1->second.first)return false;
	if (vertex1->second.first >= vertex2->second.first)return false;
	if (vertex0->second.first >= vertex2->second.first)return false;

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
	return true;
}
//use
void decide_path_del_add(std::multimap<int, int>& path_next_id, std::multimap<int, int>& path_prev_id, std::set<std::pair<int, int>>& del_path, std::set < std::pair<int, int>>& add_path, std::map<int, std::pair<int, int>>& index_to_pair) {
	bool flg = true;
	while (flg) {
		//ここ無限ループ-->最大のgapをカットして解消
		flg = false;
		//printf("input id\n");
		//for (auto itr = path_next_id.begin(); itr != path_next_id.end(); itr++) {
		//	//printf("%d->%d\n", itr->first, itr->second);
		//	printf("path_next.insert(std::make_pair(%d,%d))\n", itr->first, itr->second);
		//}
		std::vector<std::vector<std::pair<int, int>>> all_cycle = cycle_enumerate(path_next_id);
		//printf("cycle_enumerate fin\n");

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
			/*
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
			*/

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
					//printf("target %d\n", target);
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
					bool path_del_flg0 = decide_path_del_add(path_next_id, path_prev_id, del_path, add_path, index_to_pair, path_prev01, path_prev02, path_prev12);
					bool path_del_flg1 = decide_path_del_add(path_next_id, path_prev_id, del_path, add_path, index_to_pair, path_next01, path_next02, path_next12);
					if (!path_del_flg0 && !path_del_flg1) {
						//gapが最大になる部分を切断
						int max_gap = 0;
						std::pair<int, int> path_del_wide_gap;
						for (auto itr2 = itr->second.begin(); itr2 != itr->second.end(); itr2++) {
							//printf("%d %d - %d %d\n"
							//	, index_to_pair.at(itr2->first).first
							//	, index_to_pair.at(itr2->first).second
							//	, index_to_pair.at(itr2->second).first
							//	, index_to_pair.at(itr2->second).second
							//);
							if (max_gap < abs(index_to_pair.at(itr2->first).first - index_to_pair.at(itr2->second).first) / 10) {
								max_gap = abs(index_to_pair.at(itr2->first).first - index_to_pair.at(itr2->second).first) / 10;
								path_del_wide_gap.first = itr2->first;
								path_del_wide_gap.second = itr2->second;
								if (index_to_pair.at(itr2->first).first > index_to_pair.at(itr2->second).first) {
									std::swap(path_del_wide_gap.first, path_del_wide_gap.second);
								}
							}
						}


						//printf("gap %d (%d,%d)\n", max_gap, path_del_wide_gap.first, path_del_wide_gap.second);

						del_path.insert(path_del_wide_gap);
						for (auto itr3 = path_next_id.begin(); itr3 != path_next_id.end(); ) {
							if (itr3->first == path_del_wide_gap.first && itr3->second == path_del_wide_gap.second) {
								itr3 = path_next_id.erase(itr3);
							}
							else {
								itr3++;
							}
						}
						for (auto itr3 = path_prev_id.begin(); itr3 != path_prev_id.end(); ) {
							if (itr3->first == path_del_wide_gap.second && itr3->second == path_del_wide_gap.first) {
								itr3 = path_prev_id.erase(itr3);
							}
							else {
								itr3++;
							}
						}

					}
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
	////}
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
		//printf("%d %d\n", itr->first, itr->second);
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



//use
void decide_path_del_add_debug(std::multimap<int, int>& path_next_id, std::multimap<int, int>& path_prev_id, std::set<std::pair<int, int>>& del_path, std::set < std::pair<int, int>>& add_path, std::map<int, std::pair<int, int>>& index_to_pair) {
	printf("start\n");
	bool flg = true;
	while (flg) {
		flg = false;
		printf("input id\n");

		std::set<int>tmp_id;
		for (auto itr = path_next_id.begin(); itr != path_next_id.end(); itr++) {
			printf("%d->%d\n", itr->first, itr->second);
			//printf("path_next.insert(std::make_pair(%d, %d));\n", itr->first, itr->second);
			tmp_id.insert(itr->first);
			tmp_id.insert(itr->second);

		}
		for (auto itr = tmp_id.begin(); itr != tmp_id.end(); itr++) {
			//printf("id:%d\n", *itr);
			int count_tmp = 0;
			if (path_next_id.count(*itr) != 0) {
				auto range = path_next_id.equal_range(*itr);
				for (auto res = range.first; res != range.second; res++) {
					if (count_tmp == 0) {
						printf("%d", res->second);
						count_tmp++;
					}
					else printf(",%d", res->second);
				}
			}
			if (path_prev_id.count(*itr) != 0) {
				auto range = path_prev_id.equal_range(*itr);
				for (auto res = range.first; res != range.second; res++) {
					if (count_tmp == 0) {
						printf("%d", res->second);
						count_tmp++;
					}
					else printf(",%d", res->second);
				}
			}
			printf("\n");
		}

		std::vector<std::vector<std::pair<int, int>>> all_cycle = cycle_enumerate(path_next_id);
		printf("enum fin\n");
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
	printf("end\n");

	//printf("loop fin\n");
	//for (auto itr = add_path.begin(); itr != add_path.end(); itr++) {
	//	printf("add:%d %d\n", itr->first, itr->second);
	//}
	//for (auto itr = del_path.begin(); itr != del_path.end(); itr++) {
	//	printf("del:%d %d\n", itr->first, itr->second);
	//}


}
void result_no_bipartite_graph_debug(output_format& group, std::vector< std::pair<std::pair<int, int>, std::pair<int, int>>>& path, int& path_id) {
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
		decide_path_del_add_debug(path_next_id, path_prev_id, del_path, add_path, index_to_pair);
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

Chain_baselist result_uniqiuebase_cycle(output_format& group, int& path_id) {
	std::multimap <std::pair<int, int>, std::pair<int, int>> path_next;
	std::multimap <std::pair<int, int>, std::pair<int, int>> path_prev;
	std::set<std::pair<int, int>> all_vertex;
	std::set<std::pair<int, int>> start, end;
	for (auto itr = group.comfirmed_path.begin(); itr != group.comfirmed_path.end(); itr++) {
		auto base0 = std::make_pair(std::get<0>(itr->second), std::get<2>(itr->second));
		auto base1 = std::make_pair(std::get<1>(itr->second), std::get<3>(itr->second));
		path_next.insert(std::make_pair(base0, base1));
		path_prev.insert(std::make_pair(base1, base0));

		all_vertex.insert(base0);
		all_vertex.insert(base1);
	}
	for (auto itr = group.cut_path.begin(); itr != group.cut_path.end(); itr++) {
		auto base0 = std::make_pair(std::get<0>(itr->second), std::get<2>(itr->second));
		auto base1 = std::make_pair(std::get<1>(itr->second), std::get<3>(itr->second));
		path_next.insert(std::make_pair(base0, base1));
		path_prev.insert(std::make_pair(base1, base0));

		all_vertex.insert(base0);
		all_vertex.insert(base1);
	}
	for (auto itr = path_next.begin(); itr != path_next.end(); itr++) {
		int count = path_next.count(itr->first);
		if (count < 2) {
			itr = std::next(itr, count - 1);
			continue;
		}
		auto range = path_next.equal_range(itr->first);
		for (auto res = range.first; res != range.second; res++) {
			if (std::next(res, 1) != range.second && res->second.first != std::next(res, 1)->second.first) {
				start.insert(res->first);
				break;
			}
		}

		itr = std::next(itr, count - 1);
	}
	for (auto itr = path_prev.begin(); itr != path_prev.end(); itr++) {
		int count = path_prev.count(itr->first);
		if (count < 2) {
			itr = std::next(itr, count - 1);
			continue;
		}
		auto range = path_prev.equal_range(itr->first);
		for (auto res = range.first; res != range.second; res++) {
			if (std::next(res, 1) != range.second && res->second.first != std::next(res, 1)->second.first) {
				end.insert(res->first);
				break;
			}
		}

		itr = std::next(itr, count - 1);
	}
	//if (group.groupid == 265) {
	//	printf("\n");
	//	printf("start size=%d\n", start.size());
	//	printf("end size=%d\n", end.size());
	//	for (auto itr = start.begin(); itr != start.end(); itr++) {
	//		printf("start:%d,%d\n", itr->first, itr->second);
	//	}
	//	for (auto itr = end.begin(); itr != end.end(); itr++) {
	//		printf("end:%d,%d\n", itr->first, itr->second);
	//	}
	//}
	std::vector < std::pair<std::set<std::pair<int, int>>, std::set<std::pair<int, int>>>>all_unique_path;


	for (int gap = 4; gap < 10; gap++) {

		for (auto itr_s = start.begin(); itr_s != start.end(); itr_s++) {
			for (auto itr_e = end.begin(); itr_e != end.end(); itr_e++) {
				if (itr_e->first / 10 - itr_s->first / 10 != gap)
					continue;
				// //4PL以上離れていること
				// if (itr_e->first / 10 - itr_s->first / 10 < 3)continue;
				/// //距離は最大10PL
				// if (itr_e->first / 10 - itr_s->first / 10 >= 10)continue;
				std::multimap<int, int>path_next_id;
				std::multimap<int, int>path_prev_id;
				std::map<int, std::pair<int, int>>index_to_pair;
				std::map< std::pair<int, int>, int>pair_to_index;
				int count = 0;
				std::set<std::pair<int, int>> add, now;
				add.insert(*itr_s);
				index_to_pair.insert(std::make_pair(count, *itr_s));
				pair_to_index.insert(std::make_pair(*itr_s, count));
				count++;
				while (add.size() != 0) {
					now = add;
					add.clear();
					for (auto itr = now.begin(); itr != now.end(); itr++) {
						if (path_next.count(*itr) == 0)continue;
						auto range = path_next.equal_range(*itr);
						for (auto res = range.first; res != range.second; res++) {
							if (res->second.first > itr_e->first)continue;
							auto res_in = pair_to_index.insert(std::make_pair(res->second, count));
							//insertに失敗-->探索済み
							if (res_in.second) {
								index_to_pair.insert(std::make_pair(count, res->second));
								add.insert(res->second);
								count++;
							}
							auto id0 = pair_to_index.find(res->first)->second;
							auto id1 = pair_to_index.find(res->second)->second;
							path_next_id.insert(std::make_pair(id0, id1));
						}
					}
				}

				auto start_id = pair_to_index.find(*itr_s);
				auto end_id = pair_to_index.find(*itr_e);
				//到達可能なpathは無い
				if (end_id == pair_to_index.end()) continue;
				//printf("(%d,%d)->(%d,%d)\n", itr_s->first, itr_s->second, itr_e->first, itr_e->second);

				std::vector<int >hist;
				std::vector<std::vector<std::pair<int, int>>>all_path_id;
				Enumerate_Path(path_next_id, hist, start_id->second, end_id->second, all_path_id);
				if (all_path_id.size() < 2)continue;
				std::vector<std::set<std::pair<int, int>>>all_path;
				for (int i = 0; i < all_path_id.size(); i++) {
					std::set<std::pair<int, int>> path;
					bool flg = true;
					for (int j = 0; j < all_path_id[i].size(); j++) {
						if (!flg)break;
						auto edge0 = index_to_pair.find(all_path_id[i][j].first);
						auto edge1 = index_to_pair.find(all_path_id[i][j].second);
						if (edge0->second.first >= edge1->second.first)flg = false;
						path.insert(edge0->second);
						path.insert(edge1->second);
					}
					if (flg) {
						all_path.push_back(path);
					}
				}
				for (int i = 0; i < all_path.size(); i++) {
					std::set <int> seg;
					for (auto itr = all_path[i].begin(); itr != all_path[i].end(); itr++) {
						seg.insert(itr->first);
					}
					for (int j = i + 1; j < all_path.size(); j++) {
						bool flg = true;
						for (auto itr = all_path[j].begin(); itr != all_path[j].end(); itr++) {

							if (!flg)break;
							if (*itr == *itr_e)continue;
							if (*itr == *itr_s)continue;
							if (seg.count(itr->first) == 1)flg = false;
						}
						if (flg) {
							all_unique_path.push_back(std::make_pair(all_path[i], all_path[j]));
						}
					}
				}
			}
		}
	}

	std::set<std::pair<std::pair<int, int>, std::pair<int, int>>>add_path;
	for (int i = 0; i < all_unique_path.size(); i++) {
		std::set<std::pair<int, int>>all_base;
		for (auto itr = all_unique_path[i].first.begin(); itr != all_unique_path[i].first.end(); itr++) {
			all_base.insert(*itr);
		}
		for (auto itr = all_unique_path[i].second.begin(); itr != all_unique_path[i].second.end(); itr++) {
			all_base.insert(*itr);
		}

		for (auto itr = all_base.begin(); itr != all_base.end(); itr++) {
			if (std::next(itr, 1) != all_base.end()) {
				std::pair<std::pair<int, int>, std::pair<int, int>> path_tmp;
				path_tmp.first = *itr;
				path_tmp.second = *(std::next(itr, 1));
				if (path_next.count(path_tmp.first) == 0)add_path.insert(path_tmp);
				else {
					auto range = path_next.equal_range(path_tmp.first);
					bool add_flg = true;
					for (auto res = range.first; res != range.second; res++) {
						if (res->first == path_tmp.first && res->second == path_tmp.second) {
							add_flg = false;
							break;
						}
					}
					if (add_flg)add_path.insert(path_tmp);
				}
			}
		}
	}


	//for (auto itr = add_path.begin(); itr != add_path.end(); itr++) {
	//	printf("add path:(%d,%d)-->(%d,%d)\n", itr->first.first, itr->first.second, itr->second.first, itr->second.second);
	//}

	Chain_baselist chains;
	chains.groupid = group.groupid;
	chains.trackid = group.trackid;
	for (auto itr = path_next.begin(); itr != path_next.end(); itr++) {
		chains.btset.insert(std::make_pair(itr->first.first, itr->first.second));
		chains.btset.insert(std::make_pair(itr->second.first, itr->second.second));
		chains.ltlist.insert(std::make_tuple(itr->first.first, itr->second.first, itr->first.second, itr->second.second));
	}
	for (auto itr = add_path.begin(); itr != add_path.end(); itr++) {
		chains.btset.insert(std::make_pair(itr->first.first, itr->first.second));
		chains.btset.insert(std::make_pair(itr->second.first, itr->second.second));
		chains.ltlist.insert(std::make_tuple(itr->first.first, itr->second.first, itr->first.second, itr->second.second));
	}

	chains.set_usepos();
	std::vector<Linklet> link = chains.make_link_list();
	l2c::Cdat cdat = l2c_x(chains.btset, link, chains.usepos, 1, false);
	chains = link_convolution(chains, cdat);
	return chains;

}

Chain_baselist remove_branch_hangnail(Chain_baselist b) {
	boost::unordered_multimap<Segment, Segment> link_id_up, link_id_down;
	Segment seg0, seg1;
	std::set<Segment> seg_all, down_seg_edge, up_seg_edge;
	for (auto itr = b.ltlist.begin(); itr != b.ltlist.end(); itr++) {
		seg0.pos = std::get<0>(*itr);
		seg0.rawid = std::get<2>(*itr);
		seg1.pos = std::get<1>(*itr);
		seg1.rawid = std::get<3>(*itr);
		seg_all.insert(seg0);
		seg_all.insert(seg1);
		link_id_up.insert(std::make_pair(seg0, seg1));
		link_id_down.insert(std::make_pair(seg1, seg0));
	}
	for (auto itr = seg_all.begin(); itr != seg_all.end(); itr++) {
		if (link_id_down.count(*itr) == 0) {
			down_seg_edge.insert(*itr);
		}
		if (link_id_up.count(*itr) == 0) {
			up_seg_edge.insert(*itr);
		}
	}
	std::set<Segment> remove_vertex;
	std::set<std::tuple<int, int, int, int>>remove_path;
	detect_branch_hangnail(link_id_up, link_id_down, up_seg_edge, remove_vertex, remove_path);
	//printf("\n");
	//for (auto itr = remove_vertex.begin(); itr != remove_vertex.end(); itr++) {
	//	printf("remove:%3d %10d\n", itr->pos, itr->rawid);
	//}
	//for (auto itr = remove_path.begin(); itr != remove_path.end(); itr++) {
	//	printf("remove path:%3d %10d\n", std::get<0>(*itr), std::get<1>(*itr), std::get<2>(*itr), std::get<3>(*itr));
	//}

	detect_branch_hangnail(link_id_down, link_id_up, down_seg_edge, remove_vertex, remove_path);
	//for (auto itr = remove_vertex.begin(); itr != remove_vertex.end(); itr++) {
	//	printf("remove:%3d %10d\n", itr->pos, itr->rawid);
	//}
	//for (auto itr = remove_path.begin(); itr != remove_path.end(); itr++) {
	//	printf("remove path:%3d %10d\n", std::get<0>(*itr), std::get<1>(*itr), std::get<2>(*itr), std::get<3>(*itr));
	//}

	Chain_baselist ret;
	for (auto itr = b.btset.begin(); itr != b.btset.end(); itr++) {
		seg0.pos = itr->first;
		seg0.rawid = itr->second;
		if (remove_vertex.count(seg0) == 1)continue;
		ret.btset.insert(*itr);
	}
	for (auto itr = b.ltlist.begin(); itr != b.ltlist.end(); itr++) {
		if (remove_path.count(*itr) == 1)continue;
		ret.ltlist.insert(*itr);
	}
	ret.groupid = b.groupid;
	ret.trackid = b.trackid;
	ret.target_track = b.target_track;
	ret.set_usepos();
	return ret;
}

//pathを2部グラフに分解
std::vector < std::pair<std::set<std::pair<int, int>>, std::set<std::pair<int, int>>>> Divide_bipartite_graph(output_format& g, std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>& path) {
	std::vector < std::pair<std::set<std::pair<int, int>>, std::set<std::pair<int, int>>>>bipartite_graph_v;
	boost::unordered_multimap <int, int> path_next;
	boost::unordered_multimap <int, int> path_prev;
	std::map<int, std::pair<int, int>>index_base;
	std::map<std::pair<int, int>, int>base_index;
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
			std::pair<std::set<std::pair<int, int>>, std::set<std::pair<int, int>>> bipartite_graph;
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

			printf("%d %d no bipartite graph\n", res->second.first, res->second.second);
			finished.insert(target);
		}


	}


	return bipartite_graph_v;
}
//2部グラフ一通りの解決
void result_bipartite_graph_0(
	std::vector < std::pair<std::set<std::pair<int, int>>, std::set<std::pair<int, int>>>>& bipartite_graph_v,
	std::list<std::pair<int, std::tuple<int, int, int, int>>>& branch_path,
	std::set<std::tuple<int, int, int, int>>& del_path) {

	//bipartite_graph_v　2部グラフの頂点
	//single_path&branch_path 2部グラフを作るpath
	for (auto& bg : bipartite_graph_v) {
		std::map<int, std::pair<int, int>> index_to_vertex;
		std::map< std::pair<int, int>, int> vertex_to_index;
		std::set<int>up, down;
		int count = 0;
		//頂点にiDを振り分ける
		for (auto itr = bg.first.begin(); itr != bg.first.end(); itr++) {
			std::pair<int, int> base0 = *itr;
			auto res = vertex_to_index.insert(std::make_pair(*itr, count));
			if (!res.second)continue;
			else {
				index_to_vertex.insert(std::make_pair(count, *itr));
				up.insert(count);
				count++;
			}
		}
		for (auto itr = bg.second.begin(); itr != bg.second.end(); itr++) {
			std::pair<int, int> base0 = *itr;
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
			std::pair<int, int> base0, base1;
			base0.first = std::get<0>(itr->second);
			base1.first = std::get<1>(itr->second);
			base0.second = std::get<2>(itr->second);
			base1.second = std::get<3>(itr->second);

			auto index0 = vertex_to_index.find(base0);
			auto index1 = vertex_to_index.find(base1);
			if (index0 == vertex_to_index.end() || index1 == vertex_to_index.end())continue;
			if (up.count(index0->second) == 0)continue;
			if (down.count(index1->second) == 0)continue;

			all_path.push_back(std::make_pair(index0->second, index1->second));
		}
		//printf("all_path\n");
		for (auto itr = all_path.begin(); itr != all_path.end(); itr++) {
			//printf("%d-%d\n", itr->first, itr->second);
			//printf("path_next.insert(std::make_pair(%d, %d));\n", itr->first, itr->second);
		}


		//とれる組み合わせの列挙
		std::vector<std::vector<std::pair<int, int>>> enumerat_all_path = Enumeration(all_path);

		//printf("enumerate end\n");
		//組が1通りの場合のみ解決
		if (enumerat_all_path.size() == 1) {
			std::set <std::pair<int, int>>selected_path(enumerat_all_path[0].begin(), enumerat_all_path[0].end());
			for (auto itr = all_path.begin(); itr != all_path.end(); itr++) {
				if (selected_path.count(*itr) == 0) {
					auto base0 = index_to_vertex.find(itr->first);
					auto base1 = index_to_vertex.find(itr->second);
					del_path.insert(std::make_tuple(base0->second.first, base1->second.first, base0->second.second, base1->second.second));
				}
			}
		}
	}
	return;
}




bool Judge_complete_bipartite_graph(boost::unordered_multimap <int, int>& path_prev, boost::unordered_multimap <int, int>& path_next, std::set<int>& up, std::set<int>& down) {
	bool flg = true;
	for (auto itr = up.begin(); itr != up.end(); itr++) {
		if (path_next.count(*itr) != down.size())return false;
	}
	for (auto itr = down.begin(); itr != down.end(); itr++) {
		if (path_prev.count(*itr) != up.size())return false;
	}
	return true;
}
void gragh_cut(output_format& group, std::set<std::pair<int32_t, int64_t>>& btset, std::vector<Linklet>& ltlist, std::vector<int32_t>& usepos) {

	std::map<int, std::pair<int, int>> index_vertex;
	std::map<std::pair<int, int>, int> vertex_index;
	std::set<int> all_edge;
	boost::unordered_multimap <int, int> path_prev, path_next;
	std::pair<int, int> seg0, seg1;
	int count = 0;
	for (auto itr = group.cut_path.begin(); itr != group.cut_path.end(); itr++) {
		seg0.first = std::get<0>(itr->second);
		seg0.second = std::get<2>(itr->second);
		seg1.first = std::get<1>(itr->second);
		seg1.second = std::get<3>(itr->second);

		auto res = index_vertex.insert(std::make_pair(count, seg0));
		if (res.second) {
			vertex_index.insert(std::make_pair(seg0, count));
			count++;
		}
		res = index_vertex.insert(std::make_pair(count, seg1));
		if (res.second) {
			vertex_index.insert(std::make_pair(seg1, count));
			count++;
		}

		auto id0 = vertex_index.find(seg0)->second;
		auto id1 = vertex_index.find(seg1)->second;
		all_edge.insert(id0);
		all_edge.insert(id1);
		path_next.insert(std::make_pair(id0, id1));
		path_prev.insert(std::make_pair(id1, id0));

	}

	std::unordered_set<int> fin_edge;
	std::set <std::pair<int, int>> delete_path;
	for (auto itr = all_edge.begin(); itr != all_edge.end(); itr++) {
		if (fin_edge.count(*itr) == 1)continue;
		//2部グラフ判定
		std::set<int>up, down;
		//2部グラフではない
		if (!Judge_bipartite_graph(*itr, path_prev, path_next, up, down)) {
			if (path_next.count(*itr) != 0) {
				auto range = path_next.equal_range(*itr);
				for (auto res = range.first; res != range.second; res++) {
					delete_path.insert(*res);
				}
			}
			fin_edge.insert(*itr);
		}
		//2部グラフである
		else {
			//完全2部グラフでない
			//完全2部グラフK (n,n) ではない
			if (!Judge_complete_bipartite_graph(path_prev, path_next, up, down)
				||
				up.size() != down.size()) {
				for (auto itr2 = up.begin(); itr2 != up.end(); itr2++) {
					auto range = path_next.equal_range(*itr2);
					for (auto res = range.first; res != range.second; res++) {
						delete_path.insert(*res);
					}
					fin_edge.insert(*itr2);
				}
			}
			else {
				//完全2部グラフK (n,n) である
				for (auto itr2 = up.begin(); itr2 != up.end(); itr2++) {
					fin_edge.insert(*itr2);
				}
			}

		}
	}

	std::set <std::tuple<int, int, int, int>> delete_path2;
	for (auto itr = delete_path.begin(); itr != delete_path.end(); itr++) {
		auto base0 = index_vertex.find(itr->first);
		auto base1 = index_vertex.find(itr->second);
		if (base0 == index_vertex.end() || base1 == index_vertex.end()) {
			fprintf(stderr, "error function[gragh_cut]\n");
			exit(1);
		}
		delete_path2.insert(std::make_tuple(base0->second.first, base1->second.first, base0->second.second, base1->second.second));
	}

	ltlist.clear();
	btset.clear();
	ltlist.reserve(group.comfirmed_path.size() + group.cut_path.size());
	std::set<int>pos;
	for (auto itr = group.comfirmed_path.begin(); itr != group.comfirmed_path.end(); itr++) {
		ltlist.emplace_back(std::get<0>(itr->second), std::get<1>(itr->second), std::get<2>(itr->second), std::get<3>(itr->second));
		btset.insert(std::make_pair(std::get<0>(itr->second), std::get<2>(itr->second)));
		btset.insert(std::make_pair(std::get<1>(itr->second), std::get<3>(itr->second)));
		pos.insert(std::get<0>(itr->second));
		pos.insert(std::get<1>(itr->second));
	}
	for (auto itr = group.cut_path.begin(); itr != group.cut_path.end(); itr++) {
		if (delete_path2.count(itr->second) == 1)continue;
		ltlist.emplace_back(std::get<0>(itr->second), std::get<1>(itr->second), std::get<2>(itr->second), std::get<3>(itr->second));
		btset.insert(std::make_pair(std::get<0>(itr->second), std::get<2>(itr->second)));
		btset.insert(std::make_pair(std::get<1>(itr->second), std::get<3>(itr->second)));
		pos.insert(std::get<0>(itr->second));
		pos.insert(std::get<1>(itr->second));
	}
	usepos.clear();
	for (auto itr = pos.begin(); itr != pos.end(); itr++) {
		usepos.push_back(*itr);
	}

	group.cut_path.clear();
	group.comfirmed_path.clear();
	for (auto itr = delete_path2.begin(); itr != delete_path2.end(); itr++) {
		group.cut_path.push_back(std::make_pair(0, *itr));
	}

}

void path_id_reroll(output_format& g, l2c::Cdat& cdat, std::vector<int32_t>& usepos) {

	int path_id = 0;
	size_t grsize = cdat.GetNumOfGroups();
	size_t possize = usepos.size();
	for (size_t grid = 0; grid < grsize; ++grid)
	{
		const Group& gr = cdat.GetGroup(grid);
		if (grid != gr.GetID()) throw std::exception("grid != gr.GetID().");//gridはgr.GetID()の戻り値と基本的に等しいはずである。
		int32_t spl = gr.GetStartPL();
		int32_t epl = gr.GetEndPL();
		//overflow==2chain以上
		if (gr.IsOverUpperLim()) {
			auto link_list = gr.GetLinklets();
			for (auto itr = link_list.begin(); itr != link_list.end(); itr++) {
				int32_t btpl0 = usepos[itr->first.GetPL()];
				int64_t btid0 = itr->first.GetRawID();
				int32_t btpl1 = usepos[itr->second.GetPL()];
				int64_t btid1 = itr->second.GetRawID();
				g.select_path.push_back(std::make_pair(path_id, std::make_tuple(btpl0, btpl1, btid0, btid1)));
			}
			path_id++;
		}
		else {
			auto link_list = gr.GetLinklets();
			for (auto itr = link_list.begin(); itr != link_list.end(); itr++) {
				int32_t btpl0 = usepos[itr->first.GetPL()];
				int64_t btid0 = itr->first.GetRawID();
				int32_t btpl1 = usepos[itr->second.GetPL()];
				int64_t btid1 = itr->second.GetRawID();
				g.comfirmed_path.push_back(std::make_pair(path_id, std::make_tuple(btpl0, btpl1, btid0, btid1)));
			}
			path_id++;
		}
	}
	for (auto itr = g.cut_path.begin(); itr != g.cut_path.end(); itr++) {
		itr->first = path_id;
		path_id++;
	}

	g.num_comfirmed_path = g.comfirmed_path.size();
	g.num_cut_path = g.cut_path.size();
	g.num_select_path = g.select_path.size();

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

l2c::Cdat l2c_x_one_chain(std::set<std::pair<int32_t, int64_t>>& btset, std::vector<Linklet>& ltlist, std::vector<int32_t>& usepos) {
	try
	{

		//std::vector<int32_t> usepos = { 310, 320, 330, 340, 350, 360, 370, 380, 390, 400, 410, 420, 430, 440, 450, 460, 470, 480, 490, 500 };
		size_t possize = usepos.size();

		//l2c関数本体を呼び出す。
		//出力されるCdatは基本的にl2c-xと同等のものになっているはずだが、若干の違いがある。
		//1. upperlimを超過したGroupは、Chainは持たず、全BaseTrackの情報だけを持った状態で出力される。
		//2. upperlim超過のGroupも含め、属す全BaseTrackの情報を持っている。
		//またplate枚数は最大64枚。変更するにはLibL2c-xのリビルドを要する。
		//デフォルトでoutput_isolated_linkletはfalse、upperlimは無効になっている。useposは与えないとエラーになる。
		l2c::Cdat cdat = l2c::MakeCdat_no_out(ltlist, opt::usepos = usepos, opt::upperlim = 1, opt::output_isolated_linklet = true);

		return cdat;

		//中身を出力する。
		auto grsize = cdat.GetNumOfGroups();
		for (size_t grid = 0; grid < grsize; ++grid)
		{
			const Group& gr = cdat.GetGroup(grid);
			size_t chsize = gr.GetNumOfChains();
			if (grid != gr.GetID()) throw std::exception("grid != gr.GetID().");//gridはgr.GetID()の戻り値と基本的に等しいはずである。
			int32_t spl = gr.GetStartPL();
			int32_t epl = gr.GetEndPL();

			//l2c-xのPLは一般に使われるPLとは意味が異なり、useposの中で"何番目のposであるか"を意味する。
			//LibL2c-xもl2c-xの仕様に倣うことにした（ただしl2c-xが1始まりなのに対してこちらは0始まり）。
			//例えば今回は、pos==310～500の20枚でChain Groupを作成しているが、
			//GetStartPL()の戻り値が"3"であったとすると、これはこのGroupがusepos[3]すなわち380から始まっているという意味になる。

			fprintf(stdout, "Group ID:%-9llu nchain:%-7llu start:%-4d end:%-4d\n", gr.GetID(), chsize, spl, epl);
			if (gr.IsOverUpperLim())
			{
				//upperlimを超過している場合、chainの情報はない。
				//ただしchainの本数はGetNumOfChainsで正しく取得できる。
				//また属すBaseTrackの情報は保持しているので、GetBaseTracksで全BaseTrackを取得できる。
				fprintf(stdout, "over upperlim. nchain:%-9lld\n", chsize);
				continue;
			}
			for (size_t ich = 0; ich < chsize; ++ich)
			{
				Chain ch = gr.GetChain(ich);//仕様上、Chainオブジェクトは一時変数として戻ってくる。
				int64_t chid = ch.GetID();
				int32_t nseg = ch.GetNSeg();
				int32_t spl = ch.GetStartPL();
				int32_t epl = ch.GetEndPL();
				fprintf(stdout, "    Chain ID:%-9lld nseg:%-3d start:%-4d end:%-4d\n", chid, nseg, spl, epl);

				for (size_t pl = 0; pl < possize; ++pl)
				{
					//ここで用いるplも上と同様に、useposの中で何番目のposかを意味する。
					//pl == 5だとすると、usepos中で（0から数えて）5番目、つまりpos==400のBaseTrack情報が返ってくる。
					BaseTrackID bt = ch.GetBaseTrack(pl);
					if (bt.IsEmpty()) continue;//IsEmptyがtrueのときは、このBaseTrackは歯抜け。
					int32_t btpl = bt.GetPL();
					int64_t btid = bt.GetRawID();
					//btplとplは厳密に一致しなければおかしい。
					if (btpl != pl) throw std::exception("BaseTrackID::PL != pl.");
					if (btset.find(std::make_pair(usepos[btpl], btid)) == btset.end())
					{
						//エラーチェック。
						//Chain内の全てのBaseTrackはbtset（ファイルから読み込まれた全BaseTrack）に必ず含まれているはず。
						//含まれていなければエラー。
						throw std::exception("BaseTrack is not found in btset.");
					}
					fprintf(stdout, "        pl:%-3d pos:%-4d rawid:%-9lld\n", btpl, usepos[pl], btid);
				}
			}

		}
	}
	catch (const std::exception& e)
	{
		fprintf(stderr, "%s\n", e.what());
	}

}

std::vector<Linklet> Chain_baselist::make_link_list() {
	std::vector<Linklet> ret;
	ret.reserve(ltlist.size());
	for (auto itr = ltlist.begin(); itr != ltlist.end(); itr++) {
		ret.emplace_back(std::get<0>(*itr), std::get<1>(*itr), std::get<2>(*itr), std::get<3>(*itr));
	}
	return ret;
}
std::vector<Linklet> Chain_baselist_compress::make_comp_link_list() {
	std::vector<Linklet> ret;
	ret.reserve(comp_ltlist.size());
	for (auto itr = comp_ltlist.begin(); itr != comp_ltlist.end(); itr++) {
		ret.emplace_back(std::get<0>(itr->first), std::get<1>(itr->first), std::get<2>(itr->first), std::get<3>(itr->first));
	}
	return ret;
}
std::vector<Linklet> Chain_baselist_compress::make_comp_link_list(std::set<std::tuple<int, int, int, int>>& cut_list) {
	std::vector<Linklet> ret;
	ret.reserve(comp_ltlist.size());
	for (auto itr = comp_ltlist.begin(); itr != comp_ltlist.end(); itr++) {
		if (cut_list.count(itr->first))continue;
		ret.emplace_back(std::get<0>(itr->first), std::get<1>(itr->first), std::get<2>(itr->first), std::get<3>(itr->first));
	}
	return ret;

}

void Chain_baselist::set_usepos() {
	std::set<int> pos_set;
	for (auto itr = btset.begin(); itr != btset.end(); itr++) {
		pos_set.insert(itr->first);
	}
	usepos.clear();
	for (auto itr = pos_set.begin(); itr != pos_set.end(); itr++) {
		usepos.push_back(*itr);
	}
}
void Chain_baselist_compress::set_comp_usepos() {
	std::set<int> pos_set;
	for (auto itr = comp_btset.begin(); itr != comp_btset.end(); itr++) {
		pos_set.insert(itr->first);
	}
	usepos.clear();
	for (auto itr = pos_set.begin(); itr != pos_set.end(); itr++) {
		usepos.push_back(*itr);
	}
}




bool judeg_overflow(l2c::Cdat& cdat) {
	size_t grsize = cdat.GetNumOfGroups();
	for (size_t grid = 0; grid < grsize; ++grid)
	{
		const Group& gr = cdat.GetGroup(grid);
		if (gr.IsOverUpperLim())return true;
	}
	return false;
}

// l2c-xで生成されたgroup情報から、basetrackとlinkletのstd::setを生成して返す。
Chain_baselist link_convolution(Chain_baselist b, l2c::Cdat& cdat) {
	Chain_baselist ret;
	ret.target_track = b.target_track;
	ret.groupid = b.groupid;
	ret.trackid = b.trackid;

	size_t grsize = cdat.GetNumOfGroups();
	size_t possize = b.usepos.size();
	for (size_t grid = 0; grid < grsize; ++grid)
	{
		const Group& gr = cdat.GetGroup(grid);
		if (grid != gr.GetID()) throw std::exception("grid != gr.GetID().");//gridはgr.GetID()の戻り値と基本的に等しいはずである。
		int32_t spl = gr.GetStartPL();
		int32_t epl = gr.GetEndPL();
		auto link_list = gr.GetLinklets();
		for (auto itr = link_list.begin(); itr != link_list.end(); itr++) {
			int32_t btpl0 = b.usepos[itr->first.GetPL()];
			int64_t btid0 = itr->first.GetRawID();
			int32_t btpl1 = b.usepos[itr->second.GetPL()];
			int64_t btid1 = itr->second.GetRawID();
			ret.btset.insert(std::make_pair(btpl0, btid0));
			ret.btset.insert(std::make_pair(btpl1, btid1));
			ret.ltlist.insert(std::make_tuple(btpl0, btpl1, btid0, btid1));
		}
	}
	ret.set_usepos();
	return ret;
}


void dfs(const std::vector < std::vector<int >>& G, int v, int p, std::vector<bool>& seen, std::stack<int>& hist, std::vector<bool>& finished, int& pos) {
	seen[v] = true;
	hist.push(v);
	//for (int i = 0; i < hist.size(); i++) {
	//	printf("%d --> ", hist[i]);
	//}
	//printf("\n");
	for (auto nv : G[v]) {
		if (nv == p) continue; // 逆流を禁止する

		// 完全終了した頂点はスルー
		if (finished[nv]) continue;

		//printf("%d --> %d\n", v, nv);
		// サイクルを検出
		if (seen[nv] && !finished[nv]) {
			pos = nv;
			return;
		}

		// 再帰的に探索
		dfs(G, nv, v, seen, hist, finished, pos);
		if (pos != -1) return;
	}
	hist.pop();
	finished[v] = true;
}
std::vector<int> outputcycle(std::stack<int>hist, int pos) {

	// サイクルを復元
	std::vector<int> cycle;
	while (!hist.empty()) {
		int t = hist.top();
		cycle.push_back(t);
		hist.pop();
		if (t == pos) break;
	}
	//for (int i = 0; i < cycle.size(); i++) {
	//	if (i + 1 != cycle.size()) {
	//		printf("%d-->", cycle[i]);
	//	}
	//	else {
	//		printf("%d\n", cycle[i]);
	//	}
	//}
	return cycle;

}

void dfs_search_cycle(const std::vector < std::vector<int >>& G, int v, int p, std::vector<bool>& seen, std::stack<int>& hist, std::vector<bool>& finished, int& pos, std::vector<std::vector<int>>& all_cycle) {
	seen[v] = true;
	hist.push(v);
	for (auto nv : G[v]) {
		if (nv == p) continue; // 逆流を禁止する

		// 完全終了した頂点はスルー
		if (finished[nv]) continue;

		//printf("%d --> %d\n", v, nv);
		// サイクルを検出
		if (seen[nv] && !finished[nv]) {
			pos = nv;
			all_cycle.push_back(outputcycle(hist, pos));
			continue;
		}

		// 再帰的に探索
		dfs_search_cycle(G, nv, v, seen, hist, finished, pos, all_cycle);
	}
	hist.pop();
	finished[v] = true;
}

void bfs_search_cycle_different_pl(Segment& start, boost::unordered_multimap<Segment, std::set<Segment>>& all_path, std::vector<std::vector<int>>& all_cycle) {

	std::vector<std::set<Segment>> path;
	std::set<Segment> finished;
	std::set<Segment> next;



}


std::vector<Segment> cycle_pickup(std::vector<Segment>& c0, std::vector<Segment>& c1) {
	//閉路1この保証あり
	std::set<Segment>btset;
	std::set<std::pair<Segment, Segment>> all_path;
	Segment seg0, seg1;
	for (int i = 0; i < c0.size(); i++) {
		btset.insert(c0[i]);
		if (i + 1 == c0.size())continue;
		all_path.insert(std::make_pair(c0[i], c0[i + 1]));
	}
	for (int i = 0; i < c1.size(); i++) {
		btset.insert(c1[i]);
		if (i + 1 == c1.size())continue;
		all_path.insert(std::make_pair(c1[i], c1[i + 1]));
	}

	std::map<int, Segment> index_vertex;
	std::map<Segment, int> vertex_index;
	int count = 0;
	for (auto itr = btset.begin(); itr != btset.end(); itr++) {
		seg0 = *itr;
		//seg0.rawid = itr->second;
		index_vertex.insert(std::make_pair(count, seg0));
		vertex_index.insert(std::make_pair(seg0, count));
		count++;
	}
	int N = count;
	std::vector < std::vector<int >> G(count);
	for (auto itr = all_path.begin(); itr != all_path.end(); itr++) {
		seg0 = itr->first;
		seg1 = itr->second;
		int index0 = vertex_index.at(seg0);
		int index1 = vertex_index.at(seg1);
		G[index0].push_back(index1);
		G[index1].push_back(index0);
	}
	// 探索
	std::vector<bool> seen, finished;
	seen.assign(N, false), finished.assign(N, false);
	int pos = -1;
	std::stack<int> hist; // 訪問履歴
	std::vector<std::set<int>>all_cycle;
	dfs(G, 0, -1, seen, hist, finished, pos);

	// サイクルを復元
	std::set<int> restored_cycle;
	while (!hist.empty()) {
		int t = hist.top();
		restored_cycle.insert(t);
		hist.pop();
		if (t == pos) break;
	}

	std::vector<Segment> cycle;
	for (auto itr = restored_cycle.begin(); itr != restored_cycle.end(); itr++) {
		cycle.push_back(index_vertex.at(*itr));
		//printf("(%d,%d)-->", index_vertex.at(*itr).pos, index_vertex.at(*itr).rawid);
		if (std::next(itr, 1) == restored_cycle.end()) {
			//printf("(%d,%d)\n", index_vertex.at(*restored_cycle.begin()).pos, index_vertex.at(*restored_cycle.begin()).rawid);
		}
	}
	return cycle;

}
bool judeg_search(std::vector<Segment>& c0, std::vector<Segment>& c1) {
	std::set<std::pair<int, int>>base_list;
	for (int i = 0; i < c0.size(); i++) {
		auto res = base_list.insert(std::make_pair(c0[i].pos, c0[i].rawid));
	}
	for (int i = 0; i < c1.size(); i++) {
		auto res = base_list.insert(std::make_pair(c1[i].pos, c1[i].rawid));
	}

	std::set<int>pos_list;
	for (auto itr = base_list.begin(); itr != base_list.end(); itr++) {
		auto res = pos_list.insert(itr->first);
		if (!res.second)return false;
	}
	return true;
}

Chain_baselist link_convolution_closed_path2(Chain_baselist b) {

	std::map<int, Segment> index_vertex;
	std::map<Segment, int> vertex_index;
	Segment seg0, seg1;
	int count = 0;
	for (auto itr = b.btset.begin(); itr != b.btset.end(); itr++) {
		seg0.pos = itr->first;
		seg0.rawid = itr->second;
		index_vertex.insert(std::make_pair(count, seg0));
		vertex_index.insert(std::make_pair(seg0, count));
		count++;
	}
	int N = count;
	std::vector < std::vector<int >> G(count);

	boost::unordered_map<Segment, std::set<Segment>>all_path;
	//同一posへのpathの切断
	for (auto itr = b.ltlist.begin(); itr != b.ltlist.end(); itr++) {
		seg0.pos = std::get<0>(*itr);
		seg0.rawid = std::get<2>(*itr);
		seg1.pos = std::get<1>(*itr);
		seg1.rawid = std::get<3>(*itr);
		auto res2 = all_path.find(seg0);
		if (res2 == all_path.end()) {
			std::set<Segment> set_tmp;
			set_tmp.insert(seg1);
			all_path.insert(std::make_pair(seg0, set_tmp));
		}
		else {
			res2->second.insert(seg1);
		}
		res2 = all_path.find(seg1);
		if (res2 == all_path.end()) {
			std::set<Segment> set_tmp;
			set_tmp.insert(seg0);
			all_path.insert(std::make_pair(seg1, set_tmp));
		}
		else {
			res2->second.insert(seg0);
		}
	}

	for (auto itr = all_path.begin(); itr != all_path.end();) {
		std::map<int, int> pos_map;
		for (auto itr2 = itr->second.begin(); itr2 != itr->second.end(); itr2++) {
			auto res = pos_map.insert(std::make_pair(itr2->pos, 1));
			if (!res.second)res.first->second++;
		}
		for (auto itr2 = itr->second.begin(); itr2 != itr->second.end();) {
			if (pos_map.at(itr2->pos) < 2) {
				itr2++;
			}
			else {
				itr2 = itr->second.erase(itr2);
			}
		}
		if (itr->second.size() == 0) {
			itr = all_path.erase(itr);
		}
		else {
			itr++;
		}
	}

	//経路の記憶
	for (auto itr = all_path.begin(); itr != all_path.end(); itr++) {
		for (auto itr2 = itr->second.begin(); itr2 != itr->second.end(); itr2++) {
			if (itr->first.pos > itr2->pos)continue;
			seg0 = itr->first;
			seg1 = *itr2;
			int index0 = vertex_index.at(seg0);
			int index1 = vertex_index.at(seg1);
			G[index0].push_back(index1);
			G[index1].push_back(index0);
		}
	}

	// 探索
	std::vector<bool> seen, finished;
	seen.assign(N, false), finished.assign(N, false);
	std::vector<std::vector<int>>all_cycle;
	for (int i = 0; i < N; i++) {
		int pos = -1;
		std::stack<int> hist; // 訪問履歴
		dfs_search_cycle(G, i, -1, seen, hist, finished, pos, all_cycle);
	}

	//for (int i = 0; i < all_cycle.size(); i++) {
	//	printf("%d %d\n", i, all_cycle[i].size());
	//	for (int j = 0; j < all_cycle[i].size(); j++) {
	//		auto ver0 = index_vertex.at(all_cycle[i][j]);
	//		auto ver1 = index_vertex.at(all_cycle[i][j]);
	//		if (j + 1 == all_cycle[i].size())ver1 = index_vertex.at(all_cycle[i][0]);
	//		else ver1 = index_vertex.at(all_cycle[i][j + 1]);
	//		if (ver0.pos > ver1.pos) std::swap(ver0, ver1);
	//		printf("%d %d %d %d\n", ver0.pos, ver0.rawid, ver1.pos, ver1.rawid);
	//	}
	//}

	//l2cを使って全経路を列挙
	std::vector<std::vector<Segment>>all_cycles;

	std::set<std::pair<int32_t, int64_t>> calc_path_btset;
	std::set<std::tuple<int, int, int, int>> set_calc_path_ltlist;
	std::vector<Linklet> calc_path_ltlist;
	std::set<int32_t> set_calc_path_usepos;
	std::vector<int32_t> calc_path_usepos;

	for (int i = 0; i < all_cycle.size(); i++) {
		for (int j = 0; j < all_cycle[i].size(); j++) {
			auto ver0 = index_vertex.at(all_cycle[i][j]);
			auto ver1 = index_vertex.at(all_cycle[i][j]);
			if (j + 1 == all_cycle[i].size())ver1 = index_vertex.at(all_cycle[i][0]);
			else ver1 = index_vertex.at(all_cycle[i][j + 1]);
			if (ver0.pos > ver1.pos) std::swap(ver0, ver1);
			calc_path_btset.insert(std::make_pair(ver0.pos, ver0.rawid));
			calc_path_btset.insert(std::make_pair(ver1.pos, ver1.rawid));
			set_calc_path_usepos.insert(ver0.pos);
			set_calc_path_usepos.insert(ver1.pos);
			set_calc_path_ltlist.insert(std::make_tuple(ver0.pos, ver1.pos, ver0.rawid, ver1.rawid));
			//printf("%d %d %d %d\n", ver0.pos, ver0.rawid, ver1.pos, ver1.rawid);
		}
	}
	calc_path_usepos.reserve(set_calc_path_usepos.size());
	for (auto itr = set_calc_path_usepos.begin(); itr != set_calc_path_usepos.end(); itr++) {
		calc_path_usepos.push_back(*itr);
	}
	calc_path_ltlist.reserve(set_calc_path_ltlist.size());
	for (auto itr = set_calc_path_ltlist.begin(); itr != set_calc_path_ltlist.end(); itr++) {
		calc_path_ltlist.emplace_back(std::get<0>(*itr), std::get<1>(*itr), std::get<2>(*itr), std::get<3>(*itr));
	}
	l2c::Cdat cdat_all_path = l2c_x(calc_path_btset, calc_path_ltlist, calc_path_usepos);

	size_t grsize = cdat_all_path.GetNumOfGroups();
	size_t possize = calc_path_usepos.size();
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
			std::vector<std::vector<Segment>>group_path;

			for (size_t ich = 0; ich < chsize; ++ich)
			{
				std::vector<Segment> chain_path;
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
					Segment seg;
					seg.pos = calc_path_usepos[btpl];
					seg.rawid = btid;
					chain_path.push_back(seg);
				}
				group_path.push_back(chain_path);
			}
			for (int i = 0; i < group_path.size(); i++) {
				for (int j = i + 1; j < group_path.size(); j++) {
					//2このchainで同一PLの異なるbasetrackを持つ場合
					if (!judeg_search(group_path[i], group_path[j]))continue;

					std::vector<Segment> cycle_seg = cycle_pickup(group_path[i], group_path[j]);
					if (!judge_cycle_uniquebase(cycle_seg)) continue;
					insert_all_cycles(all_cycles, cycle_seg);
				}
			}
		}
	}

	//for (int i = 0; i < all_cycles.size(); i++) {
	//	for (int j = 0; j < all_cycles[i].size(); j++) {
	//		printf("(%d,%d)-->", all_cycles[i][j].pos, all_cycles[i][j].rawid);
	//		if (j + 1 == all_cycles[i].size()) {
	//			printf("(%d,%d)\n", all_cycles[i][0].pos, all_cycles[i][0].rawid);
	//		}
	//	}
	//}

	Chain_baselist ret = b;
	std::map < std::pair<Segment, Segment>, std::vector<std::pair<Segment, Segment>>>replacement_link_list;
	for (int i = 0; i < all_cycles.size(); i++) {
		std::map<int, Segment> seg_map;
		for (int j = 0; j < all_cycles[i].size(); j++) {
			seg_map.insert(std::make_pair(all_cycles[i][j].pos, all_cycles[i][j]));
		}
		for (auto itr = seg_map.begin(); itr != seg_map.end(); itr++) {
			//printf("%d %d\n", itr->second.pos, itr->second.rawid);
			if (std::next(itr, 1) == seg_map.end())continue;
			auto itr2 = std::next(itr, 1);
			ret.ltlist.insert(std::make_tuple(
				itr->second.pos,
				itr2->second.pos,
				itr->second.rawid,
				itr2->second.rawid
			));
		}
	}
	ret.set_usepos();
	return ret;
}

bool judeg_search(mfile0::M_Chain& c0, mfile0::M_Chain& c1) {
	std::set<std::pair<int, int>>base_list;
	for (int i = 0; i < c0.basetracks.size(); i++) {
		auto res = base_list.insert(std::make_pair(c0.basetracks[i].pos, c0.basetracks[i].rawid));
	}
	for (int i = 0; i < c1.basetracks.size(); i++) {
		auto res = base_list.insert(std::make_pair(c1.basetracks[i].pos, c1.basetracks[i].rawid));
	}

	std::set<int>pos_list;
	for (auto itr = base_list.begin(); itr != base_list.end(); itr++) {
		auto res = pos_list.insert(itr->first);
		if (!res.second)return false;
	}
	return true;
}
std::vector<Segment> cycle_pickup(mfile0::M_Chain& c0, mfile0::M_Chain& c1) {
	//閉路1この保証あり
	std::set<Segment>btset;
	std::set<std::pair<Segment, Segment>> all_path;
	Segment seg0, seg1;
	for (int i = 0; i < c0.basetracks.size(); i++) {
		seg0.pos = c0.basetracks[i].pos;
		seg0.rawid = c0.basetracks[i].rawid;
		btset.insert(seg0);
		if (i + 1 == c0.basetracks.size())continue;
		seg1.pos = c0.basetracks[i + 1].pos;
		seg1.rawid = c0.basetracks[i + 1].rawid;

		all_path.insert(std::make_pair(seg0, seg1));
	}
	for (int i = 0; i < c1.basetracks.size(); i++) {
		seg0.pos = c1.basetracks[i].pos;
		seg0.rawid = c1.basetracks[i].rawid;
		btset.insert(seg0);
		if (i + 1 == c1.basetracks.size())continue;
		seg1.pos = c1.basetracks[i + 1].pos;
		seg1.rawid = c1.basetracks[i + 1].rawid;

		all_path.insert(std::make_pair(seg0, seg1));
	}

	std::map<int, Segment> index_vertex;
	std::map<Segment, int> vertex_index;
	int count = 0;
	for (auto itr = btset.begin(); itr != btset.end(); itr++) {
		seg0 = *itr;
		//seg0.rawid = itr->second;
		index_vertex.insert(std::make_pair(count, seg0));
		vertex_index.insert(std::make_pair(seg0, count));
		count++;
	}
	int N = count;
	std::vector < std::vector<int >> G(count);
	for (auto itr = all_path.begin(); itr != all_path.end(); itr++) {
		seg0 = itr->first;
		seg1 = itr->second;
		int index0 = vertex_index.at(seg0);
		int index1 = vertex_index.at(seg1);
		G[index0].push_back(index1);
		G[index1].push_back(index0);
	}
	// 探索
	std::vector<bool> seen, finished;
	seen.assign(N, false), finished.assign(N, false);
	int pos = -1;
	std::stack<int> hist; // 訪問履歴
	std::vector<std::set<int>>all_cycle;
	dfs(G, 0, -1, seen, hist, finished, pos);

	// サイクルを復元
	std::set<int> restored_cycle;
	while (!hist.empty()) {
		int t = hist.top();
		restored_cycle.insert(t);
		hist.pop();
		if (t == pos) break;
	}

	std::vector<Segment> cycle;
	for (auto itr = restored_cycle.begin(); itr != restored_cycle.end(); itr++) {
		cycle.push_back(index_vertex.at(*itr));
		//printf("(%d,%d)-->", index_vertex.at(*itr).pos, index_vertex.at(*itr).rawid);
		if (std::next(itr, 1) == restored_cycle.end()) {
			//printf("(%d,%d)\n", index_vertex.at(*restored_cycle.begin()).pos, index_vertex.at(*restored_cycle.begin()).rawid);
		}
	}
	return cycle;

}
bool judge_cycle_uniquebase(std::vector<Segment>& seg) {
	std::set<std::pair<int, int>>base_list;
	for (int i = 0; i < seg.size(); i++) {
		auto res = base_list.insert(std::make_pair(seg[i].pos, seg[i].rawid));
	}

	std::set<int>pos_list;
	for (auto itr = base_list.begin(); itr != base_list.end(); itr++) {
		auto res = pos_list.insert(itr->first);
		if (!res.second)return false;
	}
	return true;



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
void insert_all_cycles(std::vector<std::vector<Segment>>& all_cycles, std::vector<Segment>& input_seg) {
	std::vector < std::set<Segment>> all_cycles_set;
	for (int i = 0; i < all_cycles.size(); i++) {
		if (input_seg.size() != all_cycles[i].size())continue;
		std::set<Segment>cycles_set;
		for (int j = 0; j < all_cycles[i].size(); j++) {
			cycles_set.insert(all_cycles[i][j]);
		}
		all_cycles_set.push_back(cycles_set);
	}

	bool flg = true;
	int set_size = 0, hit_size = 0;
	for (int i = 0; i < all_cycles_set.size(); i++) {
		set_size = all_cycles_set[i].size();
		hit_size = 0;
		for (int j = 0; j < input_seg.size(); j++) {
			hit_size += all_cycles_set[i].count(input_seg[j]);
		}
		if (hit_size == set_size)flg = false;
	}
	if (flg) {
		all_cycles.push_back(input_seg);
	}

}
void detect_branch_hangnail(boost::unordered_multimap<Segment, Segment>& link, boost::unordered_multimap<Segment, Segment>& link_inv, std::set<Segment>& seg_edge, std::set<Segment>& remove_vertex, std::set<std::tuple<int, int, int, int>>& remove_path) {
	Segment seg0, seg1;
	bool flg = false;
	for (auto itr = seg_edge.begin(); itr != seg_edge.end(); itr++) {
		auto range = link_inv.equal_range(*itr);
		flg = false;
		for (auto res = range.first; res != range.second; res++) {
			if (link.count(res->second) < 2) {
				flg = false;
				break;
			}
			auto range2 = link.equal_range(res->second);
			for (auto res2 = range2.first; res2 != range2.second; res2++) {
				//分岐のもう片方も端点だった場合は残す
				if (seg_edge.count(res2->second) == 1)continue;
				flg = true;
			}
		}
		//flg=trueでささくれ判定
		if (!flg)continue;
		seg0 = *itr;
		remove_vertex.insert(seg0);

		for (auto res = range.first; res != range.second; res++) {
			seg1 = res->second;
			if (seg0.pos < seg1.pos) {
				remove_path.insert(std::make_tuple(
					seg0.pos,
					seg1.pos,
					seg0.rawid,
					seg1.rawid
				));
			}
			else {
				remove_path.insert(std::make_tuple(
					seg1.pos,
					seg0.pos,
					seg1.rawid,
					seg0.rawid
				));
			}
		}
	}
}


std::vector<std::set<int >> straight_line(const std::vector < std::vector<int >>& G, const std::vector < std::vector<int >>& G_inv, std::stack<int>& hist, std::vector<bool>& finished) {
	std::vector < std::set<int >> ret;
	int count = 0;
	std::set<int >path;
	while (!hist.empty()) {
		int t = hist.top();
		if (G[t].size() != 1) {
			count++;
		}
		else {
			finished[t] = true;
		}
		path.insert(t);
		if (G_inv[t].size() > 1) {
			ret.push_back(path);
			path.clear();
			path.insert(t);
		}

		if (count == 2)break;
		hist.pop();
	}
	ret.push_back(path);

	//for (auto itr = ret.begin(); itr != ret.end(); itr++) {
	//	if (std::next(itr, 1) != ret.end()) {
	//		printf("%d-->", *itr);
	//	}
	//	else {
	//		printf("%d\n", *itr);
	//	}
	//}
	return ret;
}
void print_hist(std::stack<int>hist) {
	while (!hist.empty()) {
		int t = hist.top();
		if (hist.size() > 1) {
			printf("%d-->", t);
		}
		else {
			printf("%d\n", t);
		}
		hist.pop();
	}
}
void DFS_straight_path(const std::vector < std::vector<int >>& G, const std::vector < std::vector<int >>& G_inv, int v, int p, std::vector<bool>& seen, std::stack<int>& hist, std::vector<bool>& finished, std::vector<std::set<int >>& clustered_index) {
	//print_hist(hist);
	seen[v] = true;
	hist.push(v);
	bool flg = true;
	for (auto nv : G[v]) {
		// 完全終了した頂点はスルー
		if (finished[nv]) continue;

		// 再帰的に探索
		DFS_straight_path(G, G_inv, nv, v, seen, hist, finished, clustered_index);
		//探索辺がfinishedになった場合
		if (finished[nv]) {
			//戻りながら分岐があるまでfinishにしていく
			//print_hist(hist);
			auto line_v = straight_line(G, G_inv, hist, finished);
			for (auto itr = line_v.begin(); itr != line_v.end(); itr++) {
				if (itr->size() < 3)continue;
				clustered_index.push_back(*itr);

			}

		}
		//if (pos != -1) return;
	}
	//hist.pop();
	finished[v] = true;
	//v = hist.top();
}
bool set_match(std::set<int>& set0, std::set<int >& set1) {
	if (set0.size() != set1.size())return false;
	int count = 0;
	for (auto itr = set0.begin(); itr != set0.end(); itr++) {
		count += set1.count(*itr);
	}
	return set0.size() == count;
}
bool judge_push_back_cluster(std::vector<std::set<int >>& clustered_index, std::vector<std::pair<std::set<int>, std::set<int >>>& rhombus_cycle_index, std::set<int >& cluster) {
	if (cluster.size() < 3)return false;
	int count = 0;
	for (int i = 0; i < clustered_index.size(); i++) {
		count = 0;
		for (auto itr = cluster.begin(); itr != cluster.end(); itr++) {
			if (clustered_index[i].count(*itr) == 1)count++;
		}
		if (count == cluster.size())return false;
		//3点以上同一 --> 同じcluster(端点の有無)
		if (count > 2) {
			for (auto itr = cluster.begin(); itr != cluster.end(); itr++) {
				clustered_index[i].insert(*itr);
			}
			return false;
		}
		else if (count == 2) {
			//ひし形閉路
			//rhombus_cycle_indexに入れる
			int match_count = 0;
			std::set<int> add0 = clustered_index[i];
			std::set<int> add1 = cluster;

			for (int j = 0; j < rhombus_cycle_index.size(); j++) {
				match_count = 0;
				std::set<int> ori0 = rhombus_cycle_index[j].first;
				std::set<int> ori1 = rhombus_cycle_index[j].second;
				if (set_match(add0, ori0))match_count += 1;
				if (set_match(add1, ori0))match_count += 1;
				if (set_match(add0, ori1))match_count += 1;
				if (set_match(add1, ori1))match_count += 1;
				if (match_count < 2) {
					rhombus_cycle_index.push_back(std::make_pair(add0, add1));
				}
				else if (match_count > 2) {
					printf("\n rhombus cycle path error\n");
					for (auto itr = cluster.begin(); itr != cluster.end(); itr++) {
						printf("%d ", *itr);
					}
					printf("\n");
					for (auto itr = clustered_index[i].begin(); itr != clustered_index[i].end(); itr++) {
						printf("%d ", *itr);
					}
					printf("\n");
				}
			}
			//同時にclustered_indexにも入る
			return true;
		}
	}
	return true;
}

void DFS_straight_path(const std::vector < std::vector<int >>& G, const std::vector < std::vector<int >>& G_inv, int start, std::vector<bool>& finished, std::vector<int>& hist, std::vector<std::set<int >>& clustered_index, std::vector<std::pair<std::set<int>, std::set<int >>>& rhombus_cycle_index) {
	//点startに立ち寄る
	hist.push_back(start);
	for (int i = 0; i < G[start].size(); i++) {
		int next = G[start][i];
		if (finished[next]) {
			//次の辺が探索済み+後方分岐の場合次の辺まで入れて探索
			if (G_inv[next].size() != 1) {
				std::set<int > cluster;
				cluster.insert(next);
				for (auto itr = hist.rbegin(); itr != hist.rend(); itr++) {
					bool flg = true;
					cluster.insert(*itr);
					if (G[*itr].size() != 1 || G_inv[*itr].size() != 1) {
						if (judge_push_back_cluster(clustered_index, rhombus_cycle_index, cluster)) {
							clustered_index.push_back(cluster);
						}
						cluster.clear();
						cluster.insert(*itr);
					}
				}
			}
		}
		DFS_straight_path(G, G_inv, next, finished, hist, clustered_index, rhombus_cycle_index);
	}
	//ここで点stratに対しての探索が終了する
	finished[start] = true;
	//startが一本道ではない場合
	if (G[start].size() != 1) {

		std::set<int > cluster;
		for (auto itr = hist.rbegin(); itr != hist.rend(); itr++) {
			bool flg = true;
			//if (std::next(itr, 1) != hist.rend()) {
			//	printf("%d-->", *itr);
			//}
			//else {
			//	printf("%d\n", *itr);
			//}
			cluster.insert(*itr);
			if (itr != hist.rbegin()) {
				if (G[*itr].size() != 1 || G_inv[*itr].size() != 1) {
					if (judge_push_back_cluster(clustered_index, rhombus_cycle_index, cluster)) {
						clustered_index.push_back(cluster);
					}
					cluster.clear();
					cluster.insert(*itr);
				}
			}

		}
	}
	hist.pop_back();
}
/////実装
void DFS_straight_path_compress(boost::unordered_multimap <int, Line>& line_prev, boost::unordered_multimap <int, Line>& line_next, int start, std::unordered_set<int>& edge, std::vector<std::vector<int>>& compress_point_list) {

	auto range = line_next.equal_range(start);
	for (auto res = range.first; res != range.second; res++) {
		std::vector<int> hist;
		hist.push_back(start);
		int next = res->second.next_edge;

		while (edge.count(next) == 0) {
			//line_next.count(next)==1&&line_prev.count(next)==1
			//の条件と等価
			hist.push_back(next);
			next = line_next.find(next)->second.next_edge;
		}
		hist.push_back(next);
		compress_point_list.push_back(hist);
	}

}
Chain_baselist_compress path_compress(Chain_baselist b) {
	std::map<int, Segment> index_vertex;
	std::map<Segment, int> vertex_index;
	Segment seg0, seg1;
	int count = 0;
	for (auto itr = b.btset.begin(); itr != b.btset.end(); itr++) {
		seg0.pos = itr->first;
		seg0.rawid = itr->second;
		index_vertex.insert(std::make_pair(count, seg0));
		vertex_index.insert(std::make_pair(seg0, count));
		count++;
	}
	int N = count;
	std::unordered_set<int> all_edge;
	boost::unordered_multimap <int, Line> line_prev, line_next;
	for (auto itr = b.ltlist.begin(); itr != b.ltlist.end(); itr++) {
		seg0.pos = std::get<0>(*itr);
		seg0.rawid = std::get<2>(*itr);
		seg1.pos = std::get<1>(*itr);
		seg1.rawid = std::get<3>(*itr);
		int index0 = vertex_index.at(seg0);
		int index1 = vertex_index.at(seg1);
		all_edge.insert(index0);
		all_edge.insert(index1);
		Line line;
		line.prev_edge = index0;
		line.next_edge = index1;
		line_next.insert(std::make_pair(index0, line));
		line_prev.insert(std::make_pair(index1, line));
	}
	//分岐,端点を記録
	std::unordered_set<int>edge;
	for (auto itr = all_edge.begin(); itr != all_edge.end(); itr++) {
		if (line_prev.count(*itr) != 1 || line_next.count(*itr) != 1) {
			edge.insert(*itr);
			//printf("branch: %d %d\n", index_vertex.at(*itr).pos, index_vertex.at(*itr).rawid);
		}
	}
	std::vector<std::vector<int>> compress_point_list;
	for (auto itr = edge.begin(); itr != edge.end(); itr++) {
		//edge-->edge の経路を縮退させる
		if (line_next.count(*itr) == 0)continue;
		DFS_straight_path_compress(line_prev, line_next, *itr, edge, compress_point_list);
	}
	//printf("compressed line\n");
	//for (int i = 0; i < compress_point_list.size(); i++) {
	//	for (int j = 0; j < compress_point_list[i].size(); j++) {
	//		if (j + 1 != compress_point_list[i].size()) {
	//			printf("(%d,%d)-->", index_vertex.at(compress_point_list[i][j]).pos, index_vertex.at(compress_point_list[i][j]).rawid);
	//		}
	//		else {
	//			printf("(%d,%d)\n", index_vertex.at(compress_point_list[i][j]).pos, index_vertex.at(compress_point_list[i][j]).rawid);
	//		}
	//	}
	//}
	std::unordered_set<int> new_edge;
	std::vector<Line>new_line;

	for (int i = 0; i < compress_point_list.size(); i++) {
		Line l;
		for (int j = 0; j < compress_point_list[i].size(); j++) {
			l.points.push_back(compress_point_list[i][j]);
		}
		l.prev_edge = *l.points.begin();
		l.next_edge = *l.points.rbegin();

		new_line.push_back(l);
		new_edge.insert(l.next_edge);
		new_edge.insert(l.prev_edge);
	}

	Chain_baselist_compress ret;
	ret.groupid = b.groupid;
	ret.target_track = b.target_track;
	ret.btset = b.btset;
	ret.ltlist = b.ltlist;

	for (auto itr = new_edge.begin(); itr != new_edge.end(); itr++) {
		auto res = index_vertex.at(*itr);
		ret.comp_btset.insert(std::make_pair(res.pos, res.rawid));
	}
	for (auto itr = new_line.begin(); itr != new_line.end(); itr++) {
		std::vector<std::pair<int32_t, int64_t>> comp_lt;
		std::tuple<int, int, int, int> lt;
		for (auto itr2 = itr->points.begin(); itr2 != itr->points.end(); itr2++) {
			auto res = index_vertex.at(*itr2);
			comp_lt.push_back(std::make_pair(res.pos, res.rawid));
		}
		auto res0 = index_vertex.at(itr->prev_edge);
		auto res1 = index_vertex.at(itr->next_edge);
		std::get<0>(lt) = res0.pos;
		std::get<1>(lt) = res1.pos;
		std::get<2>(lt) = res0.rawid;
		std::get<3>(lt) = res1.rawid;
		ret.comp_ltlist.insert(std::make_pair(lt, comp_lt));
	}

	return ret;
}

Chain_baselist_compress gragh_cut(Chain_baselist_compress& b) {

	std::map<int, Segment> index_vertex;
	std::map<Segment, int> vertex_index;
	Segment seg0, seg1;
	int count = 0;
	for (auto itr = b.comp_btset.begin(); itr != b.comp_btset.end(); itr++) {
		seg0.pos = itr->first;
		seg0.rawid = itr->second;
		index_vertex.insert(std::make_pair(count, seg0));
		vertex_index.insert(std::make_pair(seg0, count));
		count++;
	}
	int N = count;
	std::unordered_set<int> all_edge;
	boost::unordered_multimap <int, int> path_prev, path_next;
	for (auto itr = b.comp_ltlist.begin(); itr != b.comp_ltlist.end(); itr++) {
		seg0.pos = std::get<0>(itr->first);
		seg0.rawid = std::get<2>(itr->first);
		seg1.pos = std::get<1>(itr->first);
		seg1.rawid = std::get<3>(itr->first);
		int index0 = vertex_index.at(seg0);
		int index1 = vertex_index.at(seg1);
		all_edge.insert(index0);
		all_edge.insert(index1);
		path_next.insert(std::make_pair(index0, index1));
		path_prev.insert(std::make_pair(index1, index0));
	}

	std::unordered_set<int> fin_edge;
	boost::unordered_set <std::pair<int, int>> delete_path;
	for (auto itr = all_edge.begin(); itr != all_edge.end(); itr++) {
		if (fin_edge.count(*itr) == 1)continue;
		//2部グラフ判定
		std::set<int>up, down;
		//2部グラフではない
		if (!Judge_bipartite_graph(*itr, path_prev, path_next, up, down)) {
			if (path_next.count(*itr) != 0) {
				auto range = path_next.equal_range(*itr);
				for (auto res = range.first; res != range.second; res++) {
					delete_path.insert(*res);
				}
			}
			fin_edge.insert(*itr);
		}
		//2部グラフである
		else {
			//完全2部グラフでない
			//完全2部グラフK (n,n) ではない
			if (!Judge_complete_bipartite_graph(path_prev, path_next, up, down)
				||
				up.size() != down.size()) {
				for (auto itr2 = up.begin(); itr2 != up.end(); itr2++) {
					auto range = path_next.equal_range(*itr2);
					for (auto res = range.first; res != range.second; res++) {
						delete_path.insert(*res);
					}
					fin_edge.insert(*itr2);
				}
			}
			else {
				//完全2部グラフK (n,n) である
				for (auto itr2 = up.begin(); itr2 != up.end(); itr2++) {
					fin_edge.insert(*itr2);
				}
			}

		}
	}


	std::multimap<std::tuple<int, int, int, int>, std::vector<std::pair<int32_t, int64_t>>>comp_ltlist;
	for (auto itr = b.comp_ltlist.begin(); itr != b.comp_ltlist.end(); itr++) {
		seg0.pos = std::get<0>(itr->first);
		seg0.rawid = std::get<2>(itr->first);
		seg1.pos = std::get<1>(itr->first);
		seg1.rawid = std::get<3>(itr->first);
		int index0 = vertex_index.at(seg0);
		int index1 = vertex_index.at(seg1);

		if (delete_path.count(std::make_pair(index0, index1)) == 1) {
			//削除する場合
			b.cut_ltlist.insert(itr->first);
		}
	}
	return b;

}


output_format change_format(Chain_baselist_compress& b, l2c::Cdat& cdat) {
	output_format ret;
	//中身を出力する。
	size_t possize = b.usepos.size();
	//中身を出力する。
	size_t grsize = cdat.GetNumOfGroups();
	std::set<std::tuple<int, int, int, int>> saved_path;

	int path_id = 0;
	for (size_t grid = 0; grid < grsize; ++grid)
	{
		const Group& gr = cdat.GetGroup(grid);
		size_t chsize = gr.GetNumOfChains();
		if (grid != gr.GetID()) throw std::exception("grid != gr.GetID().");//gridはgr.GetID()の戻り値と基本的に等しいはずである。
		int32_t spl = gr.GetStartPL();
		int32_t epl = gr.GetEndPL();
		//fprintf(stdout, "Group ID:%-9llu nchain:%-7llu start:%-4d end:%-4d\n", gr.GetID(), chsize, spl, epl);
		if (gr.IsOverUpperLim())
		{
			auto basetracks = gr.GetBaseTracks();
			std::set<std::pair<int, int>>base_list;
			for (auto itr = basetracks.begin(); itr != basetracks.end(); itr++) {
				int pos = b.usepos[itr->GetPL()];
				int rawid = itr->GetRawID();
				base_list.insert(std::make_pair(pos, rawid));
			}
			for (auto itr = b.comp_ltlist.begin(); itr != b.comp_ltlist.end(); itr++) {
				std::pair<int, int> bt0, bt1;
				bt0.first = std::get<0>(itr->first);
				bt0.second = std::get<2>(itr->first);
				bt1.first = std::get<1>(itr->first);
				bt1.second = std::get<3>(itr->first);
				if (base_list.count(bt0) == 1 && base_list.count(bt1) == 1) {
					for (int i = 0; i < itr->second.size() - 1; i++) {
						ret.select_path.push_back(std::make_pair(path_id,
							std::make_tuple(itr->second[i].first, itr->second[i + 1].first, itr->second[i].second, itr->second[i + 1].second)));
					}
				}
				saved_path.insert(itr->first);
			}
			path_id++;
			//upperlimを超過している場合、chainの情報はない。
			//ただしchainの本数はGetNumOfChainsで正しく取得できる。
			//また属すBaseTrackの情報は保持しているので、GetBaseTracksで全BaseTrackを取得できる。
			//fprintf(stdout, "over upperlim. nchain:%-9lld\n", chsize);
			continue;
		}
		for (size_t ich = 0; ich < chsize; ++ich)
		{
			std::vector<std::pair<int, int>> chain;
			Chain ch = gr.GetChain(ich);//仕様上、Chainオブジェクトは一時変数として戻ってくる。
			int64_t chid = ch.GetID();
			int32_t nseg = ch.GetNSeg();
			int32_t spl = ch.GetStartPL();
			int32_t epl = ch.GetEndPL();
			//fprintf(stdout, "    Chain ID:%-9lld nseg:%-3d start:%-4d end:%-4d\n", chid, nseg, spl, epl);

			for (size_t pl = 0; pl < possize; ++pl)
			{
				//ここで用いるplも上と同様に、useposの中で何番目のposかを意味する。
				//pl == 5だとすると、usepos中で（0から数えて）5番目、つまりpos==400のBaseTrack情報が返ってくる。
				BaseTrackID bt = ch.GetBaseTrack(pl);
				if (bt.IsEmpty()) continue;//IsEmptyがtrueのときは、このBaseTrackは歯抜け。
				int32_t btpl = bt.GetPL();
				int64_t btid = bt.GetRawID();
				//btplとplは厳密に一致しなければおかしい。
				if (btpl != pl) throw std::exception("BaseTrackID::PL != pl.");
				if (b.btset.find(std::make_pair(b.usepos[btpl], btid)) == b.btset.end())throw std::exception("BaseTrack is not found in btset.");
				chain.push_back(std::make_pair(b.usepos[pl], btid));
			}
			for (int i = 0; i < chain.size() - 1; i++) {
				auto path_list = b.comp_ltlist.find(std::make_tuple(chain[i].first, chain[i + 1].first, chain[i].second, chain[i + 1].second));
				if (path_list == b.comp_ltlist.end()) {
					fprintf(stderr, "exception\n");
				}
				for (int j = 0; j < path_list->second.size() - 1; j++) {
					ret.comfirmed_path.push_back(std::make_pair(path_id,
						std::make_tuple(path_list->second[j].first, path_list->second[j + 1].first, path_list->second[j].second, path_list->second[j + 1].second)));
				}
				saved_path.insert(path_list->first);
			}
			path_id++;
		}
	}

	for (auto itr = b.cut_ltlist.begin(); itr != b.cut_ltlist.end(); itr++) {
		if (b.comp_ltlist.count(*itr) == 0)continue;
		auto range = b.comp_ltlist.equal_range(*itr);
		for (auto res = range.first; res != range.second; res++) {
			for (int i = 0; i < res->second.size() - 1; i++) {
				ret.cut_path.push_back(std::make_pair(path_id,
					std::make_tuple(res->second[i].first, res->second[i + 1].first, res->second[i].second, res->second[i + 1].second)));
			}
			path_id++;
			saved_path.insert(res->first);
		}
	}
	for (auto itr = b.comp_ltlist.begin(); itr != b.comp_ltlist.end(); itr++) {
		if (saved_path.count(itr->first) == 1)continue;
		if (b.comp_ltlist.count(itr->first) == 0)continue;
		auto range = b.comp_ltlist.equal_range(itr->first);
		for (auto res = range.first; res != range.second; res++) {
			for (int i = 0; i < res->second.size() - 1; i++) {
				ret.comfirmed_path.push_back(std::make_pair(path_id,
					std::make_tuple(res->second[i].first, res->second[i + 1].first, res->second[i].second, res->second[i + 1].second)));
			}
			path_id++;
			saved_path.insert(res->first);
		}
	}

	ret.num_comfirmed_path = ret.comfirmed_path.size();
	ret.num_cut_path = ret.cut_path.size();
	ret.num_select_path = ret.select_path.size();

	ret.groupid = b.groupid;
	ret.trackid = b.trackid;
	return ret;
}
output_format cut_path_organize(output_format& out) {

	std::multimap<std::pair<int, int>, std::pair<int, std::pair<int, int>>>path_next, path_prev;

	std::set<int> select_path, comfirmed_path, cut_path;
	std::pair<int, int> seg0, seg1;
	for (auto itr = out.select_path.begin(); itr != out.select_path.end(); itr++) {
		select_path.insert(itr->first);
		seg0.first = std::get<0>(itr->second);
		seg0.second = std::get<2>(itr->second);
		seg1.first = std::get<1>(itr->second);
		seg1.second = std::get<3>(itr->second);
		path_next.insert(std::make_pair(seg0, std::make_pair(itr->first, seg1)));
		path_prev.insert(std::make_pair(seg1, std::make_pair(itr->first, seg0)));
	}
	for (auto itr = out.comfirmed_path.begin(); itr != out.comfirmed_path.end(); itr++) {
		comfirmed_path.insert(itr->first);
		seg0.first = std::get<0>(itr->second);
		seg0.second = std::get<2>(itr->second);
		seg1.first = std::get<1>(itr->second);
		seg1.second = std::get<3>(itr->second);
		path_next.insert(std::make_pair(seg0, std::make_pair(itr->first, seg1)));
		path_prev.insert(std::make_pair(seg1, std::make_pair(itr->first, seg0)));
	}
	for (auto itr = out.cut_path.begin(); itr != out.cut_path.end(); itr++) {
		seg0.first = std::get<0>(itr->second);
		seg0.second = std::get<2>(itr->second);
		seg1.first = std::get<1>(itr->second);
		seg1.second = std::get<3>(itr->second);
		path_next.insert(std::make_pair(seg0, std::make_pair(itr->first, seg1)));
		path_prev.insert(std::make_pair(seg1, std::make_pair(itr->first, seg0)));
	}

	std::vector<std::pair<int, std::tuple<int, int, int, int>>> remain_cut_path;
	//std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> change_path;
	std::set<std::pair<std::pair<int, int>, std::pair<int, int>>>change_path;

	for (auto itr = out.cut_path.begin(); itr != out.cut_path.end(); itr++) {
		seg0.first = std::get<0>(itr->second);
		seg0.second = std::get<2>(itr->second);
		seg1.first = std::get<1>(itr->second);
		seg1.second = std::get<3>(itr->second);
		if (path_next.count(seg0) > 1 || path_prev.count(seg1) > 1) {
			remain_cut_path.push_back(*itr);
			continue;
		}
		else {
			change_path.insert(std::make_pair(seg0, seg1));
		}
	}
	int path_id_num = out.num_comfirmed_path + out.num_cut_path + out.num_select_path;
	std::set<std::pair<std::pair<int, int>, std::pair<int, int>> >finished_path;
	std::vector<std::set<std::pair<std::pair<int, int>, std::pair<int, int>>>> path_group;
	bool flg_pass;
	for (auto itr = change_path.begin(); itr != change_path.end(); itr++) {
		if (finished_path.count(*itr) == 1)continue;

		std::set<std::pair<std::pair<int, int>, std::pair<int, int>>> path_cluster;
		std::set<std::pair<std::pair<int, int>, std::pair<int, int>>> path_cluster_add;
		path_cluster.insert(*itr);
		int sum_p = -1;
		while (sum_p != path_cluster.size()) {
			sum_p = path_cluster.size();
			for (auto itr2 = path_cluster.begin(); itr2 != path_cluster.end(); itr2++) {
				if (path_next.count(itr2->second) != 0) {
					auto range = path_next.equal_range(itr2->second);
					for (auto res = range.first; res != range.second; res++) {
						if (change_path.count(std::make_pair(res->first, res->second.second)) == 0)continue;
						path_cluster_add.insert(std::make_pair(res->first, res->second.second));
					}
				}
				if (path_prev.count(itr2->first) != 0) {
					auto  range = path_prev.equal_range(itr2->first);
					for (auto res = range.first; res != range.second; res++) {
						if (change_path.count(std::make_pair(res->second.second, res->first)) == 0)continue;
						path_cluster_add.insert(std::make_pair(res->second.second, res->first));
					}
				}
			}
			for (auto itr2 = path_cluster_add.begin(); itr2 != path_cluster_add.end(); itr2++) {
				path_cluster.insert(*itr2);
			}
			path_cluster_add.clear();
		}

		for (auto itr2 = path_cluster.begin(); itr2 != path_cluster.end(); itr2++) {
			finished_path.insert(*itr2);
		}
		path_group.push_back(path_cluster);
	}
	for (int i = 0; i < path_group.size(); i++) {
		for (auto itr = path_group[i].begin(); itr != path_group[i].end(); itr++) {
			if (itr == path_group[i].begin()) {
				seg0 = itr->first;
				seg1 = itr->second;
			}
			if (seg0.first > itr->first.first) {
				seg0 = itr->first;
			}
			if (seg1.first < itr->second.first) {
				seg1 = itr->second;
			}
		}
		int path_id0 = 0, path_id1 = 0, path_id;
		auto res_prev = path_prev.find(seg0);
		auto res_next = path_next.find(seg1);
		//端点はid=0
		//select path  id=1
		//comfirmed path  id=2
		//そのほか(カット辺など) id=0

		if (res_prev == path_prev.end())path_id0 = 0;
		else if (select_path.count(res_prev->second.first) == 1)path_id0 = 1;
		else if (comfirmed_path.count(res_prev->second.first) == 1)path_id0 = 2;
		else path_id0 = 0;

		if (res_next == path_next.end())path_id1 = 0;
		else if (select_path.count(res_next->second.first) == 1)path_id1 = 1;
		else if (comfirmed_path.count(res_next->second.first) == 1)path_id1 = 2;
		else path_id1 = 0;
		//各辺に割り振る
		if (path_id0 == 0 && path_id1 == 1) {
			path_id = res_next->second.first;
			for (auto itr = path_group[i].begin(); itr != path_group[i].end(); itr++) {
				out.select_path.push_back(std::make_pair(path_id, std::make_tuple(
					itr->first.first, itr->second.first, itr->first.second, itr->second.second)));
			}
		}
		else if (path_id0 == 1 && path_id1 == 0) {
			path_id = res_prev->second.first;
			for (auto itr = path_group[i].begin(); itr != path_group[i].end(); itr++) {
				out.select_path.push_back(std::make_pair(path_id, std::make_tuple(
					itr->first.first, itr->second.first, itr->first.second, itr->second.second)));
			}
		}
		else if (path_id0 == 0 && path_id1 == 2) {
			path_id = res_next->second.first;
			for (auto itr = path_group[i].begin(); itr != path_group[i].end(); itr++) {
				out.comfirmed_path.push_back(std::make_pair(path_id, std::make_tuple(
					itr->first.first, itr->second.first, itr->first.second, itr->second.second)));
			}
		}
		else if (path_id0 == 2 && path_id1 == 0) {
			path_id = res_prev->second.first;
			for (auto itr = path_group[i].begin(); itr != path_group[i].end(); itr++) {
				out.comfirmed_path.push_back(std::make_pair(path_id, std::make_tuple(
					itr->first.first, itr->second.first, itr->first.second, itr->second.second)));
			}
		}
		else {
			path_id = path_id_num;
			path_id_num++;
			for (auto itr = path_group[i].begin(); itr != path_group[i].end(); itr++) {
				out.comfirmed_path.push_back(std::make_pair(path_id, std::make_tuple(
					itr->first.first, itr->second.first, itr->first.second, itr->second.second)));
			}
		}
	}
	out.cut_path = remain_cut_path;
	out.num_comfirmed_path = out.comfirmed_path.size();
	out.num_select_path = out.select_path.size();
	out.num_cut_path = out.cut_path.size();


	std::map<int, int> id_re_roll;
	for (auto itr = out.comfirmed_path.begin(); itr != out.comfirmed_path.end(); itr++) {
		id_re_roll.insert(std::make_pair(itr->first, 0));
	}
	for (auto itr = out.select_path.begin(); itr != out.select_path.end(); itr++) {
		id_re_roll.insert(std::make_pair(itr->first, 0));
	}
	for (auto itr = out.cut_path.begin(); itr != out.cut_path.end(); itr++) {
		id_re_roll.insert(std::make_pair(itr->first, 0));
	}
	int count = 0;
	for (auto itr = id_re_roll.begin(); itr != id_re_roll.end(); itr++) {
		itr->second = count;
		count++;
	}

	for (auto itr = out.comfirmed_path.begin(); itr != out.comfirmed_path.end(); itr++) {
		itr->first = id_re_roll.at(itr->first);
	}
	for (auto itr = out.select_path.begin(); itr != out.select_path.end(); itr++) {
		itr->first = id_re_roll.at(itr->first);
	}
	for (auto itr = out.cut_path.begin(); itr != out.cut_path.end(); itr++) {
		itr->first = id_re_roll.at(itr->first);
	}



	return out;
}



// hayakawa memo
// output_formatのconfirmed_pathとcut_pathを作り直す関数。
// コードを読む限り、元々のconfirmed_pathはbt0とbt1が1対1で結びつくものを選んでいるが、
// こちらではbt0, bt1, ..., btNまでが奇麗に分岐せず直鎖上に連なる場合を選出しているらしい。
// 色々と無駄が多い気がするが、なぜこんな回りくどい探索をするのかはまだ理解できていない。
void path_id_reroll(output_format& g) {
	boost::unordered_multimap <std::pair<int, int>, std::pair<int, int>> path_next;
	boost::unordered_multimap<std::pair<int, int>, std::pair<int, int>> path_prev;
	std::set<std::pair<int, int>>all_vertex;
	std::pair<int, int>ver0, ver1;

	// hayakawa memo
	// confirmed_pathはbt0とbt1が1対1で結びつく、他のbasetrackと繋がらないlinklet。
	for (auto itr = g.comfirmed_path.begin(); itr != g.comfirmed_path.end(); itr++) {
		ver0.first = std::get<0>(itr->second);// pl0
		ver0.second = std::get<2>(itr->second);// rawid0
		ver1.first = std::get<1>(itr->second);// pl1
		ver1.second = std::get<3>(itr->second);// rawid1
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
			// あるbasetrack（*itr）から繋がる下流側basetrackが0ないし2本以上の場合。
			// つまり、*itrがchainの最下流始点となっているか、下流方向への分岐点になっているとき。
			start.insert(*itr);
			continue;
		}
		auto res = path_prev.find(*itr);
		if (path_next.count(res->second) > 1) {
			// *itr == res.firstから繋がる下流側basetrackが1本のみで、
			// かつ下流側basetrackは*itrを含む多数の上流basetrackと繋がる場合。
			// つまり、*itrの下流にあるbasetrackが上流方向への分岐点になっているとき。
			start.insert(*itr);
			continue;
		}
	}
	int count = 0;
	std::set < std::pair<std::pair<int, int>, std::pair<int, int>>>finished;


	// startから始まる分岐のない直鎖の経路を探索し、confirmed_pathへ追加する。
	// またこの経路内のlinkletについては、探索を終えた扱いとしてfinishedに格納しておく。
	for (auto itr = start.begin(); itr != start.end(); itr++) {
		auto now = *itr;
		count = 0;
		while (path_next.count(now) == 1) {
			// now（!= itrの可能性あり）から上流方向へのlinkletが1本のみの場合に入るループ。
			// nowの上流側basetrack==res->secondを取得。
			auto res = path_next.find(now);
			// res->secondから下流側（now側）へのlinkletが多数ある場合、ループは打ち切る。
			if (path_prev.count(res->second) > 1)break;
			// res->secondからnow側へのlinkletが1本しかない場合、
			// 一直線の経路とみなせるのでconfirmed_pathに追加する。
			// pl0, pl1, rawid0, rawid1、というフォーマット。
			g.comfirmed_path.push_back(std::make_pair(path_id, std::make_tuple(
				res->first.first, res->second.first, res->first.second, res->second.second
			)));
			finished.insert(std::make_pair(res->first, res->second));
			now = res->second;
			count++;
		}
		if (count > 0) path_id++;
	}

	// こちらはconfirmed_pathに格納されなかった、分岐のあるlinkletを探索し、cut_pathに追加している。
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

int count_linklet_max(std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>& connected_part) {
	std::map<std::pair<int, int>, int> link_num;
	for (auto itr = connected_part.begin(); itr != connected_part.end(); itr++) {
		//printf("(%d,%d) - (%d,%d)\n", itr->first.first, itr->first.second, itr->second.first, itr->second.second);
		int pl0 = itr->first.first / 10;
		int pl1 = itr->second.first / 10;
		if (pl1 - pl0 >= 2)continue;
		auto res = link_num.insert(std::make_pair(std::make_pair(pl0, pl1), 1));
		if (!res.second) {
			res.first->second += 1;
		}
	}
	int max_num = 0;
	for (auto itr = link_num.begin(); itr != link_num.end(); itr++) {
		max_num = std::max(max_num, itr->second);
		//printf("%d %d : %d\n", itr->first.first, itr->first.second, itr->second);
	}
	return max_num;
}