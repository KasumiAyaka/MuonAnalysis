#pragma once
#include <vector>
#include <set>
#include <map>
#include <fstream>
#include <iostream>
#include <sstream>


void output_dot(std::string filename, std::set<int> &vertex0, std::set<int> &vertex1, std::multimap<int, int> path_next, std::multimap<int, int> path_prev);
void output_dot_bipartite_graph(std::string filename, std::set<int> &vertex0, std::set<int> &vertex1, std::multimap<int, int> path_next, std::multimap<int, int> path_prev);

std::vector<std::vector<std::pair<int, int>>> Enumeration(std::vector<std::pair<int,int>>&ALL_PATH);

std::multimap<int, int> make_matching(std::multimap<int, int>&path, std::set<int> &vertex0);
void path_revers(std::multimap<int, int>&path, std::multimap<int, int>&del_list, std::multimap<int, int>&path_sel, std::multimap<int, int>&path_inv);
std::multimap<int, int> path_merge(std::multimap<int, int>&path_next, std::multimap<int, int>&path_prev);

void search_path_dfs(std::vector<int>&path, int now, std::set<int> &goal, std::set<int>&seen, std::multimap<int, int>&all_path, bool &flg);
void path_increased(std::vector<int>&path, std::multimap<int, int>&path_next, std::multimap<int, int>&path_prev);
void Set_path_endpoint(std::set<int>& start, std::set<int> &goal, std::set<int>& vertex0, std::set<int> vertex1, std::multimap<int, int>&path_next, std::multimap<int, int>&path_prev);

std::set<int> Set_all_vertex(std::set<int>&vertex0, std::set<int>&vertex1);
void dfs_sort_postorder(std::vector<int>&path, int now, std::set<int>&seen, std::multimap<int, int>&all_path, std::vector<int>&label);
std::multimap<int, int> path_reverse(std::multimap<int, int>&path);
void dfs_SCC(std::vector<int>&path, int now, std::set<int>&seen, std::multimap<int, int>&all_path, std::set<int>&hist);
void dfs_cycle(std::vector<int>&path, int now, int start, std::set<int>&seen, std::multimap<int, int>&all_path, bool &flg);
std::vector<int> Get_nouse_edge(std::set<int>&all_vertex, std::multimap<int, int>&path);
std::vector<std::vector<int>> enumarate_add_path(std::vector<int>&path_size);
std::vector<std::pair<int, int>> Get_path(int edge, std::multimap<int, int>&path);

class Maximum_matching {
public:
	std::multimap<int, int> path_next;
	std::multimap<int, int> path_prev;
	std::vector<std::pair<int, int>> path_all;

	////強連結成分分解
////閉路の列挙
////始点-終点を決める
////深さ優先探索でpathの列挙

	Maximum_matching(std::multimap<int, int>&all_path, std::multimap<int, int>&maximum_matching);
	Maximum_matching(std::multimap<int, int>&maximum_matching);
	std::set<int> Get_all_vertex();
	std::vector<int>  Divide_SCC();
	Maximum_matching Dismantling_Graph(Maximum_matching &divided_graph, bool &flg);

};
Maximum_matching::Maximum_matching(std::multimap<int, int>&all_path, std::multimap<int, int>&maximum_matching) {
	path_prev = maximum_matching;
	std::set<int>use_vector;
	for (auto itr = path_prev.begin(); itr != path_prev.end(); itr++) {
		use_vector.insert(itr->first);
		use_vector.insert(itr->second);
	}
	for (auto itr = all_path.begin(); itr != all_path.end(); itr++) {
		if (use_vector.count(itr->first) != 1)continue;
		if (use_vector.count(itr->second) != 1)continue;
		int count = path_prev.count(itr->second);
		bool flg = false;
		if (count > 0) {
			auto range = path_prev.equal_range(itr->second);
			for (auto res = range.first; res != range.second; res++) {
				if (itr->first == res->second&&itr->second == res->first) {
					flg = true;
					break;
				}

			}
		}
		if (flg)continue;
		path_next.insert(*itr);
	}
}
Maximum_matching::Maximum_matching(std::multimap<int, int>&maximum_matching) {
	path_prev = maximum_matching;
	path_next.clear();
}
std::set<int> Maximum_matching::Get_all_vertex() {
	std::set<int> ret;
	for (auto itr = path_next.begin(); itr != path_next.end(); itr++) {
		ret.insert(itr->first);
		ret.insert(itr->second);
	}
	for (auto itr = path_prev.begin(); itr != path_prev.end(); itr++) {
		ret.insert(itr->first);
		ret.insert(itr->second);
	}
	return ret;
}

std::vector<int>  Maximum_matching::Divide_SCC() {
	std::vector<int> path;
	std::set<int>seen;

	//merged_pathを使って帰りがけにsort
	std::vector<int> postorder;
	std::set<int> all_vertex = Get_all_vertex();
	seen.clear();
	path.clear();
	auto merged_path = path_merge(path_next, path_prev);

	for (auto itr = all_vertex.begin(); itr != all_vertex.end(); itr++) {
		if (seen.count(*itr) == 1)continue;
		path.clear();
		int now = *itr;
		dfs_sort_postorder(path, now, seen, merged_path, postorder);
	}
	seen.clear();
	path.clear();

	std::multimap<int, int> path_next_inv = path_reverse(path_next);
	std::multimap<int, int> path_prev_inv = path_reverse(path_prev);
	merged_path = path_merge(path_next_inv, path_prev_inv);
	//帰りがけ順の逆順でdfs
	std::set<int>hist;
	for (auto itr = postorder.rbegin(); itr != postorder.rend(); itr++) {
		if (seen.count(*itr) == 1)continue;
		path.clear();
		hist.clear();
		int now = *itr;
		dfs_SCC(path, now, seen, merged_path, hist);
		if (hist.size() < 2)continue;
		path.clear();
		std::set<int>seen2;
		bool flg = true;
		dfs_cycle(path, now, now, seen2, merged_path, flg);
		//帰りがけ順なので逆順にしておく
		std::reverse(std::begin(path), std::end(path));
		return path;
		printf("Strongly Connected Component\n");
		for (auto itr2 = hist.begin(); itr2 != hist.end(); itr2++) {
			if (std::next(itr2, 1) != hist.end()) {
				printf("%d <-->", *itr2);
			}
			else {
				printf("%d \n", *itr2);
			}
		}
		for (auto itr2 = path.begin(); itr2 != path.end(); itr2++) {
			if (std::next(itr2, 1) != path.end()) {
				printf("%d -->", *itr2);
			}
			else {
				printf("%d \n", *itr2);
			}
		}
	}
	hist.clear();
	path.clear();
	return path;
}
/*
Maximum_matching Maximum_matching::Dismantling_Graph(Maximum_matching &divided_graph, bool &flg) {
	std::vector<int> cycle = Divide_SCC();
	if (cycle.size() < 2) {
		//printf("fin");
		flg = false;
		return *this;
	}

	std::pair<int, int> selected_edge = std::make_pair(cycle[0], cycle[1]);
	Maximum_matching M1(path_next, path_prev);
	Maximum_matching M2(path_next, path_prev);

	//edgeから出る辺を消す
	for (auto itr = M1.path_next.begin(); itr != M1.path_next.end();) {
		if (itr->first == selected_edge.first ||
			itr->first == selected_edge.second ||
			itr->second == selected_edge.first ||
			itr->second == selected_edge.second) {
			itr = M1.path_next.erase(itr);
		}
		else {
			itr++;
		}
	}

	//copyのほうはpathの入れかえ
	//edgeを消す
	//cycle[even]-->cycle[odd]  : prevからnextに
	//cycle[odd] -->cycle[even] : nextからprevに
	std::set<std::pair<int, int>> next_to_prev, prev_to_next;
	//i=0はselected edge
	for (int i = 1; i < cycle.size(); i++) {
		if (i + 1 == cycle.size()) {
			next_to_prev.insert(std::make_pair(cycle[i], cycle[0]));
		}
		else if (i % 2 == 1) {
			next_to_prev.insert(std::make_pair(cycle[i], cycle[i + 1]));
		}
		else {
			prev_to_next.insert(std::make_pair(cycle[i], cycle[i + 1]));
		}
	}
	for (auto itr = M2.path_next.begin(); itr != M2.path_next.end();) {
		if (next_to_prev.count(*itr) == 1) {
			M2.path_prev.insert(std::make_pair(itr->second, itr->first));
			itr = M2.path_next.erase(itr);
		}
		else {
			itr++;
		}
	}
	for (auto itr = M2.path_prev.begin(); itr != M2.path_prev.end();) {
		if (prev_to_next.count(*itr) == 1 || (itr->first == selected_edge.first&&itr->second == selected_edge.second)) {
			if (itr->first == selected_edge.first&&itr->second == selected_edge.second) {
				itr = M2.path_prev.erase(itr);
			}
			else {
				M2.path_next.insert(std::make_pair(itr->second, itr->first));
				itr = M2.path_prev.erase(itr);
			}
		}
		else {
			itr++;
		}
	}
	divided_graph = M2;

	return M1;
}
*/
Maximum_matching Maximum_matching::Dismantling_Graph(Maximum_matching &divided_graph, bool &flg) {
	std::vector<int> cycle = Divide_SCC();
	if (cycle.size() < 2) {
		//printf("fin");
		flg = false;
		return *this;
	}
	std::pair<int, int> selected_edge;
	bool edge_select = false;
	int edge_index = 0;
	for (auto itr = path_prev.begin(); itr != path_prev.end(); itr++) {
		for (int i = 0; i < cycle.size(); i++) {
			if (i + 1 == cycle.size()) {
				selected_edge.first = cycle[i];
				selected_edge.second = cycle[0];
			}
			else {
				selected_edge.first = cycle[i];
				selected_edge.second = cycle[i + 1];
			}
			if (itr->first == selected_edge.first&&itr->second == selected_edge.second) {
				edge_select = true;
				edge_index = i;
				break;
			}
			if (itr->first == selected_edge.second&&itr->second == selected_edge.first) {
				std::swap(selected_edge.first, selected_edge.second);
				edge_select = true;
				edge_index = i;
				break;
			}
		}
		if (edge_select)break;
	}
	if (!edge_select) {
		fprintf(stderr, "deleate edge cannot select\n");
		exit(1);
	}

	Maximum_matching M1(path_next, path_prev);
	Maximum_matching M2(path_next, path_prev);

	//edgeから出る辺を消す
	for (auto itr = M1.path_next.begin(); itr != M1.path_next.end();) {
		if (itr->first == selected_edge.first ||
			itr->first == selected_edge.second ||
			itr->second == selected_edge.first ||
			itr->second == selected_edge.second) {
			itr = M1.path_next.erase(itr);
		}
		else {
			itr++;
		}
	}

	//copyのほうはpathの入れかえ
	//edgeを消す
	for (auto itr = M2.path_prev.begin(); itr != M2.path_prev.end();) {
		if (itr->first == selected_edge.first&&itr->second == selected_edge.second) {
			itr = M2.path_prev.erase(itr);
		}
		else {
			itr++;
		}
	}
	//cycle[even]-->cycle[odd]  : prevからnextに
	//cycle[odd] -->cycle[even] : nextからprevに
	std::set<std::pair<int, int>> next_to_prev, prev_to_next;
	//i=edge_indexはselected edge
	if (edge_index % 2 == 0) {
		for (int i = 0; i < cycle.size(); i++) {
			if (i + 1 == cycle.size()) {
				next_to_prev.insert(std::make_pair(cycle[i], cycle[0]));
			}
			else if (i % 2 == 1) {
				next_to_prev.insert(std::make_pair(cycle[i], cycle[i + 1]));
			}
			else {
				prev_to_next.insert(std::make_pair(cycle[i], cycle[i + 1]));
			}
		}
	}
	else {
		for (int i = 0; i < cycle.size(); i++) {
			if (i + 1 == cycle.size()) {
				prev_to_next.insert(std::make_pair(cycle[i], cycle[0]));
			}
			else if (i % 2 == 1) {
				prev_to_next.insert(std::make_pair(cycle[i], cycle[i + 1]));
			}
			else {
				next_to_prev.insert(std::make_pair(cycle[i], cycle[i + 1]));
			}
		}
	}
	for (auto itr = M2.path_next.begin(); itr != M2.path_next.end();) {
		if (next_to_prev.count(*itr) == 1) {
			M2.path_prev.insert(std::make_pair(itr->second, itr->first));
			itr = M2.path_next.erase(itr);
		}
		else {
			itr++;
		}
	}
	for (auto itr = M2.path_prev.begin(); itr != M2.path_prev.end();) {
		if (prev_to_next.count(*itr) == 1 || (itr->first == selected_edge.first&&itr->second == selected_edge.second)) {
			if (itr->first == selected_edge.first&&itr->second == selected_edge.second) {
				itr = M2.path_prev.erase(itr);
			}
			else {
				M2.path_next.insert(std::make_pair(itr->second, itr->first));
				itr = M2.path_prev.erase(itr);
			}
		}
		else {
			itr++;
		}
	}
	divided_graph = M2;

	return M1;
}



void add_path_all(int loop_num, std::vector<std::vector<std::pair<int, int>>> &all_add_path, Maximum_matching &m_add, std::vector<Maximum_matching>&m_all);
bool judge_same_matching(Maximum_matching &m0, Maximum_matching &m1);
void Set_veterx(std::set<int>&vertex0, std::set<int>&vertex1, std::multimap<int, int>&path);

std::vector<std::vector<std::pair<int, int>>> Enumeration(std::vector<std::pair<int, int>>&ALL_PATH){
	//入力
	std::multimap<int, int> path_all_next;
	std::multimap<int, int> path_next;
	std::multimap<int, int> path_prev;
	for (int i = 0; i < ALL_PATH.size(); i++) {
		path_next.insert(ALL_PATH[i]);
	}
	std::set<int> vertex0, vertex1;
	path_all_next = path_next;
	Set_veterx(vertex0, vertex1, path_all_next);

	///////最大マッチングの検出

	//matching作成
	std::multimap<int, int>matching;
	matching = make_matching(path_next, vertex0);

	//matchingした辺について上下逆転
	path_revers(path_all_next, matching, path_next, path_prev);

	//ある辺についてpathを書く
	//深さ優先探索
	//matchinの増加があるか確認
	//あれば増加
	//無ければ次のpathに
	//pathの両端を見る
	std::vector<int> path;
	int now = 0;
	std::set<int>seen, start, goal;
	bool flg = true;
	Set_path_endpoint(start, goal, vertex0, vertex1, path_next, path_prev);
	auto merged_path = path_merge(path_next, path_prev);
	for (auto itr = start.begin(); itr != start.end(); ) {
		//printf("%d start\n", *itr);
		flg = true;
		now = *itr;
		seen.clear();
		path.clear();
		search_path_dfs(path, now, goal, seen, merged_path, flg);
		if (path.size() != 0&&goal.count(*path.rbegin()) == 1) {
			path_increased(path, path_next, path_prev);
			Set_path_endpoint(start, goal, vertex0, vertex1, path_next, path_prev);
			merged_path = path_merge(path_next, path_prev);

			itr = start.begin();
		}
		else {
			itr++;
		}
	}

	std::vector<Maximum_matching> M_all;
	Maximum_matching M(path_all_next, path_prev);
	M_all.push_back(M);
	for (int i = 0; i < M_all.size(); ) {
		flg = true;
		M_all[i] = M_all[i].Dismantling_Graph(M, flg);
		if (flg) {
			M_all.push_back(M);
		}
		else {
			i++;
		}
	}
	//printf("all maximum matching num %d\n", M_all.size());


	std::set<int>all_vertex;
	for (auto itr = vertex0.begin(); itr != vertex0.end(); itr++) {
		all_vertex.insert(*itr);
	}
	for (auto itr = vertex1.begin(); itr != vertex1.end(); itr++) {
		all_vertex.insert(*itr);
	}

	std::vector<Maximum_matching> output_all;
	for (int i = 0; i < M_all.size(); i++) {
		//使っていないedgeの抽出
		std::vector<int> nouse_edge = Get_nouse_edge(all_vertex, M_all[i].path_prev);
		//nouse edgeから出るpathを全て出す
		std::vector<std::vector<std::pair<int, int>>> all_add_path;
		std::vector<int> all_add_path_size;

		for (int j = 0; j < nouse_edge.size(); j++) {
			std::vector<std::pair<int, int>>add_path = Get_path(nouse_edge[j], path_all_next);
			all_add_path.push_back(add_path);
			all_add_path_size.push_back(add_path.size());
		}
		int loop_num = 0;
		Maximum_matching m_add(M_all[i].path_prev);
		std::vector<Maximum_matching>m_all;
		add_path_all(loop_num, all_add_path, m_add, m_all);
		for (int i = 0; i < m_all.size(); i++) {
			for (auto itr = m_all[i].path_next.begin(); itr != m_all[i].path_next.end(); itr++) {
				m_all[i].path_all.push_back(std::make_pair(itr->first, itr->second));
			}
			for (auto itr = m_all[i].path_prev.begin(); itr != m_all[i].path_prev.end(); itr++) {
				m_all[i].path_all.push_back(std::make_pair(itr->second, itr->first));
			}
			std::sort(m_all[i].path_all.begin(), m_all[i].path_all.end());
			bool flg = false;
			for (int j = 0; j < output_all.size(); j++) {
				if (output_all[j].path_all == m_all[i].path_all) {
					flg = true;
					break;
				}
			}
			if (!flg) output_all.push_back(m_all[i]);
		}
	}

	std::vector<std::vector<std::pair<int, int>>> ret;

	for (int i = 0; i < output_all.size(); i++) {
		//std::vector<std::pair<int, int>> path_out;
		//for (auto itr = output_all[i].path_next.begin(); itr != output_all[i].path_next.end(); itr++) {
		//	path_out.push_back(std::make_pair(itr->first, itr->second));
		//}
		//for (auto itr = output_all[i].path_prev.begin(); itr != output_all[i].path_prev.end(); itr++) {
		//	path_out.push_back(std::make_pair(itr->second, itr->first));
		//}
		ret.push_back(output_all[i].path_all);
	}
	return ret;

}
void Set_veterx(std::set<int>&vertex0, std::set<int>&vertex1, std::multimap<int, int>&path) {
	vertex0.clear();
	vertex1.clear();
	for (auto itr = path.begin(); itr != path.end(); itr++) {
		vertex0.insert(itr->first);
		vertex1.insert(itr->second);
	}
}
std::multimap<int, int> make_matching(std::multimap<int, int>&path, std::set<int> &vertex0) {
	std::multimap<int, int> matching_path;
	std::set<int> used_edge;
	for (auto itr = vertex0.begin(); itr != vertex0.end(); itr++) {
		if (path.count(*itr) == 0)continue;
		auto range = path.equal_range(*itr);
		for (auto res = range.first; res != range.second; res++) {
			if (used_edge.count(res->second) == 0) {
				used_edge.insert(res->second);
				matching_path.insert(std::make_pair(*itr, res->second));
				break;
			}

		}
	}
	return matching_path;
}
void path_revers(std::multimap<int, int>&path, std::multimap<int, int>&del_list, std::multimap<int, int>&path_sel, std::multimap<int, int>&path_inv) {
	path_sel.clear();
	path_inv.clear();
	std::multimap<int, int> ret;
	for (auto itr = path.begin(); itr != path.end(); itr++) {
		if (del_list.count(itr->first) == 0) {
			path_sel.insert(*itr);
		}
		else {
			bool flg = false;
			auto range = del_list.equal_range(itr->first);
			for (auto itr2 = range.first; itr2 != range.second; itr2++) {
				if (*itr == *itr2)flg = true;
			}
			if (!flg) {
				path_sel.insert(*itr);
			}
		}
	}
	for (auto itr = del_list.begin(); itr != del_list.end(); itr++) {
		path_inv.insert(std::make_pair(itr->second, itr->first));
	}
}

///////最大マッチング検出///////////////
void Print_hist(std::vector<int>&path) {
	for (int i = 0; i < path.size(); i++) {
		if (path.size() != i + 1) {
			printf("%d-->", path[i]);
		}
		else {
			printf("%d\n", path[i]);
		}
	}
}
void search_path_dfs(std::vector<int>&path, int now, std::set<int> &goal, std::set<int>&seen, std::multimap<int, int>&all_path, bool &flg) {
	for (int i = 0; i < path.size(); i++) {
		//loop検出
		if (path[i] == now) {
			return;
		}
	}
	path.push_back(now);
	if (goal.count(now) == 1) {
		//printf("hit\n");
		//Print_hist(path);
		flg = false;
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
		search_path_dfs(path, res->second, goal, seen, all_path, flg);
		if (!flg)return;
	}
	seen.insert(now);
	path.pop_back();

}
std::multimap<int, int> path_merge(std::multimap<int, int>&path_next, std::multimap<int, int>&path_prev) {
	std::multimap<int, int> ret;
	for (auto itr = path_next.begin(); itr != path_next.end(); itr++) {
		ret.insert(*itr);
	}
	for (auto itr = path_prev.begin(); itr != path_prev.end(); itr++) {
		ret.insert(*itr);
	}
	return ret;
}
void Set_path_endpoint(std::set<int>& start, std::set<int>& goal, std::set<int>& vertex0, std::set<int> vertex1, std::multimap<int, int>&path_next, std::multimap<int, int>&path_prev) {
	start.clear();
	goal.clear();
	bool flg = true;
	for (auto itr = vertex0.begin(); itr != vertex0.end(); itr++) {
		flg = true;
		for (auto itr2 = path_prev.begin(); itr2 != path_prev.end(); itr2++) {
			if (itr2->second == *itr) {
				flg = false;
				break;
			}
		}
		if (flg)start.insert(*itr);
	}
	for (auto itr = vertex1.begin(); itr != vertex1.end(); itr++) {
		flg = true;
		for (auto itr2 = path_prev.begin(); itr2 != path_prev.end(); itr2++) {
			if (itr2->first == *itr) {
				flg = false;
				break;
			}
		}
		if (flg)goal.insert(*itr);
	}
}
void path_increased(std::vector<int>&path, std::multimap<int, int>&path_next, std::multimap<int, int>&path_prev) {
	std::pair<int, int> edge_pair;
	for (int i = 0; i < path.size(); i++) {
		//path[i]-->path[i+1]
		if (i + 1 == path.size())continue;
		edge_pair.first = path[i];
		edge_pair.second = path[i + 1];
		if (i % 2 == 0) {
			//非マッチング-->マッチングへ変換
			if (path_next.count(edge_pair.first) == 0) {
				fprintf(stderr, "error funtion[path_increased]\n");
				exit(1);
			}
			auto range = path_next.equal_range(edge_pair.first);
			for (auto res = range.first; res != range.second; ) {
				if (res->first == edge_pair.first&&res->second == edge_pair.second) {
					res = path_next.erase(res);
				}
				else {
					res++;
				}
			}
			path_prev.insert(std::make_pair(edge_pair.second, edge_pair.first));
		}
		else {
			//マッチング-->非マッチングへ変換
			if (path_prev.count(edge_pair.first) == 0) {
				fprintf(stderr, "error funtion[path_increased]\n");
				exit(1);
			}
			auto range = path_prev.equal_range(edge_pair.first);
			for (auto res = range.first; res != range.second; ) {
				if (res->first == edge_pair.first&&res->second == edge_pair.second) {
					res = path_prev.erase(res);
				}
				else {
					res++;
				}
			}
			path_next.insert(std::make_pair(edge_pair.second, edge_pair.first));

		}
	}
}
///////////////////////////////////////

///////強連結成分分解///////////////
void dfs_sort_postorder(std::vector<int>&path, int now, std::set<int>&seen, std::multimap<int, int>&all_path, std::vector<int>&label) {
	path.push_back(now);
	seen.insert(now);
	//Print_hist(path);
	if (all_path.count(now) == 0) {
		path.pop_back();
		label.push_back(now);
		return;
	}
	auto range = all_path.equal_range(now);
	for (auto res = range.first; res != range.second; res++) {
		if (seen.count(res->second) == 1)continue;
		dfs_sort_postorder(path, res->second, seen, all_path, label);
	}
	path.pop_back();
	label.push_back(now);
	return;

}
void dfs_SCC(std::vector<int>&path, int now, std::set<int>&seen, std::multimap<int, int>&all_path, std::set<int>&hist) {
	path.push_back(now);
	hist.insert(now);
	seen.insert(now);
	//Print_hist(path);
	if (all_path.count(now) == 0) {
		path.pop_back();
		return;
	}
	auto range = all_path.equal_range(now);
	for (auto res = range.first; res != range.second; res++) {
		if (seen.count(res->second) == 1)continue;
		dfs_SCC(path, res->second, seen, all_path, hist);
	}
	path.pop_back();
	return;

}
void dfs_cycle(std::vector<int>&path, int now, int start, std::set<int>&seen, std::multimap<int, int>&all_path, bool &flg) {
	if (path.size() > 0 && now == start) {
		//printf("hit\n");
		//Print_hist(path);
		flg = false;
		return;
	}

	path.push_back(now);
	seen.insert(now);
	if (all_path.count(now) == 0) {
		path.pop_back();
		return;
	}
	auto range = all_path.equal_range(now);
	for (auto res = range.first; res != range.second; res++) {
		if (res->second != start && seen.count(res->second) == 1)continue;
		dfs_cycle(path, res->second, start, seen, all_path, flg);
		if (!flg)return;
	}
	if (!flg)return;
	path.pop_back();

}
std::set<int> Set_all_vertex(std::set<int>&vertex0, std::set<int>&vertex1) {

	std::set<int> ret;
	for (auto itr = vertex0.begin(); itr != vertex0.end(); itr++) {
		ret.insert(*itr);
	}
	for (auto itr = vertex1.begin(); itr != vertex1.end(); itr++) {
		ret.insert(*itr);
	}
	return ret;
}

std::multimap<int, int> path_reverse(std::multimap<int, int>&path) {
	std::multimap<int, int>  ret;
	for (auto itr = path.begin(); itr != path.end(); itr++) {
		ret.insert(std::make_pair(itr->second, itr->first));
	}
	return ret;
}



/////

std::vector<int> Get_nouse_edge(std::set<int>&all_vertex, std::multimap<int, int>&path) {
	std::vector<int> ret;
	std::set<int>use_vertex;
	for (auto itr = path.begin(); itr != path.end(); itr++) {
		use_vertex.insert(itr->first);
		use_vertex.insert(itr->second);
	}

	for (auto itr = all_vertex.begin(); itr != all_vertex.end(); itr++) {
		if (use_vertex.count(*itr) == 0)ret.push_back(*itr);
	}
	return ret;
}
std::vector<std::pair<int, int>> Get_path(int edge, std::multimap<int, int>&path) {
	std::vector<std::pair<int, int>>ret;
	for (auto itr = path.begin(); itr != path.end(); itr++) {
		if (itr->first == edge || itr->second == edge) {
			ret.push_back(*itr);
		}
	}
	return ret;
}

std::vector<std::vector<int>> enumarate_add_path(std::vector<int>&path_size) {
	std::vector<std::vector<int>> ret;
	int all_path_num = 1;
	for (int i = 0; i < path_size.size(); i++) {
		all_path_num *= path_size[i];
	}
	for (int i = 0; i < all_path_num; i++) {
		std::vector<int> path;
		int now_div = 1;
		for (int j = 0; j < path_size.size(); j++) {
			now_div *= path_size[j];
			path.push_back(i%now_div);
		}
		ret.push_back(path);
	}


	//for (int i = 0; i < path_size.size(); i++) {
	//	if (i + 1 != path_size.size()) {
	//		printf("%d -->", path_size[i]);
	//	}
	//	else {
	//		printf("%d\n", path_size[i]);
	//	}
	//}

	//for (int i = 0; i < ret.size(); i++) {
	//	for (int j = 0; j < ret[i].size(); j++) {
	//		if (j + 1 != ret[i].size()) {
	//			printf("%d -->", ret[i][j]);
	//		}
	//		else {
	//			printf("%d\n", ret[i][j]);
	//		}
	//	}
	//}
	//printf("fin\n");
	return ret;
}

///////出力///////////////

void output_dot(std::string filename, std::set<int> &vertex0, std::set<int> &vertex1, std::multimap<int, int> path_next, std::multimap<int, int> path_prev) {
	std::ofstream ofs(filename);
	ofs << "digraph{" << std::endl;
	for (auto itr = vertex0.begin(); itr != vertex0.end(); itr++) {
		if (std::next(itr, 1) == vertex0.end() && vertex1.size() == 0) {
			ofs << *itr << std::endl;
		}
		else {
			ofs << *itr << ",";
		}
	}
	for (auto itr = vertex1.begin(); itr != vertex1.end(); itr++) {
		if (std::next(itr, 1) == vertex1.end()) {
			ofs << *itr << std::endl;
		}
		else {
			ofs << *itr << ",";
		}
	}
	for (auto itr = path_next.begin(); itr != path_next.end(); itr++) {
		ofs << itr->first << "->" << itr->second << std::endl;
	}
	for (auto itr = path_prev.begin(); itr != path_prev.end(); itr++) {
		ofs << itr->first << "->" << itr->second << std::endl;
	}
	ofs << "}" << std::endl;

}
void output_dot_bipartite_graph(std::string filename, std::set<int> &vertex0, std::set<int> &vertex1, std::multimap<int, int> path_next, std::multimap<int, int> path_prev) {
	std::ofstream ofs(filename);
	ofs << "digraph{" << std::endl;
	ofs << "splines = false" << std::endl;

	//ofs << "graph[style = rounded color = black shape = plaintext rankdir = LR]" << std::endl;
	ofs << "graph[style = rounded color = white shape = plaintext]" << std::endl;

	ofs << "subgraph cluster_source{" << std::endl;
	for (auto itr = vertex0.begin(); itr != vertex0.end(); itr++) {
		ofs << "\tsubgraph cluster_" << *itr << "{" << std::endl;
		//ofs << "\tlabel = \"" << *itr << "\" color = white shape = plaintext;" << std::endl;
		ofs << "\tcolor = white shape = plaintext;" << std::endl;
		ofs << "\t" << *itr << std::endl;
		ofs << "\t}" << std::endl;
	}
	ofs << "}" << std::endl;

	ofs << "subgraph cluster_sink{" << std::endl;
	for (auto itr = vertex1.begin(); itr != vertex1.end(); itr++) {
		ofs << "\tsubgraph cluster_" << *itr << "{" << std::endl;
		//ofs << "\tlabel = \"" << *itr << "\" color = white shape = plaintext;" << std::endl;
		ofs << "\tcolor = white shape = plaintext;" << std::endl;
		ofs << "\t" << *itr << std::endl;
		ofs << "\t}" << std::endl;
	}
	ofs << "}" << std::endl;

	int count;
	for (auto itr = path_next.begin(); itr != path_next.end(); itr++) {
		ofs << itr->first << "->" << itr->second << std::endl;
	}
	for (auto itr = path_prev.begin(); itr != path_prev.end(); itr++) {
		ofs << itr->first << "->" << itr->second << std::endl;
	}
	ofs << "}" << std::endl;

}

///////////////////////////////////////
void add_path_all(int loop_num, std::vector<std::vector<std::pair<int, int>>> &all_add_path, Maximum_matching &m_add, std::vector<Maximum_matching>&m_all) {
	if (all_add_path.size() == 0) {
		m_all.push_back(m_add);
	}
	if (loop_num == all_add_path.size())return;
	for (int i = 0; i < all_add_path[loop_num].size(); i++) {
		Maximum_matching m_add2(m_add.path_prev);
		m_add2.path_next = m_add.path_next;
		m_add2.path_next.insert(all_add_path[loop_num][i]);
		if (loop_num + 1 == all_add_path.size()) {
			m_all.push_back(m_add2);
		}
		add_path_all(loop_num + 1, all_add_path, m_add2, m_all);
	}
}
bool judge_same_matching(Maximum_matching &m0, Maximum_matching &m1) {
	std::set<std::pair<int, int>> path0;
	for (auto itr = m0.path_next.begin(); itr != m0.path_next.end(); itr++) {
		path0.insert(std::make_pair(itr->first, itr->second));
	}
	for (auto itr = m0.path_prev.begin(); itr != m0.path_prev.end(); itr++) {
		path0.insert(std::make_pair(itr->second, itr->first));
	}

	for (auto itr = m1.path_next.begin(); itr != m1.path_next.end(); itr++) {
		if (path0.count(std::make_pair(itr->first, itr->second)) == 0)return false;
	}
	for (auto itr = m1.path_prev.begin(); itr != m1.path_prev.end(); itr++) {
		if (path0.count(std::make_pair(itr->second, itr->first)) == 0)return false;
	}
	return true;



}
